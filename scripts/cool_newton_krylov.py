#!/usr/bin/env python3
"""Newton-Krylov solve of the cool CMFGEN closure with the validated exact Jacobian action.

Exact action (Report 20):  J.(x*z) = J_dir.z + dSTEQ_rad(z)
  * J_dir  fixed-J (direct) Jacobian. Block-diagonal in depth. Native analytic
           blocks (QS_ALI_LIN, QS_LIN_FAC=0) completed by frozen-J central
           differences for the columns the native linearization misses
           (electron density, total ion density, temperature): replay passes
           that reread the converged J of the base pass.
  * dSTEQ_rad  analytic nonlocal formal-transfer tangent (linearized, with
           blue-ray history) mapped through the native dSTEQ/dJ coefficients,
           computed in a replay pass at the same state.

Every native pass starts from the same archived restart x0 and reaches the
iterate through the read-only perturbation hook POPS = POPS(x0)*(1-s_cum), so
the residual definition (forced equation selection, lagged two-photon J,
atomic data, flow, T if fixed) is identical in every pass. Updates are
x <- x*exp(-lam*s) (positivity preserving; nothing is clipped). The merit is
||M_k^-1 R(x)||_2 over active rows with M_k = J_dir(x_k).
"""
from __future__ import annotations

import argparse, json, os, struct, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'backends' / 'cmfgen'))
sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import (BlockTridiag, read_ba_blocks, read_jv_rad,  # noqa: E402
                                read_trace_record, write_directions, write_fortran)
import cnk_trap  # noqa: E402  (Report 23: exact columns of trapped levels, coupling significance, update map)

NT, ND = 134, 84                              # reset from the first base pass
SIG_FRAC = 1e-10                              # set from --sig-frac
TRACE_WEIGHT = 1e-3                           # set from --trace-weight
ULTRA_FRAC = 0.0                              # set from --ultra-frac (0: two-tier weights as before)
ULTRA_WEIGHT = 1e-8                           # set from --ultra-weight
WEIGHT_RAMP = 0.0                             # set from --weight-ramp (decades below SIG_FRAC; 0: step function)
MERIT_MODE = 'precond'                        # set from --merit
EQ_FLOOR = 0.0                                # set from --eq-floor
PC_T_FLOOR = 0.0                              # set from --pc-t-floor
STRUCT_ACTIVE = False                         # set from --struct-active
ION_ROWS = [32] + list(range(102, 133))      # charged species entering POPION (1-based)
SPLIT_BOUNDS_84 = [5, 9, 20, 30, 60]          # depth-region edges (84-radius indices) for --split
SPLIT_REGIONS = []                            # set from ND
PTC_CONSTRAINT_ROWS = (31, 131, 132, 133)     # H number, He number, charge, T rows (0-based): no population PTC
# depth regions (0-based start, end exclusive) of the dense-window least-squares log (Report 26): outer He I, outer H,
# H front, inner H front, between, He front, inner atmosphere, deep layers
DENSE_REPORT_REGIONS = ((0, 9), (9, 18), (18, 32), (32, 47), (47, 54), (54, 61), (61, 100), (100, 167))


def trace_sweep(R, base, x0, dest, cap, k):
    """Exact linear equilibrium of the trace levels (<= SIG_FRAC of their element) at fixed
    significant populations, n_e, T and J: solve the trace-trace block of the native fixed-J
    Jacobian per depth (trace levels do not feed back on J), update x(1 - s) (exact for this
    linear subsystem; capped exp(-s) where s >= 1), and run one complete base pass there.
    Significant variables are unchanged. Trace levels far from equilibrium otherwise leave
    significant equations unbalanced (e.g. He II ground through He III recombination)."""
    x, Rv = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    _, di, _ = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', NT, ND)
    xH = x[1:32].sum(0); xHe = x[32:132].sum(0)
    frac = np.ones_like(x); frac[0:32] = x[0:32] / xH; frac[32:132] = x[32:132] / xHe
    trace = (frac <= SIG_FRAC) & (Rv != 0.0)
    s = np.zeros_like(x)
    for d in range(ND):
        t_ = np.nonzero(trace[:, d])[0]
        if t_.size:
            s[t_, d] = np.linalg.solve(di[d][np.ix_(t_, t_)], Rv[t_, d])
    f = 1.0 - s
    xn = np.where(f > 1e-4, x * f, x * np.exp(-np.minimum(s, cap)))
    R.log(f'it{k:02d} trace sweep at {base.name}: {int(trace.sum())} trace variables, max|s|={np.max(np.abs(s)):.3e}, '
          f'{int(np.sum(f <= 1e-4))} capped decreases')
    return R.run(dest, 1.0 - xn / x0, warm=base / 'QSTARS_RJ.bin')


def lm_step(lm_data, wvec, radius):
    """Levenberg-Marquardt step in the Krylov subspace: y(mu) = argmin ||A y - b||^2 + mu ||W S y||^2,
    with the smallest mu (bisection in log mu) such that max|s_T| <= radius. Returns (s, mu, predicted merit)."""
    S, A, b = lm_data
    Sm = np.stack([si.ravel() for si in S], axis=1)
    WS = Sm * wvec.ravel()[:, None]
    AtA, Atb, G = A.T @ A, A.T @ b, WS.T @ WS
    shape = S[0].shape
    tidx = np.ravel_multi_index((np.full(shape[1], 133), np.arange(shape[1])), shape)

    def solve(mu):
        y = np.linalg.lstsq(AtA + mu * G, Atb, rcond=1e-14)[0] if mu > 0 else np.linalg.lstsq(A, b, rcond=1e-12)[0]
        return y, float(np.max(np.abs(Sm[tidx] @ y)))
    y, mt = solve(0.0)
    if mt <= radius:
        return (Sm @ y).reshape(shape), 0.0, float(np.linalg.norm(A @ y - b))
    scale = max(np.trace(AtA) / max(np.trace(G), 1e-300), 1e-300)
    lo, hi = -14.0, 14.0
    while solve(scale * 10.0 ** hi)[1] > radius and hi < 40:
        hi += 4.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if solve(scale * 10.0 ** mid)[1] > radius:
            lo = mid
        else:
            hi = mid
    mu = scale * 10.0 ** hi
    y, _ = solve(mu)
    return (Sm @ y).reshape(shape), float(mu), float(np.linalg.norm(A @ y - b))


def lm_step_bound(lm_data, wvec, radius, mask):
    """Report 25 (--lm-bound): Levenberg-Marquardt step in the Krylov subspace, y(mu) = argmin ||A y - b||^2 +
    mu ||W S y||^2, with the smallest mu (bisection in log mu) such that max|s| over `mask` <= radius: the step is
    re-optimized within the bound instead of being clipped afterwards (clipping a soft-mode component, H I n=2 at
    d37-38, left the rest of the step inconsistent: full steps x2.4-5.5 worse). Returns (s, mu, y)."""
    S, A, b = lm_data
    Sm = np.stack([si.ravel() for si in S], axis=1)
    WS = Sm * wvec.ravel()[:, None]
    AtA, Atb, G = A.T @ A, A.T @ b, WS.T @ WS
    Smm = Sm[mask.ravel()]
    if np.ndim(radius) > 0:
        # per-component radius (Report 25, stage 2: T bounded tighter than the populations): bound |s| / radius <= 1
        Smm = Smm / np.asarray(radius, dtype=float).ravel()[mask.ravel()][:, None]
        radius = 1.0

    def solve(mu):
        y = np.linalg.lstsq(AtA + mu * G, Atb, rcond=1e-14)[0] if mu > 0 else np.linalg.lstsq(A, b, rcond=1e-12)[0]
        return y, float(np.max(np.abs(Smm @ y)))
    y, m0 = solve(0.0)
    if m0 <= radius:
        return (Sm @ y).reshape(S[0].shape), 0.0, y
    scale = max(np.trace(AtA) / max(np.trace(G), 1e-300), 1e-300)
    lo, hi = -14.0, 14.0
    while solve(scale * 10.0 ** hi)[1] > radius and hi < 40:
        hi += 4.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if solve(scale * 10.0 ** mid)[1] > radius:
            lo = mid
        else:
            hi = mid
    mu = scale * 10.0 ** hi
    y, _ = solve(mu)
    return (Sm @ y).reshape(S[0].shape), float(mu), y


def bounded_step(lm_data, active, radius, mask, solver):
    """V1 plan P1 (--step-solver qr|box): the bounded trial step without forming A^T A. On the stored products of free_r25
    iterates 4 and 16, cond(A) ~ 3e26, so the normal equations of lm_step_bound (cond ~ 1e53) returned a step 100 % away
    from the augmented-QR solution at the same mu, predicting 0.976 / 0.9988 of the merit where the stable LM family
    gives 0.938 / 0.969 and the explicit box 0.87-0.91 (Reports/evidence/v1_exec_20261106/p1).
      qr : y(mu) = argmin ||A y - b||^2 + mu ||S y||^2 (active rows) from the QR of [A; sqrt(mu) R_s], bisection in mu
      box: min ||A y - b|| s.t. |S y| <= radius on mask, as a second-order cone program (clarabel), column-scaled QR of A
    radius: scalar or (NT, ND). Returns (s, mu or nan, y)."""
    import scipy.linalg as sla
    S, A, b = lm_data
    shape = S[0].shape
    Sm = np.stack([si.ravel() for si in S], axis=1)
    rad = np.broadcast_to(np.asarray(radius, dtype=float), shape).ravel()
    m = mask.ravel()
    C = Sm[m] / rad[m][:, None]
    if solver == 'box':
        # the interior point can stop short on these ill-conditioned spans (v1_p2: 0.917 against 0.868 for the LM-QR step,
        # which is feasible for the box): the better feasible step of the two is returned
        s_q, mu_q, y_q = bounded_step(lm_data, active, radius, mask, 'qr')
        import scipy.sparse as sp
        import clarabel
        cn = np.linalg.norm(A, axis=0); cn[cn == 0] = 1.0
        Qa, RA = np.linalg.qr(A / cn)
        c = Qa.T @ b
        n = RA.shape[1]
        Cs = sp.csc_matrix(C / cn[None, :])
        nb = Cs.shape[0]
        z1 = lambda k: sp.csc_matrix((k, 1))  # noqa: E731
        Ac = sp.vstack([sp.hstack([sp.csc_matrix((1, n)), sp.csc_matrix(np.array([[-1.0]]))]),
                        sp.hstack([sp.csc_matrix(-RA), z1(n)]), sp.hstack([Cs, z1(nb)]), sp.hstack([-Cs, z1(nb)])]).tocsc()
        bc = np.concatenate([np.zeros(1), -c, np.ones(2 * nb)])
        q = np.zeros(n + 1); q[-1] = 1.0
        st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 400
        sol = clarabel.DefaultSolver(sp.csc_matrix((n + 1, n + 1)), q, Ac, bc,
                                     [clarabel.SecondOrderConeT(1 + n), clarabel.NonnegativeConeT(2 * nb)], st).solve()
        y = np.array(sol.x[:n]) / cn
        viol = float(np.max(np.abs(C @ y))) if nb else 0.0
        if viol > 1.0:
            y = y / viol                       # an inexact interior-point solution, scaled into the box
        if np.linalg.norm(A @ y_q - b) < np.linalg.norm(A @ y - b):
            return s_q, mu_q, y_q
        return (Sm @ y).reshape(shape), float('nan'), y
    Rs = np.linalg.qr(Sm[active.ravel()], mode='r')

    def solve(mu):
        if mu <= 0:
            return sla.lstsq(A, b, cond=1e-12, lapack_driver='gelsd')[0]
        Qm, Rm = np.linalg.qr(np.vstack([A, np.sqrt(mu) * Rs]))
        return sla.solve_triangular(Rm, Qm[:A.shape[0]].T @ b)

    def ratio(y):
        return float(np.max(np.abs(C @ y))) if C.shape[0] else 0.0
    y = solve(0.0)
    if ratio(y) <= 1.0:
        return (Sm @ y).reshape(shape), 0.0, y
    scale = max(float(np.sum(A * A)) / max(float(np.sum(Rs * Rs)), 1e-300), 1e-300)
    lo, hi = -14.0, 14.0
    while ratio(solve(scale * 10.0 ** hi)) > 1.0 and hi < 40:
        hi += 4.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if ratio(solve(scale * 10.0 ** mid)) > 1.0:
            lo = mid
        else:
            hi = mid
    mu = scale * 10.0 ** hi
    y = solve(mu)
    return (Sm @ y).reshape(shape), float(mu), y


def span_basis(S, A, tau):
    """V1 plan P2 (--span-tau): orthonormal basis of the span of the directions S (each normalized), truncated at relative
    singular value tau, with its products: U_k = Z_n V_k / s_k, A U_k = A_n V_k / s_k. The least squares over the raw
    directions of v1_p2 (238 directions, 185 independent) returned coefficients up to 4e14 cancelling to |s| = 0.4,
    which turned the products' own small errors into a model error of 0.47 of the merit (sum of stored products against a
    fresh exact product of the step); in the truncated basis the predicted fit equals the fresh-product fit to 1e-4.
    Returns (basis list, A of the basis, M with y = M c, k)."""
    shape = S[0].shape
    Zm = np.stack([z.ravel() for z in S], axis=1)
    nz = np.linalg.norm(Zm, axis=0); nz[nz == 0] = 1.0
    U, sv, Vt = np.linalg.svd(Zm / nz, full_matrices=False)
    k = int(np.sum(sv > tau * sv[0]))
    M = (Vt[:k].T / sv[:k]) / nz[:, None]
    return [U[:, i].reshape(shape) for i in range(k)], A @ M, M, k


def span_lstsq(S, A, b, tau):
    """Least-squares coefficients over the directions S: plain (tau = 0, as through Report 26) or through span_basis."""
    if tau <= 0:
        return np.linalg.lstsq(A, b, rcond=1e-12)[0]
    _, Ab, M, _ = span_basis(S, A, tau)
    return M @ np.linalg.lstsq(Ab, b, rcond=1e-12)[0]


def split_groups(free_t):
    """Species groups (0-based row slices): H, He I, He II/He III/n_e (+ T)."""
    return [(0, 32), (32, 101), (101, 133)] + ([(133, 134)] if free_t else [])


def fine_regions(regions, d0, d1, width):
    """Report 25 (--split-fine): the depth regions with depths d0..d1-1 (0-based) cut out and split into regions of
    `width` depths; the rest of the regions as before (clipped). An exact partition of the depths."""
    out = []
    for a0, a1 in regions:
        if a1 <= d0 or a0 >= d1:
            out.append((a0, a1))
            continue
        if a0 < d0:
            out.append((a0, d0))
        if a1 > d1:
            out.append((d1, a1))
    out += [(e, min(e + width, d1)) for e in range(d0, d1, width)]
    return sorted(out)


def block_pieces(v, blocks, tol):
    """Exact partition of v into the pieces v[rows, d0:d1] of disjoint blocks (rows: index array or slice); pieces
    with norm <= tol*||v|| and everything outside the blocks form one remainder piece, so sum(pieces) == v."""
    out, rest = [], v.copy()
    nv = np.linalg.norm(v)
    for rows, d0, d1 in blocks:
        p = np.zeros_like(v)
        p[rows, d0:d1] = v[rows, d0:d1]
        if np.linalg.norm(p) > tol * nv:
            out.append(p)
            rest[rows, d0:d1] = 0.0
    if np.linalg.norm(rest) > 0:
        out.append(rest)
    return out


def split_pieces(v, regions, groups, tol):
    """Exact partition of v into depth-region x species-group pieces; pieces with norm
    <= tol*||v|| are merged into one remainder piece, so sum(pieces) == v."""
    return block_pieces(v, [(slice(r0, r1), d0, d1) for d0, d1 in regions for r0, r1 in groups], tol)


def front_split_pieces(v, regions, groups, tol, window):
    """Report 25 (--front-pieces): split_pieces with the depths of the hydrogen-front window (d0, d1; 0-based,
    inclusive) cut out of the regions. Inside the window H (with H-) and n_e form one piece per depth, every other
    group one piece over the window. The slow H I n = 2 + 3 / ionization mode of the front blocks is nearly singular
    locally and coupled radiatively beyond the neighbour depths; region-wide pieces gave the least squares one
    amplitude for it over 20 depths, and it fitted the front with deep directions instead (Report 24 sec. 5)."""
    w0, w1 = window[0], window[1] + 1
    regs = []
    for d0, d1 in regions:
        if d1 <= w0 or d0 >= w1:
            regs.append((d0, d1))
            continue
        if d0 < w0:
            regs.append((d0, w0))
        if d1 > w1:
            regs.append((w1, d1))
    blocks = [(slice(r0, r1), d0, d1) for d0, d1 in regs for r0, r1 in groups]
    fine = np.r_[0:32, 132]                    # H-, H I, H II and n_e
    blocks += [(fine, d, d + 1) for d in range(w0, w1)]
    for r0, r1 in groups:
        rows = np.setdiff1d(np.arange(r0, r1), fine)
        if rows.size:
            blocks.append((rows, w0, w1))
    return block_pieces(v, blocks, tol)


def split_regions(nd):
    if nd == 84:
        edges = SPLIT_BOUNDS_84
    elif nd == 167:                           # nested grid: 84-radius depth k is 167-radius depth 2k-1
        edges = [2 * b - 1 for b in SPLIT_BOUNDS_84]
    elif nd == 333:                           # nested grid: 167-radius depth k is 333-radius depth 2k-1
        edges = [4 * b - 3 for b in SPLIT_BOUNDS_84]
    else:
        edges = [round(b * nd / 84) for b in SPLIT_BOUNDS_84]
    e = [0] + edges + [nd]
    return list(zip(e[:-1], e[1:]))


def ali_stats(path, skip=0):
    """QSTARS_ALI_STATS records written by a pass ('NCF n tol t max_it n unconverged n ... max_res r ...' per J
    computation); the first `skip` lines are inherited from the restart copy (the file is opened APPEND)."""
    out = []
    for line in Path(path).read_text().splitlines()[skip:]:
        f = line.split()
        out.append({k: float(v) for k, v in zip(f[0::2], f[1::2])})
    return out


def check_one_record(path):
    """A pass's residual trace must hold exactly ONE record: CMFGEN appends, every reader takes record 0, so a
    record inherited from the restart directory would silently replace the pass's own state (Report 23, 333 grid)."""
    data = Path(path).read_bytes()[:12]
    _, nt, nd = struct.unpack('<3i', data)
    if Path(path).stat().st_size != 12 + 16 * nt * nd:
        raise RuntimeError(f'{path}: {Path(path).stat().st_size} bytes is not one ({nt}, {nd}) record; stale or appended trace')


def release_t(dest):
    """FIX_T=F alone does not release T on a restart: READ_BA_DATA_V3 reads FIXED_T back
    from the restart's BAMATPNT ('Temperature was held fixed'), and FIXPOP_IN_BA_V3 then
    zeroes the radiative-equilibrium row at all depths. Removing the stored BA pointer
    (every pass computes its BA) leaves FIXED_T at its initial .FALSE."""
    for name in ('BAMATPNT', 'BAMAT'):
        (Path(dest) / name).unlink(missing_ok=True)
    (Path(dest) / 'QSTARS_T_RELEASED').write_text('BAMATPNT/BAMAT removed so FIXED_T is not inherited\n')


class Runner:
    def __init__(self, x0dir, exe, replace, fix_t, threads, logf, vadat=()):
        self.vadat = list(vadat)
        self.x0dir = Path(x0dir).resolve(); self.exe = Path(exe).resolve()
        self.replace = Path(replace).resolve(); self.fix_t = fix_t
        self.threads = threads; self.logf = logf
        self.extra_env = {}                 # recorded per pass in run_manifest.json (extra_env)
        self.slim = True                    # compact completed passes (scripts/slim_pass.py)
        # Report 22: skip the output-only observer-frame formal solution in solver passes (build >= 29;
        # ignored by older builds). The final/verification pass (verify_cool_state.py) does not use Runner.
        self.no_obs = True

    def log(self, msg):
        line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
        print(line, flush=True)
        with open(self.logf, 'a') as f:
            f.write(line + '\n')

    def run(self, dest, s_cum, *, replay=None, dirs=None, ba=True, store_j=True, threads=None, timeout=60000,
            warm=None, extra_env=None, keep_export=False, two_phot_j=None, lin_fac=0.0, _retried=False):
        """two_phot_j (build >= 32, Report 23): a QSTARS_RJ.bin of this radius/frequency grid. The pass then keeps
        the two-photon method TWO_METH=USE_RAD with the conjugate-frequency J interpolated from that file (fixed,
        lagged data: env QSTARS_TWO_PHOT_J) instead of the OLD_DEFAULT substitution that COMP_F=T forces. Every
        pass that must share one residual definition (base, frozen-J columns, products, trials) needs the same
        file. None (default): unchanged behaviour. A QSTARS_TWO_PHOT_J entry of extra_env (runner-wide or of this
        call) is the same request: it gets the same checks and VADAT setting and must agree with two_phot_j."""
        from formal_ali_trial import prepare
        from run_component import run
        dest = Path(dest)
        tpj = str(Path(two_phot_j).resolve()) if two_phot_j is not None else None
        # The variable given through extra_env must neither bypass the guards below nor replace the requested
        # file (an empty value is 'unset', as in the native code).
        eff = {**self.extra_env, **(extra_env or {})}.get('QSTARS_TWO_PHOT_J') or None
        if eff is not None:
            eff = str(Path(eff).resolve())
            if tpj is not None and eff != tpj:
                raise ValueError(f'two_phot_j={tpj} but extra_env sets QSTARS_TWO_PHOT_J={eff}')
            tpj = eff
        if tpj is not None:
            if any(k == 'TWO_METH' and str(v).strip() != 'USE_RAD' for k, v in self.vadat):
                raise ValueError('two_phot_j requires TWO_METH=USE_RAD; the runner VADAT settings set another method')
            if not Path(tpj).is_file():
                raise FileNotFoundError(f'two_phot_j: {tpj} does not exist')
            # an older executable ignores the variable and silently runs OLD_DEFAULT: require the recorded capability
            cap = self.exe.with_suffix('.capabilities.json')
            if not (cap.is_file() and json.loads(cap.read_text()).get('qstars_two_phot_stored_j') is True):
                raise RuntimeError(f'two_phot_j needs an executable recorded with qstars_two_phot_stored_j (build >= 32): {self.exe}')
        if (dest / 'QSTARS_RESID_TRACE.bin').exists():
            # Reuse only a pass made at exactly this state, replay source and direction set.
            ok = (dest / 'QSTARS_STEP.bin').read_bytes() == np.asarray(s_cum, '<f8').tobytes(order='F')
            man = json.loads((dest / 'run_manifest.json').read_text()).get('extra_env', {})
            ok = ok and man.get('QSTARS_JV_REPLAY') == (str(Path(replay).resolve()) if replay else None)
            # same stored two-photon J, i.e. the same residual definition (absent in passes without it)
            ok = ok and (man.get('QSTARS_TWO_PHOT_J') or None) == tpj
            if dirs is not None:
                ok = ok and (dest / 'QSTARS_DIRS.bin').exists() and \
                    (dest / 'QSTARS_DIRS.bin').read_bytes()[12:] == np.asarray(dirs, '<f8').tobytes(order='F')
            # same linearization factor (QS_LIN_FAC 0: fixed-J blocks; 1: local-ALI blocks of --ali-precond)
            vad = dest / 'VADAT'             # every real pass has one
            lf = [ln.split()[0] for ln in vad.read_text().splitlines()
                  if ln.rstrip().endswith('[QS_LIN_FAC]')] if vad.is_file() else ['0']
            ok = ok and len(lf) == 1 and float(lf[0].replace('D', 'E').replace('d', 'e')) == float(lin_fac)
            if not ok:
                raise RuntimeError(f'{dest} exists but was made for a different state/direction; refusing reuse')
            check_one_record(dest / 'QSTARS_RESID_TRACE.bin')
            return dest
        if dest.exists():
            raise FileExistsError(f'{dest} exists without a residual: inspect before rerunning')
        sets = [('QS_ALI_LIN', 'T'), ('FIX_T', 'T' if self.fix_t else 'F'), ('QS_LIN_FAC', f'{float(lin_fac):g}')] + self.vadat
        if sum(k == 'QS_LIN_FAC' for k, _ in sets) != 1:
            raise ValueError('QS_LIN_FAC is set by the lin_fac argument, not through the runner VADAT settings')
        if tpj is not None:
            sets = sets + [('TWO_METH', 'USE_RAD')]
        prepare(self.x0dir, dest, 1, 1e-10, 2000, sets)
        if not keep_export:
            (dest / 'QSTARS_EXPORT').unlink(missing_ok=True)
        if not self.fix_t:
            release_t(dest)
        write_fortran(dest / 'QSTARS_STEP.bin', s_cum)
        env = {'QSTARS_PERTURB_FILE': str((dest / 'QSTARS_STEP.bin').resolve()), 'QSTARS_PERTURB_EPS': '1.0',
               'QSTARS_JV_STOP': '1', 'QSTARS_REPLACE_FORCE': str(self.replace)}
        if ba: env['QSTARS_BA_DUMP'] = '1'
        if store_j: env['QSTARS_RJ_STORE'] = '1'
        if replay: env['QSTARS_JV_REPLAY'] = str(Path(replay).resolve())
        if warm: env['QSTARS_ALI_WARM'] = str(Path(warm).resolve())   # iteration count only
        if self.no_obs: env['QSTARS_NO_OBS_FORMAL'] = '1'
        env.update(self.extra_env)
        if extra_env: env.update(extra_env)
        env.pop('QSTARS_TWO_PHOT_J', None)                       # set last: extra_env cannot replace the checked file
        if tpj is not None: env['QSTARS_TWO_PHOT_J'] = tpj       # stored conjugate-frequency J for USE_RAD
        if dirs is not None:
            write_directions(dest / 'QSTARS_DIRS.bin', dirs)
            env['QSTARS_JV_DIRS'] = str((dest / 'QSTARS_DIRS.bin').resolve())
        ali = dest / 'QSTARS_ALI_STATS'
        n_ali = len(ali.read_text().splitlines()) if ali.exists() else 0      # lines copied from the restart
        r = run(dest, self.exe, timeout, threads or self.threads, extra_env=env)
        if (dest / 'QSTARS_PERTURB_FAILURE').exists():
            raise RuntimeError(f'{dest}: perturbation rejected: {(dest / "QSTARS_PERTURB_FAILURE").read_text()}')
        tf = dest / 'QSTARS_TWO_PHOT_J_FAILURE'                  # build 32 stops on a missing/mismatched stored J
        if tf.exists():
            raise RuntimeError(f'{dest}: stored two-photon J rejected: {tf.read_text()}')
        if not (dest / 'QSTARS_RESID_TRACE.bin').exists():
            if 'exit code -11' in str(r.get('error')) and not _retried:
                # transient native SIGSEGV (seen once under heavy load; the identical pass then ran
                # normally): keep the failed pass for provenance and rerun the same pass once
                kept = dest.with_name(dest.name + '_segv')
                dest.rename(kept)
                self.log(f'{dest.name}: native SIGSEGV; failed pass kept as {kept.name}; rerunning once')
                return self.run(dest, s_cum, replay=replay, dirs=dirs, ba=ba, store_j=store_j, threads=threads,
                                timeout=timeout, warm=warm, extra_env=extra_env, keep_export=keep_export,
                                two_phot_j=two_phot_j, lin_fac=lin_fac, _retried=True)
            raise RuntimeError(f'{dest}: no residual ({r["error"]})')
        check_one_record(dest / 'QSTARS_RESID_TRACE.bin')
        self.ali_check(dest, n_ali)
        if self.slim:
            from slim_pass import slim
            slim(dest)                       # provenance-safe compaction (QSTARS_SLIM.json)
        return dest

    def ali_check(self, dest, skip):
        # Report 23: record an incomplete ALI (QS_ALI_MAX formal solutions reached before QS_ALI_TOL at some
        # frequency; CMFGEN then continues). Warning only, never fatal.
        try:
            st = ali_stats(dest / 'QSTARS_ALI_STATS', skip)
            if not st:
                self.log(f'WARNING {dest.name}: this pass wrote no QSTARS_ALI_STATS line')
            elif any(s_['unconverged'] > 0 for s_ in st):
                # the worst J computation of the pass (not the last one, which may have converged)
                s_ = max(st, key=lambda r: (r['unconverged'], r['max_res']))
                bad = sum(r['unconverged'] > 0 for r in st)
                self.log(f'WARNING {dest.name}: ALI unconverged at {int(s_["unconverged"])} of {int(s_["NCF"])} frequencies '
                         f'(max_res {s_["max_res"]:.3e}, tol {s_["tol"]:.1e}, max_it {int(s_["max_it"])}; worst of '
                         f'{bad} unconverged J computations out of {len(st)})')
        except Exception as e:                # noqa: BLE001 (diagnostic only)
            self.log(f'WARNING {dest.name}: QSTARS_ALI_STATS not parsed ({e})')


class Iterate:
    """Base pass data plus the complete direct (fixed-J) Jacobian. Block-diagonal without
    advection; with INC_ADV/INC_AD the native off-diagonal blocks (advection d/dr, adiabatic
    dT/dR) are kept and only the diagonal blocks are completed."""

    def __init__(self, base, x0, cols, h, with_t, hom=None, rj=None):
        # rj: stored-J file of the state when `base` is a frozen-J replay of it (Report 23, --two-phot-j:
        # the residual and blocks come from the replay with the state's own conjugate-frequency J)
        self.dir = Path(base)
        self.x, self.R = read_trace_record(self.dir / 'QSTARS_RESID_TRACE.bin')
        self.hom = hom
        self.s_cum = 1.0 - self.x / x0
        lo, di, up = read_ba_blocks(self.dir / 'QSTARS_BA_TRACE.bin', NT, ND)
        lo, di, up = lo.copy(), di.copy(), up.copy()
        self.tridiagonal = bool(np.any(lo) or np.any(up))

        def offdiag(c):
            # native off-diagonal contribution of column c perturbed at ALL depths (frozen-J
            # columns perturb every depth at once); exact for the linear advection/adiabatic terms
            o = np.zeros((NT, ND))
            o[:, 1:] += lo[1:, :, c].T
            o[:, :-1] += up[:-1, :, c].T
            return o
        # Replace the ne / ion / T columns by frozen-J central differences.
        colfd = {}
        for name, (p, m) in cols.items():
            _, rp = read_trace_record(Path(p) / 'QSTARS_RESID_TRACE.bin')
            _, rm = read_trace_record(Path(m) / 'QSTARS_RESID_TRACE.bin')
            colfd[name] = (rp - rm) / (2 * h)     # dR/d(fractional change of the column set)
        fd_diag = {n: colfd[n] - offdiag({'ne': 132, 'hii': 31, 'T': 133}[n]) for n in colfd}
        for d in range(ND):
            di[d][:, 132] = fd_diag['ne'][:, d]
            if with_t:
                di[d][:, 133] = fd_diag['T'][:, d]
        # Ion levels: native column (explicit dependence) + common POPION term per unit population.
        self.ion_missing = None
        if 'hii' in colfd:
            miss = np.stack([fd_diag['hii'][:, d] - di[d][:, 31] for d in range(ND)], axis=1) / self.x[31][None, :]
            self.ion_missing = miss
            for d in range(ND):
                for r in ION_ROWS:
                    di[d][:, r - 1] += miss[:, d] * self.x[r - 1, d]
        if hom is not None and hom.get('mode') == 'ptc':
            # pseudo-transient continuation on T (physical relaxation direction): RE row
            # G_T = R_RE - delta*c*(T - T_k)/T_k with T_k this base's T and c the RE row's gross
            # (l1) scale here. Net heating raises T; the term vanishes at any fixed point, so the
            # converged state is the root of the unmodified residual. Jacobian: -delta*c on the
            # T diagonal (fractional variables, T = T_k at the base); products are unchanged.
            self.ptc_c = (np.abs(lo[:, 133, :]) + np.abs(di[:, 133, :]) + np.abs(up[:, 133, :])).sum(axis=1)
            self.ptc_Tref = self.x[133].copy()
            for d in range(ND):
                di[d][133, 133] -= hom['delta'] * self.ptc_c[d]
            self.R_native = self.R.copy()
            self.R = self.hom_resid(self.R, self.x)
        elif hom is not None and hom.get('mode', 'pin') == 'shift':
            # residual-shift (Newton) homotopy on the RE row: R_T(lam) = R_RE - (1-lam)*R_RE(start);
            # a constant shift, so the Jacobian is unchanged
            self.R_native = self.R.copy()
            self.R = self.hom_resid(self.R, self.x)
        elif hom is not None:
            # T homotopy: R_T(lam) = lam*R_RE + (1-lam)*c*(T - T_fix)/T_fix (fractional variables:
            # d/ds_T of the second term is (1-lam)*c*T/T_fix); RE row of the direct blocks scaled by lam.
            # (frontpin: lam is per depth, 0 where T is pinned)
            lam = np.broadcast_to(np.asarray(pin_lam(hom), dtype=float), (ND,))
            lo[:, 133, :] *= lam[:, None]; di[:, 133, :] *= lam[:, None]; up[:, 133, :] *= lam[:, None]
            for d in range(ND):
                di[d][133, 133] += (1.0 - lam[d]) * hom['c'][d] * self.x[133, d] / hom['Tfix'][d]
            self.R_native = self.R.copy()
            self.R = self.hom_resid(self.R, self.x)
        self.Jdir = BlockTridiag(lo, di, up)
        di_pc = di
        if PC_T_FLOOR > 0:
            # preconditioner/merit only: the fixed-J RE-row T diagonal is 1e-7..1e-2 of the row's l1 scale
            # (near-singular); add PC_T_FLOOR x row scale to it in M. Products/Jacobian stay exact.
            di_pc = di.copy()
            for d in range(ND):
                rT = np.abs(lo[d][133]).sum() + np.abs(di[d][133]).sum() + np.abs(up[d][133]).sum()
                sgn = 1.0 if di[d][133, 133] >= 0 else -1.0
                di_pc[d][133, 133] += sgn * PC_T_FLOOR * rT
        self.M = BlockTridiag(lo, di_pc, up).factor()
        self._pc = (lo, di_pc, up)          # preconditioner blocks (install_exact_columns replaces columns)
        # structural activity from the native residual (FIXPOP rows are exactly zero there)
        self.active = ((self.R_native if hom is not None else self.R) != 0.0)
        if STRUCT_ACTIVE:
            # V1 plan P1: the exact number and charge rows are satisfied to the last bit by the conservation projection
            # (charge at 143-150 of 167 depths of free_r25), so R != 0 dropped them, and with them the n_e variable: no
            # direction moved n_e there and the projection changed it afterwards, unscored. Their activity comes from
            # the direct blocks instead.
            for r in PTC_CONSTRAINT_ROWS[:3]:
                self.active[r] |= np.any(di[:, r, :] != 0.0, axis=1)
        if hom is not None:
            self.active[133] = True
        self.rj = Path(rj) if rj is not None else self.dir / 'QSTARS_RJ.bin'
        xH = self.x[1:32].sum(0); xHe = self.x[32:132].sum(0)
        frac = np.ones_like(self.x); frac[0:32] = self.x[0:32] / xH; frac[32:132] = self.x[32:132] / xHe
        self.frac = frac
        self.W = significance_weights(frac) * self.active
        self.sig = (frac > SIG_FRAC) & self.active     # variables in the convergence test (extend_significance)
        self.coupling = None
        self.ali_inc = None                 # set by install_ali_increment
        self.pop_ptc = None                 # set by add_pop_ptc (Report 25)
        rs = (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T   # (NT, ND): row l1-norm of the direct blocks
        # equation-significance floor: rows whose gross-rate scale is below EQ_FLOOR x the largest at
        # that depth (e.g. He II<->He III balance where He III is absent, scale 1e-44) are not scaled up
        self.rs_gross = rs                  # unfloored (--pop-ptc scale)
        rs = np.maximum(rs, EQ_FLOOR * rs.max(axis=0, keepdims=True))
        self.D = np.where(rs > 0, 1.0 / np.where(rs > 0, rs, 1.0), 0.0)

    def install_ali_increment(self, ali_dir, fac=1.0):
        """Report 23: add the native local-ALI radiation response to the PRECONDITIONER blocks, all columns:
        M += fac * (BA[QS_LIN_FAC=1] - BA[QS_LIN_FAC=0]) from a frozen-J replay of this state with QS_LIN_FAC=1
        (ali_dir). The fixed-J blocks hold untrapped line and continuum rates; the local response cancels them
        where the transition is optically thick (its H I n=2 column equals the exact local column to 0.1-5 %).
        Not an accurate Jacobian (Report 18) and nearly singular at the ionization fronts: preconditioner only.
        Jdir and every product are unchanged. The RE row follows the homotopy scaling of this iterate."""
        xa, _ = read_trace_record(Path(ali_dir) / 'QSTARS_RESID_TRACE.bin')
        if not np.array_equal(xa, self.x):
            raise RuntimeError(f'{ali_dir}: local-ALI pass is not at this state')
        inc = [a_ - b_ for a_, b_ in zip(read_ba_blocks(Path(ali_dir) / 'QSTARS_BA_TRACE.bin', NT, ND),
                                         read_ba_blocks(self.dir / 'QSTARS_BA_TRACE.bin', NT, ND))]
        hom = self.hom
        if hom is not None and hom.get('mode', 'pin') == 'pin':
            lam = np.broadcast_to(np.asarray(pin_lam(hom), dtype=float), (ND,))
            for c in inc:
                c[:, 133, :] *= lam[:, None]
        lo, di, up = (m_ + fac * c for m_, c in zip(self._pc, inc))
        self.ali_inc = (self._pc[1].copy(), fac * inc[1])   # diagonal blocks before / of the increment (--front-probe)
        self._pc = (lo, di, up)
        self.M = BlockTridiag(lo, di, up).factor()
        return float(max(np.abs(c).max() for c in inc))

    def floor_preconditioner(self, tau, skip=None):
        """Report 24: raise the singular values of every equilibrated diagonal block of the PRECONDITIONER below
        tau x its largest to that floor (cnk_trap.floor_blocks). The local-ALI blocks are nearly singular at the
        ionization fronts, where M^-1 inflated the measure to 1e3-3e3 and the least squares chased those modes.
        Jdir and every product are unchanged. skip: (d0, d1) 0-based depths left unfloored (Report 25,
        --front-no-floor: blocks with exact columns). Returns the number of singular values raised."""
        lo, di, up = self._pc
        act = self.active if skip is None else self.active & ~cnk_trap.window_mask(self.active, skip)
        di2, n = cnk_trap.floor_blocks(di, tau, act)
        self._pc = (lo, di2, up)
        self.M = BlockTridiag(lo, di2, up).factor()
        return n

    def add_pop_ptc(self, delta, scale='gross', jacobian=True):
        """Report 25: pseudo-transient continuation of the populations. delta x sign(M_ii) x (row l1-norm of the direct
        blocks) is added to the diagonal of every statistical-equilibrium row, in the Jacobian (so every product, made
        or reused, is that of the shifted operator) and in the preconditioner; the H and He number, charge and T rows
        get none. The residual is unchanged, so a converged state is a root of the unmodified equations. The front's
        slow H I n = 2 + 3 / ionization mode is nearly singular in the exact Jacobian (Report 25); the shift bounds the
        step along it to a pseudo-time relaxation, and SER takes delta to 0 as the residual falls."""
        lo, di, up = self._pc
        if scale == 'pc':
            # Report 25 (--pop-ptc-scale pc): the row scale of the preconditioner blocks (local-ALI increment included):
            # a trapped level relaxes at its net rate; its fixed-J gross (untrapped Lyman alpha, 1e3..3e4 x larger at the
            # front) froze H I n=2 in the preconditioner, so no direction could change the H ionization balance at d18
            rsum = (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T
        else:
            rsum = self.rs_gross            # each row's own gross scale, whatever --eq-floor does to D
        diag = np.einsum('kii->ik', di)
        Dm = np.where(diag >= 0, 1.0, -1.0) * rsum * self.active
        Dm[list(PTC_CONSTRAINT_ROWS)] = 0.0
        shift = delta * Dm
        idx = np.arange(di.shape[1])
        di2 = di.copy(); di2[:, idx, idx] += shift.T
        self._pc = (lo, di2, up)
        self.M = BlockTridiag(lo, di2, up).factor()
        if not jacobian:
            # Report 25 (--pc-ptc): the shift regularizes the PRECONDITIONER only (the directions); the Jacobian, the
            # products and the least squares stay exact. Inside a least squares a shifted operator absorbs residual
            # instead of damping (iterate 20: true-Jacobian prediction x0.92 shifted, x0.405 unshifted)
            return
        self.pop_ptc = shift
        jdi = self.Jdir.di.copy(); jdi[:, idx, idx] += self.pop_ptc.T
        self.Jdir = BlockTridiag(self.Jdir.lo, jdi, self.Jdir.up)

    def install_exact_columns(self, dirs, jzs, meta, stride):
        """Report 23: replace the block-tridiagonal columns of the probed variables in the PRECONDITIONER by
        their exact responses (cnk_trap.install_columns). The fixed-J column of a resonance-line level holds
        the untrapped line rates (H I n=2 diagonal ~2000x the exact one at the H I front), which hides that
        level's error from M^-1 R. Jdir, and with it every product, is unchanged."""
        lo, di, up = (m_.copy() for m_ in self._pc)
        n = cnk_trap.install_columns(lo, di, up, dirs, jzs, meta, stride)
        self._pc = (lo, di, up)
        self.M = BlockTridiag(lo, di, up).factor()
        return n

    def extend_significance(self, coupling, rounds=8, rows=()):
        """Report 23: variables that count as significant (least-squares weight 1 and the convergence test):
        population fraction > SIG_FRAC of the element, n_e and T as before, plus the listed 1-based `rows`
        at every active depth (--sig-rows: levels given exact preconditioner columns, e.g. H I n=2), plus
        (coupling > 0) every level whose unit relative error changes a significant equation at the same depth
        by more than `coupling` of that equation's gross scale in the current preconditioner blocks, repeated
        until no level is added. Returns the number added.
        The coupling rule is kept for the record and is NOT used: on the 167-radius cool start it promoted
        5100 levels, among them ultra-trace He II n=2 / He III (fraction 1e-33) with M^-1 R ~ 1e6, because an
        ionization equation's row belongs to the neighbouring ion, not to the variable of the same index."""
        self.sig = (self.frac > SIG_FRAC) & self.active
        self.coupling = None
        nrow = 0
        for r in rows:
            add = self.active[r - 1] & ~self.sig[r - 1]
            self.sig[r - 1] |= add
            self.W[r - 1] = np.where(add, 1.0, self.W[r - 1])
            nrow += int(add.sum())
        if coupling <= 0:
            return nrow
        lo, di, up = self._pc
        rs = (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T
        inv = np.where(rs > 0, 1.0 / np.where(rs > 0, rs, 1.0), 0.0)
        total = 0
        for _ in range(rounds):
            c = cnk_trap.coupling_significance(di, inv, self.sig, NT - 2)
            add = (c > coupling) & self.active & ~self.sig
            self.coupling = c
            if not add.any():
                break
            self.sig = self.sig | add
            self.W = np.where(add, 1.0, self.W)
            total += int(add.sum())
        return total + nrow

    def hom_resid(self, R, x):
        if self.hom is None:
            return R
        lam = pin_lam(self.hom) if self.hom.get('mode', 'pin') == 'pin' else self.hom['lam']
        R2 = R.copy()
        if self.hom.get('mode') == 'ptc':
            R2[133] = R[133] - self.hom['delta'] * self.ptc_c * (x[133] - self.ptc_Tref) / self.ptc_Tref
        elif self.hom.get('mode', 'pin') == 'shift':
            R2[133] = R[133] - (1.0 - lam) * self.hom['R0']
        else:
            R2[133] = lam * R[133] + (1.0 - lam) * self.hom['c'] * (x[133] - self.hom['Tfix']) / self.hom['Tfix']
        return R2

    def obj(self, r):
        """Weighted residual measure used by the least squares and the line-search merit.
        precond : W M^-1 r (significance-weighted frozen-J correction).
        rowscale: W D r, D = 1 / (row l1-norm of the complete direct Jacobian) per depth, i.e. each
                  equation's imbalance relative to its own gross rate sensitivities (for the RE row:
                  net/gross heating). No inverse, so near-singular fixed-J rows are not amplified."""
        if MERIT_MODE == 'rowscale':
            return self.row_weights() * self.D * r
        return self.W * self.M.solve(r)

    def row_weights(self):
        """Report 25: equation weights of the row-scaled measure. The significance weights W belong to variables; as
        row weights they gave the helium-number row the weight of the He III VARIABLE, 0 where He III is ultra-trace
        (the He front, d55-58), and the least squares broke helium conservation there (He I ground s +1.8 / -1.5 at
        d56 / d57), which the projection of every trial then undid. The exact number and charge rows weigh 1."""
        w = self.W.copy()
        w[list(PTC_CONSTRAINT_ROWS[:3])] = self.active[list(PTC_CONSTRAINT_ROWS[:3])]
        return w

    def merit(self, R):
        # Significance-weighted measure (levels below SIG_FRAC of their element weigh TRACE_WEIGHT,
        # below ULTRA_FRAC ULTRA_WEIGHT; they are numerically irrelevant to the other equations).
        y = self.M.solve(R)
        return float(np.linalg.norm(self.obj(R)[self.active])), y


def significance_weights(frac):
    """Least-squares/merit row weights: 1 above SIG_FRAC, TRACE_WEIGHT down to ULTRA_FRAC, ULTRA_WEIGHT below.
    With WEIGHT_RAMP > 0 the step at SIG_FRAC becomes continuous: log-linear from TRACE_WEIGHT at
    SIG_FRAC*10^-WEIGHT_RAMP to 1 at SIG_FRAC (Report 23: levels whose fraction crosses SIG_FRAC otherwise flip
    between weights 1 and TRACE_WEIGHT from one iterate to the next). The convergence test is unaffected."""
    w = np.where(frac > SIG_FRAC, 1.0, np.where(frac > ULTRA_FRAC, TRACE_WEIGHT, ULTRA_WEIGHT))
    if WEIGHT_RAMP > 0:
        lf = np.log10(np.maximum(frac, 1e-300)); top = np.log10(SIG_FRAC)
        t = np.clip((lf - (top - WEIGHT_RAMP)) / WEIGHT_RAMP, 0.0, 1.0)
        ramp = (frac <= SIG_FRAC) & (frac > ULTRA_FRAC) & (t > 0)
        w = np.where(ramp, TRACE_WEIGHT ** (1.0 - t), w)
    return w


def pin_lam(hom):
    """Pin-mode blend weight of the RE row: the scalar lam, or lam times the per-depth weights w (0 at depths
    whose T is pinned, 1 elsewhere) of --hom-mode frontpin."""
    return hom['lam'] * hom['w'] if hom.get('w') is not None else hom['lam']


def make_recycled_pinv(it, pairs):
    """P = M + U Q^T with P z_i = (J z_i) on the stored span (Woodbury inverse).
    Pairs may come from the previous iterate (stale J): P is only a preconditioner."""
    Zm = np.stack([p[0].ravel(order='F') for p in pairs], axis=1)
    Wm = np.stack([p[1].ravel(order='F') for p in pairs], axis=1)
    Q, Rt = np.linalg.qr(Zm)
    keep = np.abs(np.diag(Rt)) > 1e-10 * np.max(np.abs(np.diag(Rt)))
    Q, Rt, Wm = Q[:, keep], Rt[np.ix_(keep, keep)], Wm[:, keep]
    shape = pairs[0][0].shape
    MQ = np.stack([it.Jdir.matvec(Q[:, i].reshape(shape, order='F')).ravel(order='F') for i in range(Q.shape[1])], axis=1)
    U = Wm @ np.linalg.inv(Rt) - MQ
    MiU = np.stack([it.M.solve(U[:, i].reshape(shape, order='F')).ravel(order='F') for i in range(U.shape[1])], axis=1)
    cap = np.eye(U.shape[1]) + Q.T @ MiU
    def pinv(v):
        mv = it.M.solve(v).ravel(order='F')
        corr = MiU @ np.linalg.solve(cap, Q.T @ mv)
        return (mv - corr).reshape(shape, order='F')
    return pinv


def group_report(v, mask):
    """'H I 9.9e+00@d35(row 3), ...': the largest |v| over the masked variables of each species group."""
    groups = [('H-', 0, 1), ('H I', 1, 31), ('H II', 31, 32), ('He I', 32, 101), ('He II', 101, 131),
              ('He III', 131, 132), ('n_e', 132, 133), ('T', 133, 134)]
    out = []
    for name, r0, r1 in groups:
        a_ = np.where(mask[r0:r1], np.abs(v[r0:r1]), -1.0)
        if a_.max() < 0:
            continue
        r_, d_ = np.unravel_index(int(np.argmax(a_)), a_.shape)
        out.append(f'{name} {a_[r_, d_]:.2e}@d{d_ + 1}' + (f'(row {r0 + r_ + 1})' if r1 - r0 > 1 else ''))
    return ', '.join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_21/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--free-t', action='store_true')
    ap.add_argument('--max-outer', type=int, default=30)
    ap.add_argument('--krylov', type=int, default=6)
    ap.add_argument('--eta', type=float, default=0.05)
    ap.add_argument('--h', type=float, default=1e-4)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--lams', type=float, nargs='+', default=[1.0, 0.5, 0.25])
    ap.add_argument('--tol', type=float, default=1e-3,
                    help='converged when max|exact Newton step| < tol over variables holding > --sig-frac of their element (and n_e, T)')
    ap.add_argument('--sig-frac', type=float, default=1e-10)
    ap.add_argument('--trace-weight', type=float, default=1e-3)
    ap.add_argument('--ultra-frac', type=float, default=0.0,
                    help='levels holding <= this fraction of their element weigh --ultra-weight in the Krylov '
                         'objective and merit (they cannot affect any other equation); 0 disables')
    ap.add_argument('--ultra-weight', type=float, default=1e-8)
    ap.add_argument('--weight-ramp', type=float, default=0.0,
                    help='decades below --sig-frac over which the least-squares weight ramps log-linearly from the trace '
                         'weight to 1 (0: step function as before); the convergence test is unchanged')
    ap.add_argument('--col-workers', type=int, default=0, help='max concurrent frozen-J column passes (0: all)')
    ap.add_argument('--col-threads', type=int, default=1,
                    help='OpenMP threads per column pass (Report 22: a column pass does not scale; 1 is fastest)')
    ap.add_argument('--trial-threads', type=int, default=0,
                    help='OpenMP threads per line-search trial pass (0: max(1, 8 // number of trials), as Report 21)')
    ap.add_argument('--obs-formal', action='store_true',
                    help='keep the observer-frame formal solution in solver passes (default: skipped, build >= 29)')
    ap.add_argument('--full-krylov-near', type=float, default=0.0,
                    help='optional guard (default off: Report 21 behaviour): when the step\'s max significant |s| < this '
                         'x --tol, run all --krylov passes (reused same-base passes count) before the convergence test')
    ap.add_argument('--consistent-dirs', action='store_true',
                    help='V1 plan P1: every span direction is replaced by its linearized conservation projection before '
                         'its exact product (cnk_trap.project_linear); needs --struct-active')
    ap.add_argument('--step-solver', choices=['ne', 'qr', 'box'], default='ne',
                    help='V1 plan P1, bounded trial step (--lm-bound): ne = the Report 25 normal-equation LM path; qr = the '
                         'same LM family from the augmented QR system; box = the explicit component-constrained least '
                         'squares |s| <= bound (second-order cone, needs the clarabel package)')
    ap.add_argument('--span-tau', type=float, default=0.0,
                    help='V1 plan P2: least squares and bounded step in the orthonormal basis of the normalized directions, '
                         'truncated at this relative singular value (bounds coefficient amplification by 1/tau); 0: off')
    ap.add_argument('--struct-active', action='store_true',
                    help='V1 plan P1: the H number, He number and charge rows (and their variables) are active wherever '
                         'their direct-block row is nonzero, not only where the residual is nonzero')
    ap.add_argument('--pc-t-floor', type=float, default=0.0,
                    help='preconditioner only: add this fraction of the RE row l1 scale to its T diagonal in M')
    ap.add_argument('--t-homotopy', action='store_true',
                    help='T continuation: R_T(lam) = lam*R_RE + (1-lam)*c*(T-T_fix)/T_fix, lam stepped from --lam0 to 1 '
                         '(T_fix and c, the RE-row gross scale, frozen at the start state in homotopy.npz)')
    ap.add_argument('--ptc-delta0', type=float, default=0.1,
                    help='ptc mode: initial pseudo-time weight delta (in units of the RE row gross scale per unit ln T)')
    ap.add_argument('--ptc-delta-reset', type=float, default=0.0,
                    help='ptc mode: on restart, replace the stored delta by this value (logged)')
    ap.add_argument('--ptc-target-dT', type=float, default=0.05,
                    help='ptc mode with --ptc-delta0 <= 0: delta0 = max|R_RE/c| / this (largest first T step)')
    ap.add_argument('--ptc-min', type=float, default=1e-4, help='ptc mode: delta below this is set to 0 (pure Newton)')
    ap.add_argument('--ptc-switch', type=float, default=1e-2,
                    help='ptc mode: when the regularized significant step is below this, delta is set to 0')
    ap.add_argument('--pin-depths', type=int, nargs=2, default=None, metavar=('D0', 'D1'),
                    help='frontpin mode: native depths D0..D1 (1-based, inclusive) whose T is pinned')
    ap.add_argument('--pin-factor', type=float, default=1.0,
                    help='frontpin mode: pinned T = this x the start base T at the pinned depths')
    ap.add_argument('--pin-all', action='store_true',
                    help='frontpin mode: pin T at EVERY depth (start-base T, times --pin-factor on --pin-depths): '
                         'the fixed-T statistical-equilibrium problem, with every RE imbalance kept in the native residual')
    ap.add_argument('--hom-mode', choices=['pin', 'shift', 'ptc', 'frontpin'], default='pin',
                    help='pin: lam*R_RE + (1-lam)*c*(T-T_fix)/T_fix; shift: R_RE - (1-lam)*R_RE(start) (Newton homotopy, '
                         'Jacobian unchanged, path scaled by the RE residual itself)')
    ap.add_argument('--lam0', type=float, default=0.0)
    ap.add_argument('--dlam', type=float, default=0.25)
    ap.add_argument('--lam-tol', type=float, default=1e-2, help='stage tolerance on the Newton step for lam < 1')
    ap.add_argument('--stage-max-iters', type=int, default=0,
                    help='predictor-corrector: advance lam after this many corrector iterates at lam < 1 even if the '
                         'stage tolerance is not met (0: no cap); the final lam = 1 stage always uses --tol')
    ap.add_argument('--natural', action='store_true',
                    help='line search / trust region on the natural monotonicity measure ||W J_k^-1 R|| '
                         '(simplified Newton correction in the current Krylov subspace) instead of the residual merit')
    ap.add_argument('--eq-floor', type=float, default=0.0,
                    help='rowscale merit: equation scale floor relative to the largest row scale at the same depth')
    ap.add_argument('--trace-sweep', action='store_true',
                    help='DO NOT USE (kept for provenance): trace rows can hold CMFGEN\'s species-conservation '
                         'equation (He III / H II rows), so the trace-only solve is not decoupled')
    ap.add_argument('--merit', choices=['precond', 'rowscale'], default='precond',
                    help='least-squares objective and line-search merit: W M^-1 R (precond) or W D R with D the '
                         'inverse row l1-norm of the complete direct Jacobian (rowscale). Krylov directions '
                         'are W M^-1 r in both modes; the convergence criterion is unchanged.')
    ap.add_argument('--warm', action='store_true', help='ALI warm start of trial passes from the base J (build >= 23)')
    ap.add_argument('--augment', type=int, default=0,
                    help='augmented minimal-residual solve: batch the products of M^-1 R and the previous '
                         'AUGMENT Newton steps in one pass, then add preconditioned residual directions')
    ap.add_argument('--split', action='store_true',
                    help='augmented solve: batch each new weighted preconditioned residual split into depth-region x '
                         'species-group pieces (their sum is the residual) so the least squares can scale regions '
                         'independently; all products remain exact')
    ap.add_argument('--split-tol', type=float, default=1e-3,
                    help='pieces with norm below this fraction of the residual are merged into one remainder direction')
    ap.add_argument('--recycle', type=int, default=0,
                    help='use up to this many most recent exact (z, J z) pairs (current and previous '
                         'iterates) in a low-rank corrected right preconditioner; 0 disables')
    ap.add_argument('--start-base', type=Path, default=None, help='accepted base pass to continue from')
    ap.add_argument('--start-k', type=int, default=0)
    ap.add_argument('--resume', action='store_true',
                    help='V1 plan P0: continue from the accepted state of OUT/accepted_state.json (base and iterate)')
    ap.add_argument('--allow-unaccepted-base', action='store_true',
                    help='start from a trial pass that its run\'s history did not accept (refused otherwise)')
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'),
                    help='additional VADAT numerical-resolution settings applied to every pass')
    ap.add_argument('--start-step', type=Path, default=None,
                    help='QSTARS_STEP.bin of a state (relative to x0) at which to create the first base pass')
    ap.add_argument('--t-log-cap', type=float, default=None,
                    help='trust region on the temperature update: |lambda s_T| <= this per iterate (numerical '
                         'damping of the update only; convergence is judged on the unlimited Newton step)')
    ap.add_argument('--lm-trust', type=float, default=None,
                    help='Levenberg-Marquardt trust region in the Krylov subspace (augmented solve): the step '
                         'minimizes ||W M^-1 (R - J S y)||^2 + mu ||W S y||^2 with mu chosen so that the '
                         'temperature step |s_T| <= radius (initial value given here); the radius adapts '
                         '(x2 after a full step with actual/predicted reduction > 0.5, x lambda after a damped '
                         'step). Convergence is judged on the unregularized (mu=0) Newton step.')
    ap.add_argument('--t-cap-mode', choices=['clip', 'scale'], default='clip',
                    help='how --t-log-cap limits the temperature step: clip each |s_T| (as before), or scale the WHOLE '
                         'Newton step by min(1, cap / max|s_T|) so T and populations stay on the same Newton direction')
    ap.add_argument('--lm-t-only', action='store_true',
                    help='Report 23: the Levenberg-Marquardt penalty acts on the temperature components of the step only '
                         '(mu ||S_T y||^2), so populations are solved consistently with the limited T step instead of '
                         'being shrunk with it')
    ap.add_argument('--lm-trust-max', type=float, default=1.0)
    ap.add_argument('--lm-trust-min', type=float, default=0.005)
    ap.add_argument('--component-log-cap', type=float, default=np.log(1e4),
                    help='per-component log-step limit applied before the line search; components below it '
                         '(all major species so far) keep the exact Newton direction')
    ap.add_argument('--probe-rows', type=int, nargs='+', default=[],
                    help='Report 23: 1-based variables whose block-tridiagonal columns in the PRECONDITIONER are '
                         'replaced by exact responses from one batched product pass per base (levels of optically '
                         'thick resonance lines, e.g. 2 3 = H I n=1, n=2); preconditioner only')
    ap.add_argument('--probe-stride', type=int, default=4,
                    help='each probe direction perturbs every STRIDE-th depth (STRIDE directions per probed row)')
    ap.add_argument('--front-probe', action='store_true',
                    help='Report 25: exact block-tridiagonal columns in the PRECONDITIONER for every variable with a '
                         'radiation response inside the hydrogen ionization fronts (depths with H II fraction < '
                         '--front-hii, widened by --front-margin), from batched product passes per base (front*); '
                         'installed after the local-ALI increment and before the SVD floor; preconditioner only')
    ap.add_argument('--front-rows', type=int, nargs='+', default=[],
                    help='1-based variables to probe (default: H-, H I, H II, n_e and every column whose local-ALI '
                         'increment changes a significant equation in the window by more than --front-inc-tol)')
    ap.add_argument('--front-inc-tol', type=float, default=1e-4)
    ap.add_argument('--front-hii', type=float, default=0.5)
    ap.add_argument('--front-margin', type=int, default=2)
    ap.add_argument('--front-stride', type=int, default=3,
                    help='each front direction perturbs one variable at every STRIDE-th window depth (>= 3: the '
                         'neighbour blocks are exact too)')
    ap.add_argument('--front-batch', type=int, default=18, help='directions per front product pass')
    ap.add_argument('--front-workers', type=int, default=3, help='front product passes at once')
    ap.add_argument('--front-threads', type=int, default=3, help='OpenMP threads per front product pass')
    ap.add_argument('--front-no-floor', action='store_true',
                    help='leave the front-window blocks (exact columns) out of --pc-svd-floor')
    ap.add_argument('--dense-window', type=int, nargs=2, action='append', default=[], metavar=('D0', 'D1'),
                    help='Report 26 (repeatable): dense exact Jacobian on depths D0..D1 (1-based): one direction per '
                         'variable and depth (no colouring), its exact product (full nonlocal column) in the span of the '
                         'least squares; product passes dense* (--dense-batch directions each, --front-workers at once, '
                         '--front-threads each), reused after a stop')
    ap.add_argument('--dense-rows', type=int, nargs='+', default=[],
                    help='1-based variables of the dense windows (default per window as --front-probe: H-, H I, H II, '
                         'n_e, T when free, and every column whose local-ALI increment changes a significant equation in '
                         'the window by more than --front-inc-tol)')
    ap.add_argument('--dense-batch', type=int, default=60, help='directions per dense product pass')
    ap.add_argument('--geodesic', action='store_true',
                    help='Report 26: after the trials, the second-order term of the residual along the step (from the '
                         'trials with lambda <= --geodesic-max-lam) and its correction t (J t = Q in the subspace); a '
                         'second round of trials on the curved path lambda s + lambda^2 t (labels g*). At the inner H '
                         'front the residual is quadratic along the step (iterate 3: full step 0.389 of the merit^2 at '
                         'depths 33-47 against 0.203 predicted, every other region linear)')
    ap.add_argument('--geodesic-max-lam', type=float, default=0.5,
                    help='trials with lambda above this are left out of the second-order fit (unless fewer than two remain)')
    ap.add_argument('--geodesic-alpha', type=float, default=1.0,
                    help='the correction is scaled down to ||W t|| <= alpha ||W s|| (the step it corrects)')
    ap.add_argument('--recycle-dense', action='store_true',
                    help='Report 26: at a base without --dense-window, the dense columns of the latest earlier iterate '
                         'that has them, as approximate pairs (radiation part from their own base, direct part exact)')
    ap.add_argument('--front-pieces', action='store_true',
                    help='Report 25 (with --split): inside the hydrogen-front window (as --front-probe) the H/n_e part '
                         'of every residual direction is split per depth (one Krylov direction per front depth)')
    ap.add_argument('--coarse-front', type=float, default=0.0,
                    help='Report 25: two-level preconditioner with the slow modes of the hydrogen-front blocks (window as '
                         '--front-probe; singular value below this fraction of the block\'s largest, before the floor) '
                         'and their exact products (coarse* passes); 0: off')
    ap.add_argument('--coarse-weight', type=float, default=0.5,
                    help='only slow modes with at least this share of the left vector on H-, H I, H II, n_e rows')
    ap.add_argument('--coarse-mode', choices=['ls', 'galerkin', 'span'], default='ls',
                    help="span (Report 25): the slow modes and their exact products only enter the span (with "
                         "--coarse-span); the preconditioner is not changed")
    ap.add_argument('--coarse-span', action='store_true', help='also put the slow modes (exact products) in the span')
    ap.add_argument('--pop-ptc', type=float, default=0.0,
                    help='Report 25: pseudo-transient continuation of the populations, initial delta (x the row gross '
                         'scale on every statistical-equilibrium diagonal; SER on the row-scaled residual); 0: off')
    ap.add_argument('--pop-ptc-min', type=float, default=1e-5, help='delta below this is set to 0 (pure Newton)')
    ap.add_argument('--pc-ptc', type=float, default=0.0,
                    help='Report 25: the population-PTC shift (scale --pop-ptc-scale) in the preconditioner only')
    ap.add_argument('--pop-ptc-scale', choices=['gross', 'pc'], default='gross',
                    help='population-PTC row scale: the direct (fixed-J) gross rates, or the preconditioner blocks '
                         '(local-ALI increment included: trapped levels at their net rates)')
    ap.add_argument('--split-fine', type=int, nargs=3, default=None, metavar=('D0', 'D1', 'W'),
                    help='Report 25: depths D0..D1 (1-based) split into residual-piece regions of W depths (with --split)')
    ap.add_argument('--target-dirs', type=int, default=0,
                    help='Report 25: per Krylov pass, also the preconditioned correction of each of this many cells with '
                         'the largest row-scaled residual, alone')
    ap.add_argument('--lm-bound', type=float, default=None,
                    help='Report 25: the line-search step is the Levenberg-Marquardt step in the subspace with the '
                         'smallest damping that keeps every |s| (populations, n_e) <= this; the convergence test keeps '
                         'the unbounded Newton step')
    ap.add_argument('--recycle-stale', type=int, default=0,
                    help='Report 25: also use the directions of the last K iterates, radiation part of their products '
                         'from their own bases, direct part recomputed exactly')
    ap.add_argument('--lm-bound-t', type=float, default=None,
                    help='Report 25, stage 2 (T free): radius of the T component in --lm-bound (fractional T change)')
    ap.add_argument('--trace-ion-balance', action='store_true',
                    help='Report 25: in every trial step, He III where it is trace (no direction moves it) takes the '
                         'component that zeroes the linearized He II / He III balance row (the replaced He II ground row)')
    ap.add_argument('--adaptive-cap', action='store_true',
                    help='Report 25: per-component limits (start --component-log-cap); a limited component whose update '
                         'reversed direction from the previous accepted update gets its limit halved, others grow back')
    ap.add_argument('--accept-rowscale', action='store_true',
                    help='Report 25: line-search trials are accepted on the row-scaled residual of the base\'s direct '
                         'blocks (the measure population PTC reduces) instead of the natural measure')
    ap.add_argument('--pop-ptc-switch', type=float, default=1e-2,
                    help='population-PTC step below this: delta -> 0 at the same base')
    ap.add_argument('--krylov-split', type=int, default=1,
                    help='Report 25: run each Krylov product pass as this many concurrent sub-passes (jvNNp0, ...)')
    ap.add_argument('--krylov-threads', type=int, default=0, help='OpenMP threads per Krylov sub-pass (0: --threads)')
    ap.add_argument('--sig-coupling', type=float, default=0.0,
                    help='Report 23: a level also counts as significant (weight 1, convergence test) when a unit '
                         'relative error in it changes a significant equation at its depth by more than this '
                         'fraction of that equation\'s gross scale; 0 (default): population fraction only')
    ap.add_argument('--ali-precond', action='store_true',
                    help='Report 23: add the native local-ALI radiation response of every column to the '
                         'PRECONDITIONER blocks, from one frozen-J replay of each base with QS_LIN_FAC=1 '
                         '(<base>_ali; fast with build >= 33). Exact --probe-rows columns are installed on top')
    ap.add_argument('--pc-svd-floor', type=float, default=0.0,
                    help='Report 24: floor the singular values of each equilibrated diagonal preconditioner block at '
                         'this fraction of its largest (bounded amplification of nearly singular local-ALI blocks); 0: off')
    ap.add_argument('--freeze-t', action='store_true',
                    help='Report 24: the fixed-T stage (--pin-all) with T removed from the step: no direction (residual '
                         'piece or previous step) has a temperature component, so T stays exactly at T_fix and the least '
                         'squares cannot use T as an instrument against the weakly weighted pin row')
    ap.add_argument('--conserve', action='store_true',
                    help='Report 24: project every line-search trial state onto the exact linear constraints (H and He '
                         'totals per depth equal to those of x0, n_e = charge sum); the per-component limiter breaks them. '
                         'A --lams value of 0 then tries the projected base itself')
    ap.add_argument('--trial-workers', type=int, default=0,
                    help='at most this many line-search trial passes at once (0: all at once, as before)')
    ap.add_argument('--sig-rows', type=int, nargs='+', default=[],
                    help='Report 23: 1-based variables that count as significant at every active depth whatever '
                         'their population fraction (weight 1, convergence test); meant for levels with exact '
                         'preconditioner columns (--probe-rows), e.g. 3 = H I n=2')
    ap.add_argument('--update-map', choices=['exp', 'linear'], default='exp',
                    help='trial state x*f(lambda s): exp(-s) as through Report 22, or linear (1 - s for s <= 0.5, '
                         'positive exponential tail beyond): exact Newton update for large population increases')
    ap.add_argument('--two-phot-j', action='store_true',
                    help='Report 23 (build >= 32): two-photon method USE_RAD with the conjugate-frequency J of the '
                         'state itself. The residual and blocks of each base are taken from a frozen-J replay of '
                         'the base with its own stored J as the two-photon J (<base>_tp); columns, products and '
                         'line-search trials of that base use the same file. Without it COMP_F=T forces OLD_DEFAULT')
    a = ap.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    # V1 plan P0: the start state is checked before any native pass; a rejected trial never becomes a base
    if a.resume:
        ptr = cnk_trap.read_accepted_pointer(out)
        if ptr is None:
            raise SystemExit(f'--resume: no {cnk_trap.ACCEPTED_POINTER} in {out}')
        if a.start_base or a.start_step:
            raise SystemExit('--resume takes the base from the pointer: no --start-base/--start-step')
        a.start_base, a.start_k = Path(ptr['base']), ptr['it']
    why_start = cnk_trap.check_restart_base(a.start_base) if a.start_base else None
    if why_start is not None and not a.allow_unaccepted_base:
        raise SystemExit(f'refusing to start from a rejected trial: {why_start} (--allow-unaccepted-base overrides)')
    R = Runner(a.x0, a.exe, a.replace, not a.free_t, a.threads, out / 'driver.log', [tuple(v) for v in a.vadat])
    R.no_obs = not a.obs_formal
    hist_path = out / 'history.json'
    hist = json.loads(hist_path.read_text()) if hist_path.exists() else []
    # x0 from an unperturbed base pass.
    global NT, ND, SIG_FRAC, TRACE_WEIGHT, ULTRA_FRAC, ULTRA_WEIGHT
    SIG_FRAC, TRACE_WEIGHT = a.sig_frac, a.trace_weight
    ULTRA_FRAC, ULTRA_WEIGHT = a.ultra_frac, a.ultra_weight
    global WEIGHT_RAMP
    WEIGHT_RAMP = a.weight_ramp
    global MERIT_MODE, EQ_FLOOR, PC_T_FLOOR, STRUCT_ACTIVE
    MERIT_MODE, EQ_FLOOR, PC_T_FLOOR = a.merit, a.eq_floor, a.pc_t_floor
    STRUCT_ACTIVE = a.struct_active
    if a.span_tau > 0 and a.geodesic:
        raise SystemExit('--geodesic combines the raw-direction products: not with --span-tau')
    if a.consistent_dirs and not a.struct_active:
        raise SystemExit('--consistent-dirs moves n_e at every depth: it needs --struct-active')
    nd_restart = int(next(l.split()[0] for l in (a.x0 / 'MODEL').read_text().splitlines() if 'Number of depth points' in l))
    nt_restart = int(next(l.split()[0] for l in (a.x0 / 'MODEL').read_text().splitlines() if 'Total number of variables' in l))
    NT, ND = nt_restart, nd_restart
    global SPLIT_REGIONS
    SPLIT_REGIONS = split_regions(ND)
    if a.split_fine:
        d0_, d1_, w_ = a.split_fine
        SPLIT_REGIONS = fine_regions(SPLIT_REGIONS, d0_ - 1, min(d1_, ND), w_)
    base = R.run(out / 'it00_base', np.zeros((NT, ND)))
    x0, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    cons_tot = cnk_trap.conservation_totals(x0)
    if a.conserve:
        xc_, pst_ = cnk_trap.project_conservation(x0, cons_tot)
        if pst_['ne'] > 1e-6:
            raise RuntimeError(f'--conserve: x0 does not satisfy charge conservation (|factor - 1| = {pst_["ne"]:.3e})')
    if x0.shape != (NT, ND):
        raise SystemExit(f'restart dimensions {x0.shape} differ from MODEL ({NT},{ND})')
    k = 0
    if a.start_base:
        base = a.start_base.resolve(); k = a.start_k
        if why_start is not None:
            R.log(f'WARNING starting from an unaccepted trial (--allow-unaccepted-base): {why_start}')
        R.log(f'continuing from {base.name} as iterate {k}')
    if a.start_step:
        from qstars.nlte_newton import read_fortran
        st = read_fortran(a.start_step, (NT, ND))
        base = R.run(out / 'start_base', st); k = a.start_k
        R.log(f'starting from state {a.start_step} as iterate {k}')
    if (a.start_base and not a.resume) or a.start_step or not (out / cnk_trap.ACCEPTED_POINTER).exists():
        cnk_trap.write_accepted_pointer(out, k, base, float('nan'),
                                        'start of a driver run' + (' (--resume)' if a.resume else ''))
    hom = None
    if a.t_homotopy and a.hom_mode == 'ptc':
        hs = out / 'ptc_state.json'
        if hs.exists():
            st_ = json.loads(hs.read_text())
            if a.ptc_delta_reset > 0:
                R.log(f'PTC: delta reset {st_["delta"]:.3g} -> {a.ptc_delta_reset:.3g} (--ptc-delta-reset)')
                st_['history'].append(dict(delta=a.ptc_delta_reset, reason='reset on restart'))
                st_['delta'] = a.ptc_delta_reset
        else:
            d0 = a.ptc_delta0
            if d0 <= 0:
                # automatic scale: the first pseudo-step changes T by at most --ptc-target-dT where the
                # PTC term dominates (|s_T| ~ |R_RE|/(delta*c))
                _, Rs = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
                lo_, di_, up_ = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', NT, ND)
                c = (np.abs(lo_[:, 133, :]) + np.abs(di_[:, 133, :]) + np.abs(up_[:, 133, :])).sum(axis=1)
                d0 = float(np.max(np.abs(Rs[133]) / c)) / a.ptc_target_dT
                R.log(f'PTC: automatic delta0 = max|R_RE/c| / {a.ptc_target_dT} = {d0:.3g}')
            st_ = dict(delta=d0, history=[])
        st_.update(lam=1.0, lam_done=1.0, dlam=0.0, stage_iters=0)
        hom = dict(mode='ptc', lam=1.0, delta=st_['delta'])
        R.log(f'T pseudo-transient continuation: delta={st_["delta"]:.4g} (state {hs})')
    elif a.t_homotopy and a.hom_mode == 'frontpin':
        # Coupled existence test (Report 23): T held at T_fix = pin_factor x T(start base) on --pin-depths through
        # the pin row c*(T - T_fix)/T_fix (RE row and its radiation response removed there), every other equation
        # unmodified and Newton-converged with the usual test. The RE imbalance of the native residual at the
        # pinned depths is then the coupled net heating at that front temperature, populations in statistical
        # equilibrium. lam = 1 (no stages).
        if not a.pin_depths:
            raise SystemExit('--hom-mode frontpin needs --pin-depths D0 D1')
        hp = out / 'homotopy_frontpin.npz'
        if hp.exists():
            z = np.load(hp)
            hom = dict(Tfix=z['Tfix'], c=z['c'], w=z['w'], R0=None, mode='pin')
        else:
            xs, Rs = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
            lo_, di_, up_ = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', NT, ND)
            c = (np.abs(lo_[:, 133, :]) + np.abs(di_[:, 133, :]) + np.abs(up_[:, 133, :])).sum(axis=1)
            win = np.zeros(ND, dtype=bool); win[a.pin_depths[0] - 1:a.pin_depths[1]] = True
            w = np.zeros(ND) if a.pin_all else np.where(win, 0.0, 1.0)
            Tfix = xs[133].copy(); Tfix[win] *= a.pin_factor
            hom = dict(Tfix=Tfix, c=c, w=w, R0=None, mode='pin')
            np.savez(hp, Tfix=Tfix, c=c, w=w, pin_depths=np.array(a.pin_depths), pin_factor=a.pin_factor,
                     pin_all=a.pin_all, source=str(base))
        hs = out / 'homotopy_frontpin_state.json'
        st_ = dict(lam=1.0, dlam=0.0, lam_done=1.0, stage_iters=0)
        hom.update(lam=1.0)
        R.log(f'T pinned at depths {a.pin_depths[0]}-{a.pin_depths[1]}: T_fix = {a.pin_factor:g} x start T '
              f'({hp}); every other equation unmodified')
    elif a.t_homotopy:
        tag = '' if a.hom_mode == 'pin' else '_shift'
        hp = out / f'homotopy{tag}.npz'
        if hp.exists():
            z = np.load(hp)
            hom = dict(Tfix=z['Tfix'], c=z['c'], R0=z['R0'] if 'R0' in z.files else None, mode=a.hom_mode)
        else:
            xs, Rs = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
            lo_, di_, up_ = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', NT, ND)
            c = (np.abs(lo_[:, 133, :]) + np.abs(di_[:, 133, :]) + np.abs(up_[:, 133, :])).sum(axis=1)
            hom = dict(Tfix=xs[133].copy(), c=c, R0=Rs[133].copy(), mode=a.hom_mode)
            np.savez(hp, Tfix=hom['Tfix'], c=c, R0=hom['R0'], source=str(base))
        hs = out / f'homotopy{tag}_state.json'
        st_ = json.loads(hs.read_text()) if hs.exists() else dict(lam=a.lam0, dlam=a.dlam, lam_done=None, stage_iters=0)
        hom.update(lam=st_['lam'])
        R.log(f'T homotopy: lam={st_["lam"]:.4g} dlam={st_["dlam"]:.4g} (T_fix from {hp})')
    caps, prev_logf = None, None
    caps_path, prev_path = out / 'component_caps.npy', out / 'component_prev_logf.npy'
    if a.adaptive_cap:
        # Report 25: per-component limits; a limit cycle across a component's root halves its limit (cnk_trap.adapt_caps)
        caps = np.load(caps_path) if caps_path.exists() else np.full((NT, ND), a.component_log_cap)
        if prev_path.exists():
            prev_logf = np.load(prev_path)
        elif hist and (Path(hist[-1]['base']) / 'QSTARS_RESID_TRACE.bin').exists():
            xb_, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
            xp_, _ = read_trace_record(Path(hist[-1]['base']) / 'QSTARS_RESID_TRACE.bin')
            prev_logf = np.log(np.maximum(xb_, 1e-300) / np.maximum(xp_, 1e-300))
        R.log(f'adaptive component limits: {int(np.sum(caps < a.component_log_cap))} below {a.component_log_cap:g} '
              f'({caps_path.name}); previous update {"known" if prev_logf is not None else "unknown"}')
    pop, pop_path = None, out / 'pop_ptc_state.json'
    if a.pop_ptc > 0:
        pop = json.loads(pop_path.read_text()) if pop_path.exists() else dict(delta=a.pop_ptc, history=[])
        R.log(f'population PTC: delta {pop["delta"]:.4g} ({pop_path})')
    while k < a.max_outer:
        if a.trace_sweep:
            base = trace_sweep(R, base, x0, out / f'it{k:02d}_trace', a.component_log_cap, k)
        itdir = out / f'it{k:02d}_from_{base.name}'
        itdir.mkdir(exist_ok=True)
        # --- frozen-J column passes at the base state (parallel) ---
        xk, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
        s_k = 1.0 - xk / x0
        names = ['ne', 'hii'] + (['T'] if a.free_t else [])
        rows = {'ne': 133, 'hii': 32, 'T': 134}
        # frozen-J columns depend only on the state: shared per base (homotopy stages reuse them)
        coldir = (out / f'cols_{base.name}') if a.t_homotopy else itdir
        coldir.mkdir(exist_ok=True)
        jobs = []
        for n in names:
            for sg in (+1, -1):
                e = np.zeros((NT, ND)); e[rows[n] - 1, :] = sg * a.h
                st = 1.0 - (1.0 - s_k) * (1.0 + e)
                jobs.append((n, sg, coldir / f'col_{n}_{"p" if sg > 0 else "m"}', st))
        nw = min(len(jobs), a.col_workers) if a.col_workers else len(jobs)
        # Report 23, --two-phot-j: one two-photon J file (the base's own stored J) defines the residual of this
        # iterate: base residual and blocks from a frozen-J replay with it, and the same file in every column,
        # product and trial pass. The next accepted base is replayed with its own J, so a converged state
        # satisfies the equations with its own conjugate-frequency J (native USE_RAD semantics, no lag).
        tpj = (base / 'QSTARS_RJ.bin') if a.two_phot_j else None
        resdir = base
        if tpj is not None:
            resdir = base.with_name(base.name + '_tp')
            R.run(resdir, s_k, replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False, threads=a.col_threads,
                  two_phot_j=tpj)
        alidir = base.with_name(base.name + '_ali') if a.ali_precond else None
        with ThreadPoolExecutor(nw + (1 if alidir else 0)) as ex:
            fut = ex.submit(R.run, alidir, s_k, replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False,
                            threads=a.col_threads, two_phot_j=tpj, lin_fac=1.0) if alidir else None
            list(ex.map(lambda j: R.run(j[2], j[3], replay=base / 'QSTARS_RJ.bin', ba=False, store_j=False,
                                        threads=a.col_threads, two_phot_j=tpj), jobs))
            if fut is not None:
                try:
                    fut.result()
                except RuntimeError as e:
                    # The native local response diverges where 1 - sigma*L* <= 0 (scattering-dominated or inverted
                    # infrared cells; VAR_JREL_V5 stops with QSTARS_ALI_LIN_FAILURE). Preconditioner only: this base
                    # keeps the fixed-J blocks; the pass directory is kept as evidence.
                    R.log(f'it{k:02d} WARNING local-ALI pass failed ({e}); this base uses the fixed-J preconditioner')
                    alidir = None
        cols = {n: (coldir / f'col_{n}_p', coldir / f'col_{n}_m') for n in names}

        def make_iterate():
            it_ = Iterate(resdir, x0, cols, a.h, a.free_t, hom, rj=base / 'QSTARS_RJ.bin')
            if alidir is not None:
                y_ = it_.M.solve(it_.R)
                it_.install_ali_increment(alidir)
                R.log(f'it{k:02d} local-ALI increment in the preconditioner ({alidir.name}): max|M^-1 R| over '
                      f'fraction-significant variables {np.max(np.abs(y_[it_.sig])):.3e} (fixed-J blocks) -> '
                      f'{np.max(np.abs(it_.M.solve(it_.R)[it_.sig])):.3e}')
            return it_
        it = make_iterate()
        if (hom is not None and hom['lam'] < 1.0 and a.stage_max_iters > 0
                and st_['stage_iters'] >= a.stage_max_iters):
            # predictor-corrector: enough corrector iterates at this lam; advance before solving
            st_['lam_done'] = hom['lam']
            st_['lam'] = min(1.0, hom['lam'] + st_['dlam'])
            R.log(f'it{k:02d} homotopy stage lam={st_["lam_done"]:.4g} capped after {st_["stage_iters"]} corrector '
                  f'iterates (stage tolerance {a.lam_tol} not required); advancing to lam={st_["lam"]:.4g}')
            st_['stage_iters'] = 0
            hom['lam'] = st_['lam']
            hs.write_text(json.dumps(st_, indent=1))
            it = make_iterate()
        # --- linear solve with exact products (tangent replay passes) ---
        # Objective: || W M^-1 (R - J s) ||, the significance-weighted frozen-J metric
        # (every rate equation weighted by its own fixed-J rates; ultra-trace rows weigh
        # TRACE_WEIGHT), never the raw residual whose rows span ~50 decades.
        nprod = 0

        def apply_J_batch(zs, name=None, threads=None):
            # name: a pass outside the numbered Krylov sequence (the exact-column probes)
            nonlocal nprod
            if name is None and a.krylov_split > 1 and len(zs) > 1:
                # Report 25: one Krylov pass as concurrent sub-passes jvNNp0, jvNNp1, ... (the products of different
                # directions are independent; same directions, same order)
                parts = np.array_split(np.arange(len(zs)), min(a.krylov_split, len(zs)))
                pname = f'jv{nprod:02d}'
                with ThreadPoolExecutor(len(parts)) as ex:
                    res_ = list(ex.map(lambda j: apply_J_batch([zs[i] for i in parts[j]], name=f'{pname}p{j}',
                                                               threads=a.krylov_threads or None), range(len(parts))))
                nprod += 1
                return [jz for part in res_ for jz in part]
            d = itdir / (name or f'jv{nprod:02d}')
            zst = np.stack(zs, axis=2)
            if d.exists() and not ((d / 'QSTARS_DIRS.bin').exists() and (d / 'QSTARS_STEP.bin').exists() and
                                   (d / 'QSTARS_DIRS.bin').read_bytes()[12:] == np.asarray(zst, '<f8').tobytes(order='F')
                                   and (d / 'QSTARS_STEP.bin').read_bytes() ==
                                   np.asarray(it.s_cum, '<f8').tobytes(order='F')):
                # a pass left by an earlier run of this iterate with other directions (e.g. made before same-base
                # products were available for reuse): kept for provenance, not reused
                sup = out / 'superseded'; sup.mkdir(exist_ok=True)
                tgt = sup / f'{itdir.name}_{d.name}_other_directions'
                d.rename(tgt)
                with open(sup / 'README.txt', 'a') as f_:
                    f_.write(f'{time.strftime("%Y-%m-%d %H:%M")} {itdir.name}/{d.name} -> {tgt.name}: made for other '
                             f'directions; superseded by a pass with the current direction set.\n')
                R.log(f'it{k:02d} {d.name}: existing pass had other directions; moved to superseded/{tgt.name}')
            R.run(d, it.s_cum, replay=it.rj, dirs=zst, ba=False, store_j=False, two_phot_j=tpj, threads=threads)
            _, drad = read_jv_rad(d / 'QSTARS_JV_RAD.bin')
            if name is None:
                nprod += 1
            if it.hom is not None and it.hom.get('mode', 'pin') == 'pin':
                drad[133] *= np.asarray(pin_lam(it.hom), dtype=float)[..., None]   # (ND, ndirs)
            jzs = [it.Jdir.matvec(z) + drad[:, :, i] for i, z in enumerate(zs)]
            # pairs.npz holds the products of the unshifted Jacobian (without --pop-ptc's diagonal)
            jz_exact = jzs if it.pop_ptc is None else [j - it.pop_ptc * z for j, z in zip(jzs, zs)]
            np.savez(d / 'pairs.npz', z=np.stack(zs, 2), jz=np.stack(jz_exact, 2))
            return jzs

        def consistent(z):
            # V1 plan P1 (--consistent-dirs): the direction as the conservation projection will apply it, so its exact
            # product prices the n_e (and H/He total) change of every trial
            return cnk_trap.project_linear(it.x, z) * it.active if a.consistent_dirs else z

        def cdir(c_):
            return consistent(cnk_trap.cell_direction(c_, it.x.shape))

        # Report 23: exact columns of radiatively trapped levels in the preconditioner, then the significant set
        if a.probe_rows:
            pz, pmeta = cnk_trap.probe_directions(a.probe_rows, it.active, a.probe_stride)
            y_fz = it.M.solve(it.R)
            pjz = apply_J_batch(pz, name='probe')
            ncol = it.install_exact_columns(pz, pjz, pmeta, a.probe_stride)
            y_ex = it.M.solve(it.R)
            R.log(f'it{k:02d} exact columns for rows {a.probe_rows} (stride {a.probe_stride}, {len(pz)} products, '
                  f'{ncol} columns): max|M^-1 R| over the probed rows {np.max(np.abs(y_fz[[r - 1 for r in a.probe_rows]])):.3e} '
                  f'(fixed-J columns) -> {np.max(np.abs(y_ex[[r - 1 for r in a.probe_rows]])):.3e}')
        # Report 25: exact columns at the hydrogen ionization fronts. The local-ALI blocks there have a nearly
        # singular H I n = 2 / n = 3 / ionization mode (relative singular value 1e-6 at it04_trial_l1) that the fixed-J
        # blocks do not have; it turned a deep direction's small far-field action at the front into 97-99.6 % of its
        # preconditioned action (Report 24). Every column with a radiation response inside the front window gets its
        # exact block-tridiagonal response from batched product passes; the probe pairs stay out of the span.
        front_rec = None
        fwin = cnk_trap.front_window(it.x, a.front_hii, a.front_margin) if a.front_pieces else None
        if a.front_probe:
            front = cnk_trap.front_window(it.x, a.front_hii, a.front_margin)
            if front is None:
                R.log(f'it{k:02d} front probe: H II fraction >= {a.front_hii:g} at every depth; no exact columns')
            else:
                always = list(range(1, 33)) + [133]          # H-, H I, H II, n_e
                if a.front_rows:
                    frows = sorted(set(a.front_rows))
                elif it.ali_inc is not None:
                    frows = cnk_trap.front_rows(it.ali_inc[1], it.ali_inc[0], it.sig, it.active, front,
                                                a.front_inc_tol, always=always)
                else:
                    frows = always
                if a.freeze_t:
                    frows = [r for r in frows if r != 134]
                fz, fmeta = cnk_trap.probe_directions(frows, cnk_trap.window_mask(it.active, front), a.front_stride)
                chunks = [fz[i:i + a.front_batch] for i in range(0, len(fz), a.front_batch)]
                c_fz = cnk_trap.block_conditioning(it._pc[1], it.active, front)
                y_fz = it.M.solve(it.R)
                R.log(f'it{k:02d} front probe: depths {front[0] + 1}-{front[1] + 1} (H II fraction < {a.front_hii:g} '
                      f'+ {a.front_margin}), {len(frows)} rows {frows}, stride {a.front_stride}: {len(fz)} directions '
                      f'in {len(chunks)} passes')
                with ThreadPoolExecutor(max(1, min(a.front_workers, len(chunks)))) as ex:
                    fjz = [jz for part in ex.map(lambda j: apply_J_batch(chunks[j], name=f'front{j:02d}',
                                                                         threads=a.front_threads), range(len(chunks)))
                           for jz in part]
                ncol = it.install_exact_columns(fz, fjz, fmeta, a.front_stride)
                c_ex = cnk_trap.block_conditioning(it._pc[1], it.active, front)
                y_ex = it.M.solve(it.R)
                R.log(f'it{k:02d} front probe: {ncol} exact columns installed; smallest relative singular value of '
                      f'the front blocks (median / min) {np.median(c_fz):.2e} / {np.min(c_fz):.2e} -> '
                      f'{np.median(c_ex):.2e} / {np.min(c_ex):.2e}; max|M^-1 R| over fraction-significant variables '
                      f'{np.max(np.abs(y_fz[it.sig])):.3e} -> {np.max(np.abs(y_ex[it.sig])):.3e}')
                front_rec = dict(front_window=[front[0] + 1, front[1] + 1], front_rows=frows,
                                 front_stride=a.front_stride, front_dirs=len(fz), front_columns=ncol,
                                 front_no_floor=a.front_no_floor)
        # Report 26: dense exact sub-domain Jacobian. The colored front columns (above) give the preconditioner a
        # block-tridiagonal front, but the soft modes are held by nonlocal radiative coupling (trapped Lya, He I 584 A)
        # and the floor overwrites them; a least squares over the colored probes themselves cut the residual of depths
        # 33-47 from 0.214 to 0.089 of the total (Reports/evidence/report26/span_check.py). Here each window variable at
        # each window depth gets its own exact product, which enters the span directly.
        dense, dense_rec = None, None
        if a.dense_window:
            always = list(range(1, 33)) + [133] + ([134] if a.free_t and not a.freeze_t else [])
            wins = []
            for d0_, d1_ in a.dense_window:
                w_ = (d0_ - 1, d1_ - 1)
                if a.dense_rows:
                    rws = sorted(set(a.dense_rows))
                elif it.ali_inc is not None:
                    rws = cnk_trap.front_rows(it.ali_inc[1], it.ali_inc[0], it.sig, it.active, w_, a.front_inc_tol,
                                              always=always)
                else:
                    rws = always
                if a.freeze_t:
                    rws = [r for r in rws if r != 134]
                wins.append((w_, rws))
            cells = cnk_trap.dense_cells(wins, it.active)
            batches = [cells[i:i + a.dense_batch] for i in range(0, len(cells), a.dense_batch)]
            for (w_, rws) in wins:
                R.log(f'it{k:02d} dense window: depths {w_[0] + 1}-{w_[1] + 1}, {len(rws)} rows {rws}')
            R.log(f'it{k:02d} dense windows: {len(cells)} single-cell directions in {len(batches)} passes '
                  f'({a.dense_batch} each, {a.front_workers} at once, {a.front_threads} threads)')
            for j, bt in enumerate(batches):
                pd_ = itdir / f'dense{j:02d}'
                if pd_.exists() and not (pd_ / 'QSTARS_RESID_TRACE.bin').exists():
                    # a pass cut off by a stop (power loss, kill): kept for provenance, made again
                    sup = out / 'superseded'; sup.mkdir(exist_ok=True)
                    tgt = sup / f'{itdir.name}_{pd_.name}_incomplete_{time.strftime("%Y%m%d%H%M")}'
                    pd_.rename(tgt)
                    with open(sup / 'README.txt', 'a') as f_:
                        f_.write(f'{time.strftime("%Y-%m-%d %H:%M")} {itdir.name}/{pd_.name} -> {tgt.name}: '
                                 f'incomplete (no residual) after a stop; made again.\n')
                    R.log(f'it{k:02d} {pd_.name}: incomplete pass moved to superseded/{tgt.name}')
            t0_ = time.time()
            with ThreadPoolExecutor(max(1, min(a.front_workers, len(batches)))) as ex:
                djz = list(ex.map(lambda j: apply_J_batch([cdir(c_) for c_ in batches[j]],
                                                          name=f'dense{j:02d}', threads=a.front_threads),
                                  range(len(batches))))
            dS = [cdir(c_) for c_ in cells]
            dJS = [jz for part in djz for jz in part]
            dense = (dS, dJS)
            R.log(f'it{k:02d} dense windows: {len(dS)} exact columns in {(time.time() - t0_) / 60:.1f} min')
            dense_rec = dict(dense_windows=[[w_[0] + 1, w_[1] + 1] for w_, _ in wins],
                             dense_rows=[rws for _, rws in wins], dense_dirs=len(dS), dense_batch=a.dense_batch)
        elif a.recycle_dense:
            # Report 26: the latest earlier iterate with dense passes (another base): approximate pairs, the radiation part
            # of each product from its own base (QSTARS_JV_RAD) and the direct part exact at this base
            for prevdir in sorted(out.glob('it*_from_*'), key=lambda q: q.name, reverse=True):
                if prevdir.resolve() == itdir.resolve() or prevdir.name.endswith(f'_from_{base.name}'):
                    continue
                dps = sorted(prevdir.glob('dense[0-9][0-9]'))
                if not dps or not all((q / 'pairs.npz').exists() and (q / 'QSTARS_JV_RAD.bin').exists() for q in dps):
                    continue
                dS, dJS = [], []
                for q in dps:
                    zs_ = np.load(q / 'pairs.npz')['z']
                    _, drad_ = read_jv_rad(q / 'QSTARS_JV_RAD.bin')
                    for i_ in range(zs_.shape[2]):
                        if a.freeze_t and np.any(zs_[133, :, i_] != 0.0):
                            continue
                        dS.append(zs_[:, :, i_]); dJS.append(it.Jdir.matvec(zs_[:, :, i_]) + drad_[:, :, i_])
                dense = (dS, dJS)
                R.log(f'it{k:02d} recycled {len(dS)} dense columns of {prevdir.name} (radiation part from its base, '
                      f'direct part exact at this base)')
                dense_rec = dict(dense_recycled_from=prevdir.name, dense_dirs=len(dS))
                break
        if a.pc_ptc > 0:
            it.add_pop_ptc(a.pc_ptc, scale=a.pop_ptc_scale, jacobian=False)
            R.log(f'it{k:02d} preconditioner PTC shift: {a.pc_ptc:g} x row scale ({a.pop_ptc_scale}) on the '
                  f'statistical-equilibrium diagonals of M only')
        if pop is not None and pop['delta'] > 0:
            it.add_pop_ptc(pop['delta'], scale=a.pop_ptc_scale)
            R.log(f'it{k:02d} population PTC: delta {pop["delta"]:.4g} x row gross scale on the statistical-equilibrium '
                  f'diagonals (Jacobian and preconditioner)')
            if dense is not None:
                # the dense products were made before the shift: give them the shifted operator's diagonal too
                dense = (dense[0], [j + it.pop_ptc * z for z, j in zip(*dense)])
        pc_prefloor = it._pc[1]                 # slow modes (--coarse-front) are taken before the floor
        if a.pc_svd_floor > 0:
            y_ = it.M.solve(it.R)
            skip = (front if front_rec is not None and a.front_no_floor else None)
            nraised = it.floor_preconditioner(a.pc_svd_floor, skip=skip)
            R.log(f'it{k:02d} preconditioner SVD floor {a.pc_svd_floor:g}: {nraised} singular values raised; max|M^-1 R| over '
                  f'fraction-significant variables {np.max(np.abs(y_[it.sig])):.3e} -> {np.max(np.abs(it.M.solve(it.R)[it.sig])):.3e}')
        nadd = it.extend_significance(a.sig_coupling, rows=a.sig_rows)
        if a.sig_coupling > 0 or a.sig_rows:
            R.log(f'it{k:02d} significant set: {nadd} variables added (rows {a.sig_rows}'
                  f'{f", coupling > {a.sig_coupling:g}" if a.sig_coupling > 0 else ""}) to the '
                  f'{int(it.sig.sum()) - nadd} significant by population fraction')
        # Report 25: two-level preconditioner. The exact front blocks keep a nearly singular slow mode (H I n = 2 + 3
        # against ionization, s/s_max 1e-6 at depths 23-26) that only nonlocal radiative coupling stabilizes; the slow
        # modes get exact products J Z (all of their coupling) and a coarse correction (cnk_trap.CoarseCorrected).
        coarse = None
        if a.coarse_front > 0:
            cwin = cnk_trap.front_window(it.x, a.front_hii, a.front_margin)
            if cwin is not None:
                Yc, Zc, cinfo = cnk_trap.slow_modes(pc_prefloor, it.active, cwin, a.coarse_front, rows=np.r_[0:32, 132],
                                                    min_weight=a.coarse_weight, exclude=[133] if a.freeze_t else [])
                chunks = [Zc[i:i + a.front_batch] for i in range(0, len(Zc), a.front_batch)]
                R.log(f'it{k:02d} coarse space: {len(Zc)} slow modes (s/s_max < {a.coarse_front:g}, left vector >= '
                      f'{a.coarse_weight:g} on H/n_e) at depths {cwin[0] + 1}-{cwin[1] + 1}, {len(chunks)} product passes')
                with ThreadPoolExecutor(max(1, min(a.front_workers, len(chunks)))) as ex:
                    JZc = [jz for part in ex.map(lambda j: apply_J_batch(chunks[j], name=f'coarse{j:02d}',
                                                                         threads=a.front_threads), range(len(chunks)))
                           for jz in part]
                cond_ = None
                if a.coarse_mode != 'span':
                    y_ = it.M.solve(it.R)
                    it.M = cnk_trap.CoarseCorrected(it.M, Yc, Zc, JZc, W=it.W * it.active, mode=a.coarse_mode)
                    y2_ = it.M.solve(it.R)
                    cond_ = it.M.cond
                    R.log(f'it{k:02d} coarse correction ({a.coarse_mode}): condition {it.M.cond:.3e}; max|M^-1 R| over '
                          f'fraction-significant variables {np.max(np.abs(y_[it.sig])):.3e} -> {np.max(np.abs(y2_[it.sig])):.3e}')
                coarse = (Zc, JZc)
                rec_coarse = dict(coarse_modes=len(Zc), coarse_window=[cwin[0] + 1, cwin[1] + 1], coarse_tau=a.coarse_front,
                                  coarse_mode=a.coarse_mode, coarse_cond=cond_, coarse_span=a.coarse_span)
        phi, s_fz = it.merit(it.R)
        R.log(f'it{k:02d} base {base.name}: merit={phi:.4e} max|M^-1 R|={np.max(np.abs(s_fz)):.3e} '
              f'||R||={np.linalg.norm(it.R):.4e}')
        if a.probe_rows or a.sig_coupling > 0 or a.sig_rows or a.ali_precond or a.front_probe:
            R.log(f'it{k:02d} max|M^-1 R| over significant variables by group: ' + group_report(s_fz, it.sig))

        def wm(r):
            return it.W * it.M.solve(r)
        b = wm(it.R)                        # first direction: weighted preconditioned residual
        bo = it.obj(it.R)                   # least-squares objective (== b in precond mode)
        beta = np.linalg.norm(bo)
        res_hist = [1.0]
        reused = 0
        if a.augment:
            # Augmented minimal residual (GCR-type): first pass batches M^-1 R with the
            # previous accepted Newton steps; each later pass adds the preconditioned residual.
            prev, prev_src = [], []
            for kk in range(k - 1, max(-1, k - 1 - a.augment), -1):
                cands = sorted(out.glob(f'it{kk:02d}_from_*/newton_step.npy'), key=lambda q: q.stat().st_mtime)
                if cands:
                    prev.append(np.load(cands[-1]) * it.active); prev_src.append(cands[-1].parent.resolve())
                    if a.freeze_t:
                        prev[-1][133] = 0.0
            def pieces(v):
                if a.freeze_t:
                    v = v.copy(); v[133] = 0.0        # Report 24: fixed-T stage, no temperature component
                if not a.split:
                    return [v.copy()]
                if a.front_pieces and fwin is not None:
                    return front_split_pieces(v, SPLIT_REGIONS, split_groups(a.free_t and not a.freeze_t),
                                              a.split_tol, fwin)
                return split_pieces(v, SPLIT_REGIONS, split_groups(a.free_t and not a.freeze_t), a.split_tol)
            def targeted(r):
                # Report 25 (--target-dirs K): the preconditioned correction of each of the K cells (row, depth) with the
                # largest row-scaled residual, alone: the coordinated local change that cell needs (H ionization balance
                # at d18: n=1 and n=2 +21 % together), which region-wide pieces cannot express; its exact product carries
                # the nonlocal coupling
                if a.target_dirs <= 0:
                    return []
                nr_ = 134 if (a.free_t and not a.freeze_t and hom is None) else 133
                rw = np.abs((it.row_weights() * it.D * it.active)[:nr_] * r[:nr_])
                out_ = []
                for c_ in np.argsort(rw.ravel())[::-1][:a.target_dirs]:
                    i_, d_ = np.unravel_index(c_, rw.shape)
                    e_ = np.zeros_like(r); e_[i_, d_] = r[i_, d_]
                    z_ = wm(e_)
                    if a.freeze_t:
                        z_[133] = 0.0
                    if np.any(z_):
                        out_.append(z_)
                return out_
            dirs = pieces(b) + prev + targeted(it.R)   # b = W M^-1 R (significance-weighted)
            S, JS = [], []
            if coarse is not None and a.coarse_span:
                S, JS = list(coarse[0]), list(coarse[1])     # Report 25: the slow modes with their exact products
            # The probe pairs are NOT added to the span: each perturbs one level by a common factor on every
            # STRIDE-th depth of the whole model, and the least squares used them as blunt instruments
            # (results/nk167_trap/fixedT iterate 0: period-6 sawtooth in the trial states' luminosity, every
            # trial refused by the native inner-luminosity guard). They serve the preconditioner only.
            y = None
            sigmask = it.sig

            reused_passes, reused_dirs = 0, set()
            reused_logical = set()                   # (iterate dir, jvNN): sub-passes jvNNpK of one pass count once

            def full_budget_needed(yv, done):
                # optional guard: near convergence, judge the 1e-3 test on the full-budget step (done = passes so far,
                # reused same-base passes included)
                if a.full_krylov_near <= 0 or done >= a.krylov:
                    return False
                st_ = sum(yi * si for yi, si in zip(yv, S))
                return float(np.max(np.abs(st_[sigmask]))) < a.full_krylov_near * a.tol
            # Report 22: exact products already made at this SAME base (earlier iterate at the same state, e.g.
            # after a PTC delta change) are reused: drad depends only on the state, replayed J and direction,
            # so jz = Jdir(current delta) z + drad is the exact product. Byte-identical state and replay required.
            for prevdir in sorted(out.glob(f'it*_from_{base.name}')):
                if prevdir.resolve() == itdir.resolve():
                    continue
                for jvd in sorted(prevdir.glob('jv*')):
                    if not (jvd / 'pairs.npz').exists() or not (jvd / 'QSTARS_JV_RAD.bin').exists():
                        continue
                    if (jvd / 'QSTARS_STEP.bin').read_bytes() != np.asarray(it.s_cum, '<f8').tobytes(order='F'):
                        continue
                    man_ = json.loads((jvd / 'run_manifest.json').read_text()).get('extra_env', {})
                    if man_.get('QSTARS_JV_REPLAY') != str(Path(it.rj).resolve()):
                        continue
                    if man_.get('QSTARS_TWO_PHOT_J') != (str(Path(tpj).resolve()) if tpj is not None else None):
                        continue                # made with another two-photon J: another residual definition
                    zs_ = np.load(jvd / 'pairs.npz')['z']
                    _, drad_ = read_jv_rad(jvd / 'QSTARS_JV_RAD.bin')
                    if it.hom is not None and it.hom.get('mode', 'pin') == 'pin':
                        drad_[133] *= np.asarray(pin_lam(it.hom), dtype=float)[..., None]
                    for i_ in range(zs_.shape[2]):
                        if a.freeze_t and np.any(zs_[133, :, i_] != 0.0):
                            continue            # a direction with a temperature component: not in the fixed-T span
                        S.append(zs_[:, :, i_]); JS.append(it.Jdir.matvec(zs_[:, :, i_]) + drad_[:, :, i_])
                        reused += 1
                    reused_logical.add((prevdir.name, jvd.name[:4])); reused_dirs.add(prevdir.resolve())
                    reused_passes = len(reused_logical)
            if a.recycle_stale > 0:
                # Report 25 (--recycle-stale K): the directions of the last K iterates (other bases) as approximate
                # pairs, the radiation part of each product from its own base (QSTARS_JV_RAD) and the direct part
                # recomputed exactly here; the state moves little per iterate (log changes <= 0.3) and the fit stalled
                # at 0.87 with ~100 fresh directions. The native line search judges the resulting step
                nst = 0
                for kk in range(k - 1, max(-1, k - 1 - a.recycle_stale), -1):
                    for prevdir in sorted(out.glob(f'it{kk:02d}*_from_*')):
                        if prevdir.name.endswith(f'_from_{base.name}'):
                            continue            # a same-base pass is reused exactly above
                        for jvd in sorted(prevdir.glob('jv*')):
                            if not (jvd / 'pairs.npz').exists() or not (jvd / 'QSTARS_JV_RAD.bin').exists():
                                continue
                            zs_ = np.load(jvd / 'pairs.npz')['z']
                            _, drad_ = read_jv_rad(jvd / 'QSTARS_JV_RAD.bin')
                            if it.hom is not None and it.hom.get('mode', 'pin') == 'pin':
                                drad_[133] *= np.asarray(pin_lam(it.hom), dtype=float)[..., None]
                            for i_ in range(zs_.shape[2]):
                                if a.freeze_t and np.any(zs_[133, :, i_] != 0.0):
                                    continue
                                S.append(zs_[:, :, i_]); JS.append(it.Jdir.matvec(zs_[:, :, i_]) + drad_[:, :, i_])
                                nst += 1
                if nst:
                    reused += nst
                    R.log(f'it{k:02d} recycled {nst} directions of iterates {max(0, k - a.recycle_stale)}-{k - 1} '
                          f'(radiation part from their own bases, direct part exact at this base)')
            ndense = 0
            if dense is not None:
                S += dense[0]; JS += dense[1]; ndense = len(dense[0])
            if reused or ndense:
                A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
                y = span_lstsq(S, A, bo.ravel(), a.span_tau)
                res = float(np.linalg.norm(bo.ravel() - A @ y) / beta)
                res_hist.append(res)
                if ndense:
                    rr_ = (bo.ravel() - A @ y).reshape(bo.shape)
                    R.log(f'it{k:02d} least squares with {ndense} dense-window columns'
                          f'{f" and {reused} reused same-base products" if reused else ""}: relative residual '
                          f'{res:.3e}; share of the base residual^2 by depths, before -> after: ' + ', '.join(
                              f'{a_ + 1}-{b_} {np.sum(bo[:, a_:b_] ** 2) / beta ** 2:.3f}->{np.sum(rr_[:, a_:b_] ** 2) / beta ** 2:.3f}'
                              for a_, b_ in DENSE_REPORT_REGIONS))
                else:
                    R.log(f'it{k:02d} reused {reused} exact same-base products ({reused_passes} passes): '
                          f'relative residual {res:.3e}')
                # a previous Newton step made at this same base is a combination of the reused directions
                # (s = sum y_i z_i, its product sum y_i J z_i): already in the span, no new product
                prev = [v for v, src in zip(prev, prev_src) if src not in reused_dirs]
                if res < a.eta and not full_budget_needed(y, reused_passes):
                    dirs = []
                else:
                    r_raw = it.R - sum(yi * ji for yi, ji in zip(y, JS))
                    dirs = pieces(wm(r_raw)) + prev + targeted(r_raw)
            # reused same-base passes count against the --krylov budget of exact-product passes per step
            if reused and reused_passes >= a.krylov and dirs:
                R.log(f'it{k:02d} Krylov budget ({a.krylov} passes) covered by reused same-base passes: no new pass')
            for pass_ in range(max(0, a.krylov - reused_passes)):
                if not dirs:
                    break
                dirs = [consistent(d_) for d_ in dirs]
                for d_, jd in zip(dirs, apply_J_batch(dirs)):
                    S.append(d_); JS.append(jd)
                A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
                y = span_lstsq(S, A, bo.ravel(), a.span_tau)
                res = float(np.linalg.norm(bo.ravel() - A @ y) / beta)
                res_hist.append(res)
                R.log(f'it{k:02d} augmented solve: {len(S)} directions, {nprod} passes, relative residual {res:.3e}')
                if res < a.eta:
                    if full_budget_needed(y, reused_passes + pass_ + 1):
                        R.log(f'it{k:02d} step near the convergence threshold: running the full Krylov budget')
                    else:
                        break
                r_raw = it.R - sum(yi * ji for yi, ji in zip(y, JS))
                dirs = pieces(wm(r_raw)) + targeted(r_raw)   # weighted preconditioned residual
            s = sum(yi * si for yi, si in zip(y, S))
            if a.freeze_t and np.any(s[133] != 0.0):
                raise RuntimeError('--freeze-t: the Newton step has a temperature component')
            lm_data = (S, np.stack([it.obj(j).ravel() for j in JS], axis=1), bo.ravel())
            span_k = None
            if a.span_tau > 0:
                # V1 plan P2: the bounded step in the orthonormal basis of the span (span_basis), not over raw directions
                Sb, Ab, _, kb = span_basis(S, lm_data[1], a.span_tau)
                lm_data = (Sb, Ab, bo.ravel())
                R.log(f'it{k:02d} span basis: {kb} of {len(S)} normalized directions above relative singular value '
                      f'{a.span_tau:g}')
                span_k = kb
        else:
            lm_data = None
            m = a.krylov
            V = [b / beta]; Z = []; H = np.zeros((m + 1, m)); y = None
            for j in range(m):
                z = V[j] * it.active
                Z.append(z)
                w = wm(apply_J_batch([z])[0])
                for _ in range(2):
                    for i in range(j + 1):
                        hij = np.sum(V[i] * w); H[i, j] += hij; w = w - hij * V[i]
                H[j + 1, j] = np.linalg.norm(w)
                V.append(w / H[j + 1, j])
                e1 = np.zeros(j + 2); e1[0] = beta
                y, *_ = np.linalg.lstsq(H[:j + 2, :j + 1], e1, rcond=None)
                res = float(np.linalg.norm(H[:j + 2, :j + 1] @ y - e1) / beta)
                res_hist.append(res)
                R.log(f'it{k:02d} GMRES {j + 1}: relative linear residual {res:.3e}')
                if res < a.eta:
                    break
            s = sum(yi * zi for yi, zi in zip(y, Z))
        smax = float(np.max(np.abs(s)))
        sig = it.sig
        smax_sig = float(np.max(np.abs(s[sig])))
        R.log(f'it{k:02d} Newton step over levels with fraction > {a.sig_frac:g}'
              f'{f" or coupling > {a.sig_coupling:g}" if a.sig_coupling > 0 else ""}'
              f'{f" or rows {a.sig_rows}" if a.sig_rows else ""}: max|s|={smax_sig:.3e}')
        if a.probe_rows or a.sig_coupling > 0 or a.sig_rows or a.ali_precond:
            R.log(f'it{k:02d} Newton step, max|s| over significant variables by group: ' + group_report(s, sig))
        R.log(f'it{k:02d} Newton step: max|s|={smax:.3e} (products {nprod})')
        if a.natural:
            # Natural (affine-invariant) monotonicity, Deuflhard: compare the simplified Newton
            # correction J_k^-1 R(trial), estimated in the current Krylov subspace with the existing
            # exact products (no new passes), against the current Newton correction. Measured in
            # variable space with the significance weights, so ill-conditioned rows (the cool-layer
            # RE row) are neither amplified nor hidden.
            Sm_nat = np.stack([si.ravel() for si in lm_data[0]], axis=1)
            A_nat = lm_data[1]
            vw = (it.W * it.active).ravel()

            def natural(Rx):
                y_ = np.linalg.lstsq(A_nat, it.obj(Rx).ravel(), rcond=1e-12)[0]
                return float(np.linalg.norm(vw * (Sm_nat @ y_)))
            phi_nat = natural(it.R)
            R.log(f'it{k:02d} natural measure ||W s_bar|| = {phi_nat:.4e} (residual merit {phi:.4e})')
            phi_res, phi = phi, phi_nat
        nrow = 134 if (a.free_t and not a.freeze_t and hom is None) else 133   # the RE row counts when T is free
        wrow = (it.row_weights() * it.D * it.active)[:nrow]

        def rowscale_merit(Rx):
            # Report 25 (--accept-rowscale): each equation's imbalance relative to its gross rates, this base's direct
            # blocks, populations and n_e, and the RE row in stage 2 (no inverse: the soft front modes are not amplified)
            return float(np.linalg.norm(wrow * Rx[:nrow]))
        if a.accept_rowscale:
            phi = rowscale_merit(it.R)
            R.log(f'it{k:02d} line-search merit: row-scaled residual {phi:.4e}')
        rec = dict(it=k, base=str(base), merit=phi, max_frozen_correction=float(np.max(np.abs(s_fz))),
                   resid_norm=float(np.linalg.norm(it.R)), gmres=res_hist, products=nprod, max_step=smax,
                   max_step_significant=smax_sig, products_reused=reused)
        if a.augment and a.span_tau > 0:
            rec['span_tau'], rec['span_basis'] = a.span_tau, span_k
        if hom is not None:
            rec['lam'] = hom['lam']
        if a.probe_rows:
            rec.update(probe_rows=list(a.probe_rows), probe_stride=a.probe_stride)
        if front_rec is not None:
            rec.update(front_rec)
        if dense_rec is not None:
            rec.update(dense_rec)
        if fwin is not None:
            rec['front_pieces_window'] = [fwin[0] + 1, fwin[1] + 1]
        if coarse is not None:
            rec.update(rec_coarse)
        if pop is not None:
            rec['pop_ptc'] = pop['delta']
        if a.sig_coupling > 0 or a.sig_rows:
            rec.update(sig_coupling=a.sig_coupling, sig_rows=list(a.sig_rows), significant_added=nadd,
                       significant=int(it.sig.sum()))
        if a.update_map != 'exp':
            rec['update_map'] = a.update_map
        if tpj is not None:
            rec['two_phot_j'] = str(tpj)
        if a.ali_precond:
            rec['ali_precond'] = True
        if a.pc_svd_floor > 0:
            rec['pc_svd_floor'] = a.pc_svd_floor
        try:                                  # the base pass's own ALI record is the file's last line
            rec['base_ali'] = {k_: ali_stats(base / 'QSTARS_ALI_STATS')[-1][k_] for k_ in ('unconverged', 'max_res', 'max_it')}
        except Exception:                     # noqa: BLE001 (diagnostic only)
            pass
        stage_tol = a.tol if (hom is None or hom['lam'] >= 1.0) else a.lam_tol
        capped = False
        if hom is not None and hom.get('mode') == 'ptc':
            rec['ptc_delta'] = hom['delta']
            if hom['delta'] > 0 and smax_sig < max(a.ptc_switch, a.tol):
                # regularized step already small: judge convergence on the unmodified Newton step
                rec['status'] = 'ptc_switch_same_base'
                hist.append(rec); hist_path.write_text(json.dumps(hist, indent=1))
                np.save(itdir / 'newton_step_ptc.npy', s)
                R.log(f'it{k:02d} PTC step {smax_sig:.3e} < {a.ptc_switch}; delta {hom["delta"]:.3g} -> 0 (pure Newton at this base)')
                hom['delta'] = 0.0; st_['delta'] = 0.0
                st_['history'].append(dict(it=k, delta=0.0, reason='switch'))
                hs.write_text(json.dumps(st_, indent=1))
                k += 1
                continue
        if pop is not None and pop['delta'] > 0 and smax_sig < max(a.pop_ptc_switch, stage_tol):
            rec['status'] = 'ptc_switch_same_base'
            hist.append(rec); hist_path.write_text(json.dumps(hist, indent=1))
            np.save(itdir / 'newton_step_pop_ptc.npy', s)
            R.log(f'it{k:02d} population-PTC step {smax_sig:.3e} < {max(a.pop_ptc_switch, stage_tol):g}: delta '
                  f'{pop["delta"]:.3g} -> 0 (pure Newton at this base)')
            pop['delta'] = 0.0; pop['history'].append(dict(it=k, delta=0.0, reason='switch'))
            pop_path.write_text(json.dumps(pop, indent=1))
            k += 1
            continue
        if smax_sig < stage_tol:
            rec['converged'] = True
            rec['status'] = 'converged' if (hom is None or hom['lam'] >= 1.0) else 'stage_converged_same_base'
            hist.append(rec); hist_path.write_text(json.dumps(hist, indent=1))
            np.save(itdir / 'newton_step.npy', s)
            if hom is None or hom['lam'] >= 1.0:
                R.log(f'CONVERGED: max Newton step over significant levels {smax_sig:.3e} < {a.tol} (all levels {smax:.3e})')
                break
            # stage converged: advance lam (adaptive step) at the same state
            st_['dlam'] = min(2.0 * st_['dlam'], 0.5) if st_['stage_iters'] <= 2 else st_['dlam']
            st_['lam_done'] = hom['lam']
            st_['lam'] = min(1.0, hom['lam'] + st_['dlam'])
            st_['stage_iters'] = 0
            hom['lam'] = st_['lam']
            hs.write_text(json.dumps(st_, indent=1))
            R.log(f'it{k:02d} homotopy stage lam={st_["lam_done"]:.4g} '
                  f'{"converged" if rec["converged"] else "capped"} (max|s| {smax_sig:.3e}, stage tolerance {stage_tol}); '
                  f'advancing to lam={st_["lam"]:.4g}')
            k += 1
            continue
        np.save(itdir / 'newton_step.npy', s)
        # --- parallel line-search trials (each is a full base pass) ---
        # Per-component log-step limiter (numerical damping of the update, not of the state):
        # ultra-trace levels whose populations must change by many e-folds move by at most
        # component_log_cap per iterate; convergence is judged on the unlimited Newton step.
        # Trial steps. Default: lambda in a.lams along the (component-limited) Newton step.
        # With --lm-trust (T released): Levenberg-Marquardt steps in the Krylov subspace at radii
        # Delta, Delta/2, Delta/4 on |s_T| (a curved trust-region path), unless the Newton step
        # already satisfies the smallest radius.
        steps = {}
        lm_mode = a.free_t and a.lm_trust is not None and lm_data is not None
        if lm_mode:
            radius = next((h_['lm_trust_next'] for h_ in reversed(hist) if 'lm_trust_next' in h_), a.lm_trust)
            radii = [radius * f for f in a.lams]
            if np.max(np.abs(s[133])) <= radii[-1]:
                lm_mode = False
                R.log(f'it{k:02d} Newton step |s_T|={np.max(np.abs(s[133])):.3e} within radius {radii[-1]:.3g}: lambda trials')
            else:
                rec['lm_trust'] = radius
                rec['lm'] = {}
                unconstrained_seen = False
                wlm = it.W * it.active
                if a.lm_t_only:
                    wlm = np.zeros_like(wlm); wlm[133] = it.active[133]
                for r_ in radii:
                    s_r, mu, pred = lm_step(lm_data, wlm, r_)
                    if mu == 0.0 and unconstrained_seen:
                        R.log(f'it{k:02d} LM radius {r_:.4g}: same unconstrained step as a larger radius; no trial')
                        continue
                    unconstrained_seen = unconstrained_seen or mu == 0.0
                    if a.natural:          # linear-model remaining correction after the step
                        pred = float(np.linalg.norm(vw * (s - s_r).ravel()))
                    steps[f'r{r_:.4g}'] = s_r
                    rec['lm'][f'r{r_:.4g}'] = dict(mu=mu, predicted_merit=pred, max_sT=float(np.max(np.abs(s_r[133]))),
                                                   max_s=float(np.max(np.abs(s_r))))
                    R.log(f'it{k:02d} LM radius {r_:.4g}: mu={mu:.3e} max|s_T|={np.max(np.abs(s_r[133])):.3e} '
                          f'max|s|={np.max(np.abs(s_r)):.3e} predicted merit {pred:.4e}')
                np.save(itdir / 'lm_steps.npy', np.stack([steps[k_] for k_ in steps], axis=2))
        if not lm_mode:
            sc0 = s.copy()
            if a.lm_bound is not None and lm_data is not None:
                # Report 25: the trial step re-optimized within |s| <= lm_bound (the convergence test above keeps the
                # unbounded Newton step)
                bmask = it.active.copy()
                rad = np.full(s.shape, a.lm_bound)
                if a.free_t and not a.freeze_t and a.lm_bound_t is not None:
                    rad[133] = a.lm_bound_t       # stage 2: T in the bound with its own radius
                else:
                    bmask[133] = False
                if np.max(np.abs(s[bmask]) / rad[bmask]) > 1.0:
                    if a.step_solver == 'ne':
                        sc0, mu_b, _ = lm_step_bound(lm_data, it.active.astype(float), rad, bmask)
                    else:
                        sc0, mu_b, _ = bounded_step(lm_data, it.active, rad, bmask, a.step_solver)
                    pred_b = float(np.linalg.norm(lm_data[1] @ np.linalg.lstsq(np.stack([si.ravel() for si in lm_data[0]], 1),
                                                                                 sc0.ravel(), rcond=None)[0] - lm_data[2])
                                   / np.linalg.norm(lm_data[2]))
                    R.log(f'it{k:02d} trial step bounded to |s| <= {a.lm_bound:g}'
                          f'{f" (T: {a.lm_bound_t:g})" if bmask[133].any() else ""} (mu {mu_b:.2e}); its linear residual '
                          f'{pred_b:.3e} against {res_hist[-1]:.3e} unbounded')
                    rec.update(lm_bound=a.lm_bound, lm_bound_mu=mu_b, lm_bound_linear=pred_b)
                    # V1 plan (Report 26 review): the planned step and its linear-model merit at every trial lambda, so the
                    # actual/predicted reduction rho of the accepted trial is computed exactly afterwards
                    np.save(itdir / 'trial_step.npy', sc0)
                    yb_ = np.linalg.lstsq(np.stack([si.ravel() for si in lm_data[0]], 1), sc0.ravel(), rcond=None)[0]
                    nb_ = np.linalg.norm(lm_data[2])
                    rec['predicted_rel_by_lam'] = {f'l{l_:.4g}': float(np.linalg.norm(l_ * (lm_data[1] @ yb_) - lm_data[2]) / nb_)
                                                   for l_ in a.lams}
            if a.free_t and a.t_log_cap is not None and a.lm_trust is None:
                ntl = int(np.sum(np.abs(s[133]) > a.t_log_cap))
                if a.t_cap_mode == 'scale':
                    fac = min(1.0, a.t_log_cap / max(float(np.max(np.abs(s[133]))), 1e-300))
                    sc0 = fac * s
                    R.log(f'it{k:02d} temperature trust region |s_T| <= {a.t_log_cap:g}: whole step scaled by {fac:.3g} '
                          f'(max requested |s_T| = {np.max(np.abs(s[133])):.3e})')
                    rec['t_step_scale'] = fac
                else:
                    sc0[133] = np.clip(s[133], -a.t_log_cap, a.t_log_cap)
                    R.log(f'it{k:02d} temperature trust region |s_T| <= {a.t_log_cap:g} active at {ntl} depths '
                          f'(max requested |s_T| = {np.max(np.abs(s[133])):.3e})')
                rec['t_limited_depths'] = ntl
            for lam_ in a.lams:
                steps[f'l{lam_:.4g}'] = lam_ * sc0
        with np.errstate(divide='ignore'):
            nlim = max(int(np.sum((np.abs(v) > a.component_log_cap) if a.update_map == 'exp' else
                                  (np.abs(np.log(cnk_trap.update_factor(v, 'linear', np.inf))) > a.component_log_cap)))
                       for v in steps.values())
        R.log(f'it{k:02d} component limiter active on {nlim} of {s.size} variables')
        rec['limited_components'] = nlim

        def trial(item):
            lab, st_ = item
            capv = a.component_log_cap if caps is None else caps
            xt = it.x * cnk_trap.update_factor(st_, a.update_map, capv)
            if a.trace_ion_balance:
                # Report 25: He III where it is trace balances the replaced He II ground row (cnk_trap.trace_ion_balance),
                # given the REALIZED (limited) update of every other component: balancing against the planned step made
                # He III run away at the He front once the limiter truncated the rest (row 102 at d60: 0.03 -> 4.4)
                s_real = 1.0 - xt / it.x
                s_real[131] = 0.0
                Jr = it.Jdir.matvec(s_real)
                if it.pop_ptc is not None:
                    Jr = Jr - it.pop_ptc * s_real          # the unshifted direct action
                mask_, sc_ = cnk_trap.trace_ion_balance(it.R, Jr, it.Jdir.di, it.frac, 1.0)
                f_ = cnk_trap.update_factor(sc_, a.update_map, capv if np.isscalar(capv) else capv[131])
                xt[131, mask_] = it.x[131, mask_] * f_[mask_]
                R.log(f'it{k:02d} trial {lab}: trace-ion balance at {int(mask_.sum())} depths, max|s(He III)| '
                      f'{float(np.max(np.abs(sc_))) if mask_.any() else 0.0:.3g}')
            if a.conserve:
                # Report 24: back onto the exact number and charge conservation rows (broken by the limiter)
                xt, pst = cnk_trap.project_conservation(xt, cons_tot)
                R.log(f'it{k:02d} trial {lab}: conservation projection, largest |factor - 1|: H {pst["h"]:.3e}, '
                      f'He {pst["he"]:.3e}, n_e {pst["ne"]:.3e}')
            d = out / f'it{k + 1:02d}_trial_{lab}'
            st_cum = 1.0 - xt / x0
            res = d / 'run_result.json'
            if d.exists() and not (d / 'QSTARS_RESID_TRACE.bin').exists() and res.exists() and \
                    (d / 'QSTARS_STEP.bin').read_bytes() == np.asarray(st_cum, '<f8').tobytes(order='F'):
                err = json.loads(res.read_text()).get('error')
                R.log(f'it{k:02d} trial {lab}: rejected earlier at this exact state ({err})')
                return lab, d, float('inf')
            try:
                R.run(d, st_cum, threads=a.trial_threads or max(1, 8 // len(steps)), warm=it.rj if a.warm else None,
                      two_phot_j=tpj)
            except RuntimeError as e:
                # A trial state the native code itself refuses (e.g. CMFGEN's STOP_IF_MAJOR_WARNING
                # guard: inner normalized luminosity out by a factor of 2) is not an admissible
                # iterate: it is rejected like a trial that raises the merit. The pass is kept.
                R.log(f'it{k:02d} trial {lab}: rejected, native pass failed ({e})')
                return lab, d, float('inf')
            xr, Rt = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
            Rt = it.hom_resid(Rt, xr)
            if a.accept_rowscale:
                return lab, d, rowscale_merit(Rt)
            return lab, d, (natural(Rt) if a.natural else it.merit(Rt)[0])
        with ThreadPoolExecutor(min(len(steps), a.trial_workers) if a.trial_workers > 0 else len(steps)) as ex:
            trials = list(ex.map(trial, list(steps.items())))
        for lab, d, ph in trials:
            R.log(f'it{k:02d} trial {lab}: merit {ph:.4e} (current {phi:.4e})')
        if a.geodesic and lm_data is not None and not lm_mode:
            # Report 26: geodesic acceleration (Transtrum & Sethna). The trials give the remainder of the linear model,
            # E_lam = R(x_lam) - (R - lam J s - Jdir d_lam) (d_lam: realized minus planned step, from the projection,
            # He III balance and limiter, by the direct blocks), and from it the second-order term Q; the correction
            # t with J t = Q in the same subspace makes R(x(1 - lam s - lam^2 t)) = R - lam J s + O(lam^3).
            Sm_g = np.stack([si.ravel() for si in lm_data[0]], axis=1)
            yb_ = np.linalg.lstsq(Sm_g, sc0.ravel(), rcond=None)[0]
            Jsb_ = sum(yi * ji for yi, ji in zip(yb_, JS))
            fin = [(float(lab_[1:]), d_) for lab_, d_, ph_ in trials if np.isfinite(ph_) and lab_.startswith('l')]
            use = [t_ for t_ in fin if t_[0] <= a.geodesic_max_lam]
            if len(use) < 2:
                use = fin
            if use:
                Es_ = []
                for lam_, d_ in use:
                    xt_, Rt_ = read_trace_record(d_ / 'QSTARS_RESID_TRACE.bin')
                    Rt_ = it.hom_resid(Rt_, xt_)
                    dl_ = np.where(it.active, (1.0 - xt_ / it.x) - lam_ * sc0, 0.0)
                    Es_.append(Rt_ - (it.R - lam_ * Jsb_ - it.Jdir.matvec(dl_)))
                Q_ = cnk_trap.quadratic_term([u_[0] for u_ in use], Es_)
                q_ = it.obj(Q_).ravel()
                yt_ = np.linalg.lstsq(lm_data[1], q_, rcond=1e-12)[0]
                fit_t = float(np.linalg.norm(q_ - lm_data[1] @ yt_) / max(np.linalg.norm(q_), 1e-300))
                t_g = (Sm_g @ yt_).reshape(sc0.shape)
                if a.freeze_t:
                    t_g[133] = 0.0
                wv_ = it.W * it.active
                ratio = float(np.linalg.norm(wv_ * t_g) / max(np.linalg.norm(wv_ * sc0), 1e-300))
                gsc = min(1.0, a.geodesic_alpha / max(ratio, 1e-300))
                t_g = gsc * t_g
                qreg = (it.obj(Q_) ** 2)
                R.log(f'it{k:02d} geodesic: second-order term from trials {[u_[0] for u_ in use]}, |obj(Q)| / merit '
                      f'{np.linalg.norm(q_) / phi:.3f} (by depths: ' + ', '.join(
                          f'{a_ + 1}-{b_} {qreg[:, a_:b_].sum() / max(qreg.sum(), 1e-300):.2f}' for a_, b_ in DENSE_REPORT_REGIONS)
                      + f'); correction fit {fit_t:.3f}, ||W t|| / ||W s|| {ratio:.3f}, scaled by {gsc:.3g}; '
                      f'max|t| {np.max(np.abs(t_g)):.3f}')
                np.save(itdir / 'geodesic_t.npy', t_g)
                gsteps = {f'g{lam_:.4g}': lam_ * sc0 + lam_ ** 2 * t_g for lam_ in a.lams}
                with ThreadPoolExecutor(min(len(gsteps), a.trial_workers) if a.trial_workers > 0 else len(gsteps)) as ex:
                    gtr = list(ex.map(trial, list(gsteps.items())))
                for lab, d, ph in gtr:
                    R.log(f'it{k:02d} geodesic trial {lab}: merit {ph:.4e} (current {phi:.4e})')
                trials += gtr
                rec.update(geodesic=dict(lams=[u_[0] for u_ in use], q_over_merit=float(np.linalg.norm(q_) / phi),
                                         fit=fit_t, ratio=ratio, scale=gsc))
        if min(t_[2] for t_ in trials) >= phi:
            small = {f'{lab}_d16': v / 16.0 for lab, v in steps.items()}
            with ThreadPoolExecutor(min(len(small), a.trial_workers) if a.trial_workers > 0 else len(small)) as ex:
                more = list(ex.map(trial, list(small.items())))
            for lab, d, ph in more:
                R.log(f'it{k:02d} fallback trial {lab}: merit {ph:.4e} (current {phi:.4e})')
            trials += more
        # V1 plan P0: the best trial, the accepted one (None when no trial lowers the merit) and the outcome are kept
        # apart; only an accepted trial becomes the next base and the restart pointer
        dec = cnk_trap.trial_decision(trials, phi)
        lab, d, ph = dec['best_trial'], dec['best_dir'], dec['best_merit']
        rejected = dec['accepted'] is None
        rec.update(trials={t_[0]: t_[2] for t_ in trials}, best_trial=lab, best_merit=ph, accepted=dec['accepted'])
        if not rejected:
            rec['merit_new'] = ph
            rec['status'] = 'accepted'
            if lab.startswith('l') and '_d16' not in lab:
                rec['accepted_lam'] = float(lab[1:])
        elif (hom is not None and hom.get('mode') == 'ptc') or pop is not None or \
                (hom is not None and st_['lam_done'] is not None and st_['dlam'] > 1e-3):
            rec['status'] = 'rejected_retry_same_base'
        else:
            rec['status'] = 'stopped_no_decrease'
        if 'lm_trust' in rec:
            top = f'r{radii[0]:.4g}'
            full = dict((t_[0], t_[2]) for t_ in trials)[top]
            rho = (phi - full) / max(phi - rec['lm'][top]['predicted_merit'], 1e-300)
            if lab == top and rho > 0.5:
                nxt = min(2.0 * rec['lm_trust'], a.lm_trust_max)
            elif lab == top:
                nxt = rec['lm_trust']
            elif lab.startswith('r') and '_d16' not in lab:
                nxt = max(float(lab[1:]), a.lm_trust_min)
            else:
                nxt = max(rec['lm_trust'] / 16.0, a.lm_trust_min)
            rec.update(lm_rho_top=rho, lm_trust_next=nxt)
            R.log(f'it{k:02d} LM actual/predicted reduction at radius {radii[0]:.4g}: {rho:.3f}; next radius {nxt:.3g}')
        hist.append(rec); hist_path.write_text(json.dumps(hist, indent=1))
        if not rejected:
            cnk_trap.write_accepted_pointer(out, k + 1, d, ph, f'iterate {k} accepted trial {lab}')
        if rejected and hom is not None and hom.get('mode') == 'ptc':
            old = hom['delta']
            hom['delta'] = max(4.0 * old, a.ptc_delta0 if old == 0 else 0.0)
            st_['delta'] = hom['delta']; st_['history'].append(dict(it=k, delta=hom['delta'], reason='no decrease'))
            hs.write_text(json.dumps(st_, indent=1))
            R.log(f'it{k:02d}: no decrease; PTC delta {old:.3g} -> {hom["delta"]:.3g} at the same base')
            k += 1
            continue
        if rejected and pop is not None:
            old = pop['delta']
            pop['delta'] = max(4.0 * old, a.pop_ptc if old == 0 else 0.0)
            pop['history'].append(dict(it=k, delta=pop['delta'], reason='no decrease'))
            pop_path.write_text(json.dumps(pop, indent=1))
            R.log(f'it{k:02d}: no decrease; population PTC delta {old:.3g} -> {pop["delta"]:.3g} at the same base')
            k += 1
            continue
        if rejected:
            if hom is not None and st_['lam_done'] is not None and st_['dlam'] > 1e-3:
                st_['dlam'] /= 2.0
                st_['lam'] = min(1.0, st_['lam_done'] + st_['dlam'])
                hom['lam'] = st_['lam']; st_['stage_iters'] = 0
                hs.write_text(json.dumps(st_, indent=1))
                R.log(f'it{k:02d}: no decrease at this lam; retreating to lam={st_["lam"]:.4g} (dlam {st_["dlam"]:.4g})')
                k += 1
                continue
            R.log(f'it{k:02d}: no merit decrease along the exact Newton direction; stopping for inspection')
            break
        if hom is not None and hom.get('mode') == 'ptc' and hom['delta'] > 0:
            # switched evolution relaxation: delta scales with the RE-row imbalance (net/gross)
            _, Rn = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
            r_old = float(np.linalg.norm(it.R_native[133] / it.ptc_c))
            r_new = float(np.linalg.norm(Rn[133] / it.ptc_c))
            # delta only decreases here (larger pseudo-time steps as the RE imbalance falls); while the
            # imbalance grows (crossing a thermally unstable stretch toward the stable root) it is held,
            # and it rises only when a trial is rejected (x4, below)
            fac = min(max(r_new / max(r_old, 1e-300), 0.1), 1.0)
            old = hom['delta']
            hom['delta'] = old * fac if old * fac >= a.ptc_min else 0.0
            st_['delta'] = hom['delta']
            st_['history'].append(dict(it=k, delta=hom['delta'], re_net_gross_old=r_old, re_net_gross_new=r_new))
            hs.write_text(json.dumps(st_, indent=1))
            R.log(f'it{k:02d} PTC: RE net/gross {r_old:.3e} -> {r_new:.3e}; delta {old:.3g} -> {hom["delta"]:.3g}')
        elif hom is not None:
            st_['stage_iters'] += 1
            hs.write_text(json.dumps(st_, indent=1))
        if caps is not None:
            xn_, _ = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
            new_logf = np.log(np.maximum(xn_, 1e-300) / np.maximum(it.x, 1e-300))
            if prev_logf is not None:
                caps, nflip = cnk_trap.adapt_caps(caps, prev_logf, new_logf, a.component_log_cap)
                rec['cap_flips'] = nflip; rec['caps_below_default'] = int(np.sum(caps < a.component_log_cap))
                hist_path.write_text(json.dumps(hist, indent=1))
                R.log(f'it{k:02d} adaptive limits: {nflip} limited components reversed direction (limit halved); '
                      f'{rec["caps_below_default"]} below {a.component_log_cap:g}, smallest {float(caps.min()):.3g}')
            prev_logf = new_logf
            np.save(caps_path, caps); np.save(prev_path, prev_logf)
        if pop is not None and pop['delta'] > 0:
            # SER on the row-scaled residual (no inverse) of this base's direct blocks: delta follows ||R_new|| / ||R||
            xn_, Rn_ = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
            Rn_ = it.hom_resid(Rn_, xn_)
            nr_ = 134 if (a.free_t and not a.freeze_t and hom is None) else 133
            w_ = (it.row_weights() * it.D * it.active)[:nr_]
            rho = float(np.linalg.norm(w_ * Rn_[:nr_]) / np.linalg.norm(w_ * it.R[:nr_]))
            old = pop['delta']
            # a shortened step (lambda < 1 accepted) means the PTC step was too long for the nonlinearity: raise delta
            lam_acc = rec.get('accepted_lam', 1.0 / 16.0 if '_d16' in lab else 1.0)
            new = old * min(max(rho / np.sqrt(max(lam_acc, 1e-3)), 0.25), 4.0)
            pop['delta'] = new if new >= a.pop_ptc_min else 0.0
            pop['history'].append(dict(it=k, delta=pop['delta'], rowscale_ratio=rho, accepted_lam=lam_acc))
            pop_path.write_text(json.dumps(pop, indent=1))
            rec['pop_ptc_rowscale_ratio'] = rho
            hist_path.write_text(json.dumps(hist, indent=1))
            R.log(f'it{k:02d} population PTC: row-scaled residual ratio {rho:.3f}, accepted lambda {lam_acc:g}; delta '
                  f'{old:.3g} -> {pop["delta"]:.3g}')
        base = d
        k += 1


if __name__ == '__main__':
    sys.exit(main())
