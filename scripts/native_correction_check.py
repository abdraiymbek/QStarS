#!/usr/bin/env python3
"""Read-only prediction of CMFGEN's native correction at a base pass (the audit's gate quantity).

The verify pass (scripts/verify_cool_state.py) runs one complete native iteration with QS_ALI_LIN=T,
QS_LIN_FAC=0 and the run's forced equation table; its native correction is the solution of the
assembled fixed-J system, i.e. exactly the blocks a base pass dumps in QSTARS_BA_TRACE.bin (after
REPLACE/FIXPOP and the column scaling by POPS) applied to its residual QSTARS_RESID_TRACE.bin.
Here: y = M^-1 R with the FULL native block-tridiagonal M (native n_e/T columns, no pc floor, no PTC
term, no completion by frozen-J columns), y the fractional correction (POPS -> POPS(1 - y)).
It is NOT the exact Newton step (cool_newton_krylov.py convergence test); both are reported by the
closing phase (benchmarks/cmfgen_cool_validation_plan.json closing_phase).

Classes (active rows only; FIXPOP rows have R = 0 exactly):
  significant          fraction of its element > --sig-frac (n_e and T included, as in the NK test)
  ne, T                rows 133, 134
  trace_l1_kept        fraction <= --sig-frac whose row l1 weight on trace columns exceeds that on
  trace_l1_demoted     significant columns (scripts/trace_cleanup.py default, order-independent) / the rest
  trace_paired_own     structural pairing (structural_pairing below): own statistical-equilibrium row
  trace_paired_ionization   paired with an ionization-equation row (C_ION, replaced ground row)
  trace_unpaired       no row of the assembled matrix carries its determining equation
  all                  every active variable (the audit's last_nonlinear_correction counts all levels)

Per class the gate value is max audit_normalized(y) and its location (gate_*; criterion C of the closing phase);
the largest |y| and its location (max_abs_*) can belong to another variable when a decrease dominates.

Floating-point floor (closing_phase floating_point_floor): --noise-pair P1 P2 ... MEASURES it at this state as
max over identical-state passes P (BASE's POPS bit for bit, its equation selection) of |M^-1 (R_BASE - R_P)|
with BASE's own M (measured_floor); only this floor is used for FAIL-at-floor counts. --noise-eps is a random-sign
MODEL (noise_floor), reported only (model_floor_max), never used for classification or stopping: at it14 it gives
T 5e-4 and 418 variables above 1e-3 with EPS 3.2e-12 where the measured floor is T <= 1e-11, all <= 1.8e-4.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import BlockTridiag, read_ba_blocks, read_trace_record  # noqa: E402

# H/He variable layout of this model (0-based rows; NT = 134), CMFGEN cmfgen.f ion numbering:
# ion order H0 (H-), HI, HII (top), HeI, He2, HeIII (top) = columns of the replace table.
NT_HHE = 134
SPECIES = [dict(ions=[(0, 0, 1), (1, 1, 30)], top=31),            # (replace column, ground row, levels)
           dict(ions=[(3, 32, 69), (4, 101, 30)], top=131)]
NE_ROW, T_ROW = 132, 133


def fractions(x):
    """Fraction of its element (H: H I + H II; He: all He); n_e and T rows are 1 (as the driver)."""
    xH = x[1:32].sum(0); xHe = x[32:132].sum(0)
    frac = np.ones_like(x); frac[0:32] = x[0:32] / xH; frac[32:132] = x[32:132] / xHe
    return frac


def one_record(path, nt, nd):
    """A pass trace must hold exactly one record (a stale record copied from a restart is read first)."""
    n = Path(path).stat().st_size
    if n != 12 + 16 * nt * nd:
        raise ValueError(f'{path}: {n} bytes, expected one record ({12 + 16 * nt * nd}); stale or appended trace')


def l1_split(frac, active, di, sig_frac):
    """trace_cleanup default classification. Each trace row is compared against the ORIGINAL trace set of
    its depth (copy; the earlier in-loop view mutation made it depend on the row order)."""
    trace = (frac <= sig_frac) & active
    kept = trace.copy()
    for d in range(trace.shape[1]):
        tc = trace[:, d]
        for i in np.nonzero(tc)[0]:
            row = np.abs(di[d][i])
            if row[~tc].sum() >= row[tc].sum():
                kept[i, d] = False
    return kept, trace & ~kept


def read_replace(path, nd):
    """Equation selection actually used: rows 'kd f1..f6' (QSTARS_REPLACE_TRACE or the forced table)."""
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip()]
    rep = np.array([[int(v) for v in r[1:]] for r in rows], dtype=bool)
    if rep.shape != (nd, 6) or [int(r[0]) for r in rows] != list(range(1, nd + 1)):
        raise ValueError(f'{path}: not a {nd}-depth H/He equation-selection table')
    return rep


def structural_pairing(frac, active, rep, sig_frac):
    """Equation-variable pairing of the trace subsystem from CMFGEN's assembly (generate_full_matrix_v3.f).

    Row i is the statistical-equilibrium (SE) equation of level i, except: the top-ion row of each species
    holds the number-conservation equation (EQ_SPECIES: rows 32, 132); the ground row of ion ID holds the
    ionization equation C_ION(ID) (ion ID <-> ion ID+1 balance) where REPLACE(ID) is set at that depth.
    A trace variable is solved with its own row only if that row is its SE equation. A trace ion
    representative (top ion; ground of a replaced ion) is solved with the C_ION row linking it to an
    adjacent ion that is significant or already paired (upward chain, then downward). Conservation rows
    are never used. Whatever remains has no row of its own in the assembled matrix (its equation is
    implicit, e.g. a sum of significant rows): 'unpaired', held fixed. Returns (rows, cols, kind) with
    kind 1 = own row, 2 = ionization row, and the unpaired mask; rows/cols are (nt, nd) index arrays."""
    nt, nd = frac.shape
    trace = (frac <= sig_frac) & active
    row_of = -np.ones((nt, nd), int); kind = np.zeros((nt, nd), int)
    for d in range(nd):
        own = trace[:, d].copy()
        for sp in SPECIES:
            own[sp['top']] = False
            for c, g, n in sp['ions']:
                if rep[d, c]:
                    own[g] = False
        row_of[own, d] = np.nonzero(own)[0]; kind[own, d] = 1
        for sp in SPECIES:
            ions = sp['ions']
            reps = [g for _, g, _ in ions] + [sp['top']]          # ion representatives, last = top
            # a representative needs an ionization row if it is trace and has no SE row of its own
            need = [bool(trace[g, d]) and (k == len(ions) or bool(rep[d, ions[k][0]])) for k, g in enumerate(reps)]
            done = [not n for n in need]                          # significant, or solved with its own row
            used = set()
            for k in range(len(ions)):                            # upward: C_ION(k) determines ion k+1
                c, g, _ = ions[k]
                if rep[d, c] and done[k] and not done[k + 1]:
                    row_of[reps[k + 1], d] = g; kind[reps[k + 1], d] = 2; done[k + 1] = True; used.add(g)
            for k in range(len(ions) - 1, -1, -1):                # downward: C_ION(k) determines ion k
                c, g, _ = ions[k]
                if rep[d, c] and g not in used and done[k + 1] and not done[k]:
                    row_of[g, d] = g; kind[g, d] = 2; done[k] = True; used.add(g)
    unpaired = trace & (row_of < 0)
    return row_of, kind, unpaired


def classes(x, R, di, rep, sig_frac):
    frac = fractions(x)
    active = R != 0.0
    kept, demoted = l1_split(frac, active, di, sig_frac)
    row_of, kind, unpaired = structural_pairing(frac, active, rep, sig_frac)
    rows = np.zeros_like(active)
    ne = rows.copy(); ne[NE_ROW] = True
    tt = rows.copy(); tt[T_ROW] = True
    return frac, dict(significant=(frac > sig_frac) & active, ne=ne & active, T=tt & active,
                      trace_l1_kept=kept, trace_l1_demoted=demoted,
                      trace_paired_own=(kind == 1) & active, trace_paired_ionization=(kind == 2) & active,
                      trace_unpaired=unpaired, all=active)


def audit_normalized(y):
    """CMFGEN prints % increase and % decrease relative to the old value; the audit renormalizes the
    decrease to the new value (cmfgen_backend.nonlinear_corrections): y > 0 is a decrease."""
    return np.where(y < 0, -y, np.where(y < 0.99999, y / np.maximum(1.0 - y, 1e-300), 1e5))


def row_l1(lo, di, up):
    """(nt, nd) row l1 of the native blocks at each depth (the unit of eps_row_l1)."""
    return (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T


def noise_floor(M, lo, di, up, active, eps, samples, seed=0):
    """REPORT-ONLY model: max over `samples` sign realizations of |M^-1 r|, r = eps x (row l1 of the native blocks)
    x (+-1) on active rows; eps in row-l1 units (measured_floor records it per pair; Report 22's 3.2e-12 is
    normalized by max_depth|R_row| instead, bench_pass.rowwise, and is NOT this unit). Linear in eps and blind to
    the structure of the real noise, it overstates the floor by orders of magnitude (module docstring); never
    used for FAIL-at-floor or for stopping."""
    gross = row_l1(lo, di, up)
    rng = np.random.default_rng(seed)
    fl = np.zeros_like(gross)
    for _ in range(samples):
        fl = np.maximum(fl, np.abs(M.solve(eps * gross * rng.choice([-1.0, 1.0], size=gross.shape) * active)))
    return fl


def manifest(base):
    """Executable sha256 and thread count of a pass (run_manifest.json; None when absent)."""
    p = Path(base) / 'run_manifest.json'
    m = json.loads(p.read_text()) if p.exists() else {}
    return dict(executable_sha256=m.get('executable_sha256'), threads=m.get('threads'))


def measured_floor(M, x, R, rep, base, pairs, gross):
    """MEASURED floating-point floor at this state (closing_phase floating_point_floor): per variable, max over the
    identical-state passes P of |M^-1 (R - R_P)| with BASE's own native blocks M (the BA noise between the passes,
    <= 1.4e-9 of row scale in Report 22, enters at second order only). Each P must hold BASE's POPS bit for bit,
    its active (non-FIXPOP) rows and its equation selection. Per pair: build/thread relation to BASE and
    eps_row_l1 = max over active (row, depth) of |R - R_P| / row l1 (the unit of the report-only noise_floor).
    The declared composition: >= 2 pairs, one the same executable at another thread count, one an exact repeat."""
    nt, nd = x.shape
    act = R != 0.0
    mb = manifest(base)
    fl = np.zeros_like(x); recs, fls = [], []
    for p in map(Path, pairs):
        one_record(p / 'QSTARS_RESID_TRACE.bin', nt, nd)
        xp, Rp = read_trace_record(p / 'QSTARS_RESID_TRACE.bin')
        if not np.array_equal(xp, x):
            raise SystemExit(f'{p}: POPS differ from {base}; not an identical-state pass')
        if not np.array_equal(Rp != 0.0, act):
            raise SystemExit(f'{p}: active (non-FIXPOP) rows differ from {base}')
        rp = p / 'QSTARS_REPLACE_TRACE'
        if rp.exists() and not np.array_equal(read_replace(rp, nd), rep):
            raise SystemExit(f'{rp}: equation selection differs from {base}')
        f = np.abs(M.solve(R - Rp))
        fl = np.maximum(fl, f); fls.append(f)
        mp = manifest(p)
        recs.append(dict(path=str(p.resolve()), **mp,
                         same_build=mp['executable_sha256'] is not None and mp['executable_sha256'] == mb['executable_sha256'],
                         same_threads=mp['threads'] is not None and mp['threads'] == mb['threads'],
                         residual_identical=bool(np.array_equal(R, Rp)),
                         eps_row_l1=float(np.max(np.abs(R - Rp)[act] / gross[act]))))
    comp = dict(n_pairs=len(recs), same_build_other_threads=any(r['same_build'] and not r['same_threads'] for r in recs),
                exact_repeat=any(r['same_build'] and r['same_threads'] for r in recs))
    comp['as_declared'] = comp['n_pairs'] >= 2 and comp['same_build_other_threads'] and comp['exact_repeat']
    return fl, dict(base_manifest=mb, pairs=recs, composition=comp), fls


def structural_step(x, R, lo, di, up, rep, sig_frac):
    """Trace equilibrium at the fixed significant state: solve M[E_T, T] s_T = R[E_T] for the structural
    pairing as ONE sparse system over all depths with the native lo/up blocks (advection couples He III /
    He II rows across depths), row/column equilibrated. Returns (s, solved, unpaired, componentwise
    backward error). Exact for the linear trace subsystem (fixed J, n_e, T and significant populations)."""
    import scipy.sparse as sp, scipy.sparse.linalg as spl
    nt, nd = x.shape
    row_of, kind, unpaired = structural_pairing(fractions(x), R != 0.0, rep, sig_frac)
    T = row_of >= 0
    idx = -np.ones((nt, nd), int); idx[T] = np.arange(int(T.sum())); n = int(T.sum())
    s = np.zeros_like(x)
    if n == 0:
        return s, T, unpaired, 0.0
    I, J, V, b = [], [], [], np.zeros(n)
    for v, d in zip(*np.nonzero(T)):
        r = row_of[v, d]; e = idx[v, d]; b[e] = R[r, d]
        for blk, dd in ((lo, d - 1), (di, d), (up, d + 1)):
            if 0 <= dd < nd:
                cols = np.nonzero(T[:, dd])[0]; vals = blk[d][r, cols]; nz = vals != 0
                I += [e] * int(nz.sum()); J += idx[cols[nz], dd].tolist(); V += vals[nz].tolist()
    A = sp.csr_matrix((V, (I, J)), shape=(n, n))
    rs = 1.0 / np.maximum(abs(A).max(axis=1).toarray().ravel(), 1e-300)
    cs = 1.0 / np.maximum(abs(sp.diags(rs) @ A).max(axis=0).toarray().ravel(), 1e-300)
    z = cs * spl.spsolve((sp.diags(rs) @ A @ sp.diags(cs)).tocsc(), rs * b)
    bw = float(np.max(np.abs(A @ z - b) / (abs(A) @ np.abs(z) + np.abs(b) + 1e-300)))
    s[T] = z[idx[T]]
    return s, T, unpaired, bw


def summarize(y, cls, frac, fl=None, threshold=1e-3, model=None):
    """Per class: the gate value max audit_normalized(y) with its location (gate_*), and the largest |y| with its
    location (max_abs_*). With the MEASURED floor fl a variable above the threshold whose floor, audit-normalized
    the same way (fl/(1-fl), the decrease side), also exceeds it is FAIL-at-floor: counted, listed for class 'all',
    still a FAIL. The random-sign model is reported only (model_floor_max)."""
    out = {}
    ya = audit_normalized(y)
    fa = audit_normalized(fl) if fl is not None else None
    for name, m in cls.items():
        if not m.any():
            out[name] = dict(n=0)
            continue
        g = np.unravel_index(int(np.argmax(np.where(m, ya, -1.0))), ya.shape)
        a = np.unravel_index(int(np.argmax(np.where(m, np.abs(y), -1.0))), y.shape)
        out[name] = dict(n=int(m.sum()), audit_normalized=float(ya[g]), gate_row=int(g[0]) + 1, gate_depth=int(g[1]) + 1,
                         gate_frac=float(frac[g]), gate_y=float(y[g]), max_abs=float(abs(y[a])), max_abs_row=int(a[0]) + 1,
                         max_abs_depth=int(a[1]) + 1, max_abs_frac=float(frac[a]),
                         n_above_threshold=int(np.sum(ya[m] > threshold)))
        if fa is not None:
            above = m & (ya > threshold)
            at = above & (fa > threshold)
            out[name].update(floor_max=float(fl[m].max()), floor_audit_normalized_max=float(fa[m].max()),
                             n_fail_at_floor=int(at.sum()), n_fail_above_floor=int(np.sum(above & ~at)))
            if name == 'all':
                out[name]['fail_at_floor'] = [[int(i) + 1, int(d) + 1, float(ya[i, d]), float(fa[i, d])]
                                              for i, d in zip(*np.nonzero(at))]      # row, depth, gate value, floor
        if model is not None:
            out[name]['model_floor_max'] = float(model[m].max())
    return out


def check(base, sig_frac=1e-10, replace=None, noise_eps=None, noise_samples=8, threshold=1e-3, cleanup=False,
          noise_pairs=()):
    """Classes of y = M^-1 R at a base pass; with cleanup=True also the linear prediction after the
    structural trace cleanup (y - (0, s_T) = M^-1 (R - M s): exact for the trace subsystem, J feedback of H-
    ignored). noise_pairs: identical-state passes for the measured floor (measured_floor)."""
    base = Path(base)
    x, R = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    nt, nd = x.shape
    if nt != NT_HHE:
        raise SystemExit(f'{base}: NT={nt}; the class layout here is the 134-variable H/He model')
    one_record(base / 'QSTARS_RESID_TRACE.bin', nt, nd)
    lo, di, up = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', nt, nd)
    rp = base / 'QSTARS_REPLACE_TRACE'
    rep = read_replace(rp if rp.exists() else replace, nd)
    if replace is not None and rp.exists() and not np.array_equal(rep, read_replace(replace, nd)):
        raise SystemExit(f'{rp} differs from the forced table {replace}')
    M = BlockTridiag(lo, di, up).factor()
    y = M.solve(R)
    frac, cls = classes(x, R, di, rep, sig_frac)
    fl = model = None
    if noise_pairs:
        fl, floor, fls = measured_floor(M, x, R, rep, base, noise_pairs, row_l1(lo, di, up))
        for r, f in zip(floor['pairs'], fls):
            r['class_max'] = {n: float(f[m].max()) if m.any() else 0.0 for n, m in cls.items()}
        floor['use'] = ('declared floor (closing_phase floating_point_floor)' if floor['composition']['as_declared'] else
                        'diagnostic only: pairs are not the declared composition (>= 2: same build at another thread '
                        'count, and an exact repeat)')
    if noise_eps:
        model = noise_floor(M, lo, di, up, R != 0.0, noise_eps, noise_samples)
    out = dict(base=str(base.resolve()), nt=nt, nd=nd, sig_frac=sig_frac, threshold=threshold,
               classes=summarize(y, cls, frac, fl, threshold, model))
    if fl is not None:
        out['floor'] = floor
    if model is not None:
        out['noise_model'] = dict(eps_row_l1=noise_eps, samples=noise_samples, seed=0,
                                  use='report only; never used for FAIL-at-floor or stopping')
    if cleanup:
        s, solved, unpaired, bw = structural_step(x, R, lo, di, up, rep, sig_frac)
        out['after_structural_cleanup'] = dict(trace_solved=int(solved.sum()), unpaired=int(unpaired.sum()),
                                               max_trace_step=float(np.max(np.abs(s))), backward_error=bw,
                                               classes=summarize(y - s, cls, frac, fl, threshold, model))
    return out, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('base', type=Path, help='base pass (QSTARS_RESID_TRACE.bin, QSTARS_BA_TRACE.bin, QSTARS_REPLACE_TRACE)')
    ap.add_argument('--sig-frac', type=float, default=1e-10)
    ap.add_argument('--replace', type=Path, default=None, help='forced table (checked against the pass record)')
    ap.add_argument('--out', type=Path, default=None, help='write the JSON report here')
    ap.add_argument('--noise-pair', type=Path, nargs='+', default=[], metavar='PASS',
                    help='identical-state passes (same build at another thread count, an exact repeat): the MEASURED '
                         'floor used for FAIL-at-floor (measured_floor)')
    ap.add_argument('--noise-eps', type=float, default=None,
                    help='report-only random-sign noise model with this EPS in row-l1 units (noise_floor); never used '
                         'for FAIL-at-floor or stopping')
    ap.add_argument('--noise-samples', type=int, default=8)
    ap.add_argument('--predict-cleanup', action='store_true',
                    help='also predict the classes after the structural trace cleanup (closing_phase criterion)')
    a = ap.parse_args()
    thr = json.loads((ROOT / 'benchmarks/cmfgen_validation_plan.json').read_text())['thresholds']['iteration_max_fractional_correction']
    out, _ = check(a.base, a.sig_frac, a.replace, a.noise_eps, a.noise_samples, thr, a.predict_cleanup, a.noise_pair)
    out['prediction_only'] = ('native fixed-J correction of the next complete iteration at this state; the gate is '
                              'judged on the verify pass, never on this prediction')
    for key, blk in [('below_threshold', out)] + ([('below_threshold_after_cleanup', out['after_structural_cleanup'])]
                                                  if a.predict_cleanup else []):
        out[key] = {k: v.get('audit_normalized', 0.0) <= thr for k, v in blk['classes'].items()}
    txt = json.dumps(out, indent=1)
    if a.out:
        a.out.write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
