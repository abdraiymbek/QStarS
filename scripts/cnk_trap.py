#!/usr/bin/env python3
"""Report 23: radiatively trapped levels in the cool Newton-Krylov solve (helpers for cool_newton_krylov.py).

The direct (fixed-J) Jacobian blocks hold the full spontaneous decay of a resonance-line upper level and the full
absorption out of its lower level (H I n = 2 -> 1 at A21: -1.0e4 on the n = 2 diagonal at depth 33 of the
167-radius cool model). In an optically thick line the radiation response cancels both (exact n = 2 diagonal
-4.4 there, measured with exact products). A preconditioner made from the fixed-J blocks therefore understates
the error of such a level by the trapping factor (~2000 for Lya at the H I front): M^-1 R looks converged while
the level is far from statistical equilibrium, and no Krylov direction built from M^-1 R contains its mode.

Three independent pieces, each off unless the driver option is given:

1. probe_directions / install_columns (--probe-rows, --probe-stride): exact block-tridiagonal columns of chosen
   levels from exact products. One direction perturbs the level at every stride-th depth; the response at the
   perturbed depth and its two neighbours replaces that level's column in the PRECONDITIONER only (the products
   and the Jacobian action stay exact). The response of other perturbed depths contaminates a column by the
   nonlocal part of the response (measured for H I n = 2: 5 % median, 13 % worst at stride 4).
2. coupling_significance (--sig-coupling): a level is significant when a unit relative error in it changes a
   significant equation by more than the given fraction of that equation's gross scale, whatever its population
   fraction. Such levels get least-squares weight 1 and enter the convergence test (a stricter test).
3. update_factor (--update-map linear): the trial state x * f(s). The exponential map exp(-s) overshoots a
   large increase (s << 0) by exp(|s|)/(1+|s|); rate equations are linear in the populations, so the Newton
   update of an increase is x (1 - s). Decreases beyond s = 0.5 continue as 0.5 exp(-(s - 0.5)/0.5), which is
   positive and C1.
"""
import numpy as np


def probe_directions(rows, active, stride):
    """Directions for exact columns: for each 1-based row r and colour c, a unit relative perturbation of that
    variable at every active depth d with d % stride == c. Returns (directions, meta) with meta = (r, c)."""
    if stride < 1:
        raise ValueError('probe stride must be >= 1')
    nt, nd = active.shape
    depth = np.arange(nd)
    dirs, meta = [], []
    for r in rows:
        if not 1 <= r <= nt:
            raise ValueError(f'probe row {r} outside 1..{nt}')
        for c in range(stride):
            sel = (depth % stride == c) & active[r - 1]
            if sel.any():
                z = np.zeros((nt, nd)); z[r - 1, sel] = 1.0
                dirs.append(z); meta.append((r, c))
    return dirs, meta


def install_columns(lo, di, up, dirs, jzs, meta, stride):
    """Replace, in place, the columns of the probed variables in the block-tridiagonal matrix
    (y[:, k] = lo[k] s[:, k-1] + di[k] s[:, k] + up[k] s[:, k+1]) by the exact responses jz = J z.
    The neighbour blocks are replaced only when stride >= 3 (the neighbours of a perturbed depth are then
    unperturbed in the same direction). Returns the number of columns replaced."""
    nd = di.shape[0]
    n = 0
    for z, jz, (r, _c) in zip(dirs, jzs, meta):
        c = r - 1
        for d in np.nonzero(z[c])[0]:
            di[d][:, c] = jz[:, d]
            if stride >= 3:
                if d + 1 < nd:
                    lo[d + 1][:, c] = jz[:, d + 1]
                if d - 1 >= 0:
                    up[d - 1][:, c] = jz[:, d - 1]
            n += 1
    return n


def coupling_significance(di, scale_inv, sig, npop):
    """c[l, d] = max over significant equations i != l at depth d of |di[d][i, l]| * scale_inv[i, d]: the change
    of the most affected significant equation, relative to its gross scale, per unit relative change of
    variable l at the same depth. Population variables l < npop only (others return 0)."""
    nd, nt, _ = di.shape
    a = np.abs(di) * (scale_inv * sig).T[:, :, None]          # (nd, row i, column l)
    idx = np.arange(nt)
    a[:, idx, idx] = 0.0
    c = a.max(axis=1).T                                       # (column l, depth)
    c[npop:] = 0.0
    return c


def update_factor(s, mode, cap):
    """Multiplicative update x_new = x * f(s) for the Newton step s (native convention x <- x (1 - s)).
    exp   : f = exp(-clip(s, -cap, cap)) (the map used through Report 22).
    linear: f = 1 - s for s <= 0.5, 0.5 exp(-(s - 0.5)/0.5) beyond; limited to [exp(-cap), exp(cap)]."""
    s = np.asarray(s, dtype=float)
    if mode == 'exp':
        return np.exp(-np.clip(s, -cap, cap))
    if mode != 'linear':
        raise ValueError(f'unknown update map {mode!r}')
    f = np.where(s <= 0.5, 1.0 - s, 0.5 * np.exp(-(np.maximum(s, 0.5) - 0.5) / 0.5))
    return np.clip(f, np.exp(-cap), np.exp(cap))


def equilibrate(A):
    """Rows scaled by their largest entry, then columns (as BlockTridiag.factor does): (scaled A, r, c)."""
    r = 1.0 / np.maximum(np.max(np.abs(A), axis=1), 1e-300)
    Ar = A * r[:, None]
    c = 1.0 / np.maximum(np.max(np.abs(Ar), axis=0), 1e-300)
    return Ar * c[None, :], r, c


def block_conditioning(di, active, window):
    """Smallest / largest singular value of each equilibrated diagonal block (active rows and columns) at the
    depths of window (d0, d1), 0-based inclusive."""
    out = []
    for d in range(window[0], window[1] + 1):
        idx = np.nonzero(active[:, d])[0]
        S = np.linalg.svd(equilibrate(di[d][np.ix_(idx, idx)])[0], compute_uv=False)
        out.append(S[-1] / S[0])
    return np.array(out)


def floor_blocks(di, tau, active=None):
    """Report 24: bound the amplification of nearly singular diagonal blocks of a preconditioner.
    Each block is equilibrated as BlockTridiag.factor does (rows by their largest entry, then columns), its
    singular values below tau * (largest singular value) are raised to that floor with the singular vectors kept,
    and the equilibration is undone. Well-conditioned directions are unchanged; a nearly singular direction (the
    local-ALI blocks at the ionization fronts) is amplified by at most 1/tau. Rows/columns of inactive variables
    (active[:, d] False) are left out of the decomposition and kept as they are. Returns (new blocks, number of
    singular values raised)."""
    out = di.copy()
    nraised = 0
    nd, nt, _ = di.shape
    for d in range(nd):
        idx = np.arange(nt) if active is None else np.nonzero(active[:, d])[0]
        if idx.size == 0:
            continue
        As, r, c = equilibrate(di[d][np.ix_(idx, idx)])
        U, S, Vt = np.linalg.svd(As)
        floor = tau * S[0]
        low = S < floor
        if low.any():
            nraised += int(low.sum())
            As2 = (U * np.maximum(S, floor)) @ Vt
            out[d][np.ix_(idx, idx)] = As2 / c[None, :] / r[:, None]
    return out, nraised


# Report 24: exact linear constraints of the H/He model (0-based rows of the 134-row state).
H_ROWS = slice(1, 32)          # H I levels and H II (H- in row 0 is a separate trace species, left as it is)
HE_ROWS = slice(32, 132)       # He I, He II levels and He III
HII_ROW, HEII_ROWS, HEIII_ROW, HMI_ROW, NE_ROW = 31, slice(101, 131), 131, 0, 132


def project_conservation(x, totals):
    """Report 24: project a state onto the exact linear constraints of the model: at every depth the H and He
    populations are scaled by one common factor per species so that their totals equal `totals` (H, He; the
    totals of a state that satisfies the conservation rows, e.g. x0), and n_e is set to the charge sum
    n(H II) + sum n(He II) + 2 n(He III) - n(H-). Population ratios within a species are unchanged; T and H- are
    not touched. The component limiter clips each population's change separately, so a Newton step that keeps
    the (linear) number and charge conservation rows satisfied no longer does after clipping; repeated clipped
    steps drift off the constraint surface. Returns (projected state, dict of the largest |factor - 1|)."""
    xp = np.array(x, dtype=float, copy=True)
    h_tot, he_tot = totals
    fh = h_tot / xp[H_ROWS].sum(0)
    fhe = he_tot / xp[HE_ROWS].sum(0)
    xp[H_ROWS] *= fh
    xp[HE_ROWS] *= fhe
    ne = xp[HII_ROW] + xp[HEII_ROWS].sum(0) + 2.0 * xp[HEIII_ROW] - xp[HMI_ROW]
    fne = ne / xp[NE_ROW]
    xp[NE_ROW] = ne
    stats = dict(h=float(np.max(np.abs(fh - 1.0))), he=float(np.max(np.abs(fhe - 1.0))),
                 ne=float(np.max(np.abs(fne - 1.0))))
    return xp, stats


def conservation_totals(x):
    """(H total, He total) per depth of a state, the targets of project_conservation."""
    return x[H_ROWS].sum(0), x[HE_ROWS].sum(0)


def front_window(x, hii_frac=0.5, margin=2):
    """Report 25: depths (0-based, inclusive) of the hydrogen ionization fronts of a state: from the first to the last
    depth where H II holds less than hii_frac of the hydrogen, widened by margin depths on each side. None if
    hydrogen is ionized beyond hii_frac everywhere."""
    f = x[HII_ROW] / x[H_ROWS].sum(0)
    neutral = np.nonzero(f < hii_frac)[0]
    if neutral.size == 0:
        return None
    return max(0, int(neutral[0]) - margin), min(x.shape[1] - 1, int(neutral[-1]) + margin)


def front_rows(inc, base, sig, active, window, tol, always=()):
    """Report 25: 1-based variables whose preconditioner columns need exact responses inside `window` (d0, d1): the
    `always` rows plus every column whose local-ALI increment inc (diagonal blocks, (ND, NT, NT)) changes a
    significant equation at some window depth by more than tol of that equation's gross direct scale (row l1-norm
    of `base`). These are the columns with a radiation response: H I and H-, and ultra-trace upper levels of
    optically thick resonance lines (He I n 1P at the H front), whose fixed-J columns hold untrapped rates.
    Rows that are inactive at every window depth are left out."""
    d0, d1 = window
    rows = set(always)
    for d in range(d0, d1 + 1):
        rs = np.abs(base[d]).sum(axis=1)
        a = np.abs(inc[d]) / np.maximum(rs, 1e-300)[:, None]
        colinc = (a * sig[:, d][:, None]).max(axis=0)
        rows.update(int(r) + 1 for r in np.nonzero(active[:, d] & (colinc > tol))[0])
    return sorted(r for r in rows if active[r - 1, d0:d1 + 1].any())


def window_mask(active, window):
    """active restricted to the depths of window (d0, d1), 0-based inclusive."""
    m = np.zeros_like(active, dtype=bool)
    m[:, window[0]:window[1] + 1] = active[:, window[0]:window[1] + 1]
    return m


def dense_cells(windows, active):
    """Report 26: the cells of a dense exact sub-domain Jacobian. windows: [((d0, d1), rows), ...] with 0-based inclusive
    depths and 1-based rows; every active (row - 1, depth) of each window, each cell once (overlapping windows share
    cells), in window, row, depth order. One direction per cell (cell_direction) gives that variable's full nonlocal
    column, with no colouring and so no aliasing of the response of other perturbed depths."""
    cells, seen = [], set()
    for (d0, d1), rows in windows:
        for r in rows:
            if not 1 <= r <= active.shape[0]:
                raise ValueError(f'dense row {r} outside 1..{active.shape[0]}')
            for d in range(max(0, d0), min(active.shape[1] - 1, d1) + 1):
                if active[r - 1, d] and (r - 1, d) not in seen:
                    seen.add((r - 1, d)); cells.append((r - 1, d))
    return cells


def cell_direction(cell, shape):
    """Unit relative perturbation of one variable at one depth (a dense_cells entry)."""
    z = np.zeros(shape); z[cell] = 1.0
    return z


def quadratic_term(lams, E):
    """Report 26 (--geodesic): the second-order term Q of the residual along a trial direction, from the trials'
    remainders E_lam = R(x_lam) - (R - lam J s) (the actual residual minus its linear model) at step fractions lams.
    E_lam = lam^2 Q + lam^3 C + ...: per cell, the intercept at lam = 0 of the straight-line fit of E_lam / lam^2
    against lam (one trial: E / lam^2)."""
    lams = np.asarray(lams, dtype=float)
    if lams.size == 0:
        raise ValueError('quadratic_term needs at least one trial')
    e = [Ei / li ** 2 for Ei, li in zip(E, lams)]
    if lams.size == 1 or np.ptp(lams) == 0:
        return sum(e) / len(e)
    c = np.linalg.pinv(np.stack([np.ones_like(lams), lams], axis=1))[0]    # intercept weights
    return sum(ci * ei for ci, ei in zip(c, e))


def slow_modes(di, active, window, tau, rows=None, min_weight=0.5, exclude=()):
    """Report 25: the slow modes of the diagonal blocks inside window (d0, d1), 0-based inclusive. Per depth, every
    singular pair of the equilibrated block (active rows and columns except `exclude`) with singular value below tau x
    the largest, and (rows given) with at least min_weight of its left vector's squared norm on those rows. In physical
    scaling, as (NT, ND) arrays nonzero at that depth only: Y = r*u (the block inverse amplifies the residual component
    Y.x by 1/s) and Z = c*v (the variable direction it amplifies), each scaled to max |entry| = 1. Returns
    (Y, Z, [(depth, s/s_max), ...])."""
    nd, nt, _ = di.shape
    Y, Z, info = [], [], []
    for d in range(window[0], window[1] + 1):
        idx = np.nonzero(active[:, d])[0]
        idx = idx[~np.isin(idx, exclude)]
        As, r, c = equilibrate(di[d][np.ix_(idx, idx)])
        U, S, Vt = np.linalg.svd(As)
        for j in np.nonzero(S < tau * S[0])[0]:
            if rows is not None and np.sum(U[np.isin(idx, rows), j] ** 2) < min_weight:
                continue
            y = np.zeros((nt, nd)); z = np.zeros((nt, nd))
            yy, zz = r * U[:, j], c * Vt[j]
            y[idx, d] = yy / np.max(np.abs(yy)); z[idx, d] = zz / np.max(np.abs(zz))
            Y.append(y); Z.append(z); info.append((d, float(S[j] / S[0])))
    return Y, Z, info


class CoarseCorrected:
    """Report 25: two-level preconditioner P^-1 r = Z c + M^-1 (r - (J Z) c) with the EXACT products J Z of the slow
    modes (slow_modes). P^-1 J z = z for every z in span(Z), whatever the nonlocal coupling of the modes.
      mode 'ls'       : c = argmin ||W M^-1 (r - J Z c)|| (the fine-level preconditioned residual is minimal); never
                        singular while the columns W M^-1 J z are independent;
      mode 'galerkin' : c = E^+ Y^T r, E = Y^T J Z: the fine level only sees residuals with Y^T r = 0, the components
                        its nearly singular blocks would amplify, but E is singular when slow modes couple only through
                        their neighbours.
    M: anything with .solve of an (NT, ND) array; W: weights of the least squares (mode 'ls')."""

    def __init__(self, M, Y, Z, JZ, W=None, mode='ls', rcond=1e-12):
        self.fine, self.mode = M, mode
        self.Z = np.stack([z.ravel() for z in Z])
        self.JZ = np.stack([j.ravel() for j in JZ])
        if mode == 'galerkin':
            self.Y = np.stack([y.ravel() for y in Y])
            E = self.Y @ self.JZ.T
        elif mode == 'ls':
            self.W = np.ones(JZ[0].shape) if W is None else W
            E = np.stack([(self.W * M.solve(j)).ravel() for j in JZ], axis=1)
        else:
            raise ValueError(f'unknown coarse mode {mode!r}')
        sv = np.linalg.svd(E, compute_uv=False)
        self.cond = float(sv[0] / sv[-1])
        self.E_pinv = np.linalg.pinv(E, rcond=rcond)

    def coefficients(self, r):
        if self.mode == 'galerkin':
            return self.E_pinv @ (self.Y @ r.ravel())
        return self.E_pinv @ (self.W * self.fine.solve(r)).ravel()

    def solve(self, r):
        c = self.coefficients(r)
        return (self.Z.T @ c).reshape(r.shape) + self.fine.solve(r - (self.JZ.T @ c).reshape(r.shape))


def trace_ion_balance(R, Js, Jdi, frac, lam, row=101, col=131, sig_frac=1e-10, share=0.1):
    """Report 25: the step component of a trace upper ion that balances a replaced ground-state row. Where the He II
    ground row is CMFGEN's He II / He III ionization balance (the He III column carries at least `share` of the row's
    direct l1-norm) and He III is not significant (fraction <= sig_frac of He, least-squares weight ~0, so no Krylov
    direction moves it), the row cannot be satisfied by the solve, and the preconditioner turns its imbalance into
    He II ground corrections (s 17-61 at depths 42-45). The balance is linear in the trace ion: its component is set so
    that the linearized row vanishes for the trial step lam*s, s_col = (R_row - lam (J s)_row) / J[row, col].
    R, Js: (NT, ND) residual and exact action of the unit step (without any PTC shift); Jdi: (ND, NT, NT) direct
    diagonal blocks; frac: (NT, ND) fraction of the element. Returns (mask over depths, s_col)."""
    a = Jdi[:, row, col]
    l1 = np.abs(Jdi[:, row, :]).sum(axis=1)
    mask = (frac[col] <= sig_frac) & (np.abs(a) >= share * np.maximum(l1, 1e-300)) & (R[row] != 0.0)
    s_col = np.zeros(R.shape[1])
    s_col[mask] = (R[row, mask] - lam * Js[row, mask]) / a[mask]
    return mask, s_col


def adapt_caps(caps, prev_logf, new_logf, default, shrink=0.5, grow=1.5, floor=0.01, hit=0.9):
    """Report 25: per-component limits |log f| <= caps for the next update. A component whose limited update reversed
    direction from the previous accepted update (a limit cycle across its root: He II ground at the inner H front
    jumped between 1.5e-6 and 3.6e-6 of He for six iterates with steps +32, -3, +54, -5, ...) gets its limit
    multiplied by `shrink` (not below `floor`); every other component's limit grows back by `grow` up to `default`.
    prev_logf/new_logf: log(x_new / x_base) of the last two accepted updates. Returns the new caps."""
    flip = (np.sign(prev_logf) * np.sign(new_logf) < 0) & (np.abs(new_logf) >= hit * caps)
    return np.where(flip, np.maximum(caps * shrink, floor), np.minimum(caps * grow, default)), int(flip.sum())


ACCEPTED_POINTER = 'accepted_state.json'


def trial_decision(trials, phi):
    """V1 plan P0: the outcome of one iterate's line search. trials: (label, pass dir, merit) of every trial pass, phi:
    the base's merit in the same measure. The best trial is accepted only when it lowers the merit strictly; the history
    record keeps the best trial apart from the accepted one (iterate 16 of free_r25 recorded its best trial l0.125_d16,
    merit 2.83717e-2 against 2.83713e-2, as `accepted` although the driver rejected every trial and stopped)."""
    lab, d, ph = min(trials, key=lambda t_: t_[2])
    ok = bool(np.isfinite(ph) and ph < phi)
    return dict(best_trial=lab, best_merit=float(ph), best_dir=d, accepted=lab if ok else None,
                status='accepted' if ok else 'rejected')


def history_accepted(rec):
    """The accepted trial label of a history record, None when every trial was rejected. Records written before the
    `status` field (V1 plan P0) named their best trial `accepted` whatever its merit: judged here on merit_new < merit."""
    if 'status' in rec:
        return rec.get('accepted')
    lab = rec.get('accepted')
    if lab is None or not (rec.get('merit_new', np.inf) < rec.get('merit', -np.inf)):
        return None
    return lab


def write_accepted_pointer(out, it, base, merit, source):
    """V1 plan P0: the explicit restart pointer, OUT/accepted_state.json: the accepted state the next iterate starts
    from (written atomically). source: how it became the accepted state."""
    import json, os, time
    from pathlib import Path
    p = Path(out) / ACCEPTED_POINTER
    tmp = p.with_suffix('.tmp')
    tmp.write_text(json.dumps(dict(it=int(it), base=str(Path(base).resolve()), merit=float(merit), source=source,
                                   written=time.strftime('%Y-%m-%d %H:%M:%S')), indent=1))
    os.replace(tmp, p)


def read_accepted_pointer(out):
    import json
    from pathlib import Path
    p = Path(out) / ACCEPTED_POINTER
    return json.loads(p.read_text()) if p.exists() else None


def check_restart_base(base):
    """V1 plan P0: None if `base` may start a native solve, else the reason it may not. A pass named itNN_trial_LAB is a
    trial of iterate NN-1 of its own run directory; it is admissible only if that run's history accepted LAB there
    (history_accepted). Bases outside such a run (it00_base, start_base, another run's pass) are not judged here."""
    import json, re
    from pathlib import Path
    base = Path(base).resolve()
    m = re.fullmatch(r'it(\d+)_trial_(.+)', base.name)
    hp = base.parent / 'history.json'
    if m is None or not hp.exists():
        return None
    k, lab = int(m.group(1)) - 1, m.group(2)
    recs = [r for r in json.loads(hp.read_text()) if r.get('it') == k and 'trials' in r]
    if not recs:
        return f'{base.name}: no history record of iterate {k} in {hp}'
    acc = history_accepted(recs[-1])
    if acc != lab:
        return (f'{base.name}: iterate {k} of {hp.parent.name} did not accept trial {lab} '
                f'(accepted: {acc}; best trial merit {recs[-1].get("merit_new")} against base {recs[-1].get("merit")})')
    return None


def project_linear(x, z):
    """V1 plan P1: the linearization of project_conservation at a state x that satisfies it, for a fractional step z
    (x -> x(1 - z)): H and He rows shifted by their total's relative change, n_e set to the charge-sum change. A
    direction replaced by project_linear(x, z) changes the number and charge rows only at second order, so the
    projection of a trial no longer applies an unscored n_e change (at free_r25 iterate 16: an n_e response 16x the
    merit in the H ionization-balance row)."""
    zp = np.array(z, dtype=float, copy=True)
    zp[H_ROWS] -= (x[H_ROWS] * zp[H_ROWS]).sum(0) / x[H_ROWS].sum(0)
    zp[HE_ROWS] -= (x[HE_ROWS] * zp[HE_ROWS]).sum(0) / x[HE_ROWS].sum(0)
    zp[NE_ROW] = (x[HII_ROW] * zp[HII_ROW] + (x[HEII_ROWS] * zp[HEII_ROWS]).sum(0) + 2.0 * x[HEIII_ROW] * zp[HEIII_ROW]
                  - x[HMI_ROW] * zp[HMI_ROW]) / x[NE_ROW]
    return zp


def orthonormalize_batch(S, dirs, tol=1e-8):
    """V1 plan, council block (--ortho-dirs): new directions made orthonormal (Euclidean, on the full state vector) to
    the span of the existing directions S and to each other before their product pass, so the least squares never
    needs cancelling coefficients (raw v1_p2 spans: 238 directions, 185 independent, coefficients up to 4e14). Two
    rounds of classical Gram-Schmidt against an orthonormal basis of S, then a QR of the batch; a direction whose
    remainder is below tol of its own norm is dropped. Linear combinations of conservation-consistent directions stay
    consistent (cnk_trap.project_linear is linear and idempotent). Returns the list of kept directions."""
    if not dirs:
        return []
    shape = dirs[0].shape
    V = np.stack([d.ravel() for d in dirs], axis=1)
    n0 = np.linalg.norm(V, axis=0)
    if S:
        Q, _ = np.linalg.qr(np.stack([s.ravel() for s in S], axis=1))
        for _ in range(2):
            V = V - Q @ (Q.T @ V)
    out = []
    for j in range(V.shape[1]):
        v = V[:, j].copy()
        for _ in range(2):
            for u in out:
                v -= u * (u @ v)
        nv = np.linalg.norm(v)
        if n0[j] > 0 and nv > tol * n0[j]:
            out.append(v / nv)
    return [u.reshape(shape) for u in out]
