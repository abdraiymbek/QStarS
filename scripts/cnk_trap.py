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
        A = di[d][np.ix_(idx, idx)]
        r = 1.0 / np.maximum(np.max(np.abs(A), axis=1), 1e-300)
        Ar = A * r[:, None]
        c = 1.0 / np.maximum(np.max(np.abs(Ar), axis=0), 1e-300)
        As = Ar * c[None, :]
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
