"""Report 23: exact trapped-level columns, coupling significance and the linear update map (scripts/cnk_trap.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import cnk_trap  # noqa: E402
from qstars.nlte_newton import BlockTridiag  # noqa: E402


def banded_system(nt=5, nd=12, seed=3):
    rng = np.random.default_rng(seed)
    lo = 0.1 * rng.normal(size=(nd, nt, nt)); up = 0.1 * rng.normal(size=(nd, nt, nt))
    di = rng.normal(size=(nd, nt, nt)) + 4.0 * np.eye(nt)
    lo[0] = 0.0; up[-1] = 0.0
    return lo, di, up


def test_probe_columns_recover_a_block_tridiagonal_operator_exactly():
    lo, di, up = banded_system()
    exact = BlockTridiag(lo, di, up)
    nd, nt = exact.shape
    active = np.ones((nt, nd), dtype=bool); active[2, 5] = False          # one inactive equation: never perturbed
    # the "fixed-J" preconditioner: columns 2 and 3 (1-based) wrong by a large factor, as for a trapped level
    lo_m, di_m, up_m = lo.copy(), di.copy(), up.copy()
    for c in (1, 2):
        di_m[:, :, c] *= 2000.0; lo_m[:, :, c] = 0.0; up_m[:, :, c] = 0.0
    for stride in (3, 4):
        dirs, meta = cnk_trap.probe_directions([2, 3], active, stride)
        assert len(dirs) == 2 * stride and all(z.sum() == np.count_nonzero(z) for z in dirs)
        assert not any(z[2, 5] for z in dirs)
        jzs = [exact.matvec(z) for z in dirs]
        l2, d2, u2 = lo_m.copy(), di_m.copy(), up_m.copy()
        n = cnk_trap.install_columns(l2, d2, u2, dirs, jzs, meta, stride)
        assert n == 2 * nd - 1
        # a tridiagonal operator has no response beyond the neighbours: the probed columns are exact
        keep = np.ones(nd, dtype=bool); keep[5] = False
        np.testing.assert_allclose(d2[:, :, 1], di[:, :, 1], rtol=0, atol=1e-13)
        np.testing.assert_allclose(d2[keep][:, :, 2], di[keep][:, :, 2], rtol=0, atol=1e-13)
        np.testing.assert_allclose(l2[:, :, 1], lo[:, :, 1], rtol=0, atol=1e-13)
        np.testing.assert_allclose(u2[:, :, 1], up[:, :, 1], rtol=0, atol=1e-13)
        assert d2[5, 0, 2] == di_m[5, 0, 2]                                # the inactive column is left alone
        # untouched columns keep the preconditioner's values
        np.testing.assert_array_equal(d2[:, :, 0], di_m[:, :, 0])


def test_stride_below_three_leaves_neighbour_blocks_alone():
    lo, di, up = banded_system()
    exact = BlockTridiag(lo, di, up)
    nd, nt = exact.shape
    active = np.ones((nt, nd), dtype=bool)
    dirs, meta = cnk_trap.probe_directions([2], active, 1)
    l2, d2, u2 = np.zeros_like(lo), np.zeros_like(di), np.zeros_like(up)
    cnk_trap.install_columns(l2, d2, u2, dirs, [exact.matvec(z) for z in dirs], meta, 1)
    assert not l2.any() and not u2.any()
    # all depths perturbed at once: the diagonal block column holds the summed (contaminated) response
    np.testing.assert_allclose(d2[3, :, 1], lo[3][:, 1] + di[3][:, 1] + up[3][:, 1], atol=1e-13)
    with pytest.raises(ValueError):
        cnk_trap.probe_directions([2], active, 0)
    with pytest.raises(ValueError):
        cnk_trap.probe_directions([9], active, 2)


def test_coupling_significance_flags_a_trace_level_that_drives_a_significant_equation():
    nt, nd, npop = 4, 3, 3                      # variables: major level, trace level a, trace level b, electron density
    di = np.zeros((nd, nt, nt))
    for d in range(nd):
        di[d] = np.diag([10.0, 5.0, 5.0, 1.0])
        di[d][0, 1] = 0.5                       # trace level a changes the major equation by 5 % of its scale
        di[d][0, 2] = 1e-6                      # trace level b does not
        di[d][2, 1] = 99.0                      # coupling into a non-significant equation does not count
        di[d][3, 0] = 7.0                       # column of a non-population variable is ignored (c = 0 beyond npop)
    scale_inv = np.tile(np.array([0.1, 0.2, 0.2, 1.0])[:, None], (1, nd))
    sig = np.zeros((nt, nd), dtype=bool); sig[0] = True; sig[3] = True
    c = cnk_trap.coupling_significance(di, scale_inv, sig, npop)
    assert c.shape == (nt, nd)
    np.testing.assert_allclose(c[1], 0.05)
    np.testing.assert_allclose(c[2], 1e-7)
    np.testing.assert_allclose(c[0], 7.0)       # the major level's own coupling into the electron equation
    assert not c[3].any()
    flagged = (c > 1e-3) & ~sig
    assert flagged[1].all() and not flagged[2].any()


def test_update_maps_agree_to_first_order_and_the_linear_map_is_exact_for_a_linear_equation():
    cap = np.log(1e4)
    s = np.array([-1e-4, 1e-4])
    np.testing.assert_allclose(cnk_trap.update_factor(s, 'linear', cap), cnk_trap.update_factor(s, 'exp', cap), rtol=1e-8)
    np.testing.assert_array_equal(cnk_trap.update_factor(np.array([-3.0, 0.2, 20.0]), 'exp', cap),
                                  np.exp(-np.clip(np.array([-3.0, 0.2, 20.0]), -cap, cap)))
    # rate equation a - b x = 0 at x = x*/35: Newton step s = 1 - x*/x = -34; linear map lands on x*, exp overshoots
    x, xs = 1.0, 35.0
    s = 1.0 - xs / x
    assert cnk_trap.update_factor(s, 'linear', cap) * x == pytest.approx(xs)
    assert cnk_trap.update_factor(s, 'exp', cap) * x == pytest.approx(1e4)
    # decreases: exact up to s = 0.5, then positive, monotone and C1 across the joint
    assert cnk_trap.update_factor(0.3, 'linear', cap) == pytest.approx(0.7)
    f = cnk_trap.update_factor(np.array([0.5 - 1e-7, 0.5, 0.5 + 1e-7, 1.0, 3.0, 50.0]), 'linear', cap)
    assert np.all(np.diff(f) <= 0) and np.all(f > 0)
    assert (f[2] - f[0]) / 2e-7 == pytest.approx(-1.0, rel=1e-5)
    assert f[3] == pytest.approx(0.5 * np.exp(-1.0)) and f[-1] == pytest.approx(1e-4)
    assert cnk_trap.update_factor(-1e9, 'linear', cap) == pytest.approx(1e4)
    with pytest.raises(ValueError):
        cnk_trap.update_factor(0.0, 'log', cap)


def _driver():
    import importlib
    return importlib.import_module('cool_newton_krylov')


def test_iterate_exact_columns_change_only_the_preconditioner_and_extend_the_significant_set():
    """A trapped level: its fixed-J diagonal is 2000x the exact one. M^-1 R then hides its error; with the exact
    column installed the preconditioned residual shows it, the level becomes coupling-significant (it drives the
    major level's equation), and a level that only drives the newly significant one follows in the closure."""
    cnk = _driver()
    old = (cnk.NT, cnk.ND, cnk.SIG_FRAC)
    try:
        nt, nd = 5, 6                     # variables: major, trapped trace, ladder trace, inert trace, n_e (T absent)
        cnk.NT, cnk.ND, cnk.SIG_FRAC = nt + 1, nd, 1e-10      # NT - 2 population variables = 4
        nt = cnk.NT
        exact = np.zeros((nd, nt, nt))
        for d in range(nd):
            exact[d] = np.diag([10.0, 4.0, 6.0, 3.0, 1.0, 1.0])
            exact[d][0, 1] = 2.0          # trapped level feeds the major equation (20 % of its scale)
            exact[d][1, 2] = 3.0          # ladder level feeds the trapped level's equation
            exact[d][1, 0] = 0.5
        zero = np.zeros_like(exact)
        fixed = exact.copy(); fixed[:, 1, 1] *= 2000.0
        it = cnk.Iterate.__new__(cnk.Iterate)
        it.Jdir = cnk.BlockTridiag(zero.copy(), fixed.copy(), zero.copy())
        it._pc = (zero.copy(), fixed.copy(), zero.copy())
        it.M = cnk.BlockTridiag(*[m.copy() for m in it._pc]).factor()
        it.active = np.ones((nt, nd), dtype=bool)
        frac = np.ones((nt, nd)); frac[1] = 1e-13; frac[2] = 1e-14; frac[3] = 1e-15
        it.frac = frac
        it.W = np.where(frac > 1e-10, 1.0, 1e-3)
        R = np.zeros((nt, nd)); R[1] = 40.0                 # the trapped level is 10x out of balance
        assert np.max(np.abs(it.M.solve(R)[1])) < 6e-3      # hidden by the fixed-J diagonal
        dirs, meta = cnk_trap.probe_directions([2], it.active, 3)
        ex = cnk.BlockTridiag(zero, exact, zero)
        jdir_before = it.Jdir.di.copy()
        assert it.install_exact_columns(dirs, [ex.matvec(z) for z in dirs], meta, 3) == nd
        np.testing.assert_array_equal(it.Jdir.di, jdir_before)        # products untouched
        np.testing.assert_allclose(it.M.solve(R), cnk.BlockTridiag(zero, exact, zero).factor().solve(R))
        assert np.all(np.abs(it.M.solve(R)[1]) > 10.0)                 # the true correction is visible
        assert it.extend_significance(0.0) == 0 and it.sig.sum() == 3 * nd and it.coupling is None
        # one round: only the trapped level (coupling 2/12.5 into the major equation)
        assert it.extend_significance(1e-2, rounds=1) == nd and it.sig[1].all() and not it.sig[2].any()
        np.testing.assert_allclose(it.W[1], 1.0)
        # closure: the ladder level drives the now-significant trapped equation; the inert level stays out
        assert it.extend_significance(1e-2) == 2 * nd and it.sig[1].all() and it.sig[2].all() and not it.sig[3].any()
        np.testing.assert_allclose(it.W[2], 1.0); np.testing.assert_allclose(it.W[3], 1e-3)
        assert 'T 0.00e+00@d1' in cnk.group_report(np.zeros((134, 2)), np.ones((134, 2), dtype=bool))
        # explicit rows (--sig-rows): promoted at every active depth, without the coupling rule
        it.W = np.where(frac > 1e-10, 1.0, 1e-3); it.active[3, 0] = False
        assert it.extend_significance(0.0, rows=[4, 1]) == nd - 1 and it.sig[3, 1:].all() and not it.sig[3, 0]
        assert not it.sig[1].any() and it.W[3, 1] == 1.0 and it.W[3, 0] == 1e-3 and it.W[1, 0] == 1e-3
    finally:
        cnk.NT, cnk.ND, cnk.SIG_FRAC = old


def test_lm_step_with_temperature_only_penalty_keeps_population_components():
    """--lm-t-only: penalising only the T row limits |s_T| while the other components still fit the residual."""
    cnk = _driver()
    rng = np.random.default_rng(5)
    nt, nd, k = 134, 4, 6
    S = [rng.normal(size=(nt, nd)) for _ in range(k)]
    for z in S[:3]:
        z[133] = 0.0                              # three directions without a T component
    A = rng.normal(size=(50, k)); b = A @ np.array([1.0, 2.0, -1.0, 3.0, 2.0, -2.0])
    w_all = np.ones((nt, nd)); w_t = np.zeros((nt, nd)); w_t[133] = 1.0
    s_all, mu_all, _ = cnk.lm_step((S, A, b), w_all, 0.05)
    s_t, mu_t, pred_t = cnk.lm_step((S, A, b), w_t, 0.05)
    assert np.max(np.abs(s_all[133])) <= 0.05 + 1e-9 and np.max(np.abs(s_t[133])) <= 0.05 + 1e-9
    # with the T-only penalty the T-free directions keep (nearly) their unregularized size
    nonT = slice(0, 133)
    assert np.linalg.norm(s_t[nonT]) > 3 * np.linalg.norm(s_all[nonT])


def test_floor_blocks_bounds_amplification_and_keeps_well_conditioned_blocks():
    rng = np.random.default_rng(7)
    nd, nt = 3, 6
    di = np.stack([rng.normal(size=(nt, nt)) + 6 * np.eye(nt) for _ in range(nd)])
    U, _, Vt = np.linalg.svd(rng.normal(size=(nt, nt)))
    di[1] = (U * np.array([3.0, 2.0, 1.0, 0.5, 1e-9, 1e-12])) @ Vt * 1e4     # nearly singular, badly scaled block
    out, n = cnk_trap.floor_blocks(di, 1e-2)
    np.testing.assert_array_equal(out[0], di[0]); np.testing.assert_array_equal(out[2], di[2])
    assert n == 2
    A = out[1]; r = 1 / np.abs(A).max(1); As = A * r[:, None]; c = 1 / np.abs(As).max(0); As = As * c
    s = np.linalg.svd(As, compute_uv=False)
    assert s[-1] / s[0] > 1e-3                         # bounded amplification after the floor
    active = np.ones((nt, nd), dtype=bool); active[5, 1] = False
    out2, _ = cnk_trap.floor_blocks(di, 1e-2, active)
    np.testing.assert_array_equal(out2[1][5], di[1][5]); np.testing.assert_array_equal(out2[1][:, 5], di[1][:, 5])


def test_project_conservation_restores_totals_and_charge_and_keeps_ratios():
    rng = np.random.default_rng(11)
    nt, nd = 134, 5
    x0 = rng.uniform(1.0, 10.0, size=(nt, nd))
    x0[cnk_trap.NE_ROW] = (x0[cnk_trap.HII_ROW] + x0[cnk_trap.HEII_ROWS].sum(0) + 2 * x0[cnk_trap.HEIII_ROW]
                           - x0[cnk_trap.HMI_ROW])
    tot = cnk_trap.conservation_totals(x0)
    same, st0 = cnk_trap.project_conservation(x0, tot)          # a state on the constraints is unchanged
    np.testing.assert_allclose(same, x0, rtol=1e-14)
    assert max(st0.values()) < 1e-13
    x = x0 * np.exp(rng.normal(scale=0.7, size=(nt, nd)))       # clipped-step-like drift off the constraints
    xp, st = cnk_trap.project_conservation(x, tot)
    np.testing.assert_allclose(xp[cnk_trap.H_ROWS].sum(0), tot[0], rtol=1e-13)
    np.testing.assert_allclose(xp[cnk_trap.HE_ROWS].sum(0), tot[1], rtol=1e-13)
    ch = xp[cnk_trap.HII_ROW] + xp[cnk_trap.HEII_ROWS].sum(0) + 2 * xp[cnk_trap.HEIII_ROW] - xp[cnk_trap.HMI_ROW]
    np.testing.assert_allclose(xp[cnk_trap.NE_ROW], ch, rtol=1e-14)
    # ratios within a species, H-, and T are unchanged
    np.testing.assert_allclose(xp[2:32] / xp[1], x[2:32] / x[1], rtol=1e-13)
    np.testing.assert_allclose(xp[33:132] / xp[32], x[33:132] / x[32], rtol=1e-13)
    np.testing.assert_array_equal(xp[0], x[0]); np.testing.assert_array_equal(xp[133], x[133])
    assert st['h'] > 0.01 and st['he'] > 0.01 and st['ne'] > 0.01
    again, st2 = cnk_trap.project_conservation(xp, tot)          # idempotent
    np.testing.assert_allclose(again, xp, rtol=1e-14)


def test_front_window_spans_the_neutral_hydrogen_zone_with_a_margin():
    nt, nd = 134, 30
    x = np.ones((nt, nd))
    x[cnk_trap.HII_ROW] = 1e3                          # H II dominates the hydrogen ...
    x[cnk_trap.HII_ROW, 8:12] = 1e-1                   # ... except two neutral zones
    x[cnk_trap.HII_ROW, 15:17] = 1.0
    assert cnk_trap.front_window(x, 0.5, 2) == (6, 18)
    assert cnk_trap.front_window(x, 0.5, 20) == (0, nd - 1)
    x[cnk_trap.HII_ROW] = 1e3
    assert cnk_trap.front_window(x, 0.5, 2) is None


def test_front_rows_take_radiation_responsive_columns_and_front_probes_recover_window_columns():
    lo, di, up = banded_system(nt=6, nd=14)
    exact = BlockTridiag(lo, di, up)
    nd, nt = exact.shape
    active = np.ones((nt, nd), dtype=bool); active[4, :] = False         # row 5 inactive everywhere
    sig = np.zeros((nt, nd), dtype=bool); sig[0] = True                  # equation 1 significant
    inc = np.zeros_like(di)
    inc[6, 0, 3] = 0.5 * np.abs(di[6][0]).sum()                          # column 4 drives equation 1 at depth 7
    inc[12, 0, 2] = 0.5 * np.abs(di[12][0]).sum()                        # column 3: outside the window
    inc[6, 0, 4] = 0.5 * np.abs(di[6][0]).sum()                          # column 5: inactive
    window = (4, 9)
    rows = cnk_trap.front_rows(inc, di, sig, active, window, 1e-4, always=[1])
    assert rows == [1, 4]
    # exact columns in the window only, the rest of the preconditioner unchanged
    lo_m, di_m, up_m = lo.copy(), di.copy(), up.copy()
    for c in (0, 3):
        di_m[:, :, c] *= 50.0; lo_m[:, :, c] = 0.0; up_m[:, :, c] = 0.0
    dirs, meta = cnk_trap.probe_directions(rows, cnk_trap.window_mask(active, window), 3)
    assert len(dirs) == 6 and all(np.all(z[:, :window[0]] == 0) and np.all(z[:, window[1] + 1:] == 0) for z in dirs)
    l2, d2, u2 = lo_m.copy(), di_m.copy(), up_m.copy()
    n = cnk_trap.install_columns(l2, d2, u2, dirs, [exact.matvec(z) for z in dirs], meta, 3)
    assert n == 2 * (window[1] - window[0] + 1)
    w = slice(window[0], window[1] + 1)
    for c in (0, 3):
        np.testing.assert_allclose(d2[w, :, c], di[w, :, c], rtol=0, atol=1e-13)
        np.testing.assert_allclose(l2[window[0] + 1:window[1] + 2, :, c], lo[window[0] + 1:window[1] + 2, :, c],
                                   rtol=0, atol=1e-13)
        np.testing.assert_allclose(u2[window[0] - 1:window[1], :, c], up[window[0] - 1:window[1], :, c],
                                   rtol=0, atol=1e-13)
        np.testing.assert_array_equal(d2[:window[0], :, c], di_m[:window[0], :, c])
        np.testing.assert_array_equal(d2[window[1] + 1:, :, c], di_m[window[1] + 1:, :, c])
    np.testing.assert_array_equal(d2[:, :, [1, 2, 4, 5]], di_m[:, :, [1, 2, 4, 5]])
    cond = cnk_trap.block_conditioning(di, active, window)
    assert cond.shape == (window[1] - window[0] + 1,) and np.all((cond > 0) & (cond <= 1))


def test_dense_cells_give_the_window_newton_step_that_colored_probes_alias():
    """Report 26: one direction per (variable, depth) in each window, each cell once, inactive cells left out. With a
    nonlocal operator (coupling over the whole window, as trapped radiation), the least squares over the single-cell
    columns solves a residual made by window variables exactly; stride-3 colored probes of the same rows cannot (a
    colored direction moves its row at every third depth by one common factor)."""
    nt, nd = 4, 12
    active = np.ones((nt, nd), dtype=bool); active[2, 6] = False
    cells = cnk_trap.dense_cells([((2, 5), [1, 3]), ((4, 7), [3])], active)
    assert cells == [(0, 2), (0, 3), (0, 4), (0, 5), (2, 2), (2, 3), (2, 4), (2, 5), (2, 7)]
    with pytest.raises(ValueError):
        cnk_trap.dense_cells([((0, 1), [nt + 1])], active)
    z = cnk_trap.cell_direction((2, 7), (nt, nd))
    assert z[2, 7] == 1.0 and z.sum() == 1.0
    rng = np.random.default_rng(11)
    n = nt * nd
    J = np.diag(rng.uniform(1.0, 2.0, n)) + 0.1 * rng.normal(size=(n, n))   # nonlocal: every cell couples to every cell
    matvec = lambda v: (J @ v.ravel()).reshape(nt, nd)                     # noqa: E731
    win = [((2, 7), [1, 3])]
    cells = cnk_trap.dense_cells(win, active)
    t = rng.normal(size=len(cells))
    s_true = np.zeros((nt, nd))
    for c, tc in zip(cells, t):
        s_true[c] = tc
    R = matvec(s_true)
    S = [cnk_trap.cell_direction(c, (nt, nd)) for c in cells]
    A = np.stack([matvec(v).ravel() for v in S], axis=1)
    np.testing.assert_allclose(A[:, 0], J[:, np.ravel_multi_index(cells[0], (nt, nd))], rtol=0, atol=0)
    y, *_ = np.linalg.lstsq(A, R.ravel(), rcond=1e-12)
    assert np.linalg.norm(A @ y - R.ravel()) < 1e-10 * np.linalg.norm(R)
    np.testing.assert_allclose(sum(yi * v for yi, v in zip(y, S)), s_true, atol=1e-10)
    pz, _ = cnk_trap.probe_directions([1, 3], cnk_trap.window_mask(active, (2, 7)), 3)
    Ap = np.stack([matvec(v).ravel() for v in pz], axis=1)
    yp, *_ = np.linalg.lstsq(Ap, R.ravel(), rcond=1e-12)
    assert np.linalg.norm(Ap @ yp - R.ravel()) > 0.3 * np.linalg.norm(R)


def test_quadratic_term_recovers_the_second_order_term_and_the_geodesic_path_cancels_it():
    """Report 26 (--geodesic): from the remainders of three trials of R(lam) = R0 - lam J s + lam^2 Q + lam^3 C the
    intercept fit returns Q exactly, one trial gives E / lam^2, and the path lam s + lam^2 t with J t = Q leaves only
    third-order terms of a quadratic residual map."""
    rng = np.random.default_rng(4)
    Q, C = rng.normal(size=(3, 5)), rng.normal(size=(3, 5))
    lams = [0.5, 0.25, 0.125]
    E = [l ** 2 * Q + l ** 3 * C for l in lams]
    np.testing.assert_allclose(cnk_trap.quadratic_term(lams, E), Q, atol=1e-12)
    np.testing.assert_allclose(cnk_trap.quadratic_term([0.5], [E[0]]), Q + 0.5 * C, atol=1e-12)
    with pytest.raises(ValueError):
        cnk_trap.quadratic_term([], [])
    # quadratic map F(u) = R0 - J u + B(u, u) along u = lam s + lam^2 t
    n = 6
    J = np.eye(n) * 2.0 + 0.3 * rng.normal(size=(n, n))
    Bt = rng.normal(size=(n, n, n)) * 0.2
    F = lambda u: R0 - J @ u + np.einsum('ijk,j,k->i', Bt, u, u)       # noqa: E731
    R0 = rng.normal(size=n)
    s = np.linalg.solve(J, R0)
    Es = [F(l * s) - (R0 - l * J @ s) for l in lams]
    t = np.linalg.solve(J, cnk_trap.quadratic_term(lams, Es))
    plain, geo = [], []
    for l in (0.1, 0.05):
        plain.append(np.linalg.norm(F(l * s) - (R0 - l * J @ s)))           # second order: x4 per halving
        geo.append(np.linalg.norm(F(l * s + l ** 2 * t) - (R0 - l * J @ s)))  # third order: x8 per halving
    assert plain[0] / plain[1] == pytest.approx(4.0, rel=0.05)
    assert geo[0] / geo[1] == pytest.approx(8.0, rel=0.1) and geo[1] < 0.2 * plain[1]


def test_floor_blocks_skip_depths_through_an_inactive_mask():
    rng = np.random.default_rng(5)
    nd, nt = 4, 5
    U, _, Vt = np.linalg.svd(rng.normal(size=(nt, nt)))
    di = np.stack([(U * np.array([1.0, 0.5, 0.2, 0.1, 1e-8])) @ Vt for _ in range(nd)])
    active = np.ones((nt, nd), dtype=bool)
    skip = cnk_trap.window_mask(active, (1, 2))
    out, n = cnk_trap.floor_blocks(di, 1e-2, active & ~skip)
    assert n == 2
    np.testing.assert_array_equal(out[1:3], di[1:3])
    assert not np.array_equal(out[0], di[0]) and not np.array_equal(out[3], di[3])


def test_coarse_correction_is_exact_on_the_slow_modes():
    lo, di, up = banded_system(nt=6, nd=10, seed=9)
    nd, nt = di.shape[0], di.shape[1]
    rng = np.random.default_rng(2)
    for d in (3, 4, 6):                                  # nearly singular local blocks (a fast/slow split)
        U, _, Vt = np.linalg.svd(rng.normal(size=(nt, nt)))
        di[d] = (U * np.array([3.0, 2.0, 1.5, 1.0, 0.7, 1e-7])) @ Vt
    M = BlockTridiag(lo, di, up).factor()
    lo2 = lo.copy(); lo2[4] += 0.3 * rng.normal(size=(nt, nt))   # coupling the preconditioner lacks
    Jx = BlockTridiag(lo2, di, up)
    active = np.ones((nt, nd), dtype=bool)
    Y, Z, info = cnk_trap.slow_modes(di, active, (2, 7), 1e-4)
    assert sorted(d for d, _ in info) == [3, 4, 6] and all(s < 1e-6 for _, s in info)
    Y2, Z2, info2 = cnk_trap.slow_modes(di, active, (2, 7), 1e-4, exclude=[5])
    assert len(Z2) <= 3 and all(np.all(z[5] == 0) for z in Z2)
    JZ = [Jx.matvec(z) for z in Z]
    W = rng.uniform(0.5, 1.0, size=(nt, nd))
    Pls = cnk_trap.CoarseCorrected(M, Y, Z, JZ, W=W, mode='ls')
    for z in Z:                                          # exact on the coarse space
        np.testing.assert_allclose(Pls.solve(Jx.matvec(z)), z, rtol=0, atol=1e-8)
    r = rng.normal(size=(nt, nd))
    c = Pls.coefficients(r)                              # least-squares optimal coarse coefficients
    fine = lambda c_: np.linalg.norm(W * M.solve(r - sum(ci * j for ci, j in zip(c_, JZ))))
    assert fine(c) <= fine(np.zeros(len(Z))) and all(fine(c) <= fine(c + 1e-3 * e) for e in np.eye(len(Z)))
    Pg = cnk_trap.CoarseCorrected(M, Y, Z, JZ, mode='galerkin')
    for z in Z:
        np.testing.assert_allclose(Pg.solve(Jx.matvec(z)), z, rtol=0, atol=1e-6)
    assert all(abs(float(np.sum(y * (r - (Pg.JZ.T @ Pg.coefficients(r)).reshape(r.shape))))) < 1e-8 for y in Y)


def test_population_ptc_shifts_jacobian_and_preconditioner_diagonals_except_constraint_rows():
    cnk = _driver()
    nt, nd = 134, 3
    rng = np.random.default_rng(8)
    lo = 0.1 * rng.normal(size=(nd, nt, nt)); up = 0.1 * rng.normal(size=(nd, nt, nt)); lo[0] = 0; up[-1] = 0
    di = rng.normal(size=(nd, nt, nt)) + 5.0 * np.eye(nt) * np.where(rng.random(nt) < 0.5, 1.0, -1.0)[None, :, None]
    it = cnk.Iterate.__new__(cnk.Iterate)
    it.Jdir = cnk.BlockTridiag(lo.copy(), di.copy(), up.copy())
    it._pc = (lo.copy(), di.copy(), up.copy())
    it.M = cnk.BlockTridiag(*[m.copy() for m in it._pc]).factor()
    it.active = np.ones((nt, nd), dtype=bool); it.active[5, 1] = False
    rs = (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T
    it.D = 1.0 / np.maximum(rs, 1e3 * rs.max(axis=0))   # an --eq-floor D must not change the PTC scale
    it.rs_gross = rs
    z = rng.normal(size=(nt, nd))
    jz0 = it.Jdir.matvec(z)
    it.add_pop_ptc(0.1)
    shift = it.Jdir.matvec(z) - jz0
    expect = 0.1 * np.sign(np.einsum('kii->ik', di)) * rs * z * it.active
    expect[list(cnk.PTC_CONSTRAINT_ROWS)] = 0.0
    np.testing.assert_allclose(shift, expect, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(it.pop_ptc * z, expect, rtol=1e-12, atol=1e-12)
    # the same shift in the preconditioner; off-diagonal blocks untouched
    np.testing.assert_allclose(np.einsum('kii->ik', it._pc[1]) - np.einsum('kii->ik', di), it.pop_ptc, atol=1e-12)
    np.testing.assert_array_equal(it._pc[0], lo); np.testing.assert_array_equal(it.Jdir.up, up)
    # the shift adds to the magnitude of every shifted diagonal (pseudo-time damping)
    dd = np.einsum('kii->ik', di); m = it.pop_ptc != 0
    assert np.all(np.abs(dd + it.pop_ptc)[m] > np.abs(dd)[m])


def test_trace_ion_balance_zeroes_the_linearized_replaced_row_only_where_the_ion_is_trace():
    rng = np.random.default_rng(12)
    nt, nd = 134, 6
    Jdi = rng.normal(size=(nd, nt, nt)) * 1e-3
    Jdi[:, 101, 131] = 0.3; Jdi[:, 101, 132] = 0.3            # the He II / He III balance: He III and n_e columns
    Jdi[4, 101, 131] = 1e-6                                    # depth 5: an ordinary ground-state row
    frac = np.full((nt, nd), 1e-3); frac[131, :3] = 1e-40      # He III trace at depths 1-3 only ... 
    frac[131, 4] = 1e-40                                       # ... and at depth 5 (but not that row's variable)
    R = rng.normal(size=(nt, nd)); Js = rng.normal(size=(nt, nd))
    for lam in (1.0, 0.25):
        mask, sc = cnk_trap.trace_ion_balance(R, Js, Jdi, frac, lam)
        assert list(np.nonzero(mask)[0]) == [0, 1, 2]
        lin = R[101] - lam * Js[101] - Jdi[:, 101, 131] * sc
        np.testing.assert_allclose(lin[mask], 0.0, atol=1e-12)
        assert np.all(sc[~mask] == 0.0)



def test_adapt_caps_shrinks_only_limit_cycling_components():
    default = 0.693
    caps = np.full((3, 2), default)
    prev = np.array([[0.693, 0.1], [-0.693, 0.693], [0.0, -0.3]])
    new = np.array([[-0.693, 0.1], [0.693, 0.693], [0.2, 0.3]])     # flip at the cap, no flip, flip below the cap
    c, n = cnk_trap.adapt_caps(caps, prev, new, default)
    assert n == 2 and c[0, 0] == default / 2 and c[1, 0] == default / 2
    assert np.all(c[[0, 1, 2, 2], [1, 1, 0, 1]] == default)
    c2, n2 = cnk_trap.adapt_caps(np.full((3, 2), 0.1), new, new, default)                     # no flips: limits grow back
    assert n2 == 0 and np.allclose(c2, 0.15)
    x = np.array([1.0, 1.0]); f = cnk_trap.update_factor(np.array([-5.0, -5.0]), 'linear', np.array([0.693, 0.1]))
    np.testing.assert_allclose(f, [np.exp(0.693), np.exp(0.1)])                                       # elementwise limits
    assert x is not None


def test_population_ptc_pc_scale_uses_the_preconditioner_row_scale():
    cnk = _driver()
    nt, nd = 134, 2
    rng = np.random.default_rng(6)
    di = rng.normal(size=(nd, nt, nt)) * 0.01 + np.eye(nt)[None]
    zero = np.zeros_like(di)
    it = cnk.Iterate.__new__(cnk.Iterate)
    it.Jdir = cnk.BlockTridiag(zero.copy(), di.copy(), zero.copy())
    pc = di.copy(); pc[:, 2, :] *= 1e-3                      # a trapped level: the local ALI cancels most of its row
    it._pc = (zero.copy(), pc, zero.copy())
    it.active = np.ones((nt, nd), dtype=bool)
    it.rs_gross = np.abs(di).sum(axis=2).T
    it.add_pop_ptc(0.1, scale='pc')
    np.testing.assert_allclose(np.abs(it.pop_ptc[2]), 0.1 * np.abs(pc[:, 2, :]).sum(axis=1), rtol=1e-12)
    np.testing.assert_allclose(np.abs(it.pop_ptc[5]), 0.1 * it.rs_gross[5], rtol=1e-12)


def test_preconditioner_only_ptc_leaves_the_jacobian_exact():
    cnk = _driver()
    nt, nd = 134, 2
    rng = np.random.default_rng(10)
    di = rng.normal(size=(nd, nt, nt)) * 0.01 + np.eye(nt)[None]
    zero = np.zeros_like(di)
    it = cnk.Iterate.__new__(cnk.Iterate)
    it.Jdir = cnk.BlockTridiag(zero.copy(), di.copy(), zero.copy())
    it._pc = (zero.copy(), di.copy(), zero.copy())
    it.active = np.ones((nt, nd), dtype=bool); it.rs_gross = np.abs(di).sum(axis=2).T; it.pop_ptc = None
    it.add_pop_ptc(0.1, scale='pc', jacobian=False)
    np.testing.assert_array_equal(it.Jdir.di, di)
    assert it.pop_ptc is None and not np.array_equal(it._pc[1], di)


def test_row_scaled_measure_weighs_the_constraint_rows_whatever_their_variables_weight():
    cnk = _driver()
    nt, nd = 134, 3
    it = cnk.Iterate.__new__(cnk.Iterate)
    it.W = np.ones((nt, nd)); it.W[131] = 0.0; it.W[50] = 1e-3          # He III variable ultra-trace, a trace level
    it.active = np.ones((nt, nd), dtype=bool); it.active[131, 2] = False
    w = it.row_weights()
    assert np.all(w[131, :2] == 1.0) and w[131, 2] == 0.0 and np.all(w[50] == 1e-3) and np.all(it.W[131] == 0.0)


def test_bounded_step_takes_a_separate_radius_for_temperature():
    cnk = _driver()
    rng = np.random.default_rng(1)
    S = [rng.normal(size=(134, 3)) for _ in range(6)]
    A = rng.normal(size=(50, 6)); b = rng.normal(size=50)
    mask = np.ones((134, 3), dtype=bool)
    rad = np.full((134, 3), 0.5); rad[133] = 0.03
    s, mu, y = cnk.lm_step_bound((S, A, b), np.ones((134, 3)), rad, mask)
    assert mu > 0 and np.max(np.abs(s[133])) <= 0.03 * (1 + 1e-6) and np.max(np.abs(s[:133])) <= 0.5
    s2, mu2, _ = cnk.lm_step_bound((S, A, b), np.ones((134, 3)), 0.5, mask)       # scalar radius as before
    assert np.max(np.abs(s2)) <= 0.5 * (1 + 1e-6)


# free_r25 iterate 16 (history.json): every trial, fallbacks included, worse than the base 2.837128e-2
IT16_TRIALS = {'l1': 0.028475305384235518, 'l0.5': 0.028408308527980295, 'l0.25': 0.02838631698546318,
               'l0.125': 0.02837796207612615, 'l1_d16': 0.02837441688921275, 'l0.5_d16': 0.028372798946133797,
               'l0.25_d16': 0.028372028260444025, 'l0.125_d16': 0.028371652409520692}
IT16_PHI = 0.02837128241917269


def test_trial_decision_keeps_the_best_trial_apart_from_an_acceptance_when_all_trials_are_rejected():
    trials = [(lab, Path(f'it17_trial_{lab}'), ph) for lab, ph in IT16_TRIALS.items()]
    dec = cnk_trap.trial_decision(trials, IT16_PHI)
    assert dec['best_trial'] == 'l0.125_d16' and dec['accepted'] is None and dec['status'] == 'rejected'
    assert dec['best_merit'] == IT16_TRIALS['l0.125_d16']
    # a failed native pass (merit inf) is never accepted, a strict decrease is
    assert cnk_trap.trial_decision([('l1', None, float('inf'))], 1.0)['accepted'] is None
    assert cnk_trap.trial_decision([('l1', None, 0.5), ('l0.5', None, 0.7)], 1.0)['accepted'] == 'l1'
    assert cnk_trap.trial_decision([('l1', None, 1.0)], 1.0)['accepted'] is None


def test_restart_check_refuses_a_rejected_trial_and_admits_the_accepted_one(tmp_path):
    import json
    # legacy records (no status): iterate 15 accepted l0.5 (2.83722e-2 < 2.83790e-2); iterate 16 named its rejected best
    # trial `accepted`
    hist = [dict(it=15, base=str(tmp_path / 'it15_trial_l1'), merit=0.028379043416503255,
                 trials={'l0.5': 0.028372153681630485}, accepted='l0.5', merit_new=0.028372153681630485),
            dict(it=16, base=str(tmp_path / 'it16_trial_l0.5'), merit=IT16_PHI, trials=IT16_TRIALS,
                 accepted='l0.125_d16', merit_new=IT16_TRIALS['l0.125_d16'])]
    (tmp_path / 'history.json').write_text(json.dumps(hist))
    for n in ('it16_trial_l0.5', 'it17_trial_l0.125_d16', 'it17_trial_l1'):
        (tmp_path / n).mkdir()
    assert cnk_trap.history_accepted(hist[0]) == 'l0.5' and cnk_trap.history_accepted(hist[1]) is None
    assert cnk_trap.check_restart_base(tmp_path / 'it16_trial_l0.5') is None
    assert 'did not accept' in cnk_trap.check_restart_base(tmp_path / 'it17_trial_l0.125_d16')
    assert 'did not accept' in cnk_trap.check_restart_base(tmp_path / 'it17_trial_l1')
    # a record with the status field is taken as written
    hist[1].update(status='stopped_no_decrease', accepted=None, best_trial='l0.125_d16')
    (tmp_path / 'history.json').write_text(json.dumps(hist))
    assert 'did not accept' in cnk_trap.check_restart_base(tmp_path / 'it17_trial_l0.125_d16')
    # a base outside a numbered trial (it00_base, another run's pass) is not judged
    assert cnk_trap.check_restart_base(tmp_path / 'it00_base') is None
    cnk_trap.write_accepted_pointer(tmp_path, 16, tmp_path / 'it16_trial_l0.5', IT16_PHI, 'test')
    ptr = cnk_trap.read_accepted_pointer(tmp_path)
    assert ptr['it'] == 16 and Path(ptr['base']).name == 'it16_trial_l0.5' and ptr['merit'] == IT16_PHI
    assert not (tmp_path / 'accepted_state.tmp').exists()


def test_driver_refuses_to_start_from_the_rejected_iterate_16_trial(tmp_path):
    # the driver's start check runs before any native pass: an all-rejected iterate's best trial is never a base
    import json, subprocess
    hist = [dict(it=16, base=str(tmp_path / 'it16_trial_l0.5'), merit=IT16_PHI, trials=IT16_TRIALS,
                 accepted='l0.125_d16', merit_new=IT16_TRIALS['l0.125_d16'])]
    (tmp_path / 'history.json').write_text(json.dumps(hist))
    (tmp_path / 'it17_trial_l0.125_d16').mkdir()
    x0 = tmp_path / 'x0'; x0.mkdir()
    drv = Path(__file__).resolve().parents[1] / 'scripts' / 'cool_newton_krylov.py'
    # no executable exists: a refusal before the first native pass is the only way to this message
    r = subprocess.run([sys.executable, str(drv), '--out', str(tmp_path), '--x0', str(x0), '--exe', '/nonexistent',
                        '--replace', '/nonexistent', '--start-base', str(tmp_path / 'it17_trial_l0.125_d16'),
                        '--start-k', '17'], capture_output=True, text=True, timeout=120)
    assert r.returncode != 0 and 'refusing to start from a rejected trial' in r.stderr
    assert not (tmp_path / 'it00_base').exists() and not (tmp_path / 'accepted_state.json').exists()


IT16_DIR = Path(__file__).resolve().parents[1] / 'results/nk167_trap/free_r25/it16_from_it16_trial_l0.5'


@pytest.mark.skipif(not (IT16_DIR / 'col_ne_p').exists(), reason='archived free_r25 iterate 16 not present')
def test_struct_active_restores_the_charge_row_that_the_conservation_projection_zeroed():
    # at the projected state the charge residual is exactly 0.0 at 150 of 167 depths: R != 0 dropped the row and n_e
    import campaign_norms as cn
    cnk = _driver()
    x0 = cn.configure()
    try:
        it = cn.iterate_at(IT16_DIR, x0)
        assert int((~it.active[132]).sum()) == 150 and int((it.R[132] == 0.0).sum()) == 150
        cnk.STRUCT_ACTIVE = True
        it2 = cn.iterate_at(IT16_DIR, x0)
        assert it2.active[[31, 131, 132]].all()
        others = np.ones(it.active.shape[0], dtype=bool); others[[31, 131, 132]] = False
        np.testing.assert_array_equal(it2.active[others], it.active[others])
        np.testing.assert_array_equal(it2.R, it.R)
    finally:
        cnk.STRUCT_ACTIVE = False


def test_project_linear_is_the_first_order_projection_and_is_idempotent():
    rng = np.random.default_rng(7)
    nd = 4
    x = np.exp(rng.normal(size=(134, nd)))
    x, _ = cnk_trap.project_conservation(x, (x[cnk_trap.H_ROWS].sum(0), x[cnk_trap.HE_ROWS].sum(0)))
    tot = cnk_trap.conservation_totals(x)
    z = 0.3 * rng.normal(size=x.shape)
    zp = cnk_trap.project_linear(x, z)
    np.testing.assert_allclose(cnk_trap.project_linear(x, zp), zp, rtol=0, atol=1e-13)
    for eps in (1e-3, 1e-4):
        xp, _ = cnk_trap.project_conservation(x * (1.0 - eps * z), tot)
        # the projected trial and the linearly projected step agree to second order in eps
        assert np.max(np.abs((1.0 - xp / x) - eps * zp)) < 50 * eps ** 2
    # a projected step leaves no first-order work for the projection
    xq, st = cnk_trap.project_conservation(x * (1.0 - 1e-4 * zp), tot)
    assert max(st.values()) < 1e-7


def test_bounded_step_qr_matches_the_normal_equations_when_well_conditioned_and_box_fits_best():
    cnk = _driver()
    rng = np.random.default_rng(11)
    shape = (134, 3)
    S = [rng.normal(size=shape) for _ in range(6)]
    A = rng.normal(size=(50, 6)); b = 5.0 * rng.normal(size=50)
    active = np.ones(shape, dtype=bool); mask = active.copy(); mask[133] = False
    s_ne, mu_ne, y_ne = cnk.lm_step_bound((S, A, b), active.astype(float), 0.5, mask)
    s_qr, mu_qr, y_qr = cnk.bounded_step((S, A, b), active, 0.5, mask, 'qr')
    assert mu_ne > 0 and np.isclose(mu_qr, mu_ne, rtol=1e-6)
    np.testing.assert_allclose(y_qr, y_ne, rtol=1e-6, atol=1e-10)
    pytest.importorskip('clarabel')
    s_bx, _, y_bx = cnk.bounded_step((S, A, b), active, 0.5, mask, 'box')
    assert np.max(np.abs(s_bx[mask])) <= 0.5 * (1 + 1e-7)
    assert np.linalg.norm(A @ y_bx - b) <= np.linalg.norm(A @ y_qr - b) * (1 + 1e-7)


def test_span_basis_keeps_coefficients_bounded_on_a_collinear_span_with_noisy_products():
    # a direction repeated with a 1e-14 change and products carrying 1e-10 relative noise: the raw least squares cancels
    # huge coefficients into a fit made of noise; the truncated orthonormal basis does not
    cnk = _driver()
    rng = np.random.default_rng(5)
    shape = (134, 2)
    J = rng.normal(size=(268, 268))
    base_dirs = [rng.normal(size=shape) for _ in range(4)]
    S = base_dirs + [base_dirs[0] + 1e-14 * rng.normal(size=shape)]
    JS = [(J @ z.ravel()) * (1 + 1e-10 * rng.normal(size=268)) for z in S]
    A = np.stack(JS, 1); b = rng.normal(size=268)
    y_raw = cnk.span_lstsq(S, A, b, 0.0)
    y_tr = cnk.span_lstsq(S, A, b, 1e-8)
    nz = np.array([np.linalg.norm(z) for z in S])
    s_tr = sum(yi * zi for yi, zi in zip(y_tr, S))
    assert np.linalg.norm(y_raw * nz) > 1e6 * np.linalg.norm(sum(yi * zi for yi, zi in zip(y_raw, S)))
    assert np.linalg.norm(y_tr * nz) < 10 * np.linalg.norm(s_tr)
    # the truncated model predicts the true product of its step
    true = J @ s_tr.ravel()
    np.testing.assert_allclose(A @ y_tr, true, rtol=0, atol=1e-6 * np.linalg.norm(true))
    Sb, Ab, M, k = cnk.span_basis(S, A, 1e-8)
    assert k == 4 and np.allclose(np.stack([u.ravel() for u in Sb], 1).T @ np.stack([u.ravel() for u in Sb], 1), np.eye(4))


def test_orthonormalize_batch_drops_dependent_directions_and_keeps_conservation_consistency():
    rng = np.random.default_rng(9)
    nd = 3
    x = np.exp(rng.normal(size=(134, nd)))
    x, _ = cnk_trap.project_conservation(x, (x[cnk_trap.H_ROWS].sum(0), x[cnk_trap.HE_ROWS].sum(0)))
    S = [cnk_trap.project_linear(x, rng.normal(size=x.shape)) for _ in range(4)]
    new = [cnk_trap.project_linear(x, rng.normal(size=x.shape)) for _ in range(3)]
    new.append(2.0 * S[1] - S[3])                              # already in the span: dropped
    new.append(new[0] + 1e-3 * cnk_trap.project_linear(x, rng.normal(size=x.shape)))   # nearly new[0], independent: kept
    out = cnk_trap.orthonormalize_batch(S, new)
    assert len(out) == 4
    Q = np.stack([u.ravel() for u in out], 1)
    np.testing.assert_allclose(Q.T @ Q, np.eye(4), atol=1e-12)
    Sm = np.stack([z.ravel() for z in S], 1)
    assert np.abs(Q.T @ Sm).max() < 1e-10
    for u in out:                                              # still conservation-consistent
        np.testing.assert_allclose(cnk_trap.project_linear(x, u), u, atol=1e-12)
    assert cnk_trap.orthonormalize_batch(S, []) == []
