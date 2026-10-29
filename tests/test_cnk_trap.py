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
