"""Gate-4 frozen-source hydro step: analytic limits, parsers, guard and policy."""
import json
import os
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad, simpson, cumulative_simpson
from scipy.optimize import brentq
from qstars.constants import G, C, A, MSUN
from qstars.continuum_wind import flux_driven_reference
from qstars.atmosphere_reference import AtmosphereDomainExit
import qstars.cmfgen_wind as cw
from qstars.cmfgen_wind import (FrozenSources, hydro_step, tau_seed, read_hydro_terms, read_gencool_channels,
    momentum_residual, mach_one_guard, guard_status, ConvectiveFluxExit, classify, outer_update,
    wind_ledger, grad_input_audit, reduction_screens, rosseland_depth, critical_coefficients)

M = 1e6*MSUN; GM = G*M; A2 = 1e12; KAPPA = .34


def constant(a2=A2, gamma=.5, kappa=KAPPA):
    def f(x):
        z = np.asarray(x, float); return np.full_like(z, a2), 0*z, np.full_like(z, gamma), np.full_like(z, kappa)
    return f


def parker_case(gamma=.5):
    rc = GM*(1-gamma)/(2*A2); rb, ro = rc/4, 10*rc
    return rc, rb, ro, FrozenSources('radius', constant(gamma=gamma), (np.log(rb)-1, np.log(ro)+1))


def parker_velocity(r, rc):
    """Independent implicit isothermal Parker root u^2-ln u^2=4 ln xi+4/xi-3."""
    xi = r/rc; rhs = 4*np.log(xi)+4/xi-3
    f = lambda lu: np.exp(2*lu)-2*lu-rhs
    return np.sqrt(A2)*np.exp(brentq(f, -60., 0.) if xi < 1 else brentq(f, 0., 10.))


def parker_ln_u(xi):
    """ln(v/a) of the same implicit root, valid deep in the subsonic layer (ln u ~ -lambda)."""
    rhs = 4*np.log(xi)+4/xi-3; f = lambda lu: np.exp(2*lu)-2*lu-rhs
    return brentq(f, -rhs/2-1, 0., xtol=1e-14, rtol=1e-15) if xi < 1 else brentq(f, 0., 10., xtol=1e-15)


def passed_guard(depths=3):
    return dict(status='passed', guarded_depths=depths, opacities={'cmfgen': {}, 'mesa': {}})


def test_isothermal_parker_constant_gamma_and_independent_mdot():
    rc, rb, ro, src = parker_case()
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    assert s['status'] == 'transonic' and s['all_branches_enumerated'] and len(s['critical_points']) == 1
    assert abs(s['sonic_radius']/rc-1) < 1e-12 and abs(s['sonic_slope']-1) < 1e-8
    assert s['regularity_relative'] < 1e-6
    q = critical_coefficients(src, M, np.log(rc), .5*np.log(A2))
    assert abs(q['numerator'])/(GM/rc) < 1e-12 and abs(q['denominator'])/A2 < 1e-14 and q['dlntau_dlnr'] is None
    p = s['profile']; u = p['velocity']/np.sqrt(A2); xi = p['radius']/rc
    np.testing.assert_allclose(u*u-np.log(u*u), 4*np.log(xi)+4/xi-3, atol=5e-8)
    # tau_b-tau_o = mdot int kappa/(4 pi r^2 v) dr with v from the independent implicit root
    x = np.linspace(np.log(rb), np.log(ro), 4001); v = np.array([parker_velocity(np.exp(z), rc) for z in x])
    mdot = (10.-1e-3)/simpson(KAPPA/(4*np.pi*np.exp(x)*v), x=x)
    assert abs(s['mdot']/mdot-1) < 1e-6
    d = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_density=1e-12)
    assert abs(d['mdot']/(4*np.pi*rb*rb*1e-12*parker_velocity(rb, rc))-1) < 1e-7
    assert not s['wind_certified'] and not s['physical_impossibility_claimed']


def test_tau_mode_equals_radius_mode_for_constant_sources_and_rematches():
    rc, rb, ro, src = parker_case()
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    tsrc = FrozenSources('tau', constant(), (np.log(1e-5), np.log(1e3)))
    q = hydro_step(tsrc, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.,
                   seed=tau_seed(s)+[.01, -.02])
    assert q['status'] == 'transonic' and q['method'] == 'admissible_damped_newton'
    assert abs(q['ln_mdot']-s['ln_mdot']) < 1e-7 and abs(q['sonic_radius']/rc-1) < 1e-12
    assert abs(q['ln_sonic_tau_excess']-s['ln_sonic_tau_excess']) < 1e-7
    assert abs(q['base_log_residual']) < 1e-8 and abs(q['outer_log_residual']) < 1e-8 and abs(q['outer_excess_residual']) < 1e-8
    assert q['trials'][0]['purpose'] == 'initial' and not q['all_branches_enumerated']
    # Requested radii are never clipped onto the branch ends (radius mode rejects the same input).
    for bad in ([rb/2, rb, ro], [rb, ro, 2*ro]):
        for kw in [dict(sources=src), dict(sources=tsrc, seed=tau_seed(s))]:
            with pytest.raises(ValueError, match='outside the hydro domain'):
                hydro_step(mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10., radii=bad, **kw)
    # Base moved to another subsonic depth with its own tau (or rho): no mdot input.
    r1 = 1.7*rb
    s1 = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10., radii=[rb, r1, ro])
    tau1, rho1 = s1['profile']['tau'][1], np.exp(s1['profile']['log_density'][1])
    for kw in [dict(base_tau=tau1), dict(base_density=rho1)]:
        moved = hydro_step(src, mass=M, base_radius=r1, outer_radius=ro, outer_tau=1e-3, **kw)
        assert abs(moved['ln_mdot']-s['ln_mdot']) < 1e-5
    moved = hydro_step(tsrc, mass=M, base_radius=r1, outer_radius=ro, outer_tau=1e-3, base_tau=tau1, seed=tau_seed(s))
    assert abs(moved['ln_mdot']-s['ln_mdot']) < 1e-5
    # tau-mode requires an explicit seed; nothing is guessed silently.
    with pytest.raises(ValueError, match='seed'):
        hydro_step(tsrc, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)


@pytest.mark.parametrize('fac, rtol', [(40, 1e-10), (7000, 1e-12)])
def test_tau_mode_thin_supersonic_layer_and_underflowing_mdot(fac, rtol):
    """tau_c-tau_o << tau_o (1e-30 at r_c/40) and ln mdot ~ -1.4e4 (lambda(R_b) = 1.4e4 at r_c/7000).

    Both modes must agree with each other and with an independent quadrature of
    the implicit Parker root, tau_b-tau_o = mdot kappa/(4 pi a) int dr/(r^2 u).
    rtol=1e-12 at r_c/7000: at rtol 1e-10 the radius mode's error control on
    |ln v| ~ 1.4e4 allows ~1e-6 in ln mdot, still far inside every plan threshold.
    """
    rc = GM*.5/(2*A2); rb, ro = rc/fac, 10*rc
    s = hydro_step(FrozenSources('radius', constant(), (np.log(rb)-1, np.log(ro)+1)), mass=M, base_radius=rb,
                   outer_radius=ro, outer_tau=1e-3, base_tau=10., rtol=rtol)
    tsrc = FrozenSources('tau', constant(), (np.log(1e-6), np.log(1e3)))
    q = hydro_step(tsrc, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10., seed=tau_seed(s), rtol=rtol)
    assert s['status'] == q['status'] == 'transonic'
    assert abs(q['ln_mdot']-s['ln_mdot']) < 1e-7 and abs(q['sonic_radius']/rc-1) < 1e-12
    assert abs(q['ln_sonic_tau_excess']-s['ln_sonic_tau_excess']) < 1e-7 and q['ln_sonic_tau_excess']-np.log(1e-3) < -60
    assert max(abs(q[k]) for k in ('base_log_residual', 'outer_log_residual', 'outer_excess_residual')) < 1e-8
    lub = parker_ln_u(1/fac); g = lambda r: np.exp(lub-parker_ln_u(r/rc))*(rb/r)**2
    J = quad(g, rb, rb*(1+50/fac), epsabs=0, epsrel=1e-13, limit=400)[0]  # base-dominated: width rb^2/(2 rc)
    if fac < 1000: J += quad(g, rb*(1+50/fac), ro, epsabs=0, epsrel=1e-13, limit=400, points=[rc])[0]
    ref = np.log((10-1e-3)*4*np.pi*np.sqrt(A2)*rb*rb/(KAPPA*J))+lub
    assert abs(s['ln_mdot']-ref) < 1e-7 and abs(q['ln_mdot']-ref) < 1e-7
    if fac == 7000:
        assert q['ln_mdot'] < -1.39e4 and q['mdot'] == 0. and np.all(np.isfinite(q['profile']['log_density']))
        c = classify(q, s, mass=M, luminosity=1e43, base_radius=rb, input_checks={'nonrelativistic': True},
                     guard=passed_guard(), screens=dict(all_passed=True, passed={}))
        assert c['rule'] == 10 and c['tiring'] == 0.


def test_variable_sound_speed_radius_mode_independent_momentum_integral():
    q, gamma = .5, .2; a0 = A2; r0 = GM/(20*a0)
    def f(x):
        z = np.asarray(x, float); a2 = a0*np.exp(-q*(z-np.log(r0)))
        return a2, -q*a2, np.full_like(z, gamma), np.full_like(z, KAPPA)
    rc = (GM*(1-gamma)/((2+q)*a0*r0**q))**(1/(1-q))
    rb, ro = rc/3, 20*rc
    src = FrozenSources('radius', f, (np.log(rb), np.log(ro)))
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.,
                   radii=np.geomspace(rb, ro, 20001))
    assert s['status'] == 'transonic' and abs(s['sonic_radius']/rc-1) < 1e-12
    p = s['profile']; r, v, lr = p['radius'], p['velocity'], p['log_density']; a2 = f(np.log(r))[0]
    # int rho^-1 dP = Delta a^2 + Delta(a^2 ln rho) - int ln rho da^2; mdot cancels.
    work = (a2[-1]-a2[0])+(a2[-1]*lr[-1]-a2[0]*lr[0])-simpson(lr*(-q*a2/r), x=r)
    defect = .5*(v[-1]**2-v[0]**2)+work+GM*(1-gamma)*(1/r[0]-1/r[-1])
    assert abs(defect)/(GM/rb) < 1e-6


def test_variable_sound_speed_tau_mode_integral_and_optical_depth():
    gamma = .5; rc0 = GM*(1-gamma)/(2*A2); rb, ro = rc0/4, 10*rc0
    def ft(lt):
        t = np.exp(np.asarray(lt, float)); s = t/(1+t)
        return A2*(1+.3*s), A2*.3*s/(1+t), np.full_like(t, gamma), np.full_like(t, KAPPA)
    seed = hydro_step(parker_case()[3], mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    q = hydro_step(FrozenSources('tau', ft, (np.log(1e-5), np.log(1e3))), mass=M, base_radius=rb, outer_radius=ro,
                   outer_tau=1e-3, base_tau=10., seed=tau_seed(seed), radii=np.geomspace(rb, ro, 20001))
    assert q['status'] == 'transonic' and abs(q['base_log_residual']) < 1e-8
    p = q['profile']; r, v, lr, tau = p['radius'], p['velocity'], p['log_density'], p['tau']
    a2 = ft(np.log(tau))[0]
    work = (a2[-1]-a2[0])+(a2[-1]*lr[-1]-a2[0]*lr[0])-np.sum(.5*(lr[1:]+lr[:-1])*np.diff(a2))
    assert abs(.5*(v[-1]**2-v[0]**2)+work+GM*(1-gamma)*(1/r[0]-1/r[-1]))/(GM/rb) < 1e-6
    c = cumulative_simpson(KAPPA*np.exp(lr)*r, x=np.log(r), initial=0)
    np.testing.assert_allclose((tau-1e-3)[:-10], (c[-1]-c)[:-10], rtol=1e-6)


def test_gamma_crossing_one_critical_point_slope_and_cold_flux_driven_limit():
    rb = GM/(1000*A2); xt = np.log(3*rb); w = .2; g0, g1 = .5, 3.
    gam = lambda x: g0+(g1-g0)/(1+np.exp(-(x-xt)/w))
    dgam = lambda x: (g1-g0)*np.exp(-(x-xt)/w)/w/(1+np.exp(-(x-xt)/w))**2
    f = lambda x: (lambda z: (np.full_like(z, A2), 0*z, gam(z), np.full_like(z, KAPPA)))(np.asarray(x, float))
    ro = np.exp(xt+24*w); src = FrozenSources('radius', f, (np.log(rb), np.log(ro)))
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.,
                   radii=np.geomspace(rb, ro, 3001))
    assert s['status'] == 'transonic' and len(s['critical_points']) == 1
    rc = s['sonic_radius']; xc = np.log(rc)
    assert abs(gam(xc)-(1-2*A2*rc/GM)) < 1e-12 and gam(xc) < 1
    assert abs(s['sonic_slope']/np.sqrt(1+GM/(rc*A2)*dgam(xc)/2)-1) < 1e-7
    p = s['profile']; v, r = p['velocity'], p['radius']
    B = [v[i]**2/2-A2*np.log(v[i])-2*A2*np.log(r[i])-GM/r[i]
         -quad(lambda q: GM*gam(np.log(q))/q**2, rb, r[i], epsabs=0, epsrel=1e-13, limit=400)[0] for i in range(0, 3001, 150)]
    assert np.ptp(B)/(GM/rb) < 1e-6
    # Beyond the transition Gamma=g1: Owocki et al. eq. 24 (tiring -> 0) plus the explicit isothermal work.
    m = r >= np.exp(xt+18*w); i2 = np.flatnonzero(m)[0]; r2, v2 = r[i2], v[i2]  # Gamma=g1 to 4e-8
    wref, _ = flux_driven_reference(r[m]/r2, g1, 1e-12, v2**2/(2*GM/r2))
    pressure = (A2*np.log(v[m]/v2)+2*A2*np.log(r[m]/r2))/(GM/r2)
    assert np.max(abs(v[m]**2/(2*GM/r2)-wref-pressure)) < 1e-6
    assert np.max(abs(pressure/wref)) < .1


def test_hydrostatic_no_wind_bound_is_analytic_and_bounds_the_true_transonic_rate():
    gamma = .5; rc = GM*(1-gamma)/(2*A2); rb, ro = rc/40, rc/2
    src = FrozenSources('radius', constant(gamma=gamma), (np.log(rb), np.log(3*rc)))
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_density=1e-12)
    assert s['status'] == 'no_critical_point' and s['max_scaled_numerator'] < -1e-6
    lam = GM*(1-gamma)/A2
    exact = np.log(4*np.pi*ro*ro*1e-12*A2/np.sqrt(A2))-lam*(1/rb-1/ro)
    assert abs(s['ln_mdot_bound']-exact) < 1e-8
    true = hydro_step(src, mass=M, base_radius=rb, outer_radius=2*rc, outer_tau=1e-3, base_density=1e-12)
    assert true['status'] == 'transonic' and true['ln_mdot'] < s['ln_mdot_bound']
    t = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    J = quad(lambda r: KAPPA*np.exp(-lam*(1/rb-1/r))/A2, rb, ro, epsabs=0, epsrel=1e-12, limit=200)[0]
    exact = np.log(4*np.pi*ro*ro*(10-1e-3)/J/np.sqrt(A2))-lam*(1/rb-1/ro)
    assert abs(t['ln_mdot_bound']-exact) < 1e-7
    c = classify(t, t, mass=M, luminosity=1e43, base_radius=rb, input_checks={'nonrelativistic': True},
                 guard=passed_guard())
    assert c['classification'] == 'hydrostatic_no_wind' and c['rule'] == 5 and c['confirmation_required']
    assert not c['physical_impossibility_claimed'] and np.isfinite(c['ln_mdot_bound'])
    # Inverted optical depths are invalid input, never a NaN bound.
    with pytest.raises(ValueError, match='exceed the outer'):
        hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=10., base_tau=1e-3)


def test_non_regular_sonic_approach_is_a_rejected_crossing_not_a_failure(monkeypatch):
    """Gamma bump: N signs -,+,- (X- then O-point). Bernoulli psi(r) with dpsi/dlnr=N
    falls below psi(r_X) before r_o, so the X-point's supersonic branch folds back to
    v=a with N<0 (DOP853 stalls there) and the O-point has complex slopes: rule 4."""
    rb = GM/(1000*A2); x1, x2, w = np.log(3*rb), np.log(6*rb), .1
    lg = lambda z, c: 1/(1+np.exp(-(z-c)/w)); gam = lambda z: .5+lg(z, x1)-lg(z, x2)
    f = lambda x: (lambda z: (np.full_like(z, A2), 0*z, gam(z), np.full_like(z, KAPPA)))(np.asarray(x, float))
    ro = 40*rb; src = FrozenSources('radius', f, (np.log(rb), np.log(ro)))
    n = lambda z: 2*A2*np.exp(z)/GM-1+gam(z)
    xX = brentq(n, np.log(rb), np.log(4.5*rb)); xO = brentq(n, np.log(4.5*rb), np.log(ro))
    assert n(np.log(rb)) < 0 and n(np.log(ro)) < 0
    assert quad(lambda z: n(z)*GM/A2*np.exp(-z), xX, np.log(ro), points=[xO], limit=200)[0] < 0
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    assert s['status'] == 'no_admissible_crossing' and len(s['critical_points']) == 2
    X, O = s['critical_points']
    assert abs(np.log(X['sonic_radius'])-xX) < 1e-10 and abs(np.log(O['sonic_radius'])-xO) < 1e-10
    assert X['status'] == O['status'] == 'rejected_crossing'
    assert 'non-regularly' in X['reason'] and 'No unique accelerating' in O['reason']
    c = classify(s, s, mass=M, luminosity=1e43, base_radius=rb, input_checks={'nonrelativistic': True}, guard=passed_guard())
    assert c['rule'] == 4 and c['classification'] == 'reduction_invalid'
    # Without the stall diagnosis the candidate is unresolved: enumeration still completes,
    # every candidate is recorded, and the step is numerical_failure (rule 7), never a crossing claim.
    monkeypatch.setattr(cw, 'SONIC_STALL', 0.)
    u = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    assert u['status'] == 'numerical_failure' and [p['status'] for p in u['critical_points']] == [
        'candidate_integration_failure', 'rejected_crossing']
    assert classify(u, u, mass=M, luminosity=1e43, base_radius=rb, input_checks={'nonrelativistic': True},
                    guard=passed_guard())['rule'] == 7


def test_base_not_accelerating_and_tables_are_never_extrapolated():
    src = FrozenSources('radius', constant(gamma=1.2), (np.log(1e15), np.log(1e17)))
    assert hydro_step(src, mass=M, base_radius=1e15, outer_radius=1e17, outer_tau=1e-3, base_tau=10.)['status'] == 'base_not_accelerating'
    x = np.linspace(0, 1, 5); t = FrozenSources.tabulated('radius', x, np.exp(x), x, np.ones(5))
    np.testing.assert_allclose(t(x)[0], np.exp(x), rtol=1e-14)
    with pytest.raises(ValueError, match='no extrapolation'): t(1.0001)
    with pytest.raises(ValueError, match='outside the frozen'):
        hydro_step(t, mass=M, base_radius=np.exp(-.1), outer_radius=np.e, outer_tau=1e-3, base_tau=1.)
    with pytest.raises(ValueError, match='exactly one'):
        hydro_step(t, mass=M, base_radius=1.1, outer_radius=2., outer_tau=1e-3)


def hydro_line(R, V, err, terms, gamma, depth, tau, vs, rev):
    fr = f'{R:12.6E}' if R > 9.99e4 else f'{R:12.6f}'
    ft = f'{tau:10.2E}' if tau < .1 else (f'{tau:10.3f}' if tau < 1 else f'{tau:10.2f}')
    return ' '+fr+f'{V:17.8E}'+err+''.join(t if isinstance(t, str) else f'{t:14.4E}' for t in terms)+f'{gamma:9.2f}{depth:7d}  '+ft+f'{vs:10.2f}{rev:12.2f}'


def test_native_hydro_fixed_width_parser(tmp_path):
    head = '       R                  V     % Error        VdVdR      dPdR/ROH         g_TOT         g_RAD        g_ELEC     Gamma  Depth         Tau      Vsound   E/g_rad%'
    rows = [hydro_line(2.877415e7, 198.0001, f'{197.35:9.2f}', [1.3762e-5, -8.485e-6, -1.2587e-3, 3.4431e-4, 3.7164e-4], .21, 1, 4.08e-2, 12.64, 367.09),
            hydro_line(2.0e5, 50., '*********', [1e-4, -2e-4, -3e-3, 4e-4, '    1.0000-100'], 1234.56, 2, .5, 13.1, -1.5),
            hydro_line(2.877415e4, 1e-2, f'{-68.62:9.2f}', [6.9503e-6, -18.965, -9.2759, 6.7537, 3.7159], .42, 3, 43.29, 16.43, -143.46)]
    foot = ['  ', ' Momentum equation is: VdV/dr = - dPdR/ROH - dTPdR/ROH - g + g_RAD', '  ',
            '          Photospheric radius is:     5.8293E+05(10^10 cm)    8.3754E+04(Rsun)',
            '                   Stars mass is: ********** Msun', ' Edd. fac. (e.s.) at photosphere:  0.16706',
            ' Edd. fac. (e.s.) at inner bound:  0.23181']
    p = tmp_path/'HYDRO'; p.write_text('\n'.join([head]+rows+foot)+'\n')
    h = read_hydro_terms(p)
    np.testing.assert_allclose(h['radius'], [2.877415e14, 2.0e15, 2.877415e17])
    np.testing.assert_allclose(h['velocity'], [1e3, 5e6, 1.980001e7])
    assert np.isnan(h['percent_error'][1]) and h['overflowed_fields'] == 1
    assert h['g_elec'][1] == 1e-100 and h['gamma'][1] == 1234.56 and h['tau'][::-1].tolist() == [4.08e-2, .5, 43.29]
    np.testing.assert_allclose(h['gravity'], h['g_rad']-h['g_tot'])
    assert h['photospheric_radius'] == pytest.approx(5.8293e15, rel=1e-14) and np.isnan(h['stars_mass_Msun']) and h['gamma_es_inner'] == .23181
    assert h['native_grad_luminosity'] == 'inner_comoving_L(ND)' and not h['turbulent_pressure']
    p.write_text('\n'.join([head, rows[0][:20]+'x'+rows[0][21:]]))
    with pytest.raises(ValueError): read_hydro_terms(p)


def test_gencool_channels_close_the_native_net(tmp_path):
    def block(label, rows): return [f'   {label}']+['  '+''.join(f'{v:12.4E}' for v in r) for r in rows]+['']
    bf = np.array([[1e-10, -2e-10, 3e-10, 4e-10, -5e-10], [1e-11, 1e-11, 1e-11, 1e-11, 1e-11]])
    col = np.array([[2e-12, 3e-12, -4e-12, 5e-12, 6e-12]]); ad = np.array([[-1e-12]*5])
    net = bf.sum(0)+col.sum(0)+ad.sum(0)
    lines = []
    for sl in [slice(0, 3), slice(3, 5)]:
        lines += ['   Depth', '  '+''.join(f'{i:12d}' for i in range(sl.start+1, sl.stop+1)), '']
        lines += block('Radius [1.0E+10cm] ', [np.arange(5.)[::-1][sl]+1])
        lines += block('HI Bound-Free Cooling [ergs/cm**3/s]', bf[:, sl])
        lines += block('HI Collisional Cooling', col[:, sl])
        lines += block('Adiabatic cooling rate (V term) [ergs/cm**3/sec]', ad[:, sl])
        lines += block('Ratio of adiabatic to total cooling rate (dTdR term)[Not Incl.]', [np.ones(5)[sl]*7])
        lines += block('Net Cooling Rate [ergs/cm**3/sec]', [net[sl]])
    p = tmp_path/'GENCOOL'; p.write_text('\n'.join(lines)+'\n')
    q = read_gencool_channels(p, 5)
    assert len(q['channels']) == 3 and not any('Ratio' in k or 'Net' in k for k in q['channels'])
    np.testing.assert_allclose(q['net_cooling'], net[::-1], rtol=1e-4)
    assert np.max(q['channel_closure_relative']) < 5e-4
    np.testing.assert_allclose(q['gross_exchange_lower_bound'], .5*(abs(bf).sum(0)+abs(col).sum(0)+abs(ad).sum(0))[::-1], rtol=1e-4)
    # Adiabatic rows are flow work: kept in the closure sum, dropped from the radiative gross exchange.
    assert q['adiabatic_channels'] == ['Adiabatic cooling rate (V term) [ergs/cm**3/sec]']
    np.testing.assert_allclose(q['radiative_gross_exchange'], .5*(abs(bf).sum(0)+abs(col).sum(0))[::-1], rtol=1e-4)
    with pytest.raises(ValueError, match='depths'): read_gencool_channels(p, 6)


def test_momentum_residual_on_exact_hydrostatic_structure_and_optical_depth():
    gamma = .3; rb = GM/(200*A2); r = np.geomspace(rb, 1.2*rb, 400)
    P = 1e2*np.exp(-GM*(1-gamma)/A2*(1/rb-1/r)); rho = P/A2; g = gamma*GM/r**2
    q = momentum_residual(r, np.full_like(r, 1e-3), P, rho, g, M, sigma=np.full_like(r, -1.))
    assert np.max(abs(q['relative'])) < 1e-5 and np.allclose(q['gamma'], gamma)
    bad = momentum_residual(r, np.full_like(r, 1e-3), P, rho, 1.01*g, M, sigma=np.full_like(r, -1.))
    assert np.min(abs(bad['relative'])) > 1e-3
    np.testing.assert_allclose(bad['hydro_percent'], 200*bad['residual']/(abs(bad['dpdr_on_rho'])+abs(1.01*g-GM/r**2)))
    rr = np.geomspace(1., 10., 4001); chi = rr**-3.
    np.testing.assert_allclose(rosseland_depth(rr, chi), rr**-2/2, rtol=1e-6)
    with pytest.raises(ValueError, match='p<=1'): rosseland_depth(rr, 1/rr)


class StubEOS:
    hydrogen = .7015; radiation_constant = A
    def evaluate(self, rho, T):
        p = rho*1e12+A*T**4/3; o = np.ones_like(rho)
        return dict(p=p, chiT=o, chiRho=o, cp=1e8*o, grad_ad=.4*o, cs=np.sqrt(1.4*p/rho), free_e=o, eta=o)


class StubMesa:
    def __init__(self, mach, flux=0., opacity=1.):
        self.composition = type('X', (), {'X': .7015})(); self.radiation_constant = A
        self.mach, self.flux, self.k = np.asarray(mach), flux, opacity
    def opacity_external_eos(self, rho, T, free_e, eta): return np.full_like(rho, self.k)
    def mlt(self, x, model, alpha):
        cs = np.sqrt(1.4*x[5]/x[4]); boost = np.where(x[7] > 5, 2., 1.)
        return np.array([x[8]*(1-self.flux), self.mach*boost*cs, 0*cs, 0*cs])


def guard_inputs(n=4):
    r = np.geomspace(1e15, 2e15, n)
    return dict(state=dict(radius=r, velocity=np.full(n, 1e4), density=np.full(n, 1e-10), temperature=np.full(n, 1e4),
                           kappa_rosseland=np.full(n, .5)), sound_speed_squared=np.full(n, 1e12),
                tau_rosseland=np.geomspace(100, .1, n), moments=dict(luminosity_comoving=np.full(n, 1e43)),
                hydrogen_mass_fraction=.7015)


def test_mach_one_guard_wrapper_and_classification_order():
    inp = guard_inputs()
    ok = mach_one_guard(inp, StubEOS(), StubMesa([.5, .2, .1]), mass=M, opacity='cmfgen')
    assert ok['guarded_depths'] == 3 and ok['max_convective_mach'] == pytest.approx(.5) and ok['max_convective_flux_fraction'] == 0
    with pytest.raises(AtmosphereDomainExit) as exc:
        mach_one_guard(inp, StubEOS(), StubMesa([.5, 1.2, .1]), mass=M)
    assert exc.value.diagnostics['convective_mach'] >= 1 and not exc.value.diagnostics['physical_impossibility_claimed']
    # 'both' trips when only the MESA opacity drives Mach one.
    assert guard_status(inp, StubEOS(), StubMesa([.6]*3, opacity=10.), mass=M, opacity='cmfgen')['status'] == 'passed'
    trip = guard_status(inp, StubEOS(), StubMesa([.6]*3, opacity=10.), mass=M, opacity='both')
    assert trip['status'] == 'mach_one' and trip['diagnostics']['opacity_source'] == 'mesa'
    with pytest.raises(ConvectiveFluxExit): mach_one_guard(inp, StubEOS(), StubMesa([.5]*3, flux=.1), mass=M)
    bad = dict(inp, hydrogen_mass_fraction=.7)
    with pytest.raises(ValueError, match='composition'): mach_one_guard(bad, StubEOS(), StubMesa([.5]*3), mass=M)
    rc, rb, ro, src = parker_case()
    t = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    passed = guard_status(inp, StubEOS(), StubMesa([.5]*3), mass=M); screens = dict(all_passed=True, passed={})
    assert passed['status'] == 'passed' and set(passed['opacities']) == {'cmfgen', 'mesa'}
    kw = dict(mass=M, luminosity=1e60, base_radius=rb, input_checks={'nonrelativistic': True})
    assert classify(t, t, guard=passed, screens=screens, **kw)['classification'] == 'stationary_transonic_wind'
    # A pass counts only with both opacities (plan default while pending) and at least one guarded depth.
    one = guard_status(inp, StubEOS(), StubMesa([.5]*3), mass=M, opacity='cmfgen')
    vacuous = guard_status(dict(inp, tau_rosseland=np.full(4, .1)), StubEOS(), StubMesa([.5]*3), mass=M)
    assert one['status'] == vacuous['status'] == 'passed' and vacuous['guarded_depths'] == 0
    assert classify(t, t, guard=one, screens=screens, **kw)['rule'] == 1 and classify(t, t, guard=vacuous, screens=screens, **kw)['rule'] == 1
    assert classify(t, t, guard=trip, screens=screens, **kw)['rule'] == 2
    flux = guard_status(inp, StubEOS(), StubMesa([.5]*3, flux=.1), mass=M)
    assert classify(t, t, guard=flux, screens=screens, **kw)['rule'] == 3
    assert classify(t, t, guard=passed, screens=screens, **dict(kw, input_checks={'sync': False}))['rule'] == 1
    assert classify(t, t, guard=None, screens=screens, **kw)['rule'] == 1
    invalid = guard_status(bad, StubEOS(), StubMesa([.5]*3), mass=M)
    assert invalid['status'] == 'invalid_input' and classify(t, t, guard=invalid, screens=screens, **kw)['rule'] == 1
    assert classify(t, dict(status='multiple_admissible', reason='x'), guard=passed, screens=screens, **kw)['rule'] == 6
    assert classify(dict(status='numerical_failure', reason='x'), t, guard=passed, screens=screens, **kw)['classification'] == 'numerical_failure'
    assert classify(t, t, guard=passed, screens=screens, **dict(kw, luminosity=.5*t['mdot']*GM/rb))['classification'] == 'photon_tiring_exit'
    c = classify(t, t, guard=passed, screens=dict(all_passed=False, passed={'thermal': False}), **kw)
    assert c['rule'] == 9 and c['failed'] == ['thermal']
    assert all(not classify(t, t, guard=passed, screens=screens, **kw)[k] for k in ('wind_certified', 'physical_impossibility_claimed'))
    # Declared hydro-step tolerances are enforced (rule 7); missing values fail closed.
    tq = hydro_step(FrozenSources('tau', constant(), (np.log(1e-5), np.log(1e3))), mass=M, base_radius=rb, outer_radius=ro,
                    outer_tau=1e-3, base_tau=10., seed=tau_seed(t))
    assert classify(tq, t, guard=passed, screens=screens, **kw)['rule'] == 10
    for bad in [dict(base_log_residual=2e-8), dict(outer_log_residual=-2e-8), dict(regularity_relative=2e-6)]:
        c = classify(dict(tq, **bad), t, guard=passed, screens=screens, **kw)
        assert c['rule'] == 7 and c['unconverged'] == {'tau': list(bad)}
    c = classify({k: v for k, v in tq.items() if k != 'outer_log_residual'}, t, guard=passed, screens=screens, **kw)
    assert c['rule'] == 7 and c['unconverged'] == {'tau': ['outer_log_residual']}
    assert classify(tq, dict(t, regularity_relative=1e-5), guard=passed, screens=screens, **kw)['unconverged'] == {'radius': ['regularity_relative']}
    plan = json.loads((Path(__file__).resolve().parents[1]/cw.PLAN).read_text())['acceptance_thresholds']
    assert all(plan[k] == v for k, v in cw.THRESHOLDS.items())


def test_screens_read_the_steps_own_critical_point_and_never_clip():
    rc, rb, ro, src = parker_case()
    t = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-6, base_tau=1e-3)
    r = np.geomspace(rb, ro, 50); n = len(r); P = np.full(n, 1e-2)
    cool = dict(gross_exchange_lower_bound=np.full(n, 1e3), radiative_gross_exchange=np.full(n, 1e3),
                channel_closure_relative=np.full(n, 1e-5))
    inp = dict(state=dict(radius=r), translational_energy_density=1.5*P, cooling=cool)
    q = reduction_screens(inp, t, mass=M)
    assert not q['passed']['thermal'] and not q['all_passed'] and not q['values']['atomic_energy_included']
    q = reduction_screens(dict(inp, thermal_energy_density=1.5*P), t, mass=M)
    assert q['all_passed'] and q['values']['photon_trapping'] == pytest.approx(t['sonic_velocity']*t['sonic_tau']/C)
    vc = t['sonic_velocity']; t_flow = rc/abs(2+t['sonic_slope'])/vc
    assert q['values']['thermal_relaxation_over_flow'] == pytest.approx(1.5e-2/1e3/t_flow, rel=1e-10)
    assert q['values']['thermal_screen'] == 'gross_exchange_proxy'
    # Pending amendment: the radiative-only estimate (adiabatic flow-work rows dropped) must pass too.
    rad = dict(cool, radiative_gross_exchange=np.full(n, 1e-14))
    q = reduction_screens(dict(inp, thermal_energy_density=1.5*P, cooling=rad), t, mass=M)
    assert q['values']['thermal_relaxation_over_flow'] <= 1e-2 < q['values']['thermal_relaxation_over_flow_radiative_only']
    assert not q['passed']['thermal'] and not q['all_passed']
    noad = {k: v for k, v in cool.items() if k != 'radiative_gross_exchange'}
    assert not reduction_screens(dict(inp, thermal_energy_density=1.5*P, cooling=noad), t, mass=M)['passed']['thermal']
    with pytest.raises(ValueError, match='never clipped'):
        reduction_screens(dict(inp, state=dict(radius=np.geomspace(2*rc, ro, n))), t, mass=M)


def test_tau_mode_screens_interpolate_in_the_closure_optical_depth():
    """Before the fixed point the closure's tau(r) is not the hydro tau(r): the screen
    must read the layer at tau_c (where a^2, Gamma, kappa were frozen), not at r_c."""
    rc, rb, ro, src = parker_case()
    s = hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)
    tq = hydro_step(FrozenSources('tau', constant(), (np.log(1e-5), np.log(1e3))), mass=M, base_radius=rb,
                    outer_radius=ro, outer_tau=1e-3, base_tau=10., seed=tau_seed(s))
    r = np.geomspace(rb, ro, 60); tau = np.geomspace(1e2, 1e-4, 60); n = len(r)  # closure tau(r), not the hydro one
    E = 1.5e-2*tau; closure = 1e-6*tau**.5  # vary with depth so the lookup coordinate matters
    inp = dict(state=dict(radius=r), tau_rosseland=tau, thermal_energy_density=E, translational_energy_density=E,
               cooling=dict(gross_exchange_lower_bound=np.full(n, 1e3), radiative_gross_exchange=np.full(n, 1e3),
                            channel_closure_relative=closure))
    q = reduction_screens(inp, tq, mass=M)
    tc = tq['sonic_tau']; t_flow = tq['sonic_radius']/abs(2+tq['sonic_slope'])/tq['sonic_velocity']
    assert q['values']['coordinate'] == 'tau'
    assert q['values']['thermal_relaxation_over_flow'] == pytest.approx(1.5e-2*tc/1e3/t_flow, rel=1e-12)
    assert q['values']['gencool_channel_closure'] == pytest.approx(np.interp(np.log(tc), np.log(tau[::-1]), closure[::-1]), rel=1e-12)
    Er = np.exp(np.interp(np.log(tq['sonic_radius']), np.log(r), np.log(E)))  # what a radius lookup would read
    assert abs(np.log(Er/(1.5e-2*tc))) > 1.
    with pytest.raises(ValueError, match='never clipped'):
        reduction_screens(dict(inp, tau_rosseland=np.geomspace(1e2, 1., n)), tq, mass=M)


def test_outer_fixed_point_policy_damping_secant_trust_and_budget():
    c = -3.
    for f in [lambda x: c-1.5*(x-c), lambda x: c+1.8*(x-c)+.1*(x-c)**2]:
        h = []; x = c+.5
        for _ in range(12):
            h.append((x, f(x), 0., 0.)); u = outer_update(h, momentum_passed=True, tolerance=1e-10)
            if u['status'] == 'converged': break
            x = u['next']
        assert u['status'] == 'converged' and abs(x-c) < 1e-9
        assert outer_update(h[:1], momentum_passed=True)['rule'] == 'damped' and outer_update(h[:2], momentum_passed=True)['rule'] == 'secant'
    u = outer_update([(0., 100., 0., 0.)], momentum_passed=True); assert u['truncated'] and u['next'] == pytest.approx(np.log(10.))
    h = [(float(k), k+1., 0., 0.) for k in range(12)]
    assert outer_update(h, momentum_passed=True)['status'] == 'numerical_failure'
    # Every registered converged_when item is required; incomplete rows never pass silently.
    for bad in [[(1., 1.)], [(1., 1., 0.)]]:
        with pytest.raises(ValueError, match='rows'): outer_update(bad, momentum_passed=True)
    with pytest.raises(ValueError, match='momentum_passed'): outer_update([(1., 1., 0., 0.)], momentum_passed=None)
    assert outer_update([(1., 1., 0., 0.)], momentum_passed=True, tolerance=1e-3)['status'] == 'converged'
    u = outer_update([(1., 1., 0., 0.)], momentum_passed=False, tolerance=1e-3)
    assert u['status'] == 'mdot_converged_pending_checks' and u['failing'] == ['momentum']
    u = outer_update([(1., 1., 0., 1.)], momentum_passed=True, tolerance=1e-3)
    assert u['status'] == 'continue' and u['failing'] == ['sonic_radius']


def test_synthetic_force_density_feedback_fixed_point_matches_direct_root():
    rc0, rb, ro, src0 = parker_case()
    x0 = hydro_step(src0, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)['ln_mdot']
    def xhat(x):  # Gamma rises with mdot (density-dependent opacity), frozen within each step
        src = FrozenSources('radius', constant(gamma=.52+.01*(x-x0)), (np.log(rb)-1, np.log(ro)+1))
        return hydro_step(src, mass=M, base_radius=rb, outer_radius=ro, outer_tau=1e-3, base_tau=10.)['ln_mdot']
    root = brentq(lambda x: xhat(x)-x, x0-5, x0+5, xtol=1e-12)
    h = []; x = x0
    for _ in range(12):
        h.append((x, xhat(x), 0., 0.)); u = outer_update(h, momentum_passed=True, tolerance=1e-9)
        if u['status'] == 'converged': break
        x = u['next']
    assert u['status'] == 'converged' and abs(x-root) < 1e-8 and abs(root-x0) > 1e-2


def ledger_inputs():
    r = np.array([1e15, 2e15]); v = np.array([1e5, 4e5]); rho = 1e-10*(r[0]/r)**2*v[0]/v
    E = np.array([1., .25]); F = np.array([.2, .05])*C; P = E/3
    a2 = np.full(2, 1e12); u = 1.5*a2
    return dict(state=dict(radius=r, velocity=v, density=rho, density_clump=rho, filling_factor=np.ones(2)),
                moments=dict(E_lab=E, F_lab=F, P_lab=P), sound_speed_squared=a2, gas_pressure=rho*a2), u


def test_wind_ledger_steady_energy_and_unrealizable_moments():
    q, u = ledger_inputs(); r, v, rho = q['state']['radius'], q['state']['velocity'], q['state']['density']
    mdot = 4*np.pi*r[0]**2*rho[0]*v[0]
    base = mdot*(u[0]+v[0]**2/2-GM/r[0])+4*np.pi*r[0]**2*(q['gas_pressure'][0]*v[0]+q['moments']['F_lab'][0])
    rest = mdot*(u[1]+v[1]**2/2-GM/r[1])+4*np.pi*r[1]**2*q['gas_pressure'][1]*v[1]
    q['moments']['F_lab'][1] = (base-rest)/(4*np.pi*r[1]**2)
    E1 = max(.25, 2*abs(q['moments']['F_lab'][1])/C); q['moments']['E_lab'][1] = E1; q['moments']['P_lab'][1] = E1/3
    l = wind_ledger(q, mass=M, depths=(0, 1), specific_internal_energy=u)
    assert l['status'] == 'evaluated' and abs(l['energy_relative']) < 1e-12 and abs(l['mass_relative']) < 1e-12
    assert l['atomic_energy_included'] and not l['escaping_radiation_certified']
    q['moments']['P_lab'][1] = 2*q['moments']['E_lab'][1]
    assert wind_ledger(q, mass=M, depths=(0, 1))['status'] == 'unrealizable_interface'


def test_grad_audit_separates_inner_luminosity_definition():
    n = 6; r = np.geomspace(1e15, 3e15, n); L = 1e43*(1-.05*np.arange(n)/n); N = np.full(n, 1e10); ne = .9*N
    mu = 1.3; rho = N*mu*1.66053906660e-24; F = L/(4*np.pi*r*r); kF = .5
    g = kF*F/C; gn = g*L[0]/L
    a2 = (N+ne)*1.380649e-16*1e4/rho
    inp = dict(state=dict(radius=r, velocity=np.full(n, 1e5), density_clump=rho, electron_density=ne,
                          filling_factor=np.ones(n), chi_flux=kF*rho),
               moments=dict(radiative_acceleration=g, luminosity_comoving=L, radiative_force_density=g*rho, F_comoving=F),
               hydro=dict(g_rad=gn, sound_speed=np.sqrt(a2), percent_error=np.zeros(n), tau=np.ones(n)),
               atom_density=N, sound_speed_squared=a2, gas_pressure=a2*rho, sigma=np.full(n, -1.), tau_rosseland=np.ones(n))
    q = grad_input_audit(inp, mass=M)
    assert q['L_corrected_relative_max'] < 1e-14 and q['raw_relative_max'] == pytest.approx(.05*(n-1)/n, rel=1e-12)
    assert q['rvtj_mean_atomic_weight']['median'] == pytest.approx(mu, rel=1e-12) and q['flux_mean_opacity_relative_max'] < 1e-14
    assert q['sound_speed_relative_max'] < 1e-14 and abs(q['raw_relative_inner']) < 1e-15
    inp['hydro']['g_rad'] = gn*(1+2e-3)  # a density-convention or opacity defect is not hidden by the L scaling
    assert grad_input_audit(inp, mass=M)['L_corrected_relative_max'] == pytest.approx(2e-3/(1+2e-3), rel=1e-9)


RUN = Path(__file__).resolve().parents[1]/'results/cmfgen_ali_lambda84_long'


@pytest.mark.skipif(not (RUN/'QSTARS_SPECTRAL.bin').exists() or not (RUN/'HYDRO').exists(), reason='Archived CMFGEN export absent')
def test_archived_export_loads_and_frozen_tables_build():
    from qstars.cmfgen_wind import load_wind_inputs, frozen_sources, input_checks
    q = load_wind_inputs(RUN, hashes=False)
    assert q['hydro_sync_radius'] < 2e-6 and len(q['cooling']['channels']) > 3
    assert abs(q['hydrogen_mass_fraction']-.7015) < 1e-3
    assert all(input_checks(q).values())
    for c in ('radius', 'tau'):
        s = frozen_sources(q, 1e6*MSUN, c); assert s.coordinate == c and len(s.nodes) == len(q['state']['radius'])-1
