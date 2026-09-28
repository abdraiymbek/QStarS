import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from qstars.grey_transfer import (SphericalRays, solve_grey, solve_moving_grey,
                                  comoving_moments, TransferConvergenceError)
from qstars.continuum_wind import flux_driven_reference, solve_continuum_wind, dimensional_fields
from qstars.constants import MSUN


def test_transparent_spherical_source_all_three_moments():
    r = np.geomspace(1, 1000, 151)
    rays = SphericalRays(r, 1e-12/r**2)
    j, h, k = rays.moments(np.zeros_like(r), np.ones(rays.core_rays))
    mu = np.sqrt(1-1/r**2)
    exact_j = 1/(2*r*r*(1+mu))
    assert np.max(abs(j/exact_j-1)) < 1e-7
    assert np.max(abs(h/(.25/r**2)-1)) < 1e-7
    assert np.max(abs(k/(exact_j*(1+mu+mu*mu)/3)-1)) < 1e-7


def test_published_spherical_moments_and_independent_flux():
    r = np.unique(np.r_[np.geomspace(1/100.1, 10, 401), 10*(1-np.geomspace(1e-6, .05, 40))])
    q = solve_grey(r, r**-2, 16*np.pi**2)
    # Hummer & Rybicki Table I, final iteration, interior and escaping field.
    rr = np.array([.04, .4, 1, 4, 10])
    reference = np.array([15730, 20.12, 1.884, .07637, .01040])
    j = np.exp(PchipInterpolator(np.log(r), np.log(q.mean_intensity))(np.log(rr)))
    assert max(abs(j/reference-1)) < .005
    assert q.history[-1]['angular_luminosity_max_error'] < .005
    assert q.eddington_factor[0] < .335
    assert q.eddington_factor[-1] > .92


def test_lorentz_moment_transform_against_angular_integrals():
    mu, weight = np.polynomial.legendre.leggauss(128)
    for beta in [0, .01, .5]:
        gamma = 1/np.sqrt(1-beta**2)
        intensity = (gamma*(1-beta*mu))**-4
        j, h, k = [np.dot(weight, mu**n*intensity)/2 for n in range(3)]
        j0, h0, k0 = comoving_moments(j, h, k, beta)
        assert abs(j0-1) < 1e-12
        assert abs(h0) < 1e-12
        assert abs(k0-1/3) < 1e-12


def test_moving_radiation_does_work_not_constant_luminosity():
    r = np.unique(np.r_[np.geomspace(1, 100, 301), 1+np.geomspace(1e-7, .05, 50)])
    chi, beta = 3/r**2, np.full_like(r, .01)
    q = solve_moving_grey(r, chi, beta)
    ell = r*r*q['H']
    # Independent angular intensity, not the moment ODE's flux variable.
    _, ah0, _ = comoving_moments(q['angular_J'], q['angular_H'], q['angular_K'], beta)
    integrand = r*r*chi*beta/np.sqrt(1-beta*beta)*ah0
    work = PchipInterpolator(r, integrand).integrate(r[0], r[-1])
    assert .02 < 1-ell[-1] < .04
    assert abs(ell[-1]-ell[0]+work) < 1e-4
    assert q['history'][-1]['angular_flux_max_error'] < .005


def test_failed_closure_retains_history():
    r = np.geomspace(.01, 10, 61)
    with pytest.raises(TransferConvergenceError) as exc:
        solve_grey(r, r**-2, 1., max_iterations=1)
    assert len(exc.value.history) == 1


def test_transfer_domains_are_explicit():
    r = np.geomspace(1, 10, 20)
    with pytest.raises(ValueError):
        solve_grey(r, np.zeros_like(r), 1.)
    with pytest.raises(ValueError, match='validation domain'):
        solve_moving_grey(r, 1/r**2, np.full_like(r, .1))
    with pytest.raises(ValueError):
        solve_grey(r, 1/r**2, 1., moment_tolerance=np.nan)


def test_wind_reference_includes_finite_base_and_work():
    r = np.geomspace(1, 100, 100)
    w, ell = flux_driven_reference(r, 2, .001, .002)
    assert w[0] == .002 and ell[0] == 1
    assert max(abs(ell+.001*(w-1/r)-(1+.001*(.002-1)))) < 1e-14


def test_continuum_wind_energy_and_missing_thermodynamics():
    r = np.unique(np.r_[1, np.geomspace(1.00001, 100, 161),
                        1/(1-np.geomspace(1e-8, .2, 100))])
    q = solve_continuum_wind(r, gamma=2, tiring=.001, escape_speed_c=.0001)
    w, ell = flux_driven_reference(r, 2, .001, 1e-4)
    assert max(abs(np.sqrt(q['w']/w)-1)) < .005
    assert max(abs(q['luminosity_ratio']/ell-1)) < 1e-4
    assert np.ptp(q['total_energy_ratio']) < 1e-6
    f = dimensional_fields(q, 100*MSUN, .34)
    assert f['photospheric_radius_cm'] > f['radius_cm'][0]
    assert f['radiation_energy_erg_cm3'][-1] > 0
    assert f['gas_temperature_K'] is None and f['sonic_radius_cm'] is None
    assert q['mass_loss_is_input'] and not q['physical_envelope_match']


def test_moment_quadrature_on_irregular_knots_against_closed_form():
    from qstars.grey_transfer import _moment_propagation
    r = np.unique(np.r_[np.geomspace(.01, 10, 51), .100000001, 1.000000001])
    t = np.log(r)
    q = _moment_propagation(t, np.full_like(r, 1/3), r**-2, np.zeros_like(r), [1/3, 1.], 1e-11)
    expected_j = r**-3 + (1/100-1/1000)
    assert max(abs((3*q[0]/r**2)/expected_j-1)) < 1e-9
    assert max(abs(q[1]-1)) < 1e-13
