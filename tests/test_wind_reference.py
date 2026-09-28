import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from qstars.wind_reference import solve_wind

def test_published_wind_independent_integrals_and_asymptotes():
    # Figure 2's terminal-velocity range, not a value fitted by QStarS.
    f,s=solve_wind(1.,5.,100.,samples=4001)
    assert 1.95<s['w_infinity']<2.1
    assert abs(f['mach'][0]-.8)<1e-6
    assert s['base_optical_depth_proxy']>1
    assert 2<s['photospheric_optical_depth_proxy']<3.2
    assert s['sonic_radius_ratio']<s['photospheric_radius_ratio']
    assert abs(s['momentum_integral_residual'])<1e-6
    assert abs(s['radiation_integral_residual_relative'])<1e-6
    assert s['energy_algebraic_closure_relative']<1e-13
    assert np.all(f['luminosity']>0)

def test_advective_wind_numerical_and_outer_boundary_convergence():
    a,sa=solve_wind(7.64,5.,100.)
    b,sb=solve_wind(7.64,5.,100.,delta=1e-7,rtol=1e-10)
    r=np.geomspace(1,1e4,101)
    for k in ('velocity_ratio','luminosity','advected_luminosity','u'):
        x=PchipInterpolator(np.log(a['radius_ratio']),np.log(a[k]))(np.log(r))
        y=PchipInterpolator(np.log(b['radius_ratio']),np.log(b[k]))(np.log(r))
        assert max(abs(x-y))<1e-4
    assert a['advected_luminosity'][0]>a['luminosity'][0]*5

def test_invalid_wind_boundary_rejected():
    with pytest.raises(ValueError):solve_wind(.4,5.,100.,base_mach=1.)
