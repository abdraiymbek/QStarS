"""Optional real-library tests; retain both passing limits and counterexamples."""
import os
import numpy as np
import pytest
from qstars.freeeos_reference import FreeEOSReference
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.eos import Composition
from qstars.atmosphere_reference import integrate_atmosphere

@pytest.fixture(scope='module')
def freeeos():
    library=os.environ.get('QSTARS_FREEEOS_LIBRARY')
    if not library:pytest.skip('Set QSTARS_FREEEOS_LIBRARY for direct FreeEOS tests')
    return FreeEOSReference(library,hydrogen=.75)

@pytest.fixture(scope='module')
def mesa():
    library=os.environ.get('QSTARS_MESA_LIBRARY');directory=os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Compiled MESA atmosphere library required')
    m=MesaMicrophysics(Composition(),library,directory)
    if m.abi_version<5:pytest.skip('ABI 5 reference required')
    return m

def test_freeeos_ionization_gas_first_law_and_inverse(freeeos):
    rho=np.array([1e-16,1e-14,1e-9,1e-7]);T=np.array([6000.,8000.,5000.,1e6])
    x=freeeos.raw(rho,T);q=freeeos.evaluate(rho,T)
    assert q['h_ionized'][0]>.98 and q['free_e'][2]<1e-3
    np.testing.assert_allclose(freeeos.density_from_gas_pressure(q['pg'],T),rho,rtol=1e-11)
    h=1e-4
    for density_axis in (True,False):
        plus=freeeos.raw(rho*np.exp(h) if density_axis else rho,T if density_axis else T*np.exp(h))
        minus=freeeos.raw(rho*np.exp(-h) if density_axis else rho,T if density_axis else T*np.exp(-h))
        du=(plus[3]-minus[3])/(2*h);ds=(plus[6]-minus[6])/(2*h)
        residual=T*ds-du+(q['pg']/rho if density_axis else 0)
        assert max(abs(residual)/(abs(du)+q['pg']/rho))<1e-6
    with pytest.raises(ValueError):freeeos.density_from_gas_pressure(-1.,6000.)

def test_radiation_addition_and_dilute_pressure(freeeos):
    from qstars.constants import KB,MU
    rho,T=1e-9,1e6;q=freeeos.evaluate(rho,T)
    assert abs(q['pg']/(rho*KB*T/(MU*q['mu']))-1)<1e-4
    assert abs(q['gamma1']-4/3)<1e-4 and abs(q['grad_ad']-.25)<1e-7
    h=1e-4
    du=(freeeos.evaluate(rho,T*np.exp(h))['u']-freeeos.evaluate(rho,T*np.exp(-h))['u'])/(2*h*T)
    assert abs(du/q['cv']-1)<1e-7

def test_grey_finite_pressure_and_hydrostatic_inventory(freeeos,mesa):
    from qstars.constants import C
    q=integrate_atmosphere(freeeos,mesa,teff=8000.,gravity=3.,gas_pressure_outer=.1,
                           constant_opacity=.34,points=81,rtol=1e-9)
    F=freeeos.radiation_constant*C*8000.**4/4
    exact=.1+(q.tau[:81]-q.tau[0])*(3/.34-F/C)
    np.testing.assert_allclose(q.gas_pressure[:81],exact,rtol=1e-7)
    np.testing.assert_allclose(q.pressure-q.pressure[0],3*q.column_mass,atol=2e-7,rtol=1e-7)
    assert q.column_mass[-1]>0 and q.column_energy[-1]>0

def test_cool_supersonic_mlt_candidate_is_not_certified(mesa):
    f=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],hydrogen=mesa.composition.X,
                      radiation_constant=mesa.radiation_constant) if os.environ.get('QSTARS_FREEEOS_LIBRARY') else None
    if f is None:pytest.skip('Direct FreeEOS required')
    q=integrate_atmosphere(f,mesa,teff=6000.,gravity=1.,gas_pressure_outer=1e-6,points=121)
    mach=q.convective_velocity/q.sound_speed
    assert max(mach)>2 and q.tau[np.argmax(mach)]>10
    assert max(q.radiative_eddington_factor)>1
