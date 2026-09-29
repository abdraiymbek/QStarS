"""Independent material-energy and base-selected polytropic wind checks."""
import os
import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.integrate import simpson
from qstars.constants import G,MSUN
from qstars.material_wind import GasDerivatives,continuum_coefficients,continuum_rhs,adiabatic_eigenvalue_reference


def test_material_first_law_and_momentum_with_independent_heating():
    gamma=1.4;rho=1e-10;P=1e3;u=P/((gamma-1)*rho)
    gas=GasDerivatives(P,u,0.,P,P,0.,u)
    state=dict(radius=1e13,density=rho,velocity=8e6,mass=100*MSUN,
               gas_heating_per_volume=1e-4,radiative_acceleration=4.)
    y,t=continuum_rhs(gas,**state)
    d=-2-y;r=state['radius'];v=state['velocity']
    assert u*t-P/rho*d==pytest.approx(r*state['gas_heating_per_volume']/(rho*v),rel=1e-12)
    assert v*v*y+P/rho*(d+t)+G*state['mass']/r-r*state['radiative_acceleration']==pytest.approx(0.,abs=1.)
    q=continuum_coefficients(gas,**state)
    assert q['adiabatic_sound_speed_squared']==pytest.approx(gamma*P/rho)
    assert q['adiabatic_sound_speed_squared']!=q['isothermal_sound_speed_squared']


def test_self_selected_adiabatic_profiles_against_implicit_bernoulli_and_rematching():
    M=100*MSUN;rb=1e12;rho=1e-10;gamma=4/3;P=rho*G*M/(rb*2.5*gamma)
    q=adiabatic_eigenvalue_reference(mass=M,radius=rb,density=rho,pressure=P,rtol=1e-10)
    assert q['classification']=='stationary_wind_reference'
    assert q['base_density_relative_error']<1e-7 and q['base_pressure_relative_error']<1e-7
    # Solve independent Bernoulli algebra at each radius, choosing density
    # above/below its sonic minimum; no ODE critical derivative is reused.
    K=P/rho**gamma;B=q['terminal_speed']**2/2
    for i in [0,100,350,650,950,1000]:
        r=q['radius'][i];flux=q['mdot']/(4*np.pi*r*r)
        sonic_density=(flux**2/(gamma*K))**(1/(gamma+1))
        def energy(logrho):
            d=np.exp(logrho)
            return flux**2/(2*d*d)+gamma*K*d**(gamma-1)/(gamma-1)-G*M/r-B
        center=np.log(sonic_density)
        lo,hi=(center,center+30) if r<q['sonic_radius'] else (center-30,center)
        exact=np.exp(brentq(energy,lo,hi,xtol=1e-13))
        assert q['density'][i]==pytest.approx(exact,rel=1e-7)
    np.testing.assert_allclose(q['energy_flux'],q['energy_flux'][0],rtol=1e-8)
    work=simpson(-G*M/q['radius']**2,x=q['radius'])
    enthalpy=q['internal_energy']+q['pressure']/q['density']
    assert abs((q['velocity'][-1]**2-q['velocity'][0]**2)/2+enthalpy[-1]-enthalpy[0]-work)/(G*M/rb)<1e-6
    # Move the base to another subsonic point, supplying no previous mdot/v.
    i=200
    moved=adiabatic_eigenvalue_reference(mass=M,radius=q['radius'][i],density=q['density'][i],
                                        pressure=q['pressure'][i],gamma=gamma)
    assert moved['mdot']==pytest.approx(q['mdot'],rel=1e-7)


def test_adiabatic_existence_classification_has_a_scoped_proof():
    M=100*MSUN;r=1e12;rho=1e-10;gamma=4/3
    for lam in [3.1,4.,10.]:
        q=adiabatic_eigenvalue_reference(mass=M,radius=r,density=rho,pressure=rho*G*M/(r*lam*gamma))
        assert q['classification']=='hydrostatic_no_wind_reference' and q['mdot']==0
        surface=q['vacuum_surface_radius']
        # Independent hydrostatic enthalpy reaches zero at the reported surface.
        hbase=G*M/(r*lam)/(gamma-1)
        assert abs(hbase-G*M*(1/r-1/surface))/hbase<1e-13
    q=adiabatic_eigenvalue_reference(mass=M,radius=r,density=rho,pressure=rho*G*M/(r*1.9*gamma))
    assert q['classification']=='model_domain_exit'
    with pytest.raises(ValueError,match='domain'):
        adiabatic_eigenvalue_reference(mass=M,radius=r,density=rho,pressure=1.,gamma=5/3)


def test_freeeos_material_critical_speed_matches_native_gas_gamma1():
    from qstars.freeeos_reference import FreeEOSReference
    lib=os.environ.get('QSTARS_FREEEOS_LIBRARY')
    if not lib:pytest.skip('Native FreeEOS library required')
    eos=FreeEOSReference(lib)
    for rho,T in [(3e-11,5000.),(3e-11,10000.),(3e-11,16000.),(1e-8,2e5)]:
        raw=eos.raw(rho,T);gas=GasDerivatives.freeeos(raw)
        q=continuum_coefficients(gas,radius=1e13,density=rho,velocity=1e7,mass=100*MSUN,
                                 gas_heating_per_volume=0.,radiative_acceleration=0.)
        assert q['adiabatic_sound_speed_squared']==pytest.approx(raw[14]*np.exp(raw[0])/rho,rel=3e-7)
