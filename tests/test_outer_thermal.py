"""Independent thermal transfer and non-LTE moving-interface identities."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from qstars.constants import C, SIGMA
from qstars.thermal_transfer import ThermalBands, planck_bands, opacity_means
from qstars.grey_transfer import TransferConvergenceError
from qstars.outer_interface import InterfaceState, moving_flux, compare_faces


def example_state():
    return InterfaceState(radius=2e13, enclosed_mass=1e36, density=1e-10,
                          velocity=1e6, gas_pressure=200., gas_internal_energy=3e12,
                          radiation_energy=1000., radiation_flux=.6*C*1000.,
                          radiation_pressure_rr=500., potential=-2e14,
                          energy_gauge='test_gas_ground_state', gas_temperature=1e4)


def test_general_reynolds_flux_from_conserved_density():
    s=example_state();A=4*np.pi*s.radius**2;u=-3e5
    a,b=moving_flux(s),moving_flux(s,u)
    conserved={'mass':s.density,
               'momentum':s.density*s.velocity+s.radiation_flux/C**2,
               'total_energy':s.density*(s.gas_internal_energy+s.velocity**2/2+s.potential)+s.radiation_energy}
    for key,density in conserved.items():
        assert b[key]-a[key]==pytest.approx(-A*u*density,rel=1e-12)
    comoving=moving_flux(s,s.velocity)
    assert comoving['mass']==0
    assert comoving['total_energy']==pytest.approx(A*(s.gas_pressure*s.velocity+s.radiation_flux-s.velocity*s.radiation_energy))
    offset=2e13
    shifted=moving_flux(replace(s,gas_internal_energy=s.gas_internal_energy+offset),u)
    assert shifted['total_energy']-b['total_energy']==pytest.approx(b['mass']*offset,rel=1e-12)


def test_shock_jump_conserves_flux_without_temperature_continuity():
    # Analytic normal-shock Rankine-Hugoniot solution, ideal gamma=5/3, Mach 3.
    gamma=5/3;mach=3.;rho=1e-10;P=200.;v=mach*np.sqrt(gamma*P/rho)
    compression=(gamma+1)*mach**2/((gamma-1)*mach**2+2)
    pressure_ratio=(2*gamma*mach**2-(gamma-1))/(gamma+1)
    a=replace(example_state(),density=rho,velocity=v,gas_pressure=P,
              gas_internal_energy=P/((gamma-1)*rho))
    b=replace(a,density=rho*compression,velocity=v/compression,gas_pressure=P*pressure_ratio,
              gas_internal_energy=P*pressure_ratio/((gamma-1)*rho*compression),
              gas_temperature=a.gas_temperature*pressure_ratio/compression)
    result=compare_faces(a,b)
    assert result['conservation_passed']
    assert b.gas_temperature/a.gas_temperature>3
    assert not result['continuity_imposed'] and not result['physical_solution_certified']
    with pytest.raises(ValueError,match='gauges'):
        compare_faces(a,replace(b,energy_gauge='different'))
    with pytest.raises(ValueError,match='geometry'):
        compare_faces(a,replace(b,radius=a.radius*1.01))


def test_interface_rejects_nonrealizable_radiation_and_exposes_mdot_energy():
    a=example_state()
    with pytest.raises(ValueError,match='realizable'):
        moving_flux(replace(a,radiation_flux=1.1*C*a.radiation_energy))
    b=replace(a,velocity=a.velocity*2)
    result=compare_faces(a,b)
    assert not result['conservation_passed']
    assert result['defects']['mass']['relative']==pytest.approx(.5)
    assert abs(result['defects']['total_energy']['difference'])>0


def test_band_planck_integrals_and_derivatives_against_adaptive_quadrature():
    edges=np.array([0.,1000,1e4,1e5,np.inf]);T=np.array([3000.,1e4,1e5])
    B,dB=planck_bands(T,edges)
    for i,t in enumerate(T):
        for g in range(len(edges)-1):
            lo=edges[g]/t;hi=min(edges[g+1]/t,100)
            expected=quad(lambda x:x**3/np.expm1(x),lo,hi,epsabs=1e-12,epsrel=1e-12)[0]*15/np.pi**4*SIGMA/np.pi*t**4
            assert B[g,i]==pytest.approx(expected,rel=2e-8,abs=1e-12*SIGMA/np.pi*t**4)
    np.testing.assert_allclose(B.sum(axis=0),SIGMA/np.pi*T**4,rtol=1e-13)
    h=1e-5
    numeric=(planck_bands(T*np.exp(h),edges)[0]-planck_bands(T*np.exp(-h),edges)[0])/(2*h*T)
    np.testing.assert_allclose(dB,numeric,rtol=1e-6,atol=1e-5)


def test_thin_nongrey_equilibrium_does_not_set_gas_to_radiation_temperature():
    r=np.geomspace(1,20,81);edges=[0,20000,np.inf];star=20000.
    ka=np.array([1.,10.])[:,None]*np.ones((2,len(r)))
    model=ThermalBands(r,ka*1e-10,ka*1e-11,edges,star)
    q=model.radiative_equilibrium(tolerance=1e-9)
    W=.5*(1-np.sqrt(1-1/r**2));Bstar=planck_bands([star],edges)[0][:,0]
    exact=[]
    for w in W:
        exact.append(brentq(lambda t:np.dot([1.,10.],planck_bands([t],edges)[0][:,0]-w*Bstar),100.,star))
    np.testing.assert_allclose(q['gas_temperature'],exact,rtol=1e-7)
    assert q['gas_temperature'][-1]/q['radiation_temperature'][-1]>1.3
    assert q['thermal_balance_error']<1e-9
    means=opacity_means(ka,ka*.1,q['gas_temperature'],edges,q['J'],q['H'])
    assert max(means['flux_total']/means['rosseland_total'])>5
    assert not q['population_kinetics_solved'] and not q['wind_eigenvalue_solved']


def exact_constant_absorption(r,outer,chi,source,core):
    """Independent chord solution integrated by adaptive angular quadrature."""
    mu0=np.sqrt(1-1/r**2)
    def intensity(mu,core_ray):
        p=r*np.sqrt(max(0.,1-mu*mu));z=r*mu
        far=np.sqrt(max(0.,outer**2-p*p))
        inward=source*(-np.expm1(-chi*(far-z)))
        distance=z-np.sqrt(max(0.,1-p*p)) if core_ray else z+far
        initial=core if core_ray else 0.
        outward=source+(initial-source)*np.exp(-chi*distance)
        return outward,inward
    result=[]
    for power in range(3):
        def func(mu,core_ray):
            out,inc=intensity(mu,core_ray)
            return .5*mu**power*(out+(-1)**power*inc)
        result.append(quad(lambda mu:func(mu,False),0,mu0,epsabs=1e-12)[0]
                      +quad(lambda mu:func(mu,True),mu0,1,epsabs=1e-12)[0])
    return np.asarray(result)


def test_emission_absorption_spherical_chords_and_energy_integral():
    r=np.geomspace(1,3,1201);chi=.7
    model=ThermalBands(r,np.full((1,len(r)),chi),np.zeros((1,len(r))),[0,np.inf],1.,core_rays=64)
    q=model.fixed_temperature(np.full(len(r),.8))
    factor=SIGMA/np.pi
    for i in [0,200,600,1000,1200]:
        exact=exact_constant_absorption(r[i],r[-1],chi,.8**4,1.)*factor
        np.testing.assert_allclose(np.array([q[k][0,i] for k in ['J','H','K']]),exact,rtol=4e-4)
    assert q['energy_integral_relative_defect']<2e-4
    # Nonzero material heating is balanced by loss of radiation, not erased.
    assert abs(q['luminosity'][-1]/q['luminosity'][0]-1)>.1


def test_thermal_failures_preserve_history_and_do_not_infer_nonexistence():
    r=np.geomspace(1,10,41);a=np.array([.2,2.])[:,None]/r[None,:]**2
    model=ThermalBands(r,a,a*.2,[0,20000,np.inf],20000)
    with pytest.raises(TransferConvergenceError) as exc:
        model.radiative_equilibrium(max_iterations=1)
    assert exc.value.history and exc.value.history[0]['thermal_balance_error']>1e-3
    with pytest.raises(ValueError,match='coupling'):
        ThermalBands(r,np.zeros_like(a),a,[0,20000,np.inf],20000)
