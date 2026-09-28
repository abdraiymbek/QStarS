"""Sonic regularity, analytical spherical limits and component boundary fluxes."""
import os
import numpy as np
import pytest
from scipy.integrate import solve_ivp, simpson
from qstars.spherical_wind import (regular_sonic_slope, IdealGasOpacity, DiffusionWind,
    FreeEOSOpacity, WindBranchError, spherical_hydrostatic, match_base)
from qstars.constants import G,C,MSUN,RSUN


def test_parker_critical_slope_and_implicit_transonic_solution():
    # Isothermal heat reservoir is explicit; this is not an adiabatic wind.
    def coeff(y):return np.array([2-2*np.exp(-y[0]),np.exp(2*y[1])-1])
    slope,roots=regular_sonic_slope(coeff,0,1e-5)
    assert abs(slope-1)<1e-9
    for sign,bound in [(-1,np.log(.25)),(1,np.log(10.))]:
        start=sign*1e-6
        sol=solve_ivp(lambda x,y:[(2-2*np.exp(-x))/np.expm1(2*y[0])],
                      [start,bound],[slope*start],rtol=1e-10,atol=1e-12)
        assert sol.success
        lhs=np.exp(2*sol.y[0])-2*sol.y[0]
        rhs=4*sol.t+4*np.exp(-sol.t)-3
        np.testing.assert_allclose(lhs,rhs,atol=2e-8,rtol=2e-8)


def test_spherical_radiative_zero_and_enclosed_mass():
    micro=IdealGasOpacity();mass=1e4*MSUN;r=1e14;gamma=.95
    T=G*mass*(1-gamma)/(4*micro.rgas*r)
    rho=micro.arad/(3*micro.rgas)*(1-gamma)/gamma*T**3
    L=4*np.pi*C*G*mass*gamma/micro.kappa
    args=dict(mass=mass,radius=r,density=rho,temperature=T,luminosity=L,outer_ratio=2)
    f=spherical_hydrostatic(micro,**args)
    np.testing.assert_allclose(f['temperature'],T*r/f['radius'],rtol=1e-9)
    np.testing.assert_allclose(f['density'],rho*(r/f['radius'])**3,rtol=1e-9)
    f=spherical_hydrostatic(micro,**args,self_gravity=True,samples=2001)
    shell=simpson(4*np.pi*f['radius']**2*f['density'],x=f['radius'])
    assert abs((f['enclosed_mass'][-1]-mass)/shell-1)<1e-8


def test_wind_energy_work_and_independent_integrals():
    w=DiffusionWind(IdealGasOpacity(),100*MSUN,50*RSUN,9.71e-9,2.5e5)
    f,s=w.integrate(rtol=1e-9,endpoint='diffusion_limit',samples=1501)
    assert s['solver_success'] and s['diffusion_screen_passed']
    assert s['mass_flux_relative_defect']<1e-13
    assert s['energy_algebra_relative_defect']<1e-13
    for key in ['first_law_integral_relative_defect','momentum_integral_relative_defect','diffusion_integral_relative_defect']:
        assert s[key]<2e-7,(key,s[key])
    flux=w.boundary_flux(0.,[0.,0.])
    channels=['radiative_energy','gas_internal_energy','gas_pressure_work',
              'radiation_internal_energy','radiation_pressure_work','kinetic_energy','potential_energy']
    assert abs(sum(flux[k] for k in channels)/w.energy-1)<1e-13
    assert flux['radiation_internal_energy']/flux['radiation_pressure_work']==pytest.approx(3.)
    assert not s['escaping_wind_validated']


def test_inadmissible_sonic_point_not_no_global_root():
    with pytest.raises(WindBranchError,match='nonpositive diffusion'):
        DiffusionWind(IdealGasOpacity(),MSUN,1e16,1e-12,1e5)
    with pytest.raises(ValueError):
        DiffusionWind(IdealGasOpacity(),MSUN,-1.,1e-12,1e5)


@pytest.fixture(scope='module')
def real_micro():
    from qstars.freeeos_reference import FreeEOSReference
    from qstars.mesa_microphysics import MesaMicrophysics
    from qstars.eos import Composition
    a=os.environ.get('QSTARS_FREEEOS_LIBRARY');b=os.environ.get('QSTARS_MESA_LIBRARY');d=os.environ.get('MESA_DIR')
    if not a or not b or not d:pytest.skip('Real FreeEOS and MESA libraries required')
    mesa=MesaMicrophysics(Composition(.7,.3,0),b,d)
    return FreeEOSOpacity(FreeEOSReference(a,.7,mesa.radiation_constant),mesa)


def test_real_recombination_stalled_branch_and_validity(real_micro):
    w=DiffusionWind(real_micro,1e4*MSUN,1.152e15,3e-11,1e4)
    f,s=w.integrate(inner_ratio=.995,outer_ratio=2,rtol=1e-9,sonic_offset=1e-7)
    assert s['solver_success']
    assert s['terminations']['outward']=='approaching_second_sonic_point'
    assert f['hydrogen_ionized'][-1]<.2 and max(f['hydrogen_ionized'])>.9
    assert f['velocity'][-1]<max(f['velocity'])
    assert s['max_velocity_over_escape']<.2
    assert not s['escaping_wind_validated'] and not s['physical_impossibility_claimed']


def test_real_critical_offset_convergence(real_micro):
    w=DiffusionWind(real_micro,1e4*MSUN,1.152e15,3e-11,16000.)
    _,a=w.integrate(inner_ratio=.99,outer_ratio=2,rtol=1e-9,sonic_offset=1e-6)
    _,b=w.integrate(inner_ratio=.99,outer_ratio=2,rtol=1e-9,sonic_offset=1e-7)
    assert abs(a['endpoint_radius']/b['endpoint_radius']-1)<1e-6
    assert abs(a['endpoint_temperature']/b['endpoint_temperature']-1)<1e-5


def test_base_matching_recovers_mass_flux_without_mdot_input(real_micro):
    w=DiffusionWind(real_micro,1e4*MSUN,1.152e15,3e-11,16000.)
    eps=1e-7;x=-.01
    sol=solve_ivp(w.rhs,[-eps,x],-eps*np.array([w.slope,w.sonic_temperature_slope]),
                  method='DOP853',rtol=1e-10,atol=1e-12,max_step=.001)
    r,v,T,rho,q,L,*_=w.state(x,sol.y[:,-1])
    z,report=match_base(real_micro,mass=w.mass,radius=r,density=rho,temperature=T,
              luminosity=L+w.mdot*(v*v/2+q.h),sonic_guess=[w.rs*1.0001,w.rhos*1.001,w.ts*1.0001],rtol=1e-9)
    assert report['max_log_residual']<1e-7
    assert abs(z.mdot/w.mdot-1)<1e-6
    assert report['escape_not_certified']


def test_step_reduction_preserves_published_inputs(real_micro):
    # Some native-library cache histories produce a rejected DOP853 trial.
    # Compare accepted solutions; rejection count is not a physical invariant.
    args=(real_micro,100*MSUN,50*RSUN,9.85e-9,2.5e5)
    w=DiffusionWind(*args)
    opts=dict(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,samples=251)
    f,a=w.integrate(**opts)
    g,b=w.integrate(**opts,max_log_step=.005)
    assert a['solver_success'] and b['solver_success']
    assert abs(a['endpoint_radius']/b['endpoint_radius']-1)<1e-6
    assert not a['diffusion_screen_passed']


def test_paper_opacity_matches_table_nodes_without_mutating_default(real_micro):
    from pathlib import Path
    mesa=real_micro.mesa
    if mesa.abi_version<6:pytest.skip('Historical table audit requires bridge ABI 6')
    rho=np.array([1e-7,1e-10]);T=np.array([1e5,1e4])
    x=real_micro.eos.raw(rho,T)
    before=mesa.opacity_external_eos(rho,T,x[15],x[16]).copy()
    for mixture in ['gs98','gn93']:
        for name,logT in [(mixture,5.),('lowT_af94_gn93',3.8)]:
            path=Path(os.environ['MESA_DIR'])/'data/kap_data'/f'{name}_z0.0_x0.7.data'
            lines=path.read_text().splitlines();logR=np.array([float(t) for t in lines[5].split()])
            grid=np.loadtxt(path,skiprows=7)
            i=int(np.argmin(abs(grid[:,0]-logT)));j=int(np.argmin(abs(logR+4)))
            temp=10**grid[i,0];density=10**(logR[j]+3*grid[i,0]-18)
            actual=mesa.opacity_paper_reference(density,temp,mixture)
            assert abs(np.log10(actual)-grid[i,j+1])<1e-9
    after=mesa.opacity_external_eos(rho,T,x[15],x[16])
    np.testing.assert_array_equal(before,after)
    with pytest.raises(ValueError,match='domain'):
        mesa.opacity_paper_reference(1e-30,1e4)


def test_stall_endpoint_change_does_not_change_overlapping_solution(real_micro):
    w=DiffusionWind(real_micro,100*MSUN,50*RSUN,9.65e-9,2.5e5)
    opts=dict(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,samples=501)
    f,s=w.integrate(**opts,stalled_mach=1.5)
    g,t=w.integrate(**opts)
    assert s['terminations']['outward']=='approaching_second_sonic_point'
    assert f['mach'][-1]==pytest.approx(1.5,rel=1e-8)
    assert g['radius'][-1]>f['radius'][-1]
    from scipy.interpolate import PchipInterpolator
    for key in ['density','temperature','velocity']:
        value=PchipInterpolator(np.log(g['radius']),np.log(g[key]))(np.log(f['radius']))
        assert max(abs(value-np.log(f[key])))<1e-4
    with pytest.raises(ValueError,match='Mach one'):
        w.integrate(stalled_mach=1.)


def test_unreached_comparison_mach_cannot_disable_second_sonic_exit(real_micro):
    w=DiffusionWind(real_micro,1e4*MSUN,1.152e15,3e-11,1e4)
    opts=dict(inner_ratio=.995,outer_ratio=2,rtol=1e-9,sonic_offset=1e-7,samples=251)
    f,s=w.integrate(**opts,stalled_mach=100.)
    g,t=w.integrate(**opts)
    assert s['solver_success'] and s['terminations']['outward']=='approaching_second_sonic_point'
    assert f['mach'][-1]==pytest.approx(np.sqrt(1.01),rel=1e-8)
    assert abs(s['endpoint_radius']/t['endpoint_radius']-1)<1e-8
