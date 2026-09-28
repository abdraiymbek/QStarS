"""Numerical regressions for previously failed concentrated models."""
import copy
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from qstars.config import resolve
from qstars.evolution import physics,attempt_step,run
from qstars.structure import initialize,continuum_polytrope
from qstars import ledger,numerics,nonlinear
from qstars.constants import PI,YEAR
from qstars.transport import luminosities,inverse_transport_residual
from qstars.io import load_checkpoint


def confined(n=32):
    c=resolve({'numerics':{'zones':n,'mesh':'core_refined','solver_tolerance':1e-10},
         'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},
         'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4}})
    e,op=physics(c);s,_=initialize(c,e);ledger.initialize(s,e)
    return s,e,op,c


def test_inverse_transport_matches_forward_law():
    s,e,op,c=confined()
    for enabled in ('off','mlt_lossless'):
        c['transport']['convection']=enabled
        for amplitude in (-.005,0.,.005):
            s.T*=np.exp(amplitude*np.linspace(0,1,len(s.dm)))
            L,_=luminosities(s,e,op,c,0.)
            np.testing.assert_allclose(inverse_transport_residual(s,e,op,c,L),0.,atol=2e-14)


@pytest.mark.parametrize('mode',['newton','newton_flux'])
def test_colored_jacobian_against_full_complex_step(monkeypatch,mode):
    s,e,op,c=confined(8);c['numerics']['nonlinear_solver']=mode
    original=numerics.newton
    checked=[]
    def audit(fun,x,tolerance,**kwargs):
        jac=kwargs.get('jacobian',nonlinear.local_jacobian)(fun,x).toarray()
        dense=np.empty_like(jac)
        for k in range(len(x)):
            xx=x.astype(complex);xx[k]+=1e-25j;dense[:,k]=fun(xx).imag/1e-25
        np.testing.assert_allclose(jac,dense,rtol=2e-12,atol=1e-12)
        checked.append(True)
        return original(fun,x,tolerance,**kwargs)
    monkeypatch.setattr(numerics,'newton',audit)
    attempt_step(s,1e-4*YEAR,e,op,c)
    assert checked


def test_high_resolution_flux_step_and_energy():
    s,e,op,c=confined(256)
    out,d=attempt_step(s,1e-4*YEAR,e,op,c)
    assert d['nonlinear_residual']<1e-9
    assert abs(d['energy_residual_relative'])<1e-10
    assert max(abs(out.volume/(4*PI/3*np.diff(out.r)*(out.r[1:]**2+out.r[1:]*out.r[:-1]+out.r[:-1]**2))-1))<1e-8


def test_second_order_density_convergence_to_continuum():
    s,e,op,c=confined();b=continuum_polytrope(c,e)
    z=np.linspace(np.log(c['initial']['inner_radius_cm']/c['initial']['radius_cm']),0,30001)
    mass=b.sol(z)[0];mass[0]=0.;mass[-1]=1.
    inv=PchipInterpolator(mass,np.exp(z)*c['initial']['radius_cm'])
    errors=[]
    for n in (32,64,128):
        c['numerics'].update(zones=n,mesh='log_radius');s,_=initialize(c,e)
        q=np.r_[0,np.cumsum(s.dm)/s.dm.sum()];q[-1]=1.;r=inv(q)
        exact=s.dm/(4*PI/3*np.diff(r)*(r[1:]**2+r[1:]*r[:-1]+r[:-1]**2))
        errors.append(np.sqrt(np.dot(s.dm,np.log(s.rho/exact)**2)/s.dm.sum()))
    assert errors[0]/errors[1]>3 and errors[1]/errors[2]>3
    assert errors[-1]<.002


def test_opacity_table_complex_derivative(tmp_path):
    from qstars.opacity import OpacityTable
    p=tmp_path/'powerlaw.npz';x=np.array([-10.,-4.]);y=np.array([4.,7.])
    np.savez(p,logrho=x,logT=y,logkappa=x[:,None]+2*y[None,:],source='analytic kappa=rho*T^2')
    op=OpacityTable(p);r,T=1e-7,1e6;h=1e-25
    np.testing.assert_allclose(op(r*np.exp(1j*h),T).imag/h,op(r,T),rtol=1e-12)
    np.testing.assert_allclose(op(r,T*np.exp(1j*h)).imag/h,2*op(r,T),rtol=1e-12)


def test_boundary_remap_restart_equivalence(tmp_path):
    s,e,op,c=confined(64)
    c['numerics'].update(boundary_remesh_fraction=1e-4,initial_dt_years=1e-4,max_dt_years=1e-4)
    c['evolution']['end_time_years']=.002
    full,rf=run(c,tmp_path)
    half,rh=run(c,tmp_path,stop_after_steps=10);saved,cc=load_checkpoint(half/'checkpoint.npz')
    resumed,rr=run(cc,tmp_path,restart=saved)
    a,_=load_checkpoint(full/'checkpoint.npz');b,_=load_checkpoint(resumed/'checkpoint.npz')
    assert rf['status']==rr['status']=='completed_requested_interval'
    assert 'boundary_remapped' in (full/'events.jsonl').read_text()
    np.testing.assert_array_equal(a.T,b.T);np.testing.assert_array_equal(a.cell_volume,b.cell_volume)
    assert a.ledger==b.ledger
