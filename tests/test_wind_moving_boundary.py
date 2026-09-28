"""Reynolds transport on independently integrated stationary wind profiles."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import solve_ivp,quad
from qstars.constants import G,MSUN,RSUN
from qstars.spherical_wind import DiffusionWind,IdealGasOpacity


def wind_and_state():
    w=DiffusionWind(IdealGasOpacity(),100*MSUN,50*RSUN,9.71e-9,2.5e5)
    solutions=[]
    for sign,bound in [(-1,np.log(.7)),(1,np.log(1.3))]:
        x0=sign*1e-7
        s=solve_ivp(w.rhs,[x0,bound],x0*np.array([w.slope,w.sonic_temperature_slope]),
                    method='DOP853',rtol=1e-11,atol=1e-13,dense_output=True,max_step=.005)
        assert s.success;solutions.append(s)
    def state(r):
        x=np.log(r/w.rs);z=solutions[0 if x<0 else 1].sol(x)
        return x,z,w.state(x,z)
    return w,state


def test_reynolds_transport_for_two_moving_faces():
    w,state=wind_and_state();h=1e-4;dt=h*w.rs/w.vs
    r0=np.array([.8,1.2])*w.rs;u=np.array([.1,-.05])*w.vs;r1=r0+u*dt
    def inventory(r,kind):
        _,_,a=state(r);_,v,T,rho,q,*_=a
        specific=1. if kind=='mass' else q.gas_u+3*q.pr/rho+v*v/2-G*w.mass/r
        return 4*np.pi*r*r*rho*specific
    for kind in ['mass','total_energy']:
        # Directly integrate the swept volumes, independently of boundary code.
        delta=-quad(lambda r:inventory(r,kind),r0[0],r1[0],epsabs=0,epsrel=1e-10)[0]
        delta+=quad(lambda r:inventory(r,kind),r0[1],r1[1],epsabs=0,epsrel=1e-10)[0]
        flux=[]
        for r,speed in zip((r0+r1)/2,u):
            x,z,_=state(r);flux.append(w.boundary_flux(x,z,face_velocity=speed)[kind])
        assert abs(delta/dt-(flux[0]-flux[1]))/max(abs(flux[0]),abs(flux[1]))<2e-7


def test_comoving_face_retains_pressure_work_and_radiation():
    w,state=wind_and_state();x,z,a=state(.8*w.rs);r,v,T,rho,q,L,*_=a
    f=w.boundary_flux(x,z,face_velocity=v)
    assert abs(f['mass'])<1e-14*w.mdot
    assert f['total_energy']==pytest.approx(L+4*np.pi*r*r*v*(q.pg+q.pr),rel=1e-13)
    assert f['momentum']==pytest.approx(4*np.pi*r*r*(q.pg+q.pr),rel=1e-13)
    assert f['entropy']==pytest.approx(0.,abs=1.)
    with pytest.raises(ValueError,match='finite'):
        w.boundary_flux(x,z,face_velocity=float('nan'))


def test_energy_gauge_changes_only_relative_mass_advection():
    base=IdealGasOpacity();offset=1e14
    class Shifted:
        arad=base.arad
        def __call__(self,rho,T):
            q=base(rho,T)
            return replace(q,gas_u=q.gas_u+offset,h=q.h+offset)
    a=DiffusionWind(base,100*MSUN,50*RSUN,9.71e-9,2.5e5)
    b=DiffusionWind(Shifted(),100*MSUN,50*RSUN,9.71e-9,2.5e5)
    assert abs(a.slope/b.slope-1)<1e-7
    for speed in [0.,.2*a.vs,a.vs,1.1*a.vs]:
        f=a.boundary_flux(0.,[0.,0.],face_velocity=speed)
        g=b.boundary_flux(0.,[0.,0.],face_velocity=speed)
        assert abs((g['total_energy']-f['total_energy'])-f['mass']*offset)<1e-13*abs(a.energy)
