"""Finite-pressure grey atmosphere diagnostic with native Henyey MLT below.

This is a boundary *reference*, not an envelope evolution boundary. Pressure
at the outer optical depth is an explicit physical input, never a density
floor. No LTE, hydrostatic or thin-atmosphere validity is inferred from solver
success. Thick-layer inventories and scale heights are returned for screening.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from .constants import C

@dataclass
class AtmosphereProfile:
    tau: np.ndarray
    gas_pressure: np.ndarray
    temperature: np.ndarray
    density: np.ndarray
    opacity: np.ndarray
    pressure: np.ndarray
    inward_depth: np.ndarray
    column_mass: np.ndarray
    column_energy: np.ndarray
    radiative_gradient: np.ndarray
    temperature_gradient: np.ndarray
    convective_velocity: np.ndarray
    total_flux_eddington_factor: np.ndarray
    radiative_eddington_factor: np.ndarray
    sound_speed: np.ndarray

def integrate_atmosphere(eos,mesa,*,teff,gravity,gas_pressure_outer,
                         tau_outer=1e-6,tau_match=2/3,tau_base=100.,
                         alpha=2.,mlt_model='Henyey',rtol=1e-8,points=201,constant_opacity=None):
    """Grey radiative layer to tau_match, radiative+MLT envelope below.

    dP/dtau=g/kappa; dT/dtau=(T/P)*nabla*g/kappa in the diffusion
    region. Above the match use Eddington T(tau), subtract its radiation
    pressure gradient exactly in the gas-pressure equation. Column energy
    includes LTE radiation and gas, and is an inventory, not a flux.
    """
    if not all(np.isfinite(v) and v>0 for v in (teff,gravity,gas_pressure_outer,tau_outer,tau_match,tau_base,alpha,rtol)):
        raise ValueError('Atmosphere boundary inputs must be finite and positive')
    if not tau_outer<tau_match<tau_base or points<4:raise ValueError('Invalid atmosphere optical-depth mesh')
    if mlt_model not in ('Henyey','ML1'):raise ValueError('Unsupported native MLT model')
    if constant_opacity is not None and (not np.isfinite(constant_opacity) or constant_opacity<=0):
        raise ValueError('Opacity must be positive')
    if eos.hydrogen!=mesa.composition.X:raise ValueError('Atmosphere EOS and opacity composition differ')
    arad=eos.radiation_constant
    if abs(arad/mesa.radiation_constant-1)>1e-12:raise ValueError('Atmosphere radiation constants differ')
    flux=arad*C*teff**4/4
    def state(tau,y,thick):
        pg,T=y[:2]
        if pg<=0 or T<=0:raise ValueError('No positive-gas-pressure hydrostatic state on this branch')
        if not thick:T=teff*(.75*(tau+2/3))**.25
        rho=float(eos.density_from_gas_pressure(pg,T));q=eos.evaluate(rho,T)
        k=float(mesa.opacity_external_eos(rho,T,q['free_e'],q['eta'])) if constant_opacity is None else constant_opacity
        p=float(q['p']);gradr=3*k*p*flux/(4*arad*C*gravity*T**4)
        if thick:
            grad,velocity,_,_=mesa.mlt([q['chiT'],q['chiRho'],q['cp'],gravity,rho,p,T,k,gradr,q['grad_ad']],model=mlt_model,alpha=alpha)
            dT=T/p*float(grad)*gravity/k
        else:
            dT=3*teff**4/(16*T**3);grad=gradr;velocity=0.
        dpg=gravity/k-4*arad*T**3/3*dT
        return np.array([dpg,dT,1/(k*rho),1/k,float(q['u'])/k]),(rho,k,p,gradr,float(grad),float(velocity),k*flux/gravity/C,T,float(q['cs']))
    def rhs(logtau,y,thick):
        tau=np.exp(logtau);return tau*state(tau,y,thick)[0]
    T0=teff*(.75*(tau_outer+2/3))**.25
    y0=np.array([gas_pressure_outer,T0,0.,0.,0.])
    meshes=[np.geomspace(tau_outer,tau_match,points),np.geomspace(tau_match,tau_base,points)]
    profiles=[]
    for thick,mesh in zip((False,True),meshes):
        sol=solve_ivp(lambda x,y:rhs(x,y,thick),np.log(mesh[[0,-1]]),y0,
                      t_eval=np.log(mesh),rtol=rtol,atol=[gas_pressure_outer*rtol,teff*rtol,1.,1e-12,1.],method='DOP853')
        if not sol.success:raise RuntimeError(f'Atmosphere integration failed: {sol.message}')
        values=np.array([state(tau,y,thick)[1] for tau,y in zip(mesh,sol.y.T)])
        if not thick:sol.y[1]=values[:,7]
        profiles.append((mesh,sol.y,values));y0=sol.y[:,-1]
    tau=np.r_[profiles[0][0],profiles[1][0][1:]]
    y=np.c_[profiles[0][1],profiles[1][1][:,1:]]
    v=np.r_[profiles[0][2],profiles[1][2][1:]]
    return AtmosphereProfile(tau,y[0],y[1],v[:,0],v[:,1],v[:,2],y[2],y[3],y[4],v[:,3],v[:,4],v[:,5],v[:,6],v[:,6]*v[:,4]/v[:,3],v[:,8])
