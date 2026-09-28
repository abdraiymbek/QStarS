"""Hydrostatic initialization on the same staggered mesh as the evolution."""
import numpy as np
from scipy.optimize import least_squares
from .constants import PI,G,MSUN
from .state import State

def initialize(config,eos):
    c,n = config["initial"],config["numerics"]["zones"]
    # Cosine mass mesh resolves both boundaries more finely than the bulk.
    q = .5*(1-np.cos(np.linspace(0,np.pi,n+1)))
    if config['numerics']['mesh']=='core_refined': q=q*q
    dm = c["envelope_mass_msun"]*MSUN*np.diff(q)
    r0,R = c["inner_radius_cm"],c["radius_cm"]
    r = (r0**3+(R**3-r0**3)*q)**(1/3)
    r[0] = r0
    state = State(0.,c["bh_mass_msun"]*MSUN,dm,r,np.zeros(n+1),np.ones(n))
    vscale = 4*PI/3*(R**3-r0**3)
    pext = config["boundary"]["pressure_cgs"]
    pscale = G*state.enclosed[-1]**2/(4*PI*R**4)
    rhoscale = dm.sum()/vscale
    Kscale = pscale/rhoscale**(4/3)
    def unpack(x):
        V = vscale*np.exp(x[:n])
        rr = np.r_[r0,(r0**3+3*np.cumsum(V)/(4*PI))**(1/3)]
        rho = dm/V
        p = Kscale*np.exp(x[-1])*rho**(4/3)
        return rr,rho,p
    def residual(x):
        rr,rho,p = unpack(x)
        support = G*state.enclosed[1:]*state.weights[1:]/(4*PI*rr[1:]**4)
        hydro = np.log(p/(np.r_[p[1:],pext]+support))
        return np.r_[hydro,np.log(rr[-1]/R)]
    x0 = np.r_[np.log(np.diff(q)),0.]
    sol = least_squares(residual,x0,xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=1500)
    err = float(max(abs(residual(sol.x))))
    if not sol.success or err > 1e-8:
        raise RuntimeError(f"Hydrostatic initialization failed: residual {err:g}")
    state.r,rho,p = unpack(sol.x)
    state.T = eos.temperature(rho,p,"p")
    perturb = c["temperature_perturbation"]
    state.T *= 1+perturb*np.exp(-((q[:-1]+q[1:])/2-c["perturbation_mass_center"])**2/
                               (2*c["perturbation_mass_width"]**2))
    state.next_dt = config["numerics"]["initial_dt_years"]*365.25*86400
    return state,err

def hydrostatic_residual(state,eos,pext):
    p = eos.evaluate(state.rho,state.T)["p"]
    acceleration = 4*PI*state.r[1:]**2*(p-np.r_[p[1:],pext])/state.weights[1:]
    gravity = G*state.enclosed[1:]/state.r[1:]**2
    return (acceleration-gravity)/gravity
