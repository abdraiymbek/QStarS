"""Compatible staggered Lagrangian momentum and first-law energy solve.

Exact volume secant and inverse-radius secant make mechanical work cancel
globally to nonlinear tolerance. Heat transport is backward Euler, mechanics
implicit midpoint. No artificial damping, gravity correction, or energy reset.
"""
import numpy as np
from scipy.optimize import root
from .constants import G,PI
from .transport import luminosities,inverse_transport_residual
from .mass_transport import DomainError
from .nonlinear import newton,flux_jacobian

class StepRejected(RuntimeError):
    pass

def thermal_mechanical_step(old,dt,eos,opacity,config,luminosity):
    n = len(old.dm)
    q0 = eos.evaluate(old.rho,old.T)
    cs = np.r_[q0["cs"][0],.5*(q0["cs"][:-1]+q0["cs"][1:]),q0["cs"][-1]]
    g = G*old.enclosed/old.r**2
    pext = config["boundary"]["pressure_cgs"]
    last = {}
    solve_flux=config['numerics']['nonlinear_solver']=='newton_flux'
    L0,_=luminosities(old,eos,opacity,config,luminosity)
    Lscale=max(np.max(abs(L0)),1.)
    contrast_scale=np.maximum(abs(np.diff(np.log(q0['p'])))*.25,1e-3)
    def residual(x):
        vm = np.r_[0.,x[:n]*cs[1:]]
        state = old.copy()
        state.r = old.r+dt*vm
        state.v = 2*vm-old.v
        state.time = old.time+dt
        state.T = old.T*np.exp(np.clip(x[n:2*n],-100,100))
        if np.any(np.diff(state.r) <= 0):
            raise ValueError('Trial geometry has crossing radial faces')
        area = 4*PI/3*(state.r**2+state.r*old.r+old.r**2)
        state.cell_volume = old.volume+np.diff(area*(dt*vm))
        if np.any(state.cell_volume.real<=0):raise ValueError('Trial cell volume is nonpositive')
        q = eos.evaluate(state.rho,state.T)
        p = .5*(q["p"]+q0["p"])
        pressure_force = area[1:]*(p-np.r_[p[1:],pext])/old.weights[1:]
        gravity = G*old.enclosed[1:]/(state.r[1:]*old.r[1:])
        momentum = (2*(vm[1:]-old.v[1:])/dt-pressure_force+gravity)/g[1:]
        L,details = luminosities(state,eos,opacity,config,luminosity)
        if solve_flux:L[1:-1]=x[2*n:]*Lscale
        # Volume differences computed from secant area to avoid r^3 cancellation.
        deltaV = np.diff(area*(dt*vm))
        energy = (old.dm*(q["u"]-q0["u"])+p*deltaV+dt*np.diff(L))/(old.dm*q0["u"])
        last.update(state=state,L=L,details=details)
        if solve_flux:
            closure=inverse_transport_residual(state,eos,opacity,config,L)/contrast_scale
            return np.r_[momentum,energy,closure]
        return np.r_[momentum,energy]
    contraction=-dt*np.diff(old.v)
    shrinking=contraction>0
    predictor_scale=min(1.,float(np.min(.25*np.diff(old.r)[shrinking]/contraction[shrinking]))) if np.any(shrinking) else 1.
    initial = np.r_[predictor_scale*old.v[1:]/cs[1:],np.zeros(n)]
    if solve_flux:
        initial=np.r_[initial,L0[1:-1]/Lscale]
        sol = newton(residual,initial,config['numerics']['solver_tolerance'],jacobian=flux_jacobian)
    elif config['numerics']['nonlinear_solver']=='newton':
        sol = newton(residual,initial,config['numerics']['solver_tolerance'])
    else:
        sol = root(residual,initial,method="hybr",options={"xtol":config["numerics"]["solver_tolerance"],"maxfev":30*(2*n+1)})
    error = float(max(abs(residual(sol.x))))
    if not np.isfinite(error) or error > config["numerics"]["solver_tolerance"]*10:
        raise StepRejected(f"Coupled solve residual={error:g}; {sol.message}")
    state = last["state"]
    state.face_luminosity=last['L'].copy()
    change = float(max(np.max(abs(np.log(state.T/old.T))),np.max(abs(np.log(state.r/old.r)))))
    if change > config["numerics"]["max_fractional_change"]:
        raise StepRejected(f"Fractional state change {change:g} exceeds limit")
    state.ledger["central_generated"] += dt*luminosity
    state.ledger["radiated_energy"] += dt*last["L"][-1]
    state.ledger["boundary_work"] -= pext*4*PI/3*(state.r[-1]**3-old.r[-1]**3)
    state.steps = old.steps+1
    return state,dict(nonlinear_residual=error,nfev=int(sol.nfev),fractional_change=change,
                      luminosity=float(last["L"][-1]))
