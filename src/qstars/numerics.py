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
    midpoint=config['numerics']['mechanics']=='midpoint'
    velocity_factor=2 if midpoint else 1
    solve_flux=config['numerics']['nonlinear_solver']=='newton_flux'
    L0,_=luminosities(old,eos,opacity,config,luminosity)
    if solve_flux and old.face_luminosity is not None:
        predictor=L0.copy();predictor[1:-1]=old.face_luminosity[1:-1]
        try:
            inverse_transport_residual(old,eos,opacity,config,predictor)
        except (ValueError,FloatingPointError):pass
        else:L0=predictor
    Lscale=max(np.max(abs(L0)),1.)
    contrast_scale=np.maximum(abs(np.diff(np.log(q0['p'])))*.25,1e-3)
    def residual(x):
        vm = np.r_[0.,x[:n]*cs[1:]]
        state = old.copy()
        # Compensated face motion prevents repeated sub-ULP rounding from
        # drifting away from the independently integrated cell volumes.
        carry=old.radius_compensation if old.radius_compensation is not None else 0.
        increment=dt*vm-carry
        state.r = old.r+increment
        state.radius_compensation=(state.r-old.r)-increment
        state.v = velocity_factor*vm-(velocity_factor-1)*old.v
        state.time = old.time+dt
        state.T = old.T*np.exp(np.clip(x[n:2*n],-100,100))
        if np.any(np.diff(state.r) <= 0):
            raise ValueError('Trial geometry has crossing radial faces')
        area = 4*PI/3*(state.r**2+state.r*old.r+old.r**2)
        state.cell_volume = old.volume+np.diff(area*(dt*vm))
        if np.any(state.cell_volume.real<=0):raise ValueError('Trial cell volume is nonpositive')
        q = eos.evaluate(state.rho,state.T)
        p = .5*(q["p"]+q0["p"]) if midpoint else q['p']
        pressure_force = area[1:]*(p-np.r_[p[1:],pext])/old.weights[1:]
        gravity = G*old.enclosed[1:]/(state.r[1:]*old.r[1:])
        momentum = (velocity_factor*(vm[1:]-old.v[1:])/dt-pressure_force+gravity)/g[1:]
        L,details = luminosities(state,eos,opacity,config,luminosity)
        if solve_flux:L[1:-1]=x[2*n:]*Lscale
        # Volume differences computed from secant area to avoid r^3 cancellation.
        deltaV = np.diff(area*(dt*vm))
        # BE momentum loses 1/2 sum(w_j Delta-v_j^2) of resolved kinetic
        # energy. Deposit exactly that positive amount in neighboring cells;
        # this is derived from the discrete work identity, not a residual fix.
        dv=state.v-old.v
        heating=np.zeros(n) if midpoint else .25*old.dm*(dv[:-1]**2+dv[1:]**2)
        energy = (old.dm*(q["u"]-q0["u"])+p*deltaV+dt*np.diff(L)-heating)/(old.dm*q0["u"])
        last.update(state=state,L=L,details=details,heating=heating)
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
        sol = newton(residual,initial,config['numerics']['solver_tolerance'],
                     max_iterations=max(80,n//8),jacobian=flux_jacobian)
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
    state.ledger['numerical_kinetic_heating']+=float(np.sum(last['heating']))
    state.steps = old.steps+1
    return state,dict(nonlinear_residual=error,nfev=int(sol.nfev),fractional_change=change,
                      luminosity=float(last["L"][-1]),numerical_kinetic_heating_erg=float(np.sum(last['heating'])))
