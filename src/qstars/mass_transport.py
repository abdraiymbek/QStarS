"""Finite-volume source-cell treatment of two independent boundary streams.

Material joins/leaves the resolved domain in the first/last mass cell, at its
resolved velocity. Injected energy is enthalpy. The gravitational flux uses
the derivative of the discrete binding energy; no residual is zeroed.
"""
import numpy as np
from .constants import G,C,MSUN,YEAR,ME,KB,H,MU,PI
from .controls import value

class DomainError(RuntimeError):
    pass

def exchange(old,dt,eos,config,engine,wind_rate,terminal_speed=None):
    state = old.copy()
    feed = value(config["feeding"]["rate_msun_per_year"],old.time/YEAR)*MSUN/YEAR
    if feed < 0 or wind_rate < 0: raise DomainError("A one-way stream has negative mass flux")
    added,removed,swallowed = dt*feed,dt*wind_rate,dt*engine["rest_mass_rate"]
    if swallowed > .1*old.dm[0] or removed+added > .1*old.dm[-1]:
        from .numerics import StepRejected
        raise StepRejected("Boundary cell replacement exceeds 10%; reduce timestep")
    state.dm[0] -= swallowed
    state.dm[-1] += added-removed
    state.bh_mass += dt*engine["bh_growth_rate"]
    q = eos.evaluate(old.rho,old.T)
    f = config["feeding"]
    if f["thermodynamics"] == "surface":
        incoming_h,incoming_s = q["h"][-1],q["s"][-1]
    else:
        incoming_T=value(f['temperature_K'],old.time/YEAR)
        if added>0:
            if not config['validity']['min_temperature_K']<=incoming_T<=config['validity']['max_temperature_K']:
                raise DomainError('Supplied material lies outside the EOS temperature domain')
            ne=f['density_cgs']/MU*(eos.composition.X+eos.composition.Y/2+eos.composition.Z/2)
            phase=(2*PI*ME*KB*incoming_T/H**2)**1.5
            if config['eos']=='ideal_gas_radiation' and ne/phase>.1:raise DomainError('Supplied material is outside the nondegenerate EOS domain')
            for fraction,chi in ((eos.composition.X,13.5984),(eos.composition.Y,54.4178)):
                ratio=phase/ne*np.exp(-chi*1.602176634e-12/(KB*incoming_T))
                if config['eos']=='ideal_gas_radiation' and fraction>0 and ratio/(1+ratio)<.99:
                    raise DomainError('Supplied material fails fully ionized EOS screening')
        incoming = eos.evaluate(f["density_cgs"],incoming_T)
        incoming_h,incoming_s = incoming["h"],incoming["s"]
    U = old.dm*q["u"]
    U[0] -= swallowed*q["h"][0]
    U[-1] += added*incoming_h-removed*q["h"][-1]
    if np.any(U <= 0):
        from .numerics import StepRejected
        raise StepRejected("Boundary enthalpy flux exhausts internal energy; reduce timestep")
    state.T = eos.temperature(state.rho,U/state.dm)
    # Exact midpoint gradient of the quadratic discrete gravitational energy.
    mid = old.copy()
    mid.dm = .5*(old.dm+state.dm)
    mid.bh_mass = .5*(old.bh_mass+state.bh_mass)
    tail = np.cumsum((mid.weights/mid.r)[::-1])[::-1]
    phi = -G*(.5*(mid.enclosed[:-1]/mid.r[:-1]+mid.enclosed[1:]/mid.r[1:])+tail[1:])
    kinetic = .25*(old.v[:-1]**2+old.v[1:]**2)
    l = state.ledger
    l["supplied_mass"] += added
    l["wind_mass"] += removed
    l["accreted_rest_mass"] += swallowed
    l["mass_defect"] += config["engine"]["efficiency"]*swallowed
    l["supplied_energy"] += added*(incoming_h+phi[-1]+kinetic[-1])
    l["wind_energy"] += removed*(q["h"][-1]+phi[-1]+kinetic[-1])
    l["advected_energy"] += swallowed*(q["h"][0]+phi[0]+kinetic[0])
    l["central_potential_work"] -= G*(state.bh_mass-old.bh_mass)*np.sum(mid.weights/mid.r)
    l["supplied_entropy"] += added*incoming_s
    l["wind_entropy"] += removed*q["s"][-1]
    l["advected_entropy"] += swallowed*q["s"][0]
    delta_s = np.dot(state.dm,eos.evaluate(state.rho,state.T)["s"])-np.dot(old.dm,q["s"])
    l["mixing_entropy_change"] += delta_s-added*incoming_s+removed*q["s"][-1]+swallowed*q["s"][0]
    # Wind must reach infinity. Its binding deficit is charged to outgoing light.
    terminal = config["wind"]["terminal_speed_cm_s"] if terminal_speed is None else terminal_speed
    wind_cost = removed*max(0.,.5*terminal**2-q["h"][-1]-phi[-1]-kinetic[-1])
    return state,dict(feed_rate=feed,wind_rate=wind_rate,wind_cost=wind_cost)
