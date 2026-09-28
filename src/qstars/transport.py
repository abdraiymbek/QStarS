"""LTE radiative diffusion and optional optically thick local mixing length.

Convection follows the lossless-blob limit, not full MESA MLT. The cap is
P*c_s, a declared causal/subsonic transport limit. Saturation is diagnosed.
"""
import numpy as np
from .constants import PI,A,C,G,YEAR
from .controls import value

def luminosities(state,eos,opacity,config,inner_luminosity):
    n = len(state.dm)
    q = eos.evaluate(state.rho,state.T)
    kap = opacity(state.rho,state.T)
    rad,conv = np.zeros(n+1),np.zeros(n+1)
    rho = .5*(state.rho[1:]+state.rho[:-1])
    k = .5*(kap[1:]+kap[:-1])
    radius = state.r[1:-1]
    dr = np.diff(state.rc)
    rad[1:-1] = -4*PI*radius**2*A*C/(3*k*rho)*np.diff(state.T**4)/dr
    saturation = np.zeros(n-1)
    if config["transport"]["convection"] == "mlt_lossless":
        p = .5*(q["p"][1:]+q["p"][:-1])
        T = .5*(state.T[1:]+state.T[:-1])
        g = G*state.enclosed[1:-1]/radius**2
        hp = p/(rho*g)
        grad_ad = .5*(q["grad_ad"][1:]+q["grad_ad"][:-1])
        cp = .5*(q["cp"][1:]+q["cp"][:-1])
        delta = .5*(q["delta"][1:]+q["delta"][:-1])
        cs = .5*(q["cs"][1:]+q["cs"][:-1])
        length = config["transport"]["mixing_length_alpha"]*hp
        # Discrete entropy contrast: using an independently differenced Hp here
        # spuriously destabilizes an exact P~T^4 polytrope on a coarse mesh.
        excess = np.maximum(-hp*(np.diff(np.log(state.T))-grad_ad*np.diff(np.log(q["p"])))/dr,0)
        velocity = np.sqrt(g*delta/(8*hp))*length*np.sqrt(excess)
        flux = rho*cp*T*velocity*length/(2*hp)*excess
        saturation = flux/(p*cs)
        conv[1:-1] = 4*PI*radius**2*np.minimum(flux,p*cs)
    elif config["transport"]["convection"] != "off":
        raise ValueError("Unsupported convection closure")
    boundary = config["boundary"]
    if boundary["thermal"] == "grey":
        tau = .5*kap[-1]*state.dm[-1]/(4*PI*state.r[-1]**2) + boundary["external_optical_depth"]
        rad[-1] = 4*PI*state.r[-1]**2*A*C/3*(state.T[-1]**4-boundary["bath_temperature_K"]**4)/(tau+2/3)
    elif boundary["thermal"] == "insulated": rad[-1] = 0.
    elif boundary["thermal"] == "luminosity":
        rad[-1] = value(boundary["luminosity_erg_s"],state.time/YEAR)
    else: raise ValueError("Unknown outer thermal boundary")
    total = rad+conv
    total[0] = inner_luminosity
    return total,dict(radiative=rad,convective=conv,opacity=kap,saturation=saturation)
