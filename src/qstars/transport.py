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
    dtype = np.result_type(state.T, state.r)
    rad,conv = np.zeros(n+1,dtype=dtype),np.zeros(n+1,dtype=dtype)
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
        contrast=eos.thermal_contrast(state.rho,state.T,q)
        excess = np.maximum(-hp*contrast/dr,0)
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


def inverse_transport_residual(state,eos,opacity,config,L):
    """Temperature contrast needed by a specified face luminosity.

    Algebraic inverse of the same diffusion + capped lossless MLT law.
    Solving for flux avoids multiplying roundoff in nearly adiabatic
    gradients by an enormous convective conductance.
    """
    q=eos.evaluate(state.rho,state.T);kap=opacity(state.rho,state.T)
    rho=.5*(state.rho[1:]+state.rho[:-1]);k=.5*(kap[1:]+kap[:-1])
    radius=state.r[1:-1];dr=np.diff(state.rc)
    p=.5*(q['p'][1:]+q['p'][:-1]);T=.5*(state.T[1:]+state.T[:-1])
    g=G*state.enclosed[1:-1]/radius**2;hp=p/(rho*g)
    contrast=eos.thermal_contrast(state.rho,state.T,q)
    pref=4*PI*radius**2*A*C/(3*k*rho*dr)*state.T[:-1]**4
    dlnT=np.log1p(np.diff(state.T)/state.T[:-1])
    rad=-pref*np.expm1(4*dlnT)
    pref_ad=pref*np.exp(4*(dlnT-contrast))
    Lad=rad-pref*np.exp(4*dlnT)*np.expm1(-4*contrast)
    deltaL=L[1:-1]-Lad
    ratio=dr/hp
    # Pure radiative inversion on stable faces.
    argument=1-deltaL/pref_ad
    if config['transport']['convection']=='off':
        if np.any(argument.real<=0):raise ValueError('Radiative inverse has no positive temperature')
        return contrast+.25*(-np.log1p(-deltaL/pref_ad))
    unstable=deltaL.real>0
    D=np.zeros_like(deltaL)
    if np.any(argument[~unstable].real<=0):raise ValueError('Radiative inverse outside temperature domain')
    D[~unstable]=-.25*np.log1p(-deltaL[~unstable]/pref_ad[~unstable])/ratio[~unstable]
    if np.any(unstable):
        cp=.5*(q['cp'][1:]+q['cp'][:-1]);delta=.5*(q['delta'][1:]+q['delta'][:-1])
        cs=.5*(q['cs'][1:]+q['cs'][:-1]);length=config['transport']['mixing_length_alpha']*hp
        coefficient=4*PI*radius**2*rho*cp*T*np.sqrt(g*delta/(8*hp))*length**2/(2*hp)
        cap=4*PI*radius**2*p*cs
        m=unstable;target=deltaL[m];cc=coefficient[m];rr=ratio[m];pp=pref_ad[m];lim=cap[m]
        d=np.minimum((target/cc)**(2/3),target/(4*pp*rr))
        for _ in range(30):
            conv=cc*d**1.5;free=conv.real<lim.real
            flux=-pp*np.expm1(-4*d*rr)+np.minimum(conv,lim)
            derivative=4*pp*rr*np.exp(-4*d*rr)+np.where(free,1.5*cc*np.sqrt(d),0)
            correction=(flux-target)/derivative
            if np.max(abs(correction)/(abs(d)+1e-300))<2e-14:break
            d=np.maximum(d-correction,d*.1)
        else:raise ValueError('Inverse MLT solve failed')
        D[m]=d
    return contrast+D*ratio
