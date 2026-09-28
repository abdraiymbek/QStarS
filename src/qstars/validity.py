"""A failure to remain in the model domain is not a physical disruption."""
import numpy as np
from .constants import G,C,H,ME,KB,MU,PI
from .ledger import residuals

def check(state,eos,opacity,config):
    errors=[]
    q=eos.evaluate(state.rho,state.T)
    v=config["validity"]
    if min(state.T)<v["min_temperature_K"]: errors.append("EOS temperature below ionized-domain floor")
    if max(state.T)>v["max_temperature_K"]: errors.append("EOS temperature above configured domain")
    mach=float(max(abs(.5*(state.v[:-1]+state.v[1:]))/q["cs"]))
    if mach>v["max_mach"]: errors.append("Quasi-hydrostatic domain exceeded: radial Mach number")
    compactness=float(max(2*G*state.enclosed/(state.r*C**2)))
    if compactness>v["max_compactness"]: errors.append("Newtonian compactness limit exceeded")
    ne=state.rho/MU*(eos.composition.X+eos.composition.Y/2+eos.composition.Z/2)
    phase=(2*PI*ME*KB*state.T/H**2)**1.5
    eV=1.602176634e-12
    hydrogen_ratio=phase/ne*np.exp(-13.5984*eV/(KB*state.T))
    helium_ratio=phase/ne*np.exp(-54.4178*eV/(KB*state.T))
    if eos.composition.X>0 and min(hydrogen_ratio/(1+hydrogen_ratio))<.99:
        errors.append('Fully ionized EOS fails hydrogen Saha screening')
    if eos.composition.Y>0 and min(helium_ratio/(1+helium_ratio))<.99:
        errors.append('Fully ionized EOS fails He II/III Saha screening')
    degeneracy=ne*(H/np.sqrt(2*PI*ME*KB*state.T))**3
    if max(degeneracy)>.1: errors.append("Classical EOS fails electron nondegeneracy check")
    tau=opacity(state.rho,state.T)*state.rho*np.diff(state.r)
    if v["require_optically_thick"] and min(tau)<1: errors.append("Diffusion cell is optically thin")
    from .transport import luminosities
    _,transport=luminosities(state,eos,opacity,config,0.)
    from .constants import A
    free_stream=4*PI*state.r[1:-1]**2*C*A*.5*(state.T[:-1]**4+state.T[1:]**4)
    if np.any(abs(transport['radiative'][1:-1])>free_stream):
        errors.append('Radiative diffusion exceeds free-streaming flux')
    conservation=residuals(state,eos)
    tol=config["numerics"]["conservation_tolerance"]
    if max(abs(conservation["mass_residual_relative"]),abs(conservation["energy_residual_relative"]))>tol:
        errors.append("Conservation tolerance exceeded")
    return errors,dict(max_mach=mach,max_compactness=compactness,min_zone_optical_depth=float(min(tau)),
                       minimum_helium_ionized_fraction=float(min(helium_ratio/(1+helium_ratio))),
                       max_temperature_jump=float(max(abs(np.diff(np.log(state.T)))))),conservation

def assumptions(config):
    out=["Development solver: complete QStarS v1 validation gates are not satisfied.",
         "Fixed excision boundary; accretion luminosity laws do not imply replication of a full published boundary model.",
         "Newtonian, fixed-composition, fully ionized EOS; no recombination, nuclear evolution, rotation, magnetic fields, or GR.",
         "Inflow and outflow are separate, co-moving boundary source cells; multidimensional stream geometry is unresolved.",
         "Radial inertia is retained; no nonlinear pulsation or shock validation has been established."]
    if config["boundary"]["thermal"]=="grey":
        out.append("Grey boundary is a constant-opacity extrapolation; reported effective temperature is a flux diagnostic, not a validated cool photosphere.")
    if config["transport"]["convection"]=="mlt_lossless":
        out.append("Convection uses the optically thick lossless-blob approximation; radiative losses of blobs and time-dependent convection are omitted.")
    if config["wind"]["model"]=="dotan2011":
        out.append("Dotan launch formula assumes a sonic, porous boundary; the porous atmosphere and transonic wind are not resolved.")
    return out
