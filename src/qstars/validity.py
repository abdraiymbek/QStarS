"""A failure to remain in the model domain is not a physical disruption."""
import numpy as np
from .constants import G,C,H,ME,KB,MU,PI
from .ledger import residuals
from .numerics import StepRejected

def check(state,eos,opacity,config):
    errors=[]
    q=eos.evaluate(state.rho,state.T)
    eos_defect=float(max(eos.consistency(state.rho,state.T))) if hasattr(eos,'consistency') else 0.
    if eos_defect>.01:errors.append('MESA first-law defect exceeds 1% audit limit')
    gas_defect=float(np.max(eos.gas_pressure_limit(state.rho,state.T,q)['screened_relative_error'])) if hasattr(eos,'gas_pressure_limit') else 0.
    if gas_defect>.01:errors.append('MESA gas pressure differs by more than 1% from the dilute classical limit')
    conduction=float(np.max(q.get('conductive_fraction',0.)))
    if conduction>.01:errors.append('Omitted electron conduction exceeds 1% of diffusive heat flux')
    v=config["validity"]
    if min(state.T)<v["min_temperature_K"]: errors.append("EOS temperature below configured domain")
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
    ideal=config['eos']=='ideal_gas_radiation'
    if ideal and eos.composition.X>0 and min(hydrogen_ratio/(1+hydrogen_ratio))<.99:
        errors.append('Fully ionized EOS fails hydrogen Saha screening')
    if ideal and eos.composition.Y>0 and min(helium_ratio/(1+helium_ratio))<.99:
        errors.append('Fully ionized EOS fails He II/III Saha screening')
    degeneracy=ne*(H/np.sqrt(2*PI*ME*KB*state.T))**3
    if ideal and max(degeneracy)>.1: errors.append("Classical EOS fails electron nondegeneracy check")
    geometric=4*PI/3*np.diff(state.r)*(state.r[1:]**2+state.r[1:]*state.r[:-1]+state.r[:-1]**2)
    volume_error=float(max(abs(geometric/state.volume-1)))
    if volume_error>1e-8:raise StepRejected('Stored cell volumes disagree with radial geometry')
    kap=opacity(state.rho,state.T)
    tau=kap*state.rho*np.diff(state.r)
    # Diffusion concerns the physical gradient scale, not the numerical cell
    # width. Requiring tau_cell>1 incorrectly rejects convergent refinement.
    rho_face=.5*(state.rho[:-1]+state.rho[1:]);kap_face=.5*(kap[:-1]+kap[1:])
    kn=4*abs(np.diff(np.log(state.T)))/(np.diff(state.rc)*rho_face*kap_face)
    if v["require_optically_thick"] and max(kn)>.1:
        errors.append('Radiation gradient scale is outside the diffusion limit (Kn>0.1)')
    from .transport import luminosities
    total_luminosity,transport=luminosities(state,eos,opacity,config,0.)
    significant=abs(transport['convective'][1:-1])>.01*(abs(transport['radiative'][1:-1])+abs(transport['convective'][1:-1]))
    min_pe=float(min(transport['peclet'][significant])) if np.any(significant) else None
    if min_pe is not None and min_pe<100:
        errors.append('Lossless MLT outside efficient-convection domain (significant flux with Pe<100)')
    if max(transport['convective_mach'],default=0)>1:
        errors.append('Lossless MLT predicts supersonic convection')
    from .constants import A
    outer_kn=3*abs(total_luminosity[-1])/(4*PI*state.r[-1]**2*C*A*state.T[-1]**4)
    if v['require_optically_thick'] and outer_kn>.1:
        errors.append('Outer radiation flux is outside the LTE diffusion limit')
    free_stream=4*PI*state.r[1:-1]**2*C*A*.5*(state.T[:-1]**4+state.T[1:]**4)
    if np.any(abs(transport['radiative'][1:-1])>free_stream):
        errors.append('Radiative diffusion exceeds free-streaming flux')
    conservation=residuals(state,eos)
    tol=config["numerics"]["conservation_tolerance"]
    if max(abs(conservation["mass_residual_relative"]),abs(conservation["energy_residual_relative"]))>tol:
        raise StepRejected("Conservation tolerance exceeded")
    return errors,dict(max_mach=mach,max_compactness=compactness,min_zone_optical_depth=float(min(tau)),
                       max_eos_first_law_defect=eos_defect,
                       max_eos_dilute_gas_pressure_defect=gas_defect,
                       max_omitted_conductive_fraction=conduction,
                       max_radiation_knudsen=float(max(kn)),
                       outer_radiation_knudsen=float(outer_kn),
                       minimum_significant_convective_peclet=min_pe,
                       max_convective_mach=float(max(transport['convective_mach'],default=0)),
                       ideal_eos_helium_screen=float(min(helium_ratio/(1+helium_ratio))),
                       max_volume_geometry_error=volume_error,
                       max_temperature_jump=float(max(abs(np.diff(np.log(state.T)))))),conservation

def assumptions(config):
    out=["Development solver: complete QStarS v1 validation gates are not satisfied.",
         "Fixed excision boundary; accretion luminosity laws do not imply replication of a full published boundary model.",
         "Newtonian, fixed composition; no nuclear evolution, rotation, magnetic fields, or GR.",
         "Inflow and outflow are separate, co-moving boundary source cells; multidimensional stream geometry is unresolved.",
         "Radial inertia is retained; no nonlinear pulsation or shock validation has been established."]
    if config['eos']=='mesa':
        out.append('Direct MESA 24.08.1 EOS and opacity, fixed H/He composition; this does not supply MESA stellar evolution, convection, or atmosphere.')
    else:out.append('Fully ionized ideal gas plus radiation EOS; recombination and degeneracy excluded.')
    if config["boundary"]["thermal"]=="grey":
        out.append("Grey boundary is a constant-opacity extrapolation; reported effective temperature is a flux diagnostic, not a validated cool photosphere.")
    if config["transport"]["convection"]=="mlt_lossless":
        out.append("Convection uses the optically thick lossless-blob approximation; radiative losses of blobs and time-dependent convection are omitted.")
    if config["wind"]["model"]=="dotan2011":
        out.append("Dotan launch formula assumes a sonic, porous boundary; the porous atmosphere and transonic wind are not resolved.")
    return out
