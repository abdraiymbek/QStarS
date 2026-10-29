"""A failure to remain in the model domain is not a physical disruption.

check(..., policy=None) is the original screen set. With a declared policy
the screen registry (SCREENS) classifies every screen as enforced,
report-only, not applicable or diagnostic. 'published_assumption' is the
Ball et al. 2011 replication policy of benchmarks/ball2011_track_plan.json
(replication_mode): convective Mach > 1, radiation Knudsen and free
streaming at tau_from_surface < 10, the outer-face Eddington-relation check
(which replaces the outer-face Knudsen screen) and the STARS opacity clamp
are report-only; density inversion and the local Eddington factor are
diagnostics; everything else is enforced. 'production' enforces every
applicable screen. The Pe < 100 screen applies only to lossless closures.
With validity.require_optically_thick False the three diffusion-limit
screens are not_applicable, with that reason, rather than silently
unenforced. The opacity_clamp screen also covers the photospheric point
(rho_ph, Teff), whose kappa enters the Ball photosphere row.
"""
import numpy as np
from .constants import G,C,H,ME,KB,MU,PI
from .ledger import residuals
from .numerics import StepRejected

def check(state,eos,opacity,config,policy=None,face_luminosity=None):
    """Domain screens of an accepted state: (errors, diagnostics, conservation).

    policy=None runs the original screens unchanged (bitwise). A policy
    ('production', 'published_assumption' or a mapping, see POLICIES) runs the
    screen registry: enforced violations are errors, report-only violations
    are listed in diagnostics['screen_warnings'], and every screen's value,
    limit and action is in diagnostics['screens']. face_luminosity (the
    solved face luminosities of a flux-unknown step) replaces the forward
    transport map in the transport screens; it needs a policy.
    """
    if policy is None:
        if face_luminosity is not None:raise ValueError('Screens on solved face luminosities need a declared policy')
        return _legacy_check(state,eos,opacity,config)
    return _screened_check(state,eos,opacity,config,policy,face_luminosity)

def _legacy_check(state,eos,opacity,config):
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


# id: (kind, message). kind 'screen' may be enforced or report-only by policy,
# 'diagnostic' is never an error, 'reject' raises StepRejected.
SCREENS={
    'nonfinite_state':('screen','Nonfinite model state'),
    'eos_temperature_domain':('screen','EOS temperature outside the configured domain'),
    'mesa_first_law':('screen','MESA first-law defect exceeds 1% audit limit'),
    'mesa_ideal_backstop':('screen','MESA neutral ideal backstop is outside supported ionization physics'),
    'mesa_dilute_gas':('screen','MESA gas pressure differs by more than 1% from the dilute classical limit'),
    'omitted_conduction':('screen','Omitted electron conduction exceeds 1% of diffusive heat flux'),
    'ideal_eos_hydrogen_saha':('screen','Fully ionized EOS fails hydrogen Saha screening'),
    'ideal_eos_helium_saha':('screen','Fully ionized EOS fails He II/III Saha screening'),
    'ideal_eos_degeneracy':('screen','Classical EOS fails electron nondegeneracy check'),
    'compactness':('screen','Newtonian compactness limit exceeded'),
    'flow_mach':('screen','Quasi-hydrostatic domain exceeded: radial Mach number'),
    'photosphere_gas_pressure':('screen','Photospheric gas pressure is not positive (Gamma_ph >= 1)'),
    'radiation_knudsen_deep':('screen','Radiation gradient scale is outside the diffusion limit (Kn>0.1) at tau>=10'),
    'free_streaming_deep':('screen','Radiative diffusion exceeds free-streaming flux at tau>=10'),
    'radiation_knudsen_outer':('screen','Radiation gradient scale is outside the diffusion limit (Kn>0.1) at tau<10'),
    'free_streaming_outer':('screen','Radiative diffusion exceeds free-streaming flux at tau<10'),
    'outer_face_knudsen':('screen','Outer radiation flux is outside the LTE diffusion limit'),
    'outer_face_eddington':('screen','Outer cell departs from the Eddington grey relation'),
    'convective_mach':('screen','Mixing-length convection is supersonic'),
    'convective_peclet':('screen','Lossless MLT outside efficient-convection domain (significant flux with Pe<100)'),
    'opacity_clamp':('screen','Opacity clamped to the table edge (STARS nearest-value rule)'),
    'density_inversion':('diagnostic','Largest outward relative density increase'),
    'local_eddington_factor':('diagnostic','Largest radiative luminosity over the local Eddington luminosity'),
    'volume_geometry':('reject','Stored cell volumes disagree with radial geometry'),
    'conservation':('reject','Conservation tolerance exceeded')}

POLICIES={
    'production':dict(report_only=(),replaced=('outer_face_eddington',),
                      description='Every applicable screen is enforced.'),
    'published_assumption':dict(report_only=('convective_mach','radiation_knudsen_outer','free_streaming_outer',
                                             'outer_face_eddington','opacity_clamp'),replaced=('outer_face_knudsen',),
                                description='Ball et al. 2011 published-assumption replica: the STARS surface and '
                                            'convection assumptions are reported, not enforced (plan replication_mode).')}
EDDINGTON_RELATION_TOLERANCE=1e-8

def resolve_policy(policy):
    """A named policy or a mapping with report_only and replaced screen ids."""
    if isinstance(policy,str):
        if policy not in POLICIES:raise ValueError(f'Unknown screen policy {policy}')
        return policy,POLICIES[policy]
    rules=dict(policy);name=rules.pop('name','custom')
    unknown=set(rules.get('report_only',()))|set(rules.get('replaced',()))
    if unknown-set(SCREENS) or set(rules)-{'report_only','replaced','description'}:
        raise ValueError('Screen policy names unknown screens or keys')
    return name,dict(report_only=tuple(rules.get('report_only',())),replaced=tuple(rules.get('replaced',())),
                     description=rules.get('description',''))

def applicability(config,eos,opacity):
    """Screens that do not apply to this configuration, with the reason."""
    out={}
    model=config['transport']['convection']
    if model not in ('mlt_lossless','mlt_thesis_lossless'):
        out['convective_peclet']='Pe < 100 screens the lossless closure; '+('ML1 includes radiative losses' if model=='mlt_ml1' else 'no lossless convection')
    if config['eos']!='ideal_gas_radiation':
        for k in ('ideal_eos_hydrogen_saha','ideal_eos_helium_saha','ideal_eos_degeneracy'):out[k]='ideal-gas EOS screens only'
    if not hasattr(eos,'consistency'):out['mesa_first_law']='EOS has no first-law audit'
    if not hasattr(eos,'raw'):out['mesa_ideal_backstop']='EOS has no neutral backstop fraction'
    if not hasattr(eos,'gas_pressure_limit'):out['mesa_dilute_gas']='EOS has no dilute-gas reference'
    if not hasattr(opacity,'clamp_flags'):out['opacity_clamp']='opacity has no table-edge clamp'
    if config['boundary']['thermal']!='stars_photosphere':out['photosphere_gas_pressure']='no Ball photosphere'
    if config['boundary']['thermal'] not in ('stars_photosphere','grey'):out['outer_face_eddington']='no grey surface relation'
    if not config['validity']['require_optically_thick']:
        # The diffusion-limit screens are switched off by configuration; say so instead of
        # reporting them as enforced (the Ball plan lists radiation_knudsen_deep as always enforced).
        for k in ('radiation_knudsen_deep','radiation_knudsen_outer','outer_face_knudsen'):
            out[k]='validity.require_optically_thick is False (diffusion-limit screens disabled by configuration)'
    return out

def describe(policy,config,eos,opacity):
    """Action of every screen under a policy and configuration (validity_warnings.json)."""
    name,rules=resolve_policy(policy);skip=applicability(config,eos,opacity)
    table={}
    for key,(kind,message) in SCREENS.items():
        table[key]=dict(message=message,action=_action(key,kind,rules,skip),reason=skip.get(key))
    return dict(policy=name,description=rules['description'],screens=table)

def _action(key,kind,rules,skip):
    if kind!='screen':return kind
    if key in skip or key in rules['replaced']:return 'not_applicable'
    return 'report_only' if key in rules['report_only'] else 'enforced'

def _screened_check(state,eos,opacity,config,policy,face_luminosity):
    from .constants import A,SIGMA
    from .transport import luminosities,_mixing_length_faces
    from .mlt import MODELS,closure,conductance,radiative_gradient
    name,rules=resolve_policy(policy);skip=applicability(config,eos,opacity)
    v=config['validity'];screens={};errors=[];warnings=[]
    def record(key,value,limit,violated):
        kind,message=SCREENS[key];action=_action(key,kind,rules,skip)
        violated=bool(violated) and action not in ('not_applicable','diagnostic')
        screens[key]=dict(value=None if value is None else float(value),limit=limit,violated=violated,action=action)
        if violated and action=='enforced':errors.append(message)
        if violated and action=='report_only':warnings.append(key)
    arrays=[state.r,state.T,state.dm,state.volume,state.v]+([np.asarray(face_luminosity)] if face_luminosity is not None else [])
    nonfinite=sum(int(np.sum(~np.isfinite(a))) for a in arrays)
    record('nonfinite_state',nonfinite,0,nonfinite>0)
    if nonfinite:
        return errors,dict(policy=name,screens=screens,screen_warnings=warnings),{}
    q=eos.evaluate(state.rho,state.T)
    eos_defect=float(max(eos.consistency(state.rho,state.T))) if hasattr(eos,'consistency') else 0.
    record('mesa_first_law',eos_defect,.01,eos_defect>.01)
    backstop=float(np.max(eos.raw(state.rho,state.T)[18])) if hasattr(eos,'raw') else 0.
    record('mesa_ideal_backstop',backstop,1e-6,backstop>1e-6)
    gas_defect=float(np.max(eos.gas_pressure_limit(state.rho,state.T,q)['screened_relative_error'])) if hasattr(eos,'gas_pressure_limit') else 0.
    record('mesa_dilute_gas',gas_defect,.01,gas_defect>.01)
    conduction=float(np.max(q.get('conductive_fraction',0.)))
    record('omitted_conduction',conduction,.01,conduction>.01)
    tmin,tmax=float(min(state.T)),float(max(state.T))
    record('eos_temperature_domain',tmin if tmin<v['min_temperature_K'] else tmax,
           [v['min_temperature_K'],v['max_temperature_K']],tmin<v['min_temperature_K'] or tmax>v['max_temperature_K'])
    mach=float(max(abs(.5*(state.v[:-1]+state.v[1:]))/q['cs']))
    record('flow_mach',mach,v['max_mach'],mach>v['max_mach'])
    compactness=float(max(2*G*state.enclosed/(state.r*C**2)))
    record('compactness',compactness,v['max_compactness'],compactness>v['max_compactness'])
    ne=state.rho/MU*(eos.composition.X+eos.composition.Y/2+eos.composition.Z/2)
    phase=(2*PI*ME*KB*state.T/H**2)**1.5
    eV=1.602176634e-12
    hydrogen_ratio=phase/ne*np.exp(-13.5984*eV/(KB*state.T))
    helium_ratio=phase/ne*np.exp(-54.4178*eV/(KB*state.T))
    hydrogen=float(min(hydrogen_ratio/(1+hydrogen_ratio)));helium=float(min(helium_ratio/(1+helium_ratio)))
    degeneracy=float(max(ne*(H/np.sqrt(2*PI*ME*KB*state.T))**3))
    record('ideal_eos_hydrogen_saha',hydrogen,.99,eos.composition.X>0 and hydrogen<.99)
    record('ideal_eos_helium_saha',helium,.99,eos.composition.Y>0 and helium<.99)
    record('ideal_eos_degeneracy',degeneracy,.1,degeneracy>.1)
    geometric=4*PI/3*np.diff(state.r)*(state.r[1:]**2+state.r[1:]*state.r[:-1]+state.r[:-1]**2)
    volume_error=float(max(abs(geometric/state.volume-1)))
    screens['volume_geometry']=dict(value=volume_error,limit=1e-8,violated=volume_error>1e-8,action='reject')
    if volume_error>1e-8:raise StepRejected('Stored cell volumes disagree with radial geometry')
    kap=opacity(state.rho,state.T)
    tau=kap*state.rho*np.diff(state.r)
    bc=config['boundary']
    external=2/3 if bc['thermal']=='stars_photosphere' else bc['external_optical_depth'] if bc['thermal']=='grey' else 0.
    tau_face=external+np.cumsum(tau[::-1])[::-1][1:]
    rho_face=.5*(state.rho[:-1]+state.rho[1:]);kap_face=.5*(kap[:-1]+kap[1:])
    kn=4*abs(np.diff(np.log(state.T)))/(np.diff(state.rc)*rho_face*kap_face)
    deep=tau_face>=10
    model=config['transport']['convection'];alpha=config['transport']['mixing_length_alpha']
    if face_luminosity is not None and model in MODELS:
        L=np.asarray(face_luminosity,dtype=float)
        faces,mass=_mixing_length_faces(state,q,kap)
        out=closure(model,radiative_gradient(L[1:-1],faces['kappa'],faces['p'],faces['T'],mass),**faces,alpha=alpha)
        radiative=conductance(faces['kappa'],faces['p'],faces['T'],mass)*out['gradient']
        convective=L[1:-1]-radiative
        cs=.5*(q['cs'][1:]+q['cs'][:-1]);length=alpha*faces['p']/(faces['rho']*faces['g'])
        convective_mach=out['velocity']/cs
        peclet=3*faces['kappa']*faces['rho']**2*faces['cp']*out['velocity']*length/(4*A*C*faces['T']**3)
    else:
        L,transport=luminosities(state,eos,opacity,config,0.)
        if face_luminosity is not None:L=np.asarray(face_luminosity,dtype=float)
        radiative,convective=transport['radiative'][1:-1],transport['convective'][1:-1]
        convective_mach,peclet=transport['convective_mach'],transport['peclet']
    kn_deep=float(max(kn[deep],default=0.));kn_outer=float(max(kn[~deep],default=0.))
    record('radiation_knudsen_deep',kn_deep,.1,kn_deep>.1)
    record('radiation_knudsen_outer',kn_outer,.1,kn_outer>.1)
    free=4*PI*state.r[1:-1]**2*C*A*.5*(state.T[:-1]**4+state.T[1:]**4)
    streaming=abs(radiative)/free
    record('free_streaming_deep',max(streaming[deep],default=0.),1.,np.any(streaming[deep]>1))
    record('free_streaming_outer',max(streaming[~deep],default=0.),1.,np.any(streaming[~deep]>1))
    outer_kn=3*abs(L[-1])/(4*PI*state.r[-1]**2*C*A*state.T[-1]**4)
    record('outer_face_knudsen',outer_kn,.1,outer_kn>.1)
    if bc['thermal'] in ('stars_photosphere','grey'):
        tau_half=.5*kap[-1]*state.dm[-1]/(4*PI*state.r[-1]**2)+external
        teff4=L[-1]/(4*PI*state.r[-1]**2*SIGMA)
        bath=bc['bath_temperature_K'] if bc['thermal']=='grey' else 0.
        eddington=(.75*teff4*(tau_half+2/3)+bath**4)/state.T[-1]**4-1
    else:eddington=None
    record('outer_face_eddington',None if eddington is None else abs(eddington),EDDINGTON_RELATION_TOLERANCE,
           eddington is not None and abs(eddington)>EDDINGTON_RELATION_TOLERANCE)
    max_conv_mach=float(max(convective_mach,default=0))
    record('convective_mach',max_conv_mach,1.,max_conv_mach>1)
    significant=abs(convective)>.01*(abs(radiative)+abs(convective))
    min_pe=float(min(peclet[significant])) if np.any(significant) else None
    record('convective_peclet',min_pe,100.,min_pe is not None and min_pe<100)
    gamma_ph=None;has_photosphere=bc['thermal']=='stars_photosphere' and state.photosphere_density is not None
    if has_photosphere:
        from .photosphere import effective_temperature,eddington_factor
        teff=effective_temperature(L[-1],state.r[-1])
        kph=float(np.atleast_1d(opacity(np.atleast_1d(state.photosphere_density),np.atleast_1d(teff)))[0])
        gamma_ph=float(eddington_factor(L[-1],state.enclosed[-1],kph))
    record('photosphere_gas_pressure',gamma_ph,1.,gamma_ph is not None and gamma_ph>=1)
    clamped_mass=clamped=photosphere_clamped=None
    if hasattr(opacity,'clamp_flags'):
        flags=np.asarray(opacity.clamp_flags(state.rho,state.T))
        clamped=int(np.sum(flags!=0));clamped_mass=float(np.sum(state.dm[flags!=0]))
        # kappa_ph enters the photosphere row (f): the point (rho_ph, Teff) is screened too.
        if has_photosphere:
            photosphere_clamped=bool(np.asarray(opacity.clamp_flags(np.atleast_1d(state.photosphere_density),
                                                                    np.atleast_1d(teff)))[0]!=0)
    record('opacity_clamp',clamped_mass,0.,bool(clamped) or bool(photosphere_clamped))
    inversion=float(max(np.diff(state.rho)/state.rho[:-1],default=0.))
    record('density_inversion',inversion,None,False)
    local=abs(radiative)*kap_face/(4*PI*C*G*state.enclosed[1:-1])
    record('local_eddington_factor',float(max(local,default=0.)),None,False)
    conservation=residuals(state,eos)
    tol=config['numerics']['conservation_tolerance']
    worst=max(abs(conservation['mass_residual_relative']),abs(conservation['energy_residual_relative']))
    screens['conservation']=dict(value=float(worst),limit=tol,violated=worst>tol,action='reject')
    if worst>tol:raise StepRejected('Conservation tolerance exceeded')
    diagnostics=dict(max_mach=mach,max_compactness=compactness,min_zone_optical_depth=float(min(tau)),
                     max_eos_first_law_defect=eos_defect,max_eos_dilute_gas_pressure_defect=gas_defect,
                     max_omitted_conductive_fraction=conduction,max_radiation_knudsen=float(max(kn)),
                     outer_radiation_knudsen=float(outer_kn),minimum_significant_convective_peclet=min_pe,
                     max_convective_mach=max_conv_mach,ideal_eos_helium_screen=helium,
                     max_volume_geometry_error=volume_error,max_temperature_jump=float(max(abs(np.diff(np.log(state.T))))),
                     photosphere_eddington_factor=gamma_ph,opacity_clamped_cells=clamped,
                     opacity_clamped_mass_g=clamped_mass,opacity_photosphere_clamped=photosphere_clamped,max_density_inversion=inversion,
                     max_local_eddington_factor=float(max(local,default=0.)),
                     policy=name,screens=screens,screen_warnings=warnings)
    return errors,diagnostics,conservation
