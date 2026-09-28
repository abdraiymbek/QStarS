#!/usr/bin/env python3
"""Independent transonic components, published cases and retained failures.

This is a validation campaign, not a coupled stellar parameter sweep.
"""
import argparse,csv,hashlib,json,os,resource,time,traceback
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.interpolate import PchipInterpolator
from qstars import io
from qstars.constants import G,C,MSUN,RSUN,YEAR
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.freeeos_reference import FreeEOSReference
from qstars.spherical_wind import (FreeEOSOpacity,IdealGasOpacity,DiffusionWind,
    WindMatchError,match_base,regular_sonic_slope,spherical_hydrostatic)

p=argparse.ArgumentParser();p.add_argument('--output',default='results')
p.add_argument('--report8-profiles',type=Path,default=Path('results/20260928T044002.683023Z_atmosphere_reference'))
a=p.parse_args();start=time.perf_counter()
out=io.create_directory(a.output,'spherical_wind_validation')
io.provenance(out,dict(family='spherical_wind_components',report8_profiles=str(a.report8_profiles),
    numerical_profile_log_tolerance=.001,independent_integral_tolerance=1e-6,
    paper_comparison='unaltered printed sonic inputs; vector stroke uncertainties'))
print(out,flush=True)
(out/'driver.py').write_bytes(Path(__file__).read_bytes())
mesa=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
eos=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,mesa.radiation_constant)
micro=FreeEOSOpacity(eos,mesa)
io.json_write(out/'microphysics.json',dict(mesa=mesa.manifest,freeeos=eos.manifest,
    hydrogen=.7,helium=.3,radiation_constant=micro.arad,
    radiation='LTE Rosseland diffusion; no thin/non-LTE transfer claim',
    freeeos_energy_zero='Unmodified upstream EOS1; boundary energy uses the same gauge'))
rows=[];profiles={};wind_objects={}

def run(label,args,options=None,derivative_step=1e-4):
    t=time.perf_counter();options={} if options is None else options
    try:
        w=DiffusionWind(micro,*args,derivative_step=derivative_step)
        f,s=w.integrate(**options);np.savez_compressed(out/(label+'.npz'),**f)
        profiles[label]=f;wind_objects[label]=w
        s['mass_flux_msun_year']=w.mdot*YEAR/MSUN
        s['sonic_flux']=w.boundary_flux(0.,[0.,0.])
        s['status']='completed_component' if s['solver_success'] else 'numerical_failure'
    except Exception as e:
        s=dict(status='numerical_or_domain_failure',reason=str(e),traceback=traceback.format_exc(),
               physical_impossibility_claimed=False)
    s.update(label=label,inputs=list(args),options=options,seconds=time.perf_counter()-t)
    rows.append(s);io.json_write(out/'cases.json',rows)
    print(label,s['status'],s.get('terminations',s.get('reason')),flush=True)
    return s

# All four printed cases, not a scan adjusted to reproduce the red curve.
for rho in [9.65e-9,9.71e-9,9.85e-9,10.5e-9]:
    run('paper_'+str(rho),(100*MSUN,50*RSUN,rho,2.5e5),
        dict(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,samples=2001))

# Two prescribed sonic states in the problematic gravity/recombination regime.
# They are component controls, not evolved quasi-star solutions.
for tag,T in [('recombining',10000.),('diffusion_exit',16000.)]:
    args=(1e4*MSUN,1.152e15,3e-11,T)
    options=dict(inner_ratio=.995,outer_ratio=2.,endpoint='diffusion_limit',rtol=1e-9,
                 sonic_offset=1e-7,max_log_step=.002,samples=1001)
    run(tag,args,options)
    changes=[('loose',{'rtol':1e-6}),('tight',{'rtol':1e-10}),
             ('radial_step',{'max_log_step':.0005}),('sparse_output',{'samples':251}),
             ('sonic_offset',{'sonic_offset':1e-6}),('outer_radius',{'outer_ratio':4.}),
             ('independent_integrator',{'method':'Radau'})]
    for name,change in changes:run(tag+'_'+name,args,dict(options,**change))
    run(tag+'_critical_derivative',args,options,derivative_step=1e-5)

comparisons=[]
for label,f in profiles.items():
    tag=next((s for s in ['recombining','diffusion_exit'] if label.startswith(s+'_')),None)
    if tag is None:continue
    base=profiles[tag]
    x=np.linspace(max(np.log(f['radius'][0]),np.log(base['radius'][0])),
                  min(np.log(f['radius'][-1]),np.log(base['radius'][-1])),1001)
    errors={}
    for key in ['density','temperature','velocity','opacity','radiative_luminosity',
                'gas_pressure','radiation_enthalpy_luminosity','hydrogen_ionized']:
        u=PchipInterpolator(np.log(f['radius']),np.log(f[key]))(x)
        v=PchipInterpolator(np.log(base['radius']),np.log(base[key]))(x)
        errors[key]=float(max(abs(u-v)))
    comparisons.append(dict(case=label,max_profile_log_errors=errors,
        endpoint_radius_relative_error=float(abs(f['radius'][-1]/base['radius'][-1]-1)),
        all_profiles_below_0p1_percent=bool(max(errors.values())<.001)))
io.json_write(out/'convergence.json',comparisons)

# Base matching is tested against independently supplied slices of the same
# exact component solution; this verifies matching numerics, not a star model.
matching=[];w=wind_objects['diffusion_exit']
for x in [-.005,-.01,-.02]:
    eps=1e-8
    sol=solve_ivp(w.rhs,[-eps,x],-eps*np.array([w.slope,w.sonic_temperature_slope]),
                  method='DOP853',rtol=1e-10,atol=1e-12,max_step=.001)
    r,v,T,rho,q,L,*_=w.state(x,sol.y[:,-1])
    supplied=dict(mass=w.mass,radius=r,density=rho,temperature=T,luminosity=L+w.mdot*(v*v/2+q.h))
    try:
        recovered,s=match_base(micro,**supplied,sonic_guess=[w.rs*1.0001,w.rhos*1.001,w.ts*1.0001],rtol=1e-9)
        s.update(mass_flux_relative_error=recovered.mdot/w.mdot-1,sonic_radius_relative_error=recovered.rs/w.rs-1,
                 sonic_temperature_relative_error=recovered.ts/w.ts-1,status='matched_component')
    except Exception as e:s=dict(status='matching_failure',reason=str(e),trials=getattr(e,'trials',[]))
    s.update(log_matching_radius=x,supplied=supplied,base_gas_mach=v/np.sqrt(q.ct2),
             base_flux=w.boundary_flux(x,sol.y[:,-1]))
    matching.append(s);io.json_write(out/'matching.json',matching)
    print('match',x,s['status'],flush=True)

# Do not label failure of a root finder as physical failure of the envelope.
report8=[]
for teff,g in [(6000,1),(8000,3)]:
    path=a.report8_profiles/f'T{teff}_g{float(g)}_baseline.npz'
    f=np.load(path);rp=np.sqrt(G*1e4*MSUN/g);ip=np.argmin(abs(f['tau']-2/3))
    rb=rp-(f['inward_depth'][-1]-f['inward_depth'][ip]);rho=f['density'][-1];T=f['temperature'][-1]
    L=4*np.pi*rp**2*mesa.radiation_constant*C/4*teff**4
    inputs=dict(mass=1e4*MSUN,radius=rb,density=rho,temperature=T,luminosity=L)
    # These are boundary experiments using a failed hydrostatic atmosphere's
    # base data; they do not certify a physical stellar boundary or energy gauge.
    guess=[rb*1.003,rho*.1,12000.]
    try:
        w,s=match_base(micro,**inputs,sonic_guess=guess,rtol=1e-8,max_nfev=250)
        ff,ws=w.integrate(inner_ratio=rb/w.rs,outer_ratio=4,rtol=1e-9,sonic_offset=1e-7)
        np.savez_compressed(out/f'report8_T{teff}.npz',**ff)
        s.update(status='matched_local_component',wind=ws)
    except Exception as e:s=dict(status='matching_failure',reason=str(e),trials=getattr(e,'trials',[]))
    s.update(inputs=inputs,guess=guess,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
             no_physical_no_root_inference=True,energy_gauge_coupling_unvalidated=True)
    report8.append(s);io.json_write(out/'report8_base_attempts.json',report8)
    print('Report8',teff,s['status'],s.get('reason',''),flush=True)

# Both horizontal and vertical vector-stroke uncertainties propagated.
data=list(csv.DictReader(open('benchmarks/nakauchi2017_figure2.csv')))
paper=[]
for d in data:
    label='paper_'+str(float(d['sonic_density']));f=profiles.get(label)
    row=dict(d)
    if f is None:row['status']='no_model'
    else:
        x=np.log(float(d['radius_rsun'])*RSUN);dx=float(d['log_r_uncertainty']);y=np.log(float(d['value']));dy=float(d['log_y_uncertainty'])
        axis=np.log(f['radius']);lo=max(x-dx,axis[0]);hi=min(x+dx,axis[-1]);key=d['field']
        if lo>hi:row['status']='outside_model_coverage'
        else:
            pred=PchipInterpolator(axis,np.log(f[key]))(np.linspace(lo,hi,9))
            distance=max(0.,min(pred)-y,y-max(pred))
            row.update(status='compared',log_distance_outside_radius_band=float(distance),passed=bool(distance<=dy))
    paper.append(row)
io.json_write(out/'published_profile_comparison.json',paper)
counts={s:sum(d['status']==s for d in paper) for s in ['compared','outside_model_coverage','no_model']}
counts['passed']=sum(d.get('passed',False) for d in paper)
# Scalar targets deliberately retain approximate wording, not extra precision.
red=next((x for x in rows if x['label']=='paper_9.71e-09'),{})
io.json_write(out/'published_scalar_comparison.json',dict(
    paper='Nakauchi2017 Figure2 and Section3; different EOS/opacity realization',
    mass_loss_msun_year=dict(published_approximate=.14,computed=red.get('mass_flux_msun_year')),
    photosphere_rsun=dict(published_approximate=1400,computed=red.get('endpoint_radius',0)/RSUN),
    velocity_at_photosphere_cm_s=dict(published_approximate=2e7,computed=float(profiles['paper_9.71e-09']['velocity'][-1])),
    exact_replication_passed=False,no_retuning=True))

# Analytical spherical hydrostatic and isothermal sonic checks.
ideal=IdealGasOpacity();M=1e4*MSUN;r=1e14;gamma=.95
T=G*M*(1-gamma)/(4*ideal.rgas*r);rho=ideal.arad/(3*ideal.rgas)*(1-gamma)/gamma*T**3;L=4*np.pi*C*G*M*gamma/ideal.kappa
f=spherical_hydrostatic(ideal,mass=M,radius=r,density=rho,temperature=T,luminosity=L,outer_ratio=2)
analytic=dict(hydrostatic_temperature_max_relative=float(max(abs(f['temperature']/(T*r/f['radius'])-1))),
              hydrostatic_density_max_relative=float(max(abs(f['density']/(rho*(r/f['radius'])**3)-1))))
np.savez_compressed(out/'hydrostatic.npz',**f)
slope,_=regular_sonic_slope(lambda y:np.array([2-2*np.exp(-y[0]),np.expm1(2*y[1])]),0,1e-5)
analytic['parker_critical_slope_error']=abs(slope-1)
io.json_write(out/'analytical.json',analytic)
summary=dict(cases=len(rows),failures=sum(x['status']!='completed_component' for x in rows),
    published_profile_counts=counts,published_replication_passed=False,
    convergence_all_passed=all(c['all_profiles_below_0p1_percent'] for c in comparisons),
    matched_component_slices=sum(x['status']=='matched_component' for x in matching),
    report8_base_matches=sum(x['status']=='matched_local_component' for x in report8),
    runtime_seconds=time.perf_counter()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
    source_snapshot_verified=io.finish_provenance(out),
    outer_closure_pass=False,campaign_complete=False,v1_validated=False)
io.json_write(out/'validation.json',summary);print(json.dumps(summary,indent=2),flush=True)
