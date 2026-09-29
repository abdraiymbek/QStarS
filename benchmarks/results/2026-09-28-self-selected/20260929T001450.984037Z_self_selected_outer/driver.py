#!/usr/bin/env python3
"""Report 12: verify prerequisites; never promote them to a physical wind."""
import os,json,time,resource,traceback
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp,simpson,cumulative_trapezoid
from qstars import io
from qstars.constants import G,C,MSUN
from qstars.thermal_transfer import ThermalBands,opacity_means
from qstars.material_wind import adiabatic_eigenvalue_reference,GasDerivatives,continuum_coefficients
from qstars.outer_interface import InterfaceState,moving_flux,compare_faces
from qstars.spherical_wind import DiffusionWind,FreeEOSOpacity,match_base
from qstars.eos import Composition
from qstars.freeeos_reference import FreeEOSReference
from qstars.mesa_microphysics import MesaMicrophysics

plan=json.loads(Path('benchmarks/self_selected_plan.json').read_text())
limits=plan['thresholds'];out=io.create_directory('results','self_selected_outer');io.provenance(out,plan)
(out/'driver.py').write_bytes(Path(__file__).read_bytes());print(out,flush=True)
start=time.perf_counter();cases=[];checks=[];results={}
def check(name,value,limit):
 checks.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(np.isfinite(value) and value<=limit)))
def attempt(name,parameters,fn):
 tick=time.perf_counter();row=dict(name=name,parameters=parameters)
 try:
  result=fn();results[name]=result
  arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
  if arrays:np.savez_compressed(out/(name+'.npz'),**arrays)
  summary={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
  io.json_write(out/(name+'.json'),summary);row.update(status='completed_component',summary=summary)
 except Exception as exc:
  row.update(status='numerical_or_domain_failure',error=str(exc),traceback=traceback.format_exc(),
             history=getattr(exc,'history',getattr(exc,'trials',None)),physical_nonexistence_claimed=False)
 row['seconds']=time.perf_counter()-tick;cases.append(row);io.json_write(out/'cases.json',cases)
 print(name,row['status'],round(row['seconds'],3),flush=True)
 return results.get(name)

def thermal(n=401,outer=20.,core=32,tol=1e-7,maxiter=100):
 r=np.geomspace(1,outer,n);a=np.array([.2,2.])[:,None]/r[None,:]**2
 model=ThermalBands(r,a,a*.2,[0,20000,np.inf],20000,core_rays=core)
 q=model.radiative_equilibrium(tolerance=tol,max_iterations=maxiter)
 # Opacities below are per mass in the synthetic choice rho=r^-2, r0=1.
 rho=r**-2;means=opacity_means(a/rho,a*.2/rho,q['gas_temperature'],[0,20000,np.inf],q['J'],q['H'])
 q.update({'mean_'+k:v for k,v in means.items()});q['J_total']=q['J'].sum(axis=0);q['K_total']=q['K'].sum(axis=0)
 q['max_gas_radiation_temperature_ratio']=float(max(q['gas_temperature']/q['radiation_temperature']))
 q['max_flux_rosseland_ratio']=float(max(means['flux_total']/means['rosseland_total']))
 q['luminosity_profile_relative_defect']=float(max(abs(q['luminosity']/q['luminosity'][0]-1)))
 tau=np.array([cumulative_trapezoid(-chi[::-1],r[::-1],initial=0)[::-1] for chi in a*1.2]);q['band_extinction_optical_depth']=tau
 q['band_tau_two_thirds_radii']=[float(np.interp(2/3,t[::-1],r[::-1])) if t[0]>2/3 else None for t in tau]
 q['thermal_photosphere_certified']=False
 return q
for n in [101,201,401,801]:
 name=f'thermal_{n}';q=attempt(name,dict(n=n,outer=20,core=32),lambda n=n:thermal(n=n))
 if q:
  check(name+':thermal_balance',q['thermal_balance_error'],limits['thermal_balance_relative'])
  check(name+':energy_integral',q['energy_integral_relative_defect'],limits['independent_energy_relative'])
for name,opts in [('thermal_core64',dict(core=64)),('thermal_tight',dict(tol=1e-9)),
                  ('thermal_outer60',dict(n=601,outer=60)),('thermal_outer180',dict(n=801,outer=180)),
                  ('thermal_iteration_failure',dict(n=101,maxiter=1))]:
 attempt(name,opts,lambda opts=opts:thermal(**opts))
def compare(a,b,fields,rmin,rmax,limit,label):
 if a not in results or b not in results:
  check(label+':missing_prerequisite',1,0);return
 x=np.linspace(np.log(rmin),np.log(rmax),1001);q,p=results[a],results[b]
 for field in fields:
  av=PchipInterpolator(np.log(q['radius']),q[field])(x)
  bv=PchipInterpolator(np.log(p['radius']),p[field])(x)
  err=max(abs(av-bv)/np.maximum(abs(bv),1e-12*max(abs(bv))))
  check(label+':'+field,err,limit)
fields=['gas_temperature','J_total','radiation_flux','K_total','mean_acceleration']
for n in [101,201,401]:compare(f'thermal_{n}','thermal_801',fields,1,20,limits['radial_profile_relative'],f'radial_{n}')
compare('thermal_401','thermal_core64',fields,1,20,limits['angular_profile_relative'],'core_resolution')
compare('thermal_401','thermal_tight',fields,1,20,limits['iteration_profile_relative'],'thermal_tolerance')
compare('thermal_801','thermal_outer180',fields,1,10,limits['exterior_profile_relative'],'exterior20')
compare('thermal_outer60','thermal_outer180',fields,1,10,limits['exterior_profile_relative'],'exterior60')

# Independent frequency limit: kappa_nu proportional to nu, dilution W.
# Integral nu*B_nu is proportional to T^5: Tgas=Tstar*W^(1/5).
for bands in [16,32,64,128]:
 def run(bands=bands):
  r=np.geomspace(1,20,81);ts=20000.;e=np.r_[0,np.geomspace(ts*1e-4,ts*80,bands-1),np.inf]
  mid=np.r_[e[1]/2,np.sqrt(e[1:-2]*e[2:-1]),e[-2]*2]
  a=1e-10*(mid/ts)[:,None]/r[None,:]**2
  q=ThermalBands(r,a,a*.1,e,ts).radiative_equilibrium(tolerance=1e-9)
  W=.5*(1-np.sqrt(1-1/r**2));exact=ts*W**.2
  q.update(exact_temperature=exact,temperature_relative_error=float(max(abs(q['gas_temperature']/exact-1))))
  return q
 name=f'frequency_{bands}';q=attempt(name,dict(bands=bands),run)
 if q:check(name+':analytic_temperature',q['temperature_relative_error'],limits['thin_powerlaw_temperature_relative'])

# Thermodynamic base-selected wind, no mdot or base velocity input.
M=100*MSUN;rb=1e12;rhob=1e-10;gamma=4/3
base=dict(mass=M,radius=rb,density=rhob,pressure=rhob*G*M/(rb*2.5*gamma))
for name,opts in [('adiabatic',{}),('adiabatic_coarse',dict(samples=101,critical_offset=1e-4,rtol=1e-7)),
                  ('adiabatic_tight',dict(samples=1001,critical_offset=1e-7,rtol=1e-11,root_tolerance=1e-14)),
                  ('adiabatic_exterior',dict(samples=1001,outer_ratio=1000)),
                  ('adiabatic_cross',dict(gamma=1.4,pressure=rhob*G*M/(rb*2.3*1.4)))]:
 pars={**base,**opts};q=attempt(name,pars,lambda pars=pars:adiabatic_eigenvalue_reference(**pars))
 if q:
  check(name+':base_density',q['base_density_relative_error'],limits['adiabatic_profile_relative'])
  check(name+':base_pressure',q['base_pressure_relative_error'],limits['adiabatic_profile_relative'])
  check(name+':mass',max(abs(q['mass_flux']/q['mdot']-1)),limits['adiabatic_energy_relative'])
  check(name+':energy',max(abs(q['energy_flux']/q['energy_flux'][0]-1)),limits['adiabatic_energy_relative'])
  work=simpson(-G*M/q['radius']**2,x=q['radius']);h=q['internal_energy']+q['pressure']/q['density']
  defect=abs((q['velocity'][-1]**2-q['velocity'][0]**2)/2+h[-1]-h[0]-work)/(G*M/rb)
  check(name+':momentum_energy_integral',defect,limits['adiabatic_integral_relative'])
for name in ['adiabatic_coarse','adiabatic','adiabatic_exterior']:
 compare(name,'adiabatic_tight',['velocity','density','pressure'],rb,100*rb,limits['adiabatic_profile_relative'],name+'_profiles')
if 'adiabatic_tight' in results:
 q=results['adiabatic_tight']
 for index in [200,400,600]:
  name='adiabatic_rematch_'+str(index);pars=dict(mass=M,radius=q['radius'][index],density=q['density'][index],pressure=q['pressure'][index])
  p=attempt(name,pars,lambda pars=pars:adiabatic_eigenvalue_reference(**pars))
  if p:check(name+':mdot',abs(p['mdot']/q['mdot']-1),limits['adiabatic_profile_relative'])
for lam in [1.9,3.1,4.,10.]:
 pars={**base,'pressure':rhob*G*M/(rb*lam*gamma)}
 attempt('adiabatic_existence_'+str(lam),pars,lambda pars=pars:adiabatic_eigenvalue_reference(**pars))

# Native microphysics remains gas-only LTE; no radiation-population inference.
mesa=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
eos=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,mesa.radiation_constant);micro=FreeEOSOpacity(eos,mesa)
io.json_write(out/'microphysics.json',dict(mesa=mesa.manifest,freeeos=eos.manifest))
rows=[]
for rho in [3e-11,1e-8]:
 for T in [4000.,6000.,10000.,16000.,40000.,200000.]:
  raw=eos.raw(rho,T);gas=GasDerivatives.freeeos(raw)
  q=continuum_coefficients(gas,radius=1e13,density=rho,velocity=1e7,mass=100*MSUN,gas_heating_per_volume=0.,radiative_acceleration=0.)
  ratio=q['adiabatic_sound_speed_squared']/q['isothermal_sound_speed_squared']
  rows.append(dict(rho=rho,T=T,H_ionized=float(raw[18]),adiabatic_to_isothermal_speed_ratio=float(np.sqrt(ratio)),**q))
  check(f'native_gas_gamma1:{rho}:{T}',abs(q['adiabatic_sound_speed_squared']/(raw[14]*gas.pressure/rho)-1),3e-7)
io.json_write(out/'native_critical_states.json',rows)

w=DiffusionWind(micro,1e4*MSUN,1.152e15,3e-11,16000.)
for x in [-.003,-.01,-.03]:
 eps=1e-7;s=solve_ivp(w.rhs,[-eps,x],-eps*np.array([w.slope,w.sonic_temperature_slope]),method='DOP853',rtol=1e-11,atol=1e-13,max_step=.001)
 if not s.success:raise RuntimeError('Anchor integration failed')
 r,v,T,rho,g,L,*_=w.state(x,s.y[:,-1]);supplied=dict(mass=w.mass,radius=r,density=rho,temperature=T,luminosity=L+w.mdot*(v*v/2+g.h))
 for label,factors,budget in [('broad',[1.01,.7,.98],250),('near',[1.001,1.01,1.001],250),('exhausted',[1.01,.7,.98],1)]:
  def run(supplied=supplied,factors=factors,budget=budget):
   p,report=match_base(micro,**supplied,sonic_guess=np.array([w.rs,w.rhos,w.ts])*factors,rtol=1e-9,match_tolerance=1e-8,max_nfev=budget)
   report.update(mdot=p.mdot,mdot_relative_error=float(abs(p.mdot/w.mdot-1)),base_is_wind_slice=True,physical_envelope_test=False)
   return report
  name='diffusion_match_'+str(x)+'_'+label;q=attempt(name,dict(base=supplied,factors=factors,budget=budget),run)
  if q:
   check(name+':residual',q['max_log_residual'],limits['diffusion_match_log_residual'])
   check(name+':mass_flux',q['mdot_relative_error'],limits['diffusion_recovered_mass_flux_relative'])

# No accepted stellar envelope has been attached. Explicit staged decision.
io.json_write(out/'checks.json',checks)
summary=dict(outcome='FAIL',checks=len(checks),passed=sum(c['passed'] for c in checks),
             failed_checks=[c for c in checks if not c['passed']],cases=len(cases),
             completed_components=sum(c['status']=='completed_component' for c in cases),
             elapsed_seconds=time.perf_counter()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
             source_snapshot_verified=io.finish_provenance(out),
             physical_envelope_test=False,coupled_evolution=False,published_track=False,qsrd=False,
             remaining=['validated H/He frequency-dependent population/energy/force coupling','convective handoff',
                        'combined escaping thermodynamic eigenvalue and physical envelope matching','global branch search'])
io.json_write(out/'campaign.json',summary);print(json.dumps(summary),flush=True)
