#!/usr/bin/env python3
"""Resolve spatial/output defects without changing any physical coefficients."""
import os,json,time,resource,traceback
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from qstars import io
from qstars.constants import G,C,MSUN
from qstars.thermal_transfer import ThermalBands
from qstars.material_wind import adiabatic_eigenvalue_reference
from qstars.spherical_wind import DiffusionWind,FreeEOSOpacity,match_base
from qstars.eos import Composition
from qstars.freeeos_reference import FreeEOSReference
from qstars.mesa_microphysics import MesaMicrophysics
plan=json.loads(Path('benchmarks/self_selected_refinement_plan.json').read_text());out=io.create_directory('results','self_selected_refinement')
io.provenance(out,plan);(out/'driver.py').write_bytes(Path(__file__).read_bytes());print(out,flush=True)
start=time.perf_counter();rows=[];checks=[];results={}
def check(name,value,limit):checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
def attempt(name,fn):
 tick=time.perf_counter();row=dict(name=name)
 try:
  q=fn();results[name]=q
  np.savez_compressed(out/(name+'.npz'),**{k:v for k,v in q.items() if isinstance(v,np.ndarray)})
  summary={k:v for k,v in q.items() if not isinstance(v,np.ndarray)};io.json_write(out/(name+'.json'),summary)
  row.update(status='completed_component',summary=summary)
 except Exception as exc:row.update(status='numerical_or_domain_failure',reason=str(exc),traceback=traceback.format_exc(),history=getattr(exc,'history',getattr(exc,'trials',None)))
 row['seconds']=time.perf_counter()-tick;rows.append(row);io.json_write(out/'cases.json',rows);print(name,row['status'],round(row['seconds'],3),flush=True)
 return results.get(name)
for pars in plan['thermal_cases']:
 def run(pars=pars):
  r=np.exp(np.log(pars['outer'])*np.linspace(0,1,pars['n'])**2);a=np.array([.2,2.])[:,None]/r[None,:]**2
  q=ThermalBands(r,a,a*.2,[0,20000,np.inf],20000).radiative_equilibrium()
  q['radiative_acceleration']=4*np.pi/C*(a*1.2*q['H']).sum(axis=0)*r**2
  return q
 q=attempt(pars['name'],run)
 if q:check(pars['name']+':energy',q['energy_integral_relative_defect'],.0005)
def compare(a,b,fields,lo,hi,limit):
 if a not in results or b not in results:check(a+'_vs_'+b+':missing',1,0);return
 r=np.geomspace(lo,hi,2001)
 for key in fields:
  q,p=results[a],results[b];av=PchipInterpolator(np.log(q['radius']),q[key])(np.log(r));bv=PchipInterpolator(np.log(p['radius']),p[key])(np.log(r))
  check(a+'_vs_'+b+':'+key,max(abs(av/bv-1)),limit)
fields=['gas_temperature','radiation_energy','radiation_flux','radiation_pressure_rr','radiative_acceleration']
for name in ['cluster401','cluster801']:compare(name,'cluster1601',fields,1,20,.002)
compare('cluster_outer60','cluster_outer180',fields,1,10,.002)
M=100*MSUN;r=1e12;rho=1e-10;gamma=4/3;base=dict(mass=M,radius=r,density=rho,pressure=rho*G*M/(r*2.5*gamma))
for pars in plan['adiabatic_cases']:
 q=attempt(pars['name'],lambda pars=pars:adiabatic_eigenvalue_reference(**base,samples=pars['samples'],outer_ratio=pars['outer_ratio'],rtol=1e-11,critical_offset=1e-7,root_tolerance=1e-14))
 if q:check(pars['name']+':energy',max(abs(q['energy_flux']/q['energy_flux'][0]-1)),1e-7)
for name in ['dense2001','dense_exterior']:compare(name,'dense4001',['velocity','density','pressure'],r,100*r,1e-6)
def match():
 mesa=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
 eos=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,mesa.radiation_constant);micro=FreeEOSOpacity(eos,mesa)
 io.json_write(out/'microphysics.json',dict(mesa=mesa.manifest,freeeos=eos.manifest))
 w=DiffusionWind(micro,1e4*MSUN,1.152e15,3e-11,16000.);eps=1e-7;x=-.03
 s=solve_ivp(w.rhs,[-eps,x],-eps*np.array([w.slope,w.sonic_temperature_slope]),method='DOP853',rtol=1e-10,atol=1e-12,max_step=.001)
 if not s.success:raise RuntimeError('Reference integration failed')
 r,v,T,rho,g,L,*_=w.state(x,s.y[:,-1]);base=dict(mass=w.mass,radius=r,density=rho,temperature=T,luminosity=L+w.mdot*(v*v/2+g.h))
 p,report=match_base(micro,**base,sonic_guess=np.array([w.rs,w.rhos,w.ts])*[1.1,.1,.8],initial_trust_radius=2.,rtol=1e-9,match_tolerance=1e-8,max_nfev=250)
 report.update(mdot_relative_error=float(abs(p.mdot/w.mdot-1)),base_is_wind_slice=True,physical_envelope_test=False)
 return report
q=attempt('rejected_trial_recovery',match)
if q:
 check('recovered_mdot',q['mdot_relative_error'],1e-6)
 check('actual_inadmissible_trial_exercised',0 if q['rejected_trials']>0 else 1,0)
io.json_write(out/'checks.json',checks)
summary=dict(outcome='FAIL',checks=len(checks),passed=sum(x['passed'] for x in checks),failed_checks=[x for x in checks if not x['passed']],
             elapsed_seconds=time.perf_counter()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
             source_snapshot_verified=io.finish_provenance(out),physical_envelope_test=False,coupled_evolution=False)
io.json_write(out/'campaign.json',summary);print(json.dumps(summary),flush=True)
