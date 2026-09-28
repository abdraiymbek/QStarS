#!/usr/bin/env python3
"""Targeted resolution of failed boundary/asymptotic placements, no retuning."""
import json,time,resource,traceback
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from qstars import io
from qstars.constants import MSUN
from qstars.grey_transfer import solve_grey
from qstars.continuum_wind import solve_continuum_wind,dimensional_fields,flux_driven_reference
plan=json.loads(Path('benchmarks/escaping_followup_plan.json').read_text())
out=io.create_directory('results','escaping_extent');io.provenance(out,plan)
(out/'driver.py').write_bytes(Path(__file__).read_bytes());print(out,flush=True)
start=time.perf_counter();cases=[];checks=[]
def check(name,value,limit):checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
for pars in plan['outer_wind_domains']:
 name='outer_'+str(pars['outer']);tick=time.perf_counter();row=dict(name=name,parameters=pars)
 try:
  r=np.unique(np.r_[1,np.geomspace(1.00001,pars['outer'],pars['n']),1/(1-np.geomspace(1e-8,.2,100))])
  q=solve_continuum_wind(r,gamma=2,tiring=.001,escape_speed_c=.0001)
  fields=dimensional_fields(q,100*MSUN,.34)
  d={k:v for k,v in q.items() if isinstance(v,np.ndarray)};d['r']=r;d['velocity']=np.sqrt(q['w']);d['density']=q['extinction'];d['force']=q['H_comoving']
  np.savez_compressed(out/f'{name}.npz',**d)
  mass_integral=PchipInterpolator(fields['radius_cm'],4*np.pi*fields['radius_cm']**2*fields['density_comoving']/np.sqrt(1-q['beta']**2)).integrate(fields['radius_cm'][0],fields['radius_cm'][-1])
  w,l=flux_driven_reference(r,2,.001,1e-4)
  summary=dict(history=q['history_wind'],photosphere=q['photospheric_radius_ratio'],effective_temperature=fields['effective_temperature_K'],emergent_luminosity_ratio=float(q['luminosity_ratio'][-1]),
   emergent_angular_ratio=float(r[-1]**2*q['angular_H'][-1]),flux_factor_at_outer=float(q['H'][-1]/q['J'][-1]),eddington_at_outer=float(q['K'][-1]/q['J'][-1]),
   photospheric_tau_tail_estimate=float(.002/(q['beta'][-1]*r[-1])),
   total_energy_defect=float(max(abs(q['total_energy_ratio']/q['total_energy_ratio'][0]-1))),
   velocity_reference_error=float(max(abs(np.sqrt(q['w']/w)-1))),force_reference_error=float(max(abs(q['H_comoving']/(l/r**2)-1))),
   dimensional_wind_mass_over_central=float(mass_integral/(100*MSUN)),
   dimensional_scaling_is_physical_star=False)
  io.json_write(out/f'{name}.json',summary);row.update(status='completed_component',summary=summary)
  check(name+':energy',summary['total_energy_defect'],.001)
  check(name+':angular_flux',q['history_wind'][-1]['radiation']['angular_flux_max_error'],.005)
  check(name+':streaming',1-summary['flux_factor_at_outer'],.01)
  check(name+':velocity_reference',summary['velocity_reference_error'],.01)
  check(name+':force_reference',summary['force_reference_error'],.01)
 except Exception as exc:row.update(status='numerical_or_domain_failure',error=str(exc),traceback=traceback.format_exc(),history=getattr(exc,'history',None))
 row['seconds']=time.perf_counter()-tick;cases.append(row);io.json_write(out/'cases.json',cases);print(name,row['status'],round(row['seconds'],3),flush=True)

for n in [10000,30000]:
 if not (out/f'outer_{n}.npz').exists() or not (out/'outer_100000.npz').exists():
  check(f'outer_{n}:missing_comparison',1.,0.);continue
 a,b=dict(np.load(out/f'outer_{n}.npz')),dict(np.load(out/'outer_100000.npz'));r=np.geomspace(1,1000,1001)
 for f in ['velocity','density','J','H','K','force']:
  ia=PchipInterpolator(np.log(a['r']),np.log(a[f]));ib=PchipInterpolator(np.log(b['r']),np.log(b[f]));err=max(abs(np.exp(ia(np.log(r))-ib(np.log(r)))-1));check(f'outer_{n}:{f}',err,.002)
 sa=json.loads((out/f'outer_{n}.json').read_text());sb=json.loads((out/'outer_100000.json').read_text())
 for f in ['photosphere','effective_temperature','emergent_luminosity_ratio']:
  check(f'outer_{n}:{f}',abs(sa[f]/sb[f]-1),.002)

# The old tau>=30 counterexample remains unchanged in the previous campaign.
try:
 p=1.5;outer=10;inv=lambda tau: ((p-1)*tau+outer**(1-p))**(-1/(p-1));r=np.geomspace(inv(300),10,1601)
 q=solve_grey(r,r**-p,16*np.pi**2,core_rays=32)
 tau=(r**(1-p)-outer**(1-p))/(p-1);mask=(tau>=100)&(r>1.1*r[0]);asym=3/(p+1)*r**(-p-1)
 np.savez_compressed(out/'deep_power1p5.npz',r=r,J=q.mean_intensity,f=q.eddington_factor,tau=tau,asymptotic_J=asym,comparison_mask=mask)
 check('deep_power1p5:J',max(abs(q.mean_intensity[mask]/asym[mask]-1)),.005)
 check('deep_power1p5:f',max(abs(q.eddington_factor[mask]-1/3)),.001)
 io.json_write(out/'deep_power1p5.json',dict(history=q.history,largest_gradient_knudsen_in_comparison=float(max((p+1)*r[mask]**(p-1)))))
except Exception as exc:check('deep_power1p5:failed',1,0);io.json_write(out/'deep_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
io.json_write(out/'checks.json',checks)
summary=dict(outcome='OUTER CLOSURE FAIL',cases=cases,checks=len(checks),failed_checks=[x for x in checks if not x['passed']],elapsed_seconds=time.perf_counter()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,source_snapshot_verified=io.finish_provenance(out),coupled_evolution_run=False)
io.json_write(out/'campaign.json',summary);print(json.dumps({k:v for k,v in summary.items() if k!='cases'}),flush=True)
