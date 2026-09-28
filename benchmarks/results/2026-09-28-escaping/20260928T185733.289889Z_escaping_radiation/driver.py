#!/usr/bin/env python3
"""Report 11 independent transfer and restricted stationary control campaign.

Every attempted case is retained. A mathematical grey control cannot pass a
missing thermodynamic wind/physical matching gate.
"""
import argparse,json,resource,time,traceback
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from qstars import io
from qstars.constants import MSUN
from qstars.grey_transfer import solve_grey,solve_moving_grey,solve_scattering,SphericalRays,comoving_moments
from qstars.continuum_wind import solve_continuum_wind,flux_driven_reference,dimensional_fields

p=argparse.ArgumentParser();p.add_argument('--output',default='results');a=p.parse_args()
out=io.create_directory(a.output,'escaping_radiation');start=time.perf_counter()
limits=json.loads(Path('benchmarks/escaping_transport_tolerances.json').read_text())
io.provenance(out,dict(scope='grey radiation and constructed continuum control; no stellar coupling',grid_strategy='log radius plus optical-depth nodes; merge isolated gaps below 0.1 of both neighbors',tolerances=limits))
(out/'driver.py').write_bytes(Path(__file__).read_bytes())
for name in ['escaping_transport_tolerances.json','outer_candidates.json','hummer1971_table1.csv','grey_continuum_wind.json']:
 (out/name).write_bytes((Path('benchmarks')/name).read_bytes())
print(out,flush=True);cases=[];checks=[]
def check(name,value,limit,scope,expected_pass=True):
 row=dict(name=name,value=float(value),limit=float(limit),passed=bool(value<=limit),scope=scope,expected_pass=expected_pass)
 checks.append(row);io.json_write(out/'checks.json',checks);return row

def save(name,params,fn):
 t=time.perf_counter();row=dict(name=name,parameters=params)
 try:
  data,summary=fn();np.savez_compressed(out/f'{name}.npz',**data)
  io.json_write(out/f'{name}.json',summary);row.update(status='completed_component',summary=summary)
 except Exception as exc:
  row.update(status='numerical_or_domain_failure',error=str(exc),traceback=traceback.format_exc(),history=getattr(exc,'history',None),physical_nonexistence_claimed=False)
 row['seconds']=time.perf_counter()-t;cases.append(row);io.json_write(out/'cases.json',cases)
 print(name,row['status'],round(row['seconds'],3),flush=True)
 return row

def grid(n,outer=10,tau=100,power=2):
 inv=lambda tau: ((power-1)*tau+outer**(1-power))**(-1/(power-1))
 r=np.unique(np.r_[np.geomspace(inv(tau),outer,n),inv(np.geomspace(1e-6,min(1.,tau),40))])
 # Remove isolated nearly coincident knots from merging the two meshes; otherwise
 # tiny nodal transfer errors acquire large interpolation slopes. No physical law changes.
 while len(r)>3:
  d=np.diff(np.log(r));bad=np.flatnonzero(d[1:-1] < .1*np.minimum(d[:-2],d[2:]))
  if not len(bad):break
  r=np.delete(r,bad+2)
 return r
def static(n=801,outer=10,tau=100,power=2,tolerance=1e-7,moment_tolerance=1e-11,core_rays=32,max_iterations=50):
 r=grid(n,outer,tau,power)
 q=solve_grey(r,r**-power,16*np.pi**2,tolerance=tolerance,moment_tolerance=moment_tolerance,core_rays=core_rays,max_iterations=max_iterations)
 fields=q.fields(grey_opacity=.34)
 data=dict(r=r,J=q.mean_intensity,angular_J=q.angular_mean_intensity,f=q.eddington_factor,angular_H=q.angular_flux/(4*np.pi),tau=q.optical_depth_to_boundary,flux_factor=fields['flux_factor'])
 summary=dict(history=q.history,photospheric_radius=fields['photospheric_radius'],effective_temperature=fields['effective_temperature'],emergent_luminosity=fields['emergent_luminosity'],thermal_temperature_inferred=False)
 return data,summary
for n in [201,401,801,1601]:save(f'grey_n{n}',dict(n=n),lambda n=n:static(n=n))
for tau in [30,300]:save(f'grey_tau{tau}',dict(tau=tau),lambda tau=tau:static(tau=tau))
for outer in [30,100,300,1000]:save(f'grey_outer{outer}',dict(outer=outer),lambda outer=outer:static(outer=outer,n=1001))
for tol in [1e-5,1e-9]:save(f'grey_tol{tol}',dict(tolerance=tol),lambda tol=tol:static(tolerance=tol,moment_tolerance=1e-13 if tol==1e-9 else 1e-11))
for power in [1.5,3]:save(f'grey_power{power}',dict(power=power),lambda power=power:static(n=1201,power=power))
# Exact retained failed attempt with current frozen source.
save('grey_insufficient_iterations',dict(max_iterations=1,expected='failure'),lambda:static(n=201,max_iterations=1))

def load(name):return dict(np.load(out/f'{name}.npz')) if (out/f'{name}.npz').exists() else None
def interp(d,field,r,positive=True):
 y=np.log(d[field]) if positive else d[field]
 z=PchipInterpolator(np.log(d['r']),y)(np.log(r))
 return np.exp(z) if positive else z

def compare(name,other,lo,hi,fields,limit,label):
 d,e=load(name),load(other)
 if d is None or e is None:
  check(label+':missing_solution',1,0,'failed prerequisite for profile comparison');return None
 r=np.geomspace(lo,hi,1001);errors={} 
 for field in fields:
  z=interp(d,field,r);ref=interp(e,field,r)
  errors[field]=float(max(abs(z/ref-1)));check(f'{label}:{field}',errors[field],limit,label)
 return errors

table=np.genfromtxt('benchmarks/hummer1971_table1.csv',delimiter=',',names=True)
d=load('grey_n1601');tab=[]
if d is None:check('published_table_missing_solution',1,0,'failed published reference prerequisite')
for r,j,f in (zip(table['r'],table['J'],table['f']) if d is not None else []):
 q=dict(r=float(r),reference_J=float(j),reference_f=float(f),J=float(interp(d,'J',r)),f=float(interp(d,'f',r)))
 for field in ['J','f']:
  error=abs(q[field]/q['reference_'+field]-1);q[field+'_relative_error']=error
  check(f'published_table:{r}:{field}',error,.01 if r==10 else .005,'Hummer1971 Table I')
 tab.append(q)
io.json_write(out/'published_table_comparison.json',tab)
for n in [201,401,801]:compare(f'grey_n{n}','grey_n1601',.04,10,['J','f'],limits['grey_profile_mesh_relative'],f'grey_mesh_{n}')
for tau in [30,300]:compare(f'grey_tau{tau}','grey_n801',.04,10,['J','f'],.002,f'grey_inner_tau{tau}')
compare('grey_outer300','grey_outer1000',.04,3,['J','f'],.002,'grey_outer_placement')
for tol in [1e-5,1e-9]:compare(f'grey_tol{tol}','grey_n801',.04,10,['J','f'],1e-5,f'grey_closure_tol{tol}')
for name,power in [('grey_n1601',2),('grey_power1.5',1.5),('grey_power3',3)]:
 d=load(name)
 if d is None:
  check(name+':missing_solution',1,0,'failed diffusion prerequisite');continue
 mask=(d['tau']>=30)&(d['r']>=1.1*d['r'][0]);asym=3/(power+1)*d['r']**(-power-1)
 check(name+':diffusion_J',max(abs(d['J'][mask]/asym[mask]-1)),.005,'diffusion asymptotic')
 check(name+':diffusion_f',max(abs(d['f'][mask]-1/3)),.001,'diffusion asymptotic')
 check(name+':angular_luminosity',max(abs(d['r']**2*d['angular_H']-1)),.005,'independent angular conservation')
d=load('grey_outer1000');check('far_field_flux_factor',1-d['flux_factor'][-1] if d is not None else 1.,.01,'free streaming')
# Independent exact finite-disc vacuum solution, including all three moments.
r=np.geomspace(1,1000,151);rays=SphericalRays(r,1e-12/r**2)
j,h,k=rays.moments(np.zeros_like(r),np.ones(rays.core_rays));mu=np.sqrt(1-1/r**2);je=1/(2*r*r*(1+mu));he=.25/r**2;ke=je*(1+mu+mu*mu)/3
np.savez_compressed(out/'transparent_sphere.npz',r=r,J=j,H=h,K=k,J_exact=je,H_exact=he,K_exact=ke)
for field,actual,expected in [('J',j,je),('H',h,he),('K',k,ke)]:check('transparent:'+field,max(abs(actual/expected-1)),1e-6,'exact angular moments')

# Moving radiation with work on prescribed matter; not a wind solution.
def moving(n=401,tau=3,beta=.01,tolerance=1e-6,core_rays=32,direct_source=False):
 r=np.unique(np.r_[np.geomspace(1,100,n),1+np.geomspace(1e-7,.05,60)])
 chi=tau/r**2;b=np.full_like(r,beta)
 q=solve_scattering(r,chi,b,1.,core_rays=core_rays) if direct_source else solve_moving_grey(r,chi,b,tolerance=tolerance,core_rays=core_rays)
 data={key:val for key,val in q.items() if isinstance(val,np.ndarray)};data['r']=r
 ell=r*r*q['H'];source=r*r*chi*b/np.sqrt(1-b*b)*q['H_comoving']
 work=PchipInterpolator(r,source).integrate(r[0],r[-1])
 energy_error=abs(ell[-1]-ell[0]+work)/ell[0]
 summary={key:val for key,val in q.items() if not isinstance(val,np.ndarray)}
 summary.update(emergent_to_net_base_luminosity=float(ell[-1]/ell[0]),independent_work_defect=float(energy_error),luminosity_span=float(np.ptp(ell)/ell[0]))
 if not direct_source:
  _,ah0,_=comoving_moments(q['angular_J'],q['angular_H'],q['angular_K'],b)
  integral=PchipInterpolator(r,r*r*chi*b/np.sqrt(1-b*b)*ah0).integrate(r[0],r[-1])
  summary['angular_work_defect']=float(abs(ell[-1]-ell[0]+integral)/ell[0])
 return data,summary
for n in [201,401,801]:save(f'moving_n{n}',dict(n=n,tau=3,beta=.01),lambda n=n:moving(n=n))
for tau,beta in [(30,0),(30,.001),(30,.01),(.3,.01)]:
 name=f'moving_tau{tau}_b{beta}';row=save(name,dict(n=801,tau=tau,beta=beta),lambda tau=tau,beta=beta:moving(n=801,tau=tau,beta=beta))
 if row['status']=='completed_component':
  check(name+':angular_work',row['summary']['angular_work_defect'],.001,'moving radiation conservation')
  check(name+':angular_flux',row['summary']['history'][-1]['angular_flux_max_error'],.005,'moving angular closure')
compare('moving_n401','moving_n801',1,100,['J','H','K'],.002,'moving_mesh')
row=save('direct_source_coarse_failure',dict(tau=30,beta=0,n=301,expected='energy_accuracy_failure'),lambda:moving(n=301,tau=30,beta=0,direct_source=True))
if row['status']=='completed_component':check('direct_source_coarse_conservation',row['summary']['luminosity_span'],.005,'retained converged-but-inaccurate method',expected_pass=False)

# Stationary self-consistent force and moving radiation, supplied mass flow.
def wind(n=401,outer=100,gamma=2,tiring=.001,escape_speed_c=.0001,tolerance=1e-5,transport_tolerance=1e-6,core_rays=32):
 r=np.unique(np.r_[1,np.geomspace(1.00001,outer,n),1/(1-np.geomspace(1e-8,.2,100))])
 q=solve_continuum_wind(r,gamma=gamma,tiring=tiring,escape_speed_c=escape_speed_c,tolerance=tolerance,transport_tolerance=transport_tolerance,core_rays=core_rays)
 wref,lref=flux_driven_reference(r,gamma,tiring,1e-4)
 data={key:val for key,val in q.items() if isinstance(val,np.ndarray)};data['r']=r
 data.update(velocity=np.sqrt(q['w']),density=q['extinction'],force=q['H_comoving'],reference_w=wref,reference_L=lref,reference_force=lref/r**2)
 _,ah0,_=comoving_moments(q['angular_J'],q['angular_H'],q['angular_K'],q['beta'])
 xx=1-1/r;momentum=PchipInterpolator(xx,gamma*r*r*ah0-1).antiderivative();defect=q['w']-q['w'][0]-(momentum(xx)-momentum(0))
 work=PchipInterpolator(r,r*r*q['extinction']*q['beta']/np.sqrt(1-q['beta']**2)*ah0).antiderivative()
 workdef=q['luminosity_ratio']-1+work(r)-work(r[0])
 cg=dimensional_fields(q,100*MSUN,.34)
 # Independent analytic density/optical-depth comparison; correct the rest-current
 # Lorentz factor explicitly rather than silently treating lab and comoving rho alike.
 bref=escape_speed_c*np.sqrt(wref);cref=tiring*gamma*np.sqrt(1-bref*bref)/(bref*r*r)
 tau_ref=-PchipInterpolator(r,cref*(1-bref)/np.sqrt(1-bref*bref)).antiderivative()(r);tau_ref-=tau_ref[-1]
 rphoto_ref=float(np.interp(2/3,tau_ref[::-1],r[::-1])) if tau_ref[0]>=2/3 else None
 summary=dict(parameters=q['parameters'],history=q['history_wind'],
  velocity_reference_error=float(max(abs(np.sqrt(q['w']/wref)-1))),
  density_reference_error=float(max(abs(q['extinction']/cref-1))),
  luminosity_reference_error=float(max(abs(q['luminosity_ratio']/lref-1))),
  force_reference_error=float(max(abs(q['H_comoving']/(lref/r**2)-1))),
  total_energy_defect=float(max(abs(q['total_energy_ratio']/q['total_energy_ratio'][0]-1))),
  independent_angular_momentum_defect=float(max(abs(defect))/max(q['w'])),
  independent_angular_work_defect=float(max(abs(workdef))),
  optical_depth_base=float(q['optical_depth_radial_to_boundary'][0]),
  radiation_drag_max=float(max(q['radiation_drag_fraction'])),
  photospheric_radius_ratio=q['photospheric_radius_ratio'],analytic_photospheric_radius_ratio=rphoto_ref,
  photosphere_reference_error=abs(q['photospheric_radius_ratio']/rphoto_ref-1) if rphoto_ref and q['photospheric_radius_ratio'] else None,
  velocity_to_local_escape_at_outer=float(np.sqrt(q['w'][-1]*r[-1])),
  cg_scalars={k:v for k,v in cg.items() if not isinstance(v,np.ndarray)},
  missing=['gas EOS and T','ionization/flux mean','gas sonic regularity','mass-loss eigenvalue','physical envelope match'])
 for k,v in cg.items():
  if isinstance(v,np.ndarray):data[k]=v
 return data,summary
windcases=[('wind_n201',dict(n=201)),('wind_n401',dict(n=401)),('wind_n801',dict(n=801)),
 ('wind_tolerance',dict(n=401,tolerance=1e-7,transport_tolerance=1e-8)),
 ('wind_core64',dict(n=401,core_rays=64)),('wind_outer300',dict(n=801,outer=300)),
 ('wind_crosscase',dict(n=401,gamma=4,tiring=.0005,escape_speed_c=.0002)),
 ('wind_drag_negative_control',dict(n=401,tiring=.01,escape_speed_c=.001))]
for name,params in windcases:
 row=save(name,params,lambda params=params:wind(**params))
 if row['status']=='completed_component':
  s=row['summary'];expected=name!='wind_drag_negative_control'
  for k in ['velocity_reference_error','density_reference_error','luminosity_reference_error','force_reference_error','photosphere_reference_error']:
   if s[k] is not None:check(name+':'+k,s[k],.01,'stationary small-drag reference',expected_pass=expected)
  check(name+':energy',s['total_energy_defect'],.001,'stationary whole energy')
  check(name+':angular_momentum',s['independent_angular_momentum_defect'],.005,'stationary independent momentum')
  check(name+':angular_work',s['independent_angular_work_defect'],.001,'stationary independent work')
  check(name+':angular_flux',s['history'][-1]['radiation']['angular_flux_max_error'],.005,'stationary angular closure')
for label,a,b in [('mesh201','wind_n201','wind_n801'),('mesh401','wind_n401','wind_n801'),('tolerance','wind_tolerance','wind_n401'),('core_rays','wind_core64','wind_n401'),('outer','wind_outer300','wind_n801')]:
 if (out/f'{a}.npz').exists() and (out/f'{b}.npz').exists():compare(a,b,1,100,['velocity','density','J','H','K','force'],.002,'wind_'+label)
summary=dict(outcome='OUTER CLOSURE FAIL',cases=len(cases),completed=sum(r['status']=='completed_component' for r in cases),
 failures=[r['name'] for r in cases if r['status']!='completed_component'],checks=len(checks),
 failed_checks=[r for r in checks if not r['passed']],elapsed_seconds=time.perf_counter()-start,
 peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
 source_snapshot_verified=io.finish_provenance(out),coupled_evolution_run=False,
 physical_envelope_match=False,published_full_evolution_replicas=0,v1=False)
io.json_write(out/'campaign.json',summary);print(json.dumps(summary),flush=True)
