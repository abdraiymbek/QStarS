#!/usr/bin/env python3
"""Moving-face identities and explicit MLT exit, independent of evolution."""
import argparse,csv,dataclasses,json,os,resource,time,traceback
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
from qstars import io
from qstars.constants import G,MSUN,RSUN
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.freeeos_reference import FreeEOSReference
from qstars.spherical_wind import DiffusionWind,IdealGasOpacity,FreeEOSOpacity
from qstars.atmosphere_reference import integrate_atmosphere,AtmosphereDomainExit
p=argparse.ArgumentParser();p.add_argument('--output',default='results');a=p.parse_args()
start=time.perf_counter();out=io.create_directory(a.output,'outer_transport')
io.provenance(out,dict(scope='Stage A component accounting; no coupled stellar evolution',mass_energy_tolerance=2e-7,mlt_exit_tau_tolerance=1e-5))
(out/'driver.py').write_bytes(Path(__file__).read_bytes());print(out,flush=True)
m=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
e=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,m.radiation_constant)
io.json_write(out/'microphysics.json',dict(mesa=m.manifest,freeeos=e.manifest))
rows=[]
# Radiative + MLT candidates stop on an accepted Mach-one event, not on an
# out-of-domain RK trial; changing output sampling cannot move the event.
for model in ['Henyey','ML1']:
 for rtol,points in [(1e-7,41),(1e-9,301)]:
  try:
   q=integrate_atmosphere(e,m,teff=6000,gravity=1,gas_pressure_outer=1e-6,rtol=rtol,points=points,mlt_model=model)
   row=dict(status='unexpected_completed',max_mach=float(max(q.convective_velocity/q.sound_speed)))
  except AtmosphereDomainExit as exc:row=dict(status='validity_domain_exit',diagnostics=exc.diagnostics)
  row.update(model=model,rtol=rtol,points=points);rows.append(row);print(row,flush=True)
io.json_write(out/'mlt_exits.json',rows)
w=DiffusionWind(IdealGasOpacity(),100*MSUN,50*RSUN,9.71e-9,2.5e5)
sols=[]
for sign,bound in [(-1,np.log(.7)),(1,np.log(1.3))]:
 x0=sign*1e-7;s=solve_ivp(w.rhs,[x0,bound],x0*np.array([w.slope,w.sonic_temperature_slope]),method='DOP853',rtol=1e-11,atol=1e-13,dense_output=True,max_step=.005)
 if not s.success:raise RuntimeError(s.message)
 sols.append(s)
def state(r):
 x=np.log(r/w.rs);z=sols[0 if x<0 else 1].sol(x)
 return x,z,w.state(x,z)
transport=[]
for h in [1e-3,1e-4,1e-5]:
 dt=h*w.rs/w.vs;r0=np.array([.8,1.2])*w.rs;u=np.array([.1,-.05])*w.vs;r1=r0+u*dt
 for kind in ['mass','total_energy']:
  def inventory(r):
   _,_,a=state(r);_,v,T,rho,q,*_=a
   specific=1. if kind=='mass' else q.gas_u+3*q.pr/rho+v*v/2-G*w.mass/r
   return 4*np.pi*r*r*rho*specific
  swept=-quad(inventory,r0[0],r1[0],epsabs=0,epsrel=1e-10)[0]+quad(inventory,r0[1],r1[1],epsabs=0,epsrel=1e-10)[0]
  flux=[]
  for r,speed in zip((r0+r1)/2,u):
   x,z,_=state(r);flux.append(w.boundary_flux(x,z,face_velocity=speed)[kind])
  transport.append(dict(dimensionless_displacement=h,channel=kind,swept_inventory=swept,dt=dt,
   boundary_fluxes=flux,relative_defect=float(abs(swept/dt-flux[0]+flux[1])/max(abs(flux[0]),abs(flux[1])))))
io.json_write(out/'moving_sphere.json',transport)
# Quantify a bookkeeping-only improvement in inward plot coverage. The outer
# discrepancies stay visible. Original printed sonic states are unchanged.
coverage=[];data=list(csv.DictReader(open('benchmarks/nakauchi2017_figure2.csv')))
for rho in [9.65e-9,9.71e-9,9.85e-9,10.5e-9]:
 row=dict(sonic_density=rho)
 try:
  w=DiffusionWind(FreeEOSOpacity(e,m,'paper_gn93'),100*MSUN,50*RSUN,rho,2.5e5)
  f,s=w.integrate(inner_ratio=.05,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,stalled_mach=1.5)
  np.savez_compressed(out/f'paper_{rho}.npz',**f)
  subset=[d for d in data if float(d['sonic_density'])==rho]
  lo=f['radius'][0]/RSUN;hi=f['radius'][-1]/RSUN
  row.update(status='completed_component',summary=s,radial_range_rsun=[float(lo),float(hi)],
      total_samples=len(subset),uncovered_inner=sum(float(d['radius_rsun'])<lo for d in subset),
      uncovered_outer=sum(float(d['radius_rsun'])>hi for d in subset))
 except Exception as exc:row.update(status='numerical_or_domain_failure',error=str(exc),traceback=traceback.format_exc())
 coverage.append(row);io.json_write(out/'coverage.json',coverage);print(rho,row['status'],row.get('radial_range_rsun'),flush=True)
# Isolate the stopping convention on the stalled branch only.
w=DiffusionWind(FreeEOSOpacity(e,m,'paper_gn93'),100*MSUN,50*RSUN,9.65e-9,2.5e5)
f,s=w.integrate(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7)
np.savez_compressed(out/'stall_report9.npz',**f);io.json_write(out/'stall_report9.json',s)
summary=dict(outcome='OUTER FAIL',stage_b_run=False,physical_impossibility_claimed=False,
    moving_face_worst_relative_defect=max(t['relative_defect'] for t in transport),
    mlt_exits=len(rows),elapsed_seconds=time.perf_counter()-start,
    peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
    source_snapshot_verified=io.finish_provenance(out))
io.json_write(out/'campaign.json',summary);print(json.dumps(summary),flush=True)
