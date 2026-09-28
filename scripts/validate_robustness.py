#!/usr/bin/env python3
"""Timestep, nonlinear tolerance, remesh threshold and restart evidence."""
import argparse,copy,json
from pathlib import Path
import numpy as np
from qstars import io
from qstars.config import resolve
from qstars.evolution import physics,run
from validate_spatial import profile_difference
p=argparse.ArgumentParser();p.add_argument('--zones',type=int,default=256);p.add_argument('--end-time',type=float,default=10.)
p.add_argument('--parent',default='results');p.add_argument('--mesa-library');p.add_argument('--mesa-dir')
p.add_argument('--mechanics',choices=['midpoint','backward_euler'],default='midpoint')
p.add_argument('--mesh',choices=['log_radius','log_radius_cosine','log_radius_surface'],default='log_radius_cosine');a=p.parse_args()
out=io.create_directory(a.parent,'robustness_validation')
base=resolve({'name':'robustness','numerics':{'zones':a.zones,'mesh':a.mesh,'mechanics':a.mechanics,'solver_tolerance':1e-10,'initial_dt_years':.02,'max_dt_years':.02},
 'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4},'feeding':{'rate_msun_per_year':.01,'thermodynamics':'surface'},
 'wind':{'model':'prescribed','rate_msun_per_year':.008},'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},
 'evolution':{'end_time_years':a.end_time,'max_wall_seconds':300.},'output':{'profile_every_steps':500}})
if a.mesa_library:
 base['eos']='mesa';base['opacity']['model']='mesa';base['mesa']={'library_path':str(Path(a.mesa_library).resolve()),'mesa_dir':str(Path(a.mesa_dir).resolve())}
eos,_=physics(base);rows=[];states=[]
cases=[('timestep',v,{'max_dt_years':v,'initial_dt_years':v}) for v in (.02,.01,.005)]
cases += [('solver',v,{'solver_tolerance':v}) for v in (1e-8,1e-10,1e-12)]
cases += [('remesh',v,{'boundary_remesh_fraction':v}) for v in (.4,.2,.1)]
for axis,value,changes in cases:
 c=copy.deepcopy(base);c['name']=axis+'_'+str(value).replace('.','p');c['numerics'].update(changes)
 path,r=run(c,out);s,_=io.load_checkpoint(path/'checkpoint.npz')
 rows.append({'axis':axis,'value':value,'path':str(path),**r});states.append(s)
 io.json_write(out/'campaign.json',{'complete':False,'rows':rows})
 print(axis,value,r['status'],r['runtime_seconds'],flush=True)
for axis in ('timestep','solver','remesh'):
 indices=[i for i,row in enumerate(rows) if row['axis']==axis];ref=indices[-1]
 for i in indices[:-1]:
  if rows[i]['status']==rows[ref]['status']=='completed_requested_interval':
   rows[i]['comparison_to_strictest']=profile_difference(states[i],states[ref],eos)
# Exact checkpoint continuation, including accepted remesh history.
c=copy.deepcopy(base);c['name']='restart_part'
part,rp=run(c,out,stop_after_steps=100);s,cc=io.load_checkpoint(part/'checkpoint.npz')
continued,rc=run(cc,out,restart=s);s,_=io.load_checkpoint(continued/'checkpoint.npz')
ref=states[0]
restart_equal=all(np.array_equal(getattr(s,k),getattr(ref,k)) for k in ('r','T','dm','v','cell_volume','radius_compensation','face_luminosity')) and s.ledger==ref.ledger
summary={'complete':True,'rows':rows,'restart':{'part':str(part),'continued':str(continued),'report':rc,'bitwise_equal':restart_equal},
         'scope':'Fixed-excision hot pressure-confined two-stream model','v1_validated':False}
io.json_write(out/'campaign.json',summary);print(out)
