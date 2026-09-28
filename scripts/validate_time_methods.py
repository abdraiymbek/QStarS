#!/usr/bin/env python3
"""Compare conservative time integrators at decreasing timesteps, including profiles."""
import copy
from qstars import io
from qstars.config import load
from qstars.evolution import physics,run
from validate_spatial import profile_difference

out=io.create_directory('results','time_method_validation')
base=load('configs/hot_resolved_validation.yaml');base['evolution'].update(end_time_years=1.,max_wall_seconds=600.)
base['output'].update(profile_every_steps=1000,checkpoint_every_steps=1000)
eos,_=physics(base);rows=[];states={}
for dt in (.02,.005,.00125):
    for method in ('midpoint','backward_euler'):
        c=copy.deepcopy(base);c['name']=method+'_'+str(dt).replace('.','p')
        c['numerics'].update(mechanics=method,initial_dt_years=dt,max_dt_years=dt)
        p,r=run(c,out);s,_=io.load_checkpoint(p/'checkpoint.npz')
        states[method]=s
        rows.append(dict(method=method,dt_years=dt,path=str(p),**r))
        io.json_write(out/'campaign.json',{'complete':False,'rows':rows})
        print(method,dt,r['status'],r['runtime_seconds'],flush=True)
    if all(x['status']=='completed_requested_interval' for x in rows[-2:]):
        rows[-1]['method_difference']=profile_difference(states['midpoint'],states['backward_euler'],eos)
io.json_write(out/'campaign.json',{'complete':True,'rows':rows,'v1_validated':False})
print(out)
