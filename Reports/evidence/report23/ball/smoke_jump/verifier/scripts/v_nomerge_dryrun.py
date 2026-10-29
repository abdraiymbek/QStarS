"""DIAGNOSTIC DRYRUN (scratch only, not a registered run, not a proposed setting): continue the failed smoke run from
step 1255 with the UNPATCHED registered code and controller (events measured), but with merges switched off
(quasistatic.merge = 0) to test whether the later failures (step-1598 merge Newton failure, step-1624 singular step)
share the root cause of outer-zone coarsening. Same controller restoration as the investigator's DRYRUN."""
from v_common import *
import shutil
merge=float(sys.argv[1]) if len(sys.argv)>1 else 0.
STEPS=int(sys.argv[2]) if len(sys.argv)>2 else 500
out=OUT/f'DRYRUN_registered_controller_merge{merge:g}_{STEPS}'
if out.exists():shutil.rmtree(out)
out.mkdir()
plan=json.loads((ROOT/'benchmarks/ball2011_track_plan.json').read_text())
state,config=qs.load_checkpoint(RUN/'checkpoint.npz')
config=copy.deepcopy(config);config['quasistatic']['merge']=merge
assert 'events' not in config['stepcontrol']
state.aux['controller']=dict(dt_next=172.1438701239381*YEAR,hold=0,retries=0)
state.aux['driver']['screen_retries']=0
eos,opacity=vbt.physics(config)
options=vbt.track_options(plan,config)
t0=time.time()
_,report=qs.run(state,eos,opacity,config,out,end_bh_mass=options['landings'][0],stop_after_steps=STEPS,
                **{**options,'landings':options['landings'][:1]})
report['wall_seconds']=time.time()-t0
(out/'DRYRUN_report.json').write_text(json.dumps(qs._json(report),indent=1))
print(json.dumps({k:report[k] for k in ('status','reason','accepted_steps_this_run','attempts','rejections','final_step','final_bh_mass_msun','final_time_years','final_dt_years','wall_seconds')},indent=1))
