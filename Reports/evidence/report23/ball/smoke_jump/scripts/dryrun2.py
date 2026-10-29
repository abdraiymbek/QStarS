"""DRYRUN (scratch only, exploratory): continue the DRYRUN from step 1598 with events='counterfactual' AND the
exploratory dt-ladder Newton start (ladder_proto, not part of the proposed patch), for at most STEPS accepted steps."""
import os,sys,time,shutil,copy,json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'));sys.path.insert(0,HERE)
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
from pathlib import Path
from qstars import quasistatic as qs
import validate_ball_track as vbt,ladder_proto
ladder_proto.install(qs)
YEAR=365.25*86400
STEPS=int(sys.argv[1]) if len(sys.argv)>1 else 400
out=Path(HERE)/f'DRYRUN_continuation_counterfactual_ladder_from1598_{STEPS}'
if out.exists():shutil.rmtree(out)
out.mkdir()
plan=json.loads(Path('/Users/abraiymbek/Projects/QStarS/benchmarks/ball2011_track_plan.json').read_text())
state,config=qs.load_checkpoint(Path(HERE)/'DRYRUN_continuation_counterfactual_400'/'checkpoint.npz')
assert state.steps==1598 and config['stepcontrol']['events']=='counterfactual'
state.aux['controller']=dict(dt_next=30.899*YEAR,hold=2,retries=0);state.aux['driver']['screen_retries']=0
eos,opacity=vbt.physics(config);options=vbt.track_options(plan,config)
t0=time.time()
_,report=qs.run(state,eos,opacity,config,out,end_bh_mass=options['landings'][0],stop_after_steps=STEPS,**{**options,'landings':options['landings'][:1]})
report['wall_seconds']=time.time()-t0
report['label']='DRYRUN (scratch, exploratory): events=counterfactual plus dt-ladder Newton start, from DRYRUN step 1598'
(out/'DRYRUN_report.json').write_text(json.dumps(qs._json(report),indent=1))
print(json.dumps({k:report[k] for k in ('status','reason','accepted_steps_this_run','attempts','rejections','final_step','final_bh_mass_msun','final_time_years','final_dt_years','wall_seconds')},indent=1))
