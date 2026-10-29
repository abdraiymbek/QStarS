"""DRYRUN (scratch only, not a registered run): continue the failed smoke run from its last accepted state with the
candidate fix (stepcontrol events='counterfactual'), the smoke gate's options, for at most STEPS accepted steps."""
import os,sys,time
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'))
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import copy,json,shutil
from pathlib import Path
import qstars;assert 'patched' in qstars.__file__
from qstars import quasistatic as qs
import validate_ball_track as vbt
YEAR=365.25*86400
RUN=Path('/Users/abraiymbek/Projects/QStarS/benchmarks/results/2026-10-14-ball-track/20261014T051036.115106Z_smoke_C1_n200')
events=sys.argv[1] if len(sys.argv)>1 else 'counterfactual'
STEPS=int(sys.argv[2]) if len(sys.argv)>2 else 400
out=Path(HERE)/f'DRYRUN_continuation_{events}_{STEPS}'
if out.exists():shutil.rmtree(out)
out.mkdir()
plan=json.loads(Path('/Users/abraiymbek/Projects/QStarS/benchmarks/ball2011_track_plan.json').read_text())
state,config=qs.load_checkpoint(RUN/'checkpoint.npz')
config=copy.deepcopy(config);config['stepcontrol']['events']=events
# The controller as accepted at step 1255 (the checkpoint holds it after the 8 rejections): the first
# rejected attempt was the controller's dt_next, no hold, no retries.
state.aux['controller']=dict(dt_next=172.1438701239381*YEAR,hold=0,retries=0)
state.aux['driver']['screen_retries']=0
eos,opacity=vbt.physics(config)
options=vbt.track_options(plan,config)
t0=time.time()
_,report=qs.run(state,eos,opacity,config,out,end_bh_mass=options['landings'][0],stop_after_steps=STEPS,
                **{**options,'landings':options['landings'][:1]})
report['wall_seconds']=time.time()-t0
report['label']='DRYRUN (scratch, not registered): continuation of the smoke run from step 1255 with stepcontrol.events='+events
(out/'DRYRUN_report.json').write_text(json.dumps(qs._json(report),indent=1))
print(json.dumps({k:report[k] for k in ('status','reason','accepted_steps_this_run','attempts','rejections','final_step','final_bh_mass_msun','final_time_years','final_dt_years','wall_seconds')},indent=1))
