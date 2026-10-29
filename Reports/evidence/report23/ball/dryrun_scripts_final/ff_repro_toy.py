"""DRYRUN (analytic toy only, scratch): reproduce first-review findings e7 (relax.json lost when the post-processing of a
converged seed raises) and h1 (nonfinite start row labelled a model_domain_exit of row 0) with the pre-round tree
(argv[1] = 'before') or the current tree ('after')."""
import os,sys,json,copy,importlib.util,tempfile
from pathlib import Path
which=sys.argv[1];REPO=Path('/Users/abraiymbek/Projects/QStarS')
ROOT=Path(os.environ['FF_BEFORE']) if which=='before' else REPO
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(REPO/'tests'))
spec=importlib.util.spec_from_file_location('vbt',ROOT/'scripts/validate_ball_track.py');vbt=importlib.util.module_from_spec(spec);spec.loader.exec_module(vbt)
from test_quasistatic import toy,EOS,OPACITY
from qstars import quasistatic,validity
import qstars;assert str(ROOT) in qstars.__file__,qstars.__file__
PLAN=json.loads((REPO/'benchmarks/ball2011_track_plan.json').read_text())
_,base=toy(48,'cell')
config=copy.deepcopy(base);config['name']='toy_C0';config['relaxation']=dict(vbt.RELAXATION,mesh_cycles=40)
config['replication']=dict(run='C0',kind='convergence_campaign',plan_text='toy',policy='toy',total_mass_msun=100.,initial_bh_mass_msun=1.,
                           seed_radii_rsun=[150.,300.],landings_msun=[1.01],amendments=[])
tmp=Path(tempfile.mkdtemp(prefix=f'DRYRUN_repro_{which}_'))
registration=vbt.verify_registration()
save=quasistatic.save_profile
def bad(*a,**k):raise ValueError('injected: profile of a converged seed')
quasistatic.save_profile=bad
path=vbt.run_directory(tmp,'relax','C0',config,registration,[])
try:out=vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.],48,log=lambda *a:None);r='returned '+out['status']
except Exception as e:r=f'raised {type(e).__name__}: {e}'
print(f'[{which}] e7 save_profile error in a converged seed: relax_seeds {r}; relax.json written: {(path/"relax.json").exists()}')
quasistatic.save_profile=save
# h1: a start row flagged nonfinite
s,_=quasistatic.relax(toy(48,'cell')[0],EOS,OPACITY,config,dict(mesh_cycles=0));s=quasistatic.zero_time(s,EOS)
check=validity.check
def nonfinite(state,*a,**k):
    errors,diag,c=check(state,*a,**k)
    if state.steps==0:diag['screens']['nonfinite_state']=dict(diag['screens']['nonfinite_state'],violated=True);errors=['Nonfinite model state']+errors
    return errors,diag,c
validity.check=nonfinite
cfg=copy.deepcopy(config);cfg['stepcontrol']=dict(dt_initial_years=100.)
_,rep=quasistatic.run(s,EOS,OPACITY,cfg,tmp/'h1',policy='production',landings=[1.01*1.98847e33],max_accepted_steps=3)
print(f'[{which}] h1 nonfinite start row: status {rep["status"]} attempts {rep["attempts"]} reason {rep["reason"][:90]}')
