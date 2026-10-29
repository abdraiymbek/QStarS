"""DRYRUN attacks on (d) row-0 screens, (e) records of failing runs, (h) labels. Scratch output only."""
from rp_common import *
import glob,copy
from qstars import validity,fold,io
from qstars.numerics import StepRejected
OUT=HERE+'/DRYRUN_attack_e';os.makedirs(OUT,exist_ok=True)
registration=vbt.verify_registration()
cu=sorted(glob.glob(HERE+'/DRYRUN_out/*DRYRUN_relax_A3+A2_C1_n200'))[-1]
ROW0=cu+'/row0_6000rsun.npz';A3,A32=amend(['A3']),amend(['A3','A2'])
def files(path):return sorted(p.name for p in path.iterdir() if p.suffix in ('.json','.csv','.jsonl'))
def newest(before):
    made=[p for p in __import__('pathlib').Path(OUT).iterdir() if p not in before];return max(made,key=lambda p:p.name) if made else None
def listing():return set(__import__('pathlib').Path(OUT).iterdir())
print('=== (d) row-0 enforced screens')
check=validity.check
# d6a: a real enforced screen at the start row, via the library with a tightened applicability limit (not a registered configuration)
state,stored=quasistatic.load_checkpoint(ROW0);eos,opacity=vbt.physics(stored)
tight=copy.deepcopy(stored);tight['validity']['max_compactness']=1e-12
p=__import__('pathlib').Path(OUT)/'DRYRUN_run_tight';p.mkdir(exist_ok=True)
for f in p.iterdir():f.unlink()
_,rep=quasistatic.run(state,eos,opacity,tight,p,landings=[100*MSUN],end_bh_mass=100*MSUN,policy='published_assumption')
print('d6a run() tightened compactness:',rep['status'],'|',rep['reason'][:120],'| attempts',rep['attempts'],'| start_row_screen_errors',rep['start_row_screen_errors'])
# d6b: through the driver (smoke and track) with one enforced error injected at the start row only
def start_error(s,*a,**k):
    e,d,c=check(s,*a,**k)
    return (e+['Radiation gradient scale is outside the diffusion limit (Kn>0.1) at tau>=10'] if s.steps==0 else e),d,c
validity.check=start_error
b=listing();path,out=vbt.smoke(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,ROW0)
print('d6b smoke:',{k:out[k] for k in ('status','accepted_steps','attempts','passed','start_row_screen_errors')},'| files',files(path))
b=listing();path,summ=vbt.track(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,run_id='C1',row0=ROW0)
print('d6b track:',{k:summ.get(k) for k in ('status','label','attempts','start_row_screen_errors')},'| fold.json',(path/'fold.json').exists(),'| track.json',(path/'track.json').exists())
print('    endpoint():',{k:v for k,v in vbt.endpoint(PLAN,summ).items() if k in ('outcome','label','reason')})
print('=== (h) numerical defects at the start row')
def nonfinite(s,*a,**k):
    e,d,c=check(s,*a,**k)
    return (e+['Nonfinite model state'] if s.steps==0 else e),d,c
validity.check=nonfinite
path,summ=vbt.track(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,run_id='C1',row0=ROW0)
print('h1 start row flagged Nonfinite model state -> label:',summ.get('label'),'| status',summ.get('status'))
def backstop(s,*a,**k):
    e,d,c=check(s,*a,**k)
    return (e+['MESA first-law defect exceeds 1% audit limit'] if s.steps==0 else e),d,c
validity.check=backstop
path,summ=vbt.track(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,run_id='C1',row0=ROW0)
print('h2 start row flagged MESA first-law defect -> label:',summ.get('label'),'| status',summ.get('status'))
def conservation(s,*a,**k):
    if s.steps==0:raise StepRejected('Conservation tolerance exceeded')
    return check(s,*a,**k)
validity.check=conservation
path,summ=vbt.track(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,run_id='C1',row0=ROW0)
print('h3 start row raises Conservation tolerance exceeded -> label:',summ.get('label'),'| status',summ.get('status'),'|',summ.get('reason'))
validity.check=check
print('=== (e) records of failing runs')
run=quasistatic.run
def broken(*a,**k):raise KeyError('injected')
quasistatic.run=broken
path,summ=vbt.track(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,run_id='C1',row0=ROW0)
print('e1 track, run raises KeyError:',summ.get('status'),summ.get('label'),'| files',files(path))
path,out=vbt.smoke(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,ROW0)
print('e2 smoke, run raises KeyError:',out['status'],out['passed'],'| files',files(path))
def killed(*a,**k):raise SystemExit(137)
quasistatic.run=killed
b=listing()
try:vbt.smoke(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,ROW0)
except SystemExit:pass
path=newest(b);print('e3 smoke, run raises SystemExit:','| files',files(path),'| smoke.json status',json.load(open(path/'smoke.json'))['status'] if (path/'smoke.json').exists() else None)
quasistatic.run=run
# history_row fails with an uncaught error type at the third row (inside run): what do termination.json and smoke.json say
row=quasistatic.history_row;n=[0]
def third(*a,**k):
    n[0]+=1
    if n[0]>=1 and k.get('event')!='start':raise KeyError('injected in history_row')
    return row(*a,**k)
quasistatic.history_row=third
path,out=vbt.smoke(None,PLAN,registration,amend(['A3','A2']),LIB,MESA_DIR,OUT,ROW0)
quasistatic.history_row=row
print('e4 smoke, registered numerics (no accepted step, so history_row never fails here):',out['status'])
# fold.classify raising
classify=fold.classify
def bad(*a,**k):raise IndexError('injected in fold')
fold.classify=bad
path,summ=vbt.track(None,PLAN,registration,A32,LIB,MESA_DIR,OUT,run_id='C1',row0=ROW0)
fold.classify=classify
print('e5 track, termination trigger then fold.classify raises IndexError:',summ.get('status'),'| label',summ.get('label'),'|',summ.get('label_reason'),'| files',files(path))
# relax: an error type relax_seeds does not catch, raised inside quasistatic.relax
relax=quasistatic.relax
def keyerror(*a,**k):raise KeyError('injected in relax')
quasistatic.relax=keyerror
config=vbt.build_config(PLAN,'C1',LIB,MESA_DIR,A3);config['name']='DRYRUN_'+config['name']
path=vbt.run_directory(OUT,'DRYRUN_relax_inject',"C1",config,registration,A3)
try:
    out=vbt.relax_seeds(PLAN,config,eos,opacity,path,config['replication']['seed_radii_rsun'],200,log=lambda *a:None);print('e6 relax returned',out['status'])
except Exception as e:print('e6 relax, quasistatic.relax raises KeyError -> relax_seeds raised',type(e).__name__,'| relax.json written:',(path/'relax.json').exists(),'| files',files(path))
quasistatic.relax=relax
# relax: the post-processing of a converged seed fails (profile output)
save=quasistatic.save_profile
def noprofile(*a,**k):raise ValueError('injected in save_profile')
quasistatic.save_profile=noprofile
path=vbt.run_directory(OUT,'DRYRUN_relax_inject2',"C1",config,registration,A3)
try:
    out=vbt.relax_seeds(PLAN,config,eos,opacity,path,config['replication']['seed_radii_rsun'],200,log=lambda *a:None);print('e7 relax returned',out['status'])
except Exception as e:print('e7 relax, save_profile of a converged seed raises ValueError -> relax_seeds raised',type(e).__name__,'| relax.json written:',(path/'relax.json').exists(),'| files',files(path))
quasistatic.save_profile=save
