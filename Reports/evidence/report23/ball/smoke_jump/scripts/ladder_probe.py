import os,sys,copy
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'));sys.path.insert(0,HERE)
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
from qstars import quasistatic as qs,stepcontrol as sc
import validate_ball_track as vbt,ladder_proto
ladder_proto.install(qs)
YEAR=365.25*86400
st,config=qs.load_checkpoint(os.path.join(HERE,'DRYRUN_continuation_counterfactual_400','checkpoint.npz'))
eos,opacity=vbt.physics(config);S=sc.settings(config);before=sc.observables(st,eos,'face')
for dt in (77.25,30.9):
    new,d=qs.bondi_step(st,dt*YEAR,eos,opacity,config,mode='evolve')
    m=sc.metrics(before,sc.observables(new,eos,'face'),d,S)
    print(dt,'merges',d['merges'],'it',d['newton_iterations'],'measured',round(d['mapped_max_dlnT'],4),'event-free',d['event_free_mapped_max_dlnT'],'delta',m['delta_actual'],'rejected',sc.rejected(m['delta_actual'],S))
