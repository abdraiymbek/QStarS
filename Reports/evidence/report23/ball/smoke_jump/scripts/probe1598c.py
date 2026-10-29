"""Step 1598 (DRYRUN): continuation in dt on the event mesh, each solve started from the previous solution."""
import os,sys,copy,json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'))
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import numpy as np
from qstars import quasistatic as qs,cavity
from qstars.nonlinear import pattern_jacobian
import validate_ball_track as vbt
YEAR=365.25*86400;MSUN=1.98847e33
st,config=qs.load_checkpoint(os.path.join(HERE,'DRYRUN_continuation_counterfactual_400','checkpoint.npz'))
eos,opacity=vbt.physics(config)
s=qs.settings(config);q_=s['quasistatic'];zones=config['numerics']['zones']
def problem(dt,remesh=True):
    pred=qs.predict(st,dt*YEAR,eos,config,'evolve',None,s)
    cut,ci=cavity.slice_base(st,pred['swept'],eos)
    pre,rg=qs.regrid_base(st,cut,eos,zones,q_['base_mass_fraction'])
    post,rm=qs.remesh_base(st,pre,eos,zones,q_['split'],q_['merge'],protect=rg['region_cells']) if remesh else (pre,None)
    return qs.Problem(st,post,dt*YEAR,eos,opacity,config,'evolve',None,pred,s)
rows=[]
for remesh in (True,False):
    x=None
    for dt in np.geomspace(1e-3,77.25,25):
        P=problem(dt,remesh)
        start=P.x0 if x is None else x
        if x is not None:start=start.copy();start[3*P.n+2]=P.x0[3*P.n+2]  # growth unknown from this dt's predictor
        z,info=qs.newton(P.evaluate,start,P.pattern,q_['tolerance'],q_['max_iterations'],q_['polish_iterations'],q_['roundoff_residual'],q_['roundoff_step'])
        J=pattern_jacobian(P.evaluate,z,P.pattern).toarray();S=np.linalg.svd(J,compute_uv=False)
        n=P.n;th=z[n+1:2*n+1]-P.x0[n+1:2*n+1]
        rows.append(dict(event_mesh=remesh,dt_years=float(dt),success=bool(info['success']),iterations=info['iterations'],residual=info['residual'],sigma_min=float(S[-1]),
                         mapped_dlnT=float(np.max(abs(th))),cell=int(np.argmax(abs(th)))))
        print(json.dumps(rows[-1]))
        if not info['success']:break
        x=z
json.dump(rows,open(os.path.join(HERE,'probe1598_dt_continuation.json'),'w'),indent=1)
