import os,sys,copy
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'))
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import numpy as np
from qstars import quasistatic as qs,cavity,stepcontrol as sc
from qstars.numerics import StepRejected
import validate_ball_track as vbt
YEAR=365.25*86400;MSUN=1.98847e33
st,config=qs.load_checkpoint(os.path.join(HERE,'DRYRUN_continuation_counterfactual_400','checkpoint.npz'))
eos,opacity=vbt.physics(config)
s=qs.settings(config);q_=s['quasistatic'];zones=config['numerics']['zones']
# 1. the same steps on the pre-event mesh (no merge)
for dt in (1.,30.9,77.25):
    pred=qs.predict(st,dt*YEAR,eos,config,'evolve',None,s)
    cut,ci=cavity.slice_base(st,pred['swept'],eos)
    pre,rg=qs.regrid_base(st,cut,eos,zones,q_['base_mass_fraction'])
    P=qs.Problem(st,pre,dt*YEAR,eos,opacity,config,'evolve',None,pred,s);x,info=P.solve()
    n=P.n;th=x[n+1:2*n+1]-P.x0[n+1:2*n+1]
    print(f'pre-event mesh dt={dt}: success={info["success"]} it={info["iterations"]} residual={info["residual"]:.3g} mapped dlnT={np.max(abs(th)):.4g}')
# 2. which cells merge at step 1598
pred=qs.predict(st,30.9*YEAR,eos,config,'evolve',None,s);cut,ci=cavity.slice_base(st,pred['swept'],eos)
pre,rg=qs.regrid_base(st,cut,eos,zones,q_['base_mass_fraction'])
post,rm=qs.remesh_base(st,pre,eos,zones,q_['split'],q_['merge'],protect=rg['region_cells'])
a=np.cumsum(pre.dm);b=np.cumsum(post.dm);j=0
while abs(a[j]-b[j])<=1e-12*a[-1]:j+=1
goal=qs.targets(st,pre,eos,zones)
print(f'merge at pre-event cells {j},{j+1} of {len(pre.dm)} (x = {pre.dm[j]/goal[j]:.3f}, {pre.dm[j+1]/goal[j+1]:.3f}); merged cell centre {(b[-1]-b[j]+.5*post.dm[j])/MSUN:.4g} Msun below the surface; T = {pre.T[j]:.0f}, {pre.T[j+1]:.0f} K')
# 3. fallback F1 (registered, no amendment needed) on the event step
c=copy.deepcopy(config);c['quasistatic']['ptc']=True
for dt in (30.9,1.):
    try:
        new,d=qs.bondi_step(st,dt*YEAR,eos,opacity,c,mode='evolve')
        print(f'F1 dt={dt}: accepted solve, iterations {d["newton_iterations"]}, ptc_entered={d.get("ptc_entered")}, measured mapped dlnT={d["mapped_max_dlnT"]:.4g}, event-free={d["event_free_mapped_max_dlnT"]}')
    except StepRejected as e:print(f'F1 dt={dt}: rejected: {str(e)[:300]}')
