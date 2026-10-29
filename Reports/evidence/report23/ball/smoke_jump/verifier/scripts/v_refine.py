"""Outer-zone refinement at step 1255 (scratch, in-memory monkeypatches only). Records the Newton end state even
when the registered exit (1e-12) is missed, to estimate the converged outer T and to see whether the stall is
roundoff (tiny last step) or near-singularity (large last step)."""
from v_common import *
from qstars.mesh import split_merge
state,config,eos,opacity=load()
orig_remesh=qs.remesh_base;orig_solve=qs.Problem.solve
mode=dict(spec=None);stash={}
def custom(old,cut,eos_,zones,split=2.,merge=.5,max_passes=40,protect=0):
    spec=mode['spec']
    if spec is None:return orig_remesh(old,cut,eos_,zones,split,merge,max_passes,protect)
    n=len(cut.dm);x=np.ones(n)
    for factor,cells in spec:
        for c in cells:x[n+c]={2:2.5,4:4.5,8:9.,'m':.3}[factor]
    st,info=split_merge(cut,eos_,cut.dm/x,split,merge)
    keep={k:info[k] for k in ('splits','merges','energy_defect','entropy_change','merge_mixing_entropy','potential_reallocation','max_event_defect')}
    return st,dict(keep,passes=1)
def solve(self):
    x,info=orig_solve(self)
    if 'main' not in stash:stash['main']=(self,x,info)
    return x,info
qs.remesh_base=custom;qs.Problem.solve=solve
mbs0=state.dm.sum()-(np.cumsum(state.dm)-.5*state.dm)
idx=[-10,-9,-8,-7,-6,-5,-4,-3]
probe=mbs0[idx]
print('probe cells',[len(state.dm)+i for i in idx],'mass below surface',np.round(probe/MSUN,4))
variants=[('none',[]),('registered (merge 168+169)',None),('x2 cells 164..173',[(2,range(-10,0))]),('x2 cells 150..173',[(2,range(-24,0))]),
          ('x4 cells 160..173',[(4,range(-14,0))]),('x8 cells 160..173',[(8,range(-14,0))]),
          ('x2 160..173 + merge',[(2,range(-14,-6)),('m',[-6,-5]),(2,range(-4,0))])]
dt=1.
for name,spec in variants:
    mode['spec']=spec;stash.clear()
    try:
        new,d,m=attempt(state,dt,eos,opacity,config);ok=True;T=new.T;dm=new.dm
    except Exception as e:
        ok=False;p,x,info=stash['main'];n=p.n
        T=p.post.T*np.exp(x[n+1:2*n+1]);dm=p.post.dm
        msg=f"FAIL residual {info['residual']:.3g} last_step {info['last_step']:.3g} it {info['iterations']}"
    p,x,info=stash['main']
    mbs=dm.sum()-(np.cumsum(dm)-.5*dm)
    lnT_at=np.interp(probe[::-1],mbs[::-1],np.log(T)[::-1])[::-1]
    dp=lnT_at-np.log(state.T[idx])
    print(f"{name:28s} zones={len(dm)} {'ok it=%d'%d['newton_iterations'] if ok else msg} last_step={info['last_step']:.2e} "
          f"dlnT@probe={np.array2string(dp,precision=4,floatmode='fixed',suppress_small=True)}",flush=True)
qs.remesh_base=orig_remesh;qs.Problem.solve=orig_solve
