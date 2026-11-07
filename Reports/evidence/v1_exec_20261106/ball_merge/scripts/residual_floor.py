"""Residual floor of the quasi-static Newton system where the front is finely resolved (side check for goals 2-3).

(1) Step-1255 state, outer 14 cells refined x f (f = 1, 2, 4, 8) and (2) the merge-free diagnostic DRYRUN's final
state near 55 Msun (verifier, step 1667): Newton end state at the registered tolerance, the residual at x and under
4-ulp perturbations of x (a floor test), the worst rows and their closure-row scale max(0.25 |dlnP|, 1e-3).
Writes data/residual_floor.json."""
from bm_common import *
from qstars.mesh import split_merge
orig_remesh=qs.remesh_base;orig_solve=qs.Problem.solve;orig_effect=qs._remesh_effect
stash={};spec={}
def solve(self):
    x,info=orig_solve(self)
    if not stash.get('in_effect'):stash['last']=(self,x,info)
    return x,info
def effect(*a,**k):
    stash['in_effect']=True
    try:return orig_effect(*a,**k)
    finally:stash['in_effect']=False
def custom(old,cut,eos_,zones,split=2.,merge=.5,max_passes=40,protect=0):
    st,info=orig_remesh(old,cut,eos_,zones,split,spec.get('merge',merge),max_passes,protect)
    f=spec.get('refine',1)
    if f>1:
        n=len(st.dm);x=np.ones(n);x[n-14:]={2:2.5,4:4.5,8:9.}[f]
        st,_=split_merge(st,eos_,st.dm/x,2.,0.)
    return st,info
qs.remesh_base=custom;qs.Problem.solve=solve;qs._remesh_effect=effect
def probe(state,dt,eos,opacity,config):
    stash.clear()
    try:
        attempt(state,dt,eos,opacity,config);ok=True
    except qs.StepRejected:ok=False
    p,x,info=stash['last'];n=p.n
    f=p.evaluate(x);rng=np.random.default_rng(1)
    pert=[float(np.max(abs(p.evaluate(x*(1+4e-16*rng.standard_normal(x.size)))))) for _ in range(8)]
    top=np.argsort(-abs(f))[:3]
    def kind(k):
        if k<n:return f'momentum face {k+1}'
        if k<2*n:return f'energy cell {k-n}'
        if k<3*n-1:return f'closure face {k-2*n+1} (scale {p.contrast_scale[k-2*n]:.2e}, cells from surface {n-(k-2*n+1)})'
        return f'boundary row {k-3*n+1}'
    return dict(ok=ok,zones=n,iterations=info['iterations'],residual=info['residual'],last_step=info['last_step'],
                roundoff_exit=info.get('roundoff_exit'),perturbed_residual_min=min(pert),perturbed_residual_max=max(pert),
                worst_rows=[(kind(int(k)),float(f[k])) for k in top],
                closure_scale_at_floor_fraction=float(np.mean(p.contrast_scale<=1.0001e-3)))
res=dict(step1255=[],step1667=[])
state,config,eos,opacity=load()
for f in (1,2,4,8):
    spec.clear();spec.update(merge=0.,refine=f)
    r=probe(state,1.,eos,opacity,config);r['refine']=f;res['step1255'].append(r);print('1255',r,flush=True)
spec.clear()
ck=ROOT/'Reports/evidence/report23/ball/smoke_jump/verifier/DIAGNOSTIC_DRYRUN_registered_controller_merges_off_from_step1255/checkpoint.npz'
s2,c2=qs.load_checkpoint(ck)
c2=copy.deepcopy(c2);print('1667 checkpoint: steps',s2.steps,'M_BH',s2.bh_mass/MSUN,'zones',len(s2.dm),'merge',c2['quasistatic']['merge'])
for dt in (.25,.5,1.,2.):
    r=probe(s2,dt,eos,opacity,c2);r['dt_years']=dt;res['step1667'].append(r);print('1667',r,flush=True)
qs.remesh_base=orig_remesh;qs.Problem.solve=orig_solve;qs._remesh_effect=orig_effect
jdump(res,DATA/'residual_floor.json')
