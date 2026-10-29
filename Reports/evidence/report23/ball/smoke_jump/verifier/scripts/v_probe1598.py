"""Probe the investigator's DRYRUN end states (read-only copies in their scratch) with the UNPATCHED code.
1598: is the Newton failure tied to the step's merge? 1598/1624: outer-zone discretization error after the
accepted merges (refined x4 outer cells vs the coarse mesh)."""
from v_common import *
from qstars.mesh import split_merge
SJ=Path('/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/ball_jump')
which=sys.argv[1]
ck={'1598':SJ/'DRYRUN_continuation_counterfactual_400/checkpoint.npz','1624':SJ/'DRYRUN_continuation_counterfactual_ladder_from1598_400/checkpoint.npz'}[which]
state,config=qs.load_checkpoint(ck);config['stepcontrol'].pop('events',None)
eos,opacity=vbt.physics(config)
print('step',state.steps,'Mbh',state.bh_mass/MSUN,'zones',len(state.dm))
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
n0=len(state.dm)
idx=list(range(-16,-2))
mbs0=state.dm.sum()-(np.cumsum(state.dm)-.5*state.dm);probe=mbs0[idx]
kap=np.atleast_1d(opacity(state.rho,state.T))
print('cells',[n0+i for i in idx]);print('T',np.round(state.T[idx]));print('dlnkap/cell',np.round(np.diff(np.log(kap))[idx[:-1]],3))
variants=[('registered',None),('no event',[]),('x2 outer 20',[(2,range(-20,0))]),('x4 outer 20',[(4,range(-20,0))]),('x8 outer 20',[(8,range(-20,0))])]
ref={}
for dt in (1.,30.9):
    for name,spec in variants:
        mode['spec']=spec;stash.clear()
        try:
            new,d,m=attempt(state,dt,eos,opacity,config);ok=f"ok it={d['newton_iterations']} m/s={d['merges']}/{d['splits']} delta={m['delta_actual']:.4g}";T=new.T;dm=new.dm
        except Exception as e:
            if 'main' not in stash:print(name,dt,'FAIL before solve',e);continue
            p,x,info=stash['main'];n=p.n;T=p.post.T*np.exp(x[n+1:2*n+1]);dm=p.post.dm
            ok=f"FAIL residual {info['residual']:.3g} last_step {info['last_step']:.2e}"
        mbs=dm.sum()-(np.cumsum(dm)-.5*dm)
        lt=np.interp(probe[::-1],mbs[::-1],np.log(T)[::-1])[::-1]
        ref[(name,dt)]=lt
        print(f"{which} dt={dt:6.3g} {name:12s} zones={len(dm)} {ok}",flush=True)
    base=ref.get(('no event',dt))
    if base is not None:
        for name,_ in variants:
            if (name,dt) in ref:print(f"   {name:12s} dlnT vs no-event @probe: {np.array2string(ref[(name,dt)]-base,precision=4,floatmode='fixed',suppress_small=True)}")
qs.remesh_base=orig_remesh;qs.Problem.solve=orig_solve
