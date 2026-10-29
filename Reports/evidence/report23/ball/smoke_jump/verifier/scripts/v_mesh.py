"""Mesh-sensitivity probe at step 1255 (scratch; in-memory monkeypatch of remesh_base only, no file of the tree changed).

Each variant replaces P2 (remesh_base) by one split_merge pass with a hand-made target so that only the
named cells (indexed from the surface: -1 is the surface cell) are merged or split. T is compared with the
accepted pre-step state at the same mass below the surface."""
from v_common import *
from qstars.mesh import split_merge
state,config,eos,opacity=load()
orig=qs.remesh_base
mode=dict(spec=None)
def custom(old,cut,eos_,zones,split=2.,merge=.5,max_passes=40,protect=0):
    spec=mode['spec']
    if spec is None:
        return orig(old,cut,eos_,zones,split,merge,max_passes,protect)
    n=len(cut.dm);x=np.ones(n)
    for kind,cells in spec:
        for c in cells:
            x[n+c]=.3 if kind=='merge' else 2.5
    st,info=split_merge(cut,eos_,cut.dm/x,split,merge)
    keep={k:info[k] for k in ('splits','merges','energy_defect','entropy_change','merge_mixing_entropy','potential_reallocation','max_event_defect')}
    return st,dict(keep,passes=1)
qs.remesh_base=custom
mbs0=state.dm.sum()-(np.cumsum(state.dm)-.5*state.dm)
probe=mbs0[[-8,-7,-6,-5,-4]]  # mass below surface of original cells 166..170
variants=[('registered',None),('none',[]),
          ('merge 168+169',[('merge',[-6,-5])]),('merge 169+170',[('merge',[-5,-4])]),('merge 167+168',[('merge',[-7,-6])]),
          ('merge 150+151',[('merge',[-24,-23])]),('merge 140+141',[('merge',[-34,-33])]),
          ('split 168,169',[('split',[-6,-5])]),('split 164..173',[('split',list(range(-10,0)))]),
          ('split 150..173',[('split',list(range(-24,0)))])]
out={}
for dt in (1.,172.1438701239381):
    for name,spec in variants:
        mode['spec']=spec
        try:new,d,m=attempt(state,dt,eos,opacity,config)
        except Exception as e:
            print(f'{name:16s} dt={dt:.4g} FAIL {type(e).__name__}: {str(e)[:150]}',flush=True);continue
        mbs=new.dm.sum()-(np.cumsum(new.dm)-.5*new.dm)
        lnT_at=np.interp(probe[::-1],mbs[::-1],np.log(new.T)[::-1])[::-1]
        d_probe=lnT_at-np.log(state.T[[-8,-7,-6,-5,-4]])
        dl=np.log(new.T)-below_surface_map(state,new);k=int(np.argmax(abs(dl)))
        print(f"{name:16s} dt={dt:8.4g} zones={d['zones']} m/s={d['merges']}/{d['splits']} it={d['newton_iterations']} delta={m['delta_actual']:.4g} "
              f"profile_dlnT={m['profile_dlnT']:.4g} maxcell={k-len(new.T)} "
              f"dlnT@orig166..170={np.array2string(d_probe,precision=4,floatmode='fixed')} dlnTeff={m['dln_effective_temperature']:.2e} dlnL={m['dln_luminosity']:.2e}",flush=True)
        out[(name,dt)]=d_probe
qs.remesh_base=orig
