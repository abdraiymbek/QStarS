from v_common import *
state,config,eos,opacity=load()
n=len(state.dm)
print('steps',state.steps,'Mbh',state.bh_mass/MSUN,'zones',n,'controller',state.aux.get('controller'))
print('quasistatic split/merge',config['quasistatic'].get('split'),config['quasistatic'].get('merge'),'stepcontrol',config.get('stepcontrol'))
res={}
for merge in (None,0.45):
    c=copy.deepcopy(config)
    if merge is not None:c['quasistatic']['merge']=merge
    for dt in [1e-4,1e-2,1.,172.1438701239381]:
        t0=time.time()
        try:new,d,m=attempt(state,dt,eos,opacity,c)
        except Exception as e:
            print(merge,dt,'FAIL',type(e).__name__,str(e)[:200]);continue
        lnTo=below_surface_map(state,new);dl=np.log(new.T)-lnTo
        k=int(np.argmax(abs(dl)))
        mbs=new.dm.sum()-(np.cumsum(new.dm)-.5*new.dm)
        top=sorted(np.argsort(-abs(dl))[:5].tolist())
        print(f"merge={merge} dt={dt:.4g} zones={d['zones']} merges={d['merges']} splits={d['splits']} it={d['newton_iterations']} "
              f"delta={m['delta_actual']:.5g} comp={max((v,kk) for kk,v in m.items() if kk!='delta_actual')[1]} mapped_dlnT={d['mapped_max_dlnT']:.5g} "
              f"lag_dlnT={d['lagrangian_max_dlnT']:.5g} mapped_dlnr={d['mapped_max_dlnr']:.3g} my_maxdlnT={dl[k]:+.5f} at cell {k} "
              f"(mass below surf {mbs[k]/MSUN:.4f} Msun, T {new.T[k]:.1f}) top5 {top} effect_solved={d['remesh_effect_solved']} [{time.time()-t0:.1f}s]",flush=True)
        res[(merge,dt)]=dict(T=new.T.copy(),dm=new.dm.copy(),r=new.r.copy(),dl=dl,metrics=m,
                             diag={k:v for k,v in d.items() if isinstance(v,(int,float,str,type(None)))})
import pickle
pickle.dump(res,open(OUT/'v_repro.pkl','wb'))
