from common import *
state,config,eos,opacity=load()
print('steps',state.steps,'t_yr',state.time/YEAR,'Mbh',state.bh_mass/1.98847e33,'zones',len(state.dm),'dt_next',state.aux['controller'])
n=len(state.dm)
res={}
for dt in [1e-3,1e-2,0.1,1,10,100,172.1438701239381]:
    try:
        new,d,m=step(state,dt,eos,opacity,config)
    except Exception as e:
        print(dt,'FAIL',type(e).__name__,e);continue
    th=np.log(new.T/state.T) if len(new.T)==n else None
    # which cell carries the mapped dlnT
    print(f"dt={dt:9.3g} yr  profile_dlnT={m['profile_dlnT']:.5g} profile_dlnr={m['profile_dlnr']:.4g} lag_dlnT={d['lagrangian_max_dlnT']:.4g} "
          f"delta={m['delta_actual']:.4g} it={d['newton_iterations']} zones={d['zones']} splits={d['splits']} merges={d['merges']} region={d['base_region_cells']} "
          f"dlnL={m['dln_luminosity']:.3g} dlnR={m['dln_radius']:.3g} dlnTeff={m['dln_effective_temperature']:.3g}")
    res[dt]=dict(T=new.T,r=new.r,L=new.face_luminosity,dm=new.dm,rho_ph=new.photosphere_density)
np.savez(SCR/'repro_states.npz',**{f'{k}_{dt}':v for dt,r in res.items() for k,v in r.items() if v is not None})
