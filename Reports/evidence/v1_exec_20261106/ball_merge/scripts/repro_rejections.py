"""Goal 1: reproduce the eight registered rejections at step 1255 from the preserved checkpoint, unpatched code,
registered configuration (as stored in the checkpoint). Compares delta_actual with the run's events.jsonl."""
from bm_common import *
import re
state,config,eos,opacity=load()
print('checkpoint steps',state.steps,'M_BH',state.bh_mass/MSUN,'zones',len(state.dm),'controller',state.aux.get('controller'))
ev=[json.loads(l) for l in open(RUN/'events.jsonl')]
rej=[e for e in ev if e['event']=='rejected']
rows=[]
for e in rej:
    logged=float(re.search(r'delta_actual ([0-9.eE+-]+)',e['message']).group(1))
    t0=time.time()
    new,d,m=attempt(state,e['dt_years'],eos,opacity,config)
    k=int(np.argmax(abs(lnT_at(new,mbs(state))-np.log(state.T))))
    row=dict(dt_years=e['dt_years'],logged=logged,reproduced=m['delta_actual'],component=max((v,kk) for kk,v in m.items() if kk!='delta_actual')[1],
             rejected=bool(sc.rejected(m['delta_actual'],sc.settings(config))),zones=d['zones'],merges=d['merges'],splits=d['splits'],
             newton_iterations=d['newton_iterations'],max_change_cell=k,mass_below_surface_msun=float(mbs(state)[k]/MSUN),seconds=time.time()-t0)
    rows.append(row);print(row,flush=True)
assert all(abs(r['reproduced']-r['logged'])<=5e-6 for r in rows),'mismatch with the logged rejections'
# where the merge sits and what the front looks like
q=eos.evaluate(state.rho,state.T);kap=opacity(state.rho,state.T)
from qstars.constants import G,C,PI
Lf=state.face_luminosity;m_enc=state.enclosed
outer=[]
for k in range(len(state.dm)-14,len(state.dm)):
    gam=float(kap[k]*0.5*(Lf[k]+Lf[k+1])/(4*PI*C*G*0.5*(m_enc[k]+m_enc[k+1])))
    outer.append(dict(cell=k,mass_below_surface_msun=float(mbs(state)[k]/MSUN),T=float(state.T[k]),rho=float(state.rho[k]),kappa=float(kap[k]),
                      dlnkappa_to_next=float(np.log(kap[k+1]/kap[k])) if k+1<len(kap) else None,
                      dlnT_to_next=float(np.log(state.T[k+1]/state.T[k])) if k+1<len(kap) else None,eddington_factor=gam))
for o in outer:print({k:(round(v,5) if isinstance(v,float) else v) for k,v in o.items()})
jdump(dict(checkpoint=str(RUN/'checkpoint.npz'),steps=state.steps,bh_mass_msun=state.bh_mass/MSUN,zones=len(state.dm),
           rejections=rows,outer_cells=outer),DATA/'repro_rejections.json')
