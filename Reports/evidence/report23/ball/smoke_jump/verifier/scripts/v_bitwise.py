"""Run with VSRC=<tree copy>; dumps state arrays + diag + metrics for events measured/counterfactual at step 1255."""
from v_common import *
import pickle
state,config,eos,opacity=load()
tag=sys.argv[1]
out={}
modes=['registered'] if tag=='before' else ['registered','measured','counterfactual']
for ev in modes:
    c=copy.deepcopy(config)
    if ev!='registered':c['stepcontrol']['events']=ev
    for dt in (1.,172.1438701239381):
        new,d,m=attempt(state,dt,eos,opacity,c)
        out[(ev,dt)]=dict(T=new.T,r=new.r,dm=new.dm,L=new.face_luminosity,rc=new.radius_compensation,ledger=dict(new.ledger),
                          aux=json.dumps(new.aux,sort_keys=True),bh=new.bh_mass,cav=new.cavity_mass,rho_ph=new.photosphere_density,
                          diag={k:v for k,v in d.items() if not isinstance(v,(dict,list))},metrics=m,
                          event_free=d.get('event_free'))
        print(tag,ev,dt,'delta',m['delta_actual'],'profile_dlnT',m['profile_dlnT'],'rejected',sc.rejected(m['delta_actual'],sc.settings(c)),
              'measured mapped_dlnT',d['mapped_max_dlnT'],'event_free_mapped',d.get('event_free_mapped_max_dlnT'),flush=True)
pickle.dump(out,open(OUT/f'v_bitwise_{tag}.pkl','wb'))
