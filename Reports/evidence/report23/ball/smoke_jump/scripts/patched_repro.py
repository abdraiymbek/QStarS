import os,sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'patched','src'))
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'patched','scripts'))
import copy,json,numpy as np
import qstars;assert 'patched' in qstars.__file__,qstars.__file__
from qstars import quasistatic as qs,stepcontrol as sc
import validate_ball_track as vbt
RUN='/Users/abraiymbek/Projects/QStarS/benchmarks/results/2026-10-14-ball-track/20261014T051036.115106Z_smoke_C1_n200'
YEAR=365.25*86400
state,config=qs.load_checkpoint(RUN+'/checkpoint.npz')
eos,opacity=vbt.physics(config)
rows=[]
for events in ('measured','counterfactual'):
    c=copy.deepcopy(config);c['stepcontrol']['events']=events;S=sc.settings(c)
    before=sc.observables(state,eos,'face')
    for dt in (1e-3,1.,10.,100.,172.1438701239381,300.):
        new,d=qs.bondi_step(state,dt*YEAR,eos,opacity,c,mode='evolve')
        m=sc.metrics(before,sc.observables(new,eos,'face'),d,S)
        top=max((v,k) for k,v in m.items() if k!='delta_actual')
        rows.append(dict(events=events,dt_years=dt,delta_actual=m['delta_actual'],limiting=top[1],rejected=sc.rejected(m['delta_actual'],S),
                         measured_mapped_dlnT=d['mapped_max_dlnT'],event_free_mapped_dlnT=d['event_free_mapped_max_dlnT'],merges=d['merges'],
                         dln_luminosity=m['dln_luminosity'],dln_radius=m['dln_radius'],dln_teff=m['dln_effective_temperature']))
        print(json.dumps(rows[-1]))
json.dump(rows,open('patched_repro.json','w'),indent=1)
