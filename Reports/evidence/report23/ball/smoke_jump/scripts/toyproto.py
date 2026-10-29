import os,sys,time,copy
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/tests')
import numpy as np
import qstars;assert 'patched' in qstars.__file__
from qstars import quasistatic,stepcontrol
import test_quasistatic as T
from qstars.constants import YEAR,MSUN
t=time.time()
s,c=T.toy(128,'face',bh=5.,total=1e4,radius=3000.,eta=.1);s,d=T.relax(s,c)
c=copy.deepcopy(c);c['quasistatic']=dict(base='fixed_mass')
s,d=T.relax(s,c,dt=1e6*YEAR,direct=True)
print('relaxed in',time.time()-t,'s; zones',len(s.dm),'region',d['base_region_cells'])
z=T.zero_time(s);before=stepcontrol.observables(z,T.EOS,'face')
goal=quasistatic.targets(z,z,T.EOS,128);x=z.dm/goal
print('x range outside region',x[d['base_region_cells']:].min(),x[d['base_region_cells']:].max())
for merge in (.6,.7,.8,.9):
  for events in ('measured','counterfactual'):
    cc=copy.deepcopy(c);cc['quasistatic']['merge']=merge;cc['stepcontrol']=dict(metric='mapped',events=events);S=stepcontrol.settings(cc)
    out=[]
    for dt in (1e-6,1e-4,1e-2,1.):
        new,dd=quasistatic.bondi_step(z,dt*YEAR,T.EOS,T.OPACITY,cc,mode='evolve')
        m=stepcontrol.metrics(before,stepcontrol.observables(new,T.EOS,'face'),dd,S)
        out.append((dt,dd['merges'],dd['splits'],round(m['delta_actual'],8),round(dd['mapped_max_dlnT'],8),dd['event_free_mapped_max_dlnT']))
    print(merge,events,out)
