"""Scratch: F1 (PTC) on the analytic Ball row-0 toy (not a Ball configuration run)."""
import sys,copy,time
sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/src');sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/tests')
import numpy as np
from test_quasistatic import toy,relax,zero_time,EOS,OPACITY
from qstars import quasistatic
from qstars.constants import YEAR
from qstars.numerics import StepRejected
s,c=toy(128,'face',bh=5.,total=1e4,radius=3000.,eta=.1)
s,_=relax(s,c);s=zero_time(s)
for dt in (1.,3.,10.):
    for ptc in (False,True):
        cc=copy.deepcopy(c);cc['quasistatic']=dict(ptc=ptc)
        t=time.perf_counter()
        try:
            new,d=quasistatic.bondi_step(s,dt*YEAR,EOS,OPACITY,cc)
            print('dt',dt,'ptc',ptc,'ACCEPTED it',d['newton_iterations'],'res',d['nonlinear_residual'],'passes',d['passes'],{k:v for k,v in d.items() if k.startswith('ptc')},'%.1fs'%(time.perf_counter()-t))
        except StepRejected as e:
            print('dt',dt,'ptc',ptc,'REJECTED',str(e)[:300],'%.1fs'%(time.perf_counter()-t))
