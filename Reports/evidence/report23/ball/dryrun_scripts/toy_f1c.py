"""Scratch: PTC variants on the toy: full pairing vs differential rows only (first law + growth)."""
import sys,copy,time
sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/src');sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/tests')
import numpy as np
from test_quasistatic import toy,relax,zero_time,EOS,OPACITY
from qstars import quasistatic,cavity
from qstars.constants import YEAR
s,c=toy(128,'face',bh=5.,total=1e4,radius=3000.,eta=.1)
s,_=relax(s,c);s=zero_time(s)
S=quasistatic.settings(c)
def pairs(n,kind):
    rows,cols=quasistatic.pairing(n)
    if kind=='full':return rows,cols
    keep=np.r_[n+np.arange(n),3*n+2]   # first-law rows and the growth row
    return rows[keep],cols[keep]
for dt in (.1,1.,3.,10.):
    dtt=dt*YEAR
    pred=quasistatic.predict(s,dtt,EOS,c,'evolve',None,S)
    cut,_=cavity.slice_base(s,pred['swept'],EOS)
    post,_=quasistatic.remesh_base(s,cut,EOS,c['numerics']['zones'])
    P=quasistatic.Problem(s,post,dtt,EOS,OPACITY,c,'evolve',None,pred,S)
    xn,info=P.solve()
    print('dt',dt,'newton',info['success'],info['iterations'],'%.3g'%info['residual'])
    for kind in ('full','differential'):
        for d0 in (1.,1e-2):
            xp,ptc=quasistatic.pseudo_transient(P.evaluate,P.x0,P.pattern,pairs(P.n,kind),1e-12,delta0=d0,max_iterations=60,handover=1e6,max_shrinks=20)
            msg='   %s d0=%g: handover=%s it=%d shrinks=%d res=%.3g delta=%.3g'%(kind,d0,ptc['handover'],ptc['iterations'],ptc['shrinks'],ptc['residual'],ptc['delta'])
            if ptc['handover']:
                x2,fin=quasistatic.newton(P.evaluate,xp,P.pattern,1e-12)
                msg+=' | newton %s it=%d res=%.3g'%(fin['success'],fin['iterations'],fin['residual'])+(' |x2-xn|=%.3g'%np.max(abs(x2-xn)) if info['success'] and fin['success'] else '')
            print(msg)
