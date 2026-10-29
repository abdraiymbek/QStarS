"""Scratch: does the PTC iteration itself converge where Newton does (toy)?"""
import sys,copy,time
sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/src');sys.path.insert(0,'/Users/abraiymbek/Projects/QStarS/tests')
import numpy as np
from test_quasistatic import toy,relax,zero_time,EOS,OPACITY
from qstars import quasistatic,cavity
from qstars.constants import YEAR
s,c=toy(128,'face',bh=5.,total=1e4,radius=3000.,eta=.1)
s,_=relax(s,c);s=zero_time(s)
S=quasistatic.settings(c)
for dt in (.1,1.,3.):
    dtt=dt*YEAR
    pred=quasistatic.predict(s,dtt,EOS,c,'evolve',None,S)
    cut,_=cavity.slice_base(s,pred['swept'],EOS)
    post,_=quasistatic.remesh_base(s,cut,EOS,c['numerics']['zones'])
    P=quasistatic.Problem(s,post,dtt,EOS,OPACITY,c,'evolve',None,pred,S)
    xn,info=P.solve()
    f0=P.evaluate(P.x0)
    print('dt',dt,'newton',info['success'],info['iterations'],info['residual'],'|F0|max',np.max(abs(f0)),'argmax',int(np.argmax(abs(f0))),'n',P.n)
    for d0 in (1.,):
        hist=[]
        orig=P.evaluate
        xp,ptc=quasistatic.pseudo_transient(P.evaluate,P.x0,P.pattern,quasistatic.pairing(P.n),1e-12,**quasistatic.PTC)
        print('   ptc',{k:v for k,v in ptc.items()})
        if ptc['handover']:
            x2,fin=quasistatic.newton(P.evaluate,xp,P.pattern,1e-12)
            print('   newton after ptc',fin['success'],fin['iterations'],fin['residual'],'|x2-xn|',np.max(abs(x2-xn)) if info['success'] else None)
dtt=.1*YEAR
pred=quasistatic.predict(s,dtt,EOS,c,'evolve',None,S)
cut,_=cavity.slice_base(s,pred['swept'],EOS)
post,_=quasistatic.remesh_base(s,cut,EOS,c['numerics']['zones'])
P=quasistatic.Problem(s,post,dtt,EOS,OPACITY,c,'evolve',None,pred,S)
xn,info=P.solve()
xp,ptc=quasistatic.pseudo_transient(P.evaluate,P.x0,P.pattern,quasistatic.pairing(P.n),1e-12,**quasistatic.PTC)
x2,fin=quasistatic.newton(P.evaluate,xp,P.pattern,1e-12)
n=P.n;d=x2-xn
names=['y']*(n+1)+['th']*n+['ell']*(n-1)+['lnrho_ph','dhat','X']
idx=np.argsort(-abs(d))[:12]
print([(int(i),names[i],float(d[i]),float(xn[i]),float(x2[i])) for i in idx])
print('max |d| by block: y',np.max(abs(d[:n+1])),'th',np.max(abs(d[n+1:2*n+1])),'ell',np.max(abs(d[2*n+1:3*n])),'ph',d[3*n],'dhat',d[3*n+1],'X',d[3*n+2])
print('x0 ell range',P.x0[2*n+1:3*n].min(),P.x0[2*n+1:3*n].max(),'xn ell',xn[2*n+1:3*n].min(),xn[2*n+1:3*n].max(),'x2 ell',x2[2*n+1:3*n].min(),x2[2*n+1:3*n].max())
