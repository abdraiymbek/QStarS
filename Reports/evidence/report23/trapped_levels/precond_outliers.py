import sys, numpy as np
from pathlib import Path
S=Path('/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad')
src=(S/'precond_q.py').read_text().split("cands={")[0]
exec(src)
Mali=mk(dC,lC,uC,allc)
act=it.active; sig=it.sig
W0=it.W.copy(); W1=it.W.copy(); W1[2:31,:]=np.where(act[2:31,:],1.0,0.0)
names={0:'H-',1:'n1',2:'n2',3:'n3',31:'HII',32:'HeI1',101:'HeII1',102:'HeII2',131:'HeIII',132:'ne',133:'T'}
def grp(r):
    return 'H-' if r==0 else 'HI' if r<31 else 'HII' if r==31 else 'HeI' if r<101 else 'HeII' if r<131 else 'HeIII' if r==131 else 'ne' if r==132 else 'T'
qa=[]; 
for i,(zz,j) in enumerate(zip(Zs,JZs)):
    e=W0*(Mali.solve(j)-zz); qa.append(np.linalg.norm(e)/np.linalg.norm(W0*zz))
qa=np.array(qa); order=np.argsort(-qa)
print('q sorted (driver W):', np.round(np.sort(qa)[::-1][:12],3), '... median', np.median(qa))
for i in order[:5]:
    zz,j=Zs[i],JZs[i]; e=W0*(Mali.solve(j)-zz); e2=e**2
    top=np.dstack(np.unravel_index(np.argsort(-e2.ravel())[:6],e2.shape))[0]
    bydepth=e2.sum(0); dtop=np.argsort(-bydepth)[:5]
    print(f'pair {i}: q={qa[i]:.3f}; |Wz| {np.linalg.norm(W0*zz):.3e}; error by depth: '+', '.join(f'd{d+1}:{bydepth[d]/e2.sum():.2f}' for d in dtop)+' | top entries: '+', '.join(f'{names.get(r,grp(r)+str(r+1))}@d{c+1}:{e[r,c]:.2e}(z={zz[r,c]:.1e})' for r,c in top))
# conditioning per depth of the diagonal blocks (equilibrated), frozen vs ALI
def smin(A):
    r=1/np.maximum(np.abs(A).max(1),1e-300); A2=A*r[:,None]; c=1/np.maximum(np.abs(A2).max(0),1e-300); A2=A2*c[None,:]
    sv=np.linalg.svd(A2,compute_uv=False); return sv[-1],sv[-2]
dA=di+dC
actidx=[np.where(it.active[:,k])[0] for k in range(ND)]
print('depth: sigma_min (2nd) of equilibrated active diagonal block, fixed-J | ALI')
for d in (10,18,24,26,27,28,29,30,33,35,37,41,43,45,50,60,80,100,130):
    k=d-1; ix=np.ix_(actidx[k],actidx[k]); a=smin(di[k][ix]); b=smin(dA[k][ix])
    print(f'  d{d:3d}: {a[0]:.2e} ({a[1]:.2e}) | {b[0]:.2e} ({b[1]:.2e})   ratio {b[0]/a[0]:.2e}')
