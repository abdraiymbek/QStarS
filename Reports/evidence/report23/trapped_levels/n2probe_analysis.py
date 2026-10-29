import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record
out=ROOT/'results/nk167_pin/theta1.00'; base=out/'it05_trial_l1'
x0,_=read_trace_record(out/'it00_base'/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=0.0; ck.WEIGHT_RAMP=2.0; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
z=np.load(out/'homotopy_frontpin.npz'); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
cols={n:(out/f'cols_{base.name}'/f'col_{n}_p',out/f'cols_{base.name}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(base,x0,cols,1e-4,True,hom)
P=np.load(ROOT/'results/two_photon_test/pin100_it05_n2probe/probe.npz'); Z=P['z']; drad=P['drad'].copy(); deps=P['depths']
drad[133]*=np.asarray(ck.pin_lam(hom))[...,None]
JZ=[it.Jdir.matvec(Z[:,:,i])+drad[:,:,i] for i in range(Z.shape[2])]
print('depth  frac(n2)   J_dir[3,3]    J_exact[3,3]   ratio   (M^-1 R)_n2   (M^-1 J e)_n2')
Mb=it.M.solve(it.R)
for i,d in enumerate(deps):
    k=d-1; jd=it.Jdir.matvec(Z[:,:,i])[2,k]; je=JZ[i][2,k]; mi=it.M.solve(JZ[i])[2,k]
    if d%3==2 or d in (33,35,37): print(f'{d:4d} {it.frac[2,k]:9.2e} {jd:12.4e} {je:12.4e} {je/jd:8.2e} {Mb[2,k]:12.4e} {mi:12.4e}')
# least squares: archived it05 pairs + probe columns, different weights for the H excited rows
itd=out/f'it05_from_{base.name}'
Zs,JZs=[],[]
for j in sorted(itd.glob('jv*')):
    if (j/'pairs.npz').exists():
        Q=np.load(j/'pairs.npz'); [(Zs.append(Q['z'][:,:,i]),JZs.append(Q['jz'][:,:,i])) for i in range(Q['z'].shape[2])]
print('archived pairs',len(Zs))
act=it.active
def solve(W,S,JS,label):
    Mi=[it.M.solve(j) for j in JS]
    A=np.stack([(W*m).ravel() for m in Mi],1); b=(W*Mb).ravel()
    y=np.linalg.lstsq(A,b,rcond=1e-12)[0]; s=sum(yi*zz for yi,zz in zip(y,S)); r=(b-A@y)
    rn2=np.linalg.norm((Mb-sum(yi*m for yi,m in zip(y,Mi)))[2,25:46])/np.linalg.norm(Mb[2,25:46])
    print(f'{label:46s} LS res {np.linalg.norm(r)/np.linalg.norm(b):.3f}; n2-row M^-1 residual (d26-46) reduced to {rn2:.3f}; s_n2 d33/35/37 = {s[2,32]:.3g}/{s[2,34]:.3g}/{s[2,36]:.3g}; max|s_sig| {np.abs(s[(it.frac>1e-10)&act]).max():.3f}')
W0=it.W.copy(); W1=it.W.copy(); W1[2:31,:]=np.where(act[2:31,:],1.0,0.0)
solve(W0,Zs,JZs,'driver weights, archived pairs')
solve(W1,Zs,JZs,'H excited rows weight 1, archived pairs')
solve(W0,Zs+[Z[:,:,i] for i in range(21)],JZs+JZ,'driver weights, + n2 probe columns')
solve(W1,Zs+[Z[:,:,i] for i in range(21)],JZs+JZ,'H excited rows weight 1, + n2 probe columns')
