import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record, BlockTridiag
out=ROOT/'results/nk167_pin/theta1.00'; base=out/'it05_trial_l1'
x0,_=read_trace_record(out/'it00_base'/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=0.0; ck.WEIGHT_RAMP=2.0; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
z=np.load(out/'homotopy_frontpin.npz'); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
cols={n:(out/f'cols_{base.name}'/f'col_{n}_p',out/f'cols_{base.name}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(base,x0,cols,1e-4,True,hom)
P=np.load(ROOT/'results/two_photon_test/pin100_it05_n2probe/probe.npz'); Z=P['z']; drad=P['drad'].copy(); deps=list(P['depths'])
drad[133]*=np.asarray(ck.pin_lam(hom))[...,None]
JZ=np.stack([it.Jdir.matvec(Z[:,:,i])+drad[:,:,i] for i in range(Z.shape[2])],2)   # (NT,ND,21)
print('exact n2 column: own-depth diag, neighbours (row n2) and off-depth l1 share, by source depth')
for i,d in enumerate(deps):
    k=d-1; c=JZ[:,:,i]
    print(f'  d{d}: diag {c[2,k]:10.3e}  n2 row at d-1 {c[2,k-1]:10.3e} d+1 {c[2,k+1]:10.3e} | l1 own {np.abs(c[:,k]).sum():9.3e} d-1 {np.abs(c[:,k-1]).sum():9.3e} d+1 {np.abs(c[:,k+1]).sum():9.3e} rest {np.abs(c).sum()-np.abs(c[:,k-1:k+2]).sum():9.3e}')
for K in (1,2,3,4,6,8):
    errs=[]; errd=[]
    for i,d in enumerate(deps):
        k=d-1
        est=sum(JZ[:,k,j] for j,dd in enumerate(deps) if (dd-d)%K==0)       # response at depth d to all perturbed depths of d's colour
        ex=JZ[:,k,i]
        errs.append(np.abs(est-ex).sum()/np.abs(ex).sum()); errd.append(abs(est[2]-ex[2])/abs(ex[2]))
    print(f'stride {K}: local-column l1 error median {np.median(errs):.3f} max {np.max(errs):.3f} (depth {deps[int(np.argmax(errs))]}); diag error median {np.median(errd):.3f} max {np.max(errd):.3f}')
