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
print('Z nonzeros per dir', [(int(np.count_nonzero(Z[:,:,i])), float(Z[2,deps[i]-1,i])) for i in (0,10,20)])
print('x rows at d35: H-, HI(1..5), HII, ne, T:', it.x[0,34], it.x[1:6,34], it.x[31,34], it.x[132,34], it.x[133,34])
print('R rows at d35 (first 8 H rows, ion=row2):', it.R[0:8,34])
names={0:'H-',1:'HI eq/ion row',2:'n2',3:'n3',4:'n4',5:'n5',31:'HII',32:'HeI 1',132:'ne',133:'T'}
for i,d in enumerate(deps):
    if d not in (33,35,37): continue
    k=d-1
    jz=it.Jdir.matvec(Z[:,:,i])+drad[:,:,i]; jd=it.Jdir.matvec(Z[:,:,i])
    print(f'--- depth {d}: exact column of n2 (largest |entries|), frozen-J value alongside')
    idx=np.dstack(np.unravel_index(np.argsort(-np.abs(jz).ravel())[:14], jz.shape))[0]
    for r,c in idx:
        print(f'   row {r+1:4d} ({names.get(r,"")}) depth {c+1:4d}: exact {jz[r,c]:12.4e}   frozen {jd[r,c]:12.4e}')
    loc=np.abs(jz[:,k]).sum(); tot=np.abs(jz).sum(); nb=np.abs(jz[:,k-1:k+2]).sum()
    print(f'   l1 mass at own depth {loc/tot:.4f}, within +-1 {nb/tot:.4f}; frozen: own depth {np.abs(jd[:,k]).sum()/np.abs(jd).sum():.4f}')
    # relative to row scale (D): which rows respond most in relative terms
    rel=np.abs(jz)*it.D
    idx=np.dstack(np.unravel_index(np.argsort(-rel.ravel())[:8], rel.shape))[0]
    print('   largest row-scaled entries:', [(int(r+1),int(c+1),float(f'{rel[r,c]:.3g}')) for r,c in idx])
