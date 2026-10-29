#!/usr/bin/env python3
"""Report 21: least-squares relative residual of the driver objective at iterate 5 of results/nk_adv84 with the
driver's 6 exact products, the 10 exact T/n_e probe products (scripts/probe_columns.py), and both."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0,'src'); sys.path.insert(0,'scripts')
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record, read_jv_rad
out=Path('results/nk_adv84'); base=out/'it05_trial_l1'; itd=out/'it05_from_it05_trial_l1'
x0,_=read_trace_record(out/'it00_base'/'QSTARS_RESID_TRACE.bin')
ck.NT,ck.ND=x0.shape; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=1e-8; ck.PC_T_FLOOR=1e-3
z=np.load(out/'homotopy_shift.npz')
hom=dict(Tfix=z['Tfix'],c=z['c'],R0=z['R0'],mode='shift',lam=0.0)
cols={n:(out/f'cols_{base.name}'/f'col_{n}_p',out/f'cols_{base.name}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(base,x0,cols,1e-4,True,hom)
A=it.active
def ob(r): return (it.W*it.M.solve(r))[A]
b=ob(it.R)
Z=[];JZ=[]
for j in sorted(itd.glob('jv*')):
    P=np.load(j/'pairs.npz')
    for i in range(P['z'].shape[2]): Z.append(P['z'][:,:,i]); JZ.append(P['jz'][:,:,i])
nd=len(Z)
# probe directions
sys.path.insert(0,'scripts'); import probe_columns as pc
zs,tags=pc.directions(['T','ne'],5,ck.NT,ck.ND)
_,rad=read_jv_rad(out/'probe_it05_trial_l1'/'pass'/'QSTARS_JV_RAD.bin')
for i,zz in enumerate(zs): Z.append(zz); JZ.append(it.Jdir.matvec(zz)+rad[:,:,i])
def solve(idx):
    Amat=np.stack([ob(JZ[i]) for i in idx],1)
    y,*_=np.linalg.lstsq(Amat,b,rcond=None)
    return np.linalg.norm(b-Amat@y)/np.linalg.norm(b), y
print('driver dirs (6):', solve(range(nd))[0])
print('probe only (10):', solve(range(nd,len(Z)))[0])
print('T probes only + driver:', solve(list(range(nd))+list(range(nd,nd+5)))[0])
r,y=solve(range(len(Z)))
print('all 16:', r)
s=sum(y[k]*Z[i] for k,i in enumerate(range(len(Z))))
print('step max|s_T|', np.abs(s[133]).max(), 'max|s|', np.abs(s).max())
