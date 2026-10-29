import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck, cnk_trap
from qstars.nlte_newton import read_trace_record, BlockTridiag
out=ROOT/'results/nk167_trap/fixedT'; base=out/'it00_base'
x0,_=read_trace_record(base/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape; NT,ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=0.0; ck.WEIGHT_RAMP=2.0; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
z=np.load(out/'homotopy_frontpin.npz'); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
cols={n:(out/'cols_it00_base'/f'col_{n}_p',out/'cols_it00_base'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(out/'it00_base_tp',x0,cols,1e-4,True,hom,rj=base/'QSTARS_RJ.bin')
itold=ck.Iterate(base,x0,cols,1e-4,True,hom)
P=np.load(out/'it00_from_it00_base/probe/pairs.npz'); Z=P['z']; JZ=P['jz']
dirs=[Z[:,:,i] for i in range(Z.shape[2])]; jzs=[JZ[:,:,i] for i in range(Z.shape[2])]
_,meta=cnk_trap.probe_directions([2,3],it.active,6)
y0=it.M.solve(it.R)
W0=it.W.copy(); sig0=it.sig.copy()
lo0,di0,up0=[m.copy() for m in it._pc]
it.install_exact_columns(dirs,jzs,meta,6)
y1=it.M.solve(it.R)
# variant: exact column for n2 only
it2=ck.Iterate(out/'it00_base_tp',x0,cols,1e-4,True,hom,rj=base/'QSTARS_RJ.bin')
it2.install_exact_columns(dirs[6:],jzs[6:],meta[6:],6)
y2=it2.M.solve(it2.R)
print('two-photon effect on residual (USE_RAD replay vs OLD_DEFAULT base): n2 row d33/35/37', [(float(f'{itold.R[2,d-1]:.4g}'),float(f'{it.R[2,d-1]:.4g}')) for d in (33,35,37)])
print('T (1e4 K) profile d1..60 step 3:', ' '.join(f'{it.x[133,d]:.3f}' for d in range(0,60,3)))
print('H ionization fraction d1..60 step 3:', ' '.join(f'{it.x[31,d]/(it.x[1:32,d].sum()):.2e}' for d in range(0,60,3)))
print('depth  T    xHII    frac_n2 | M^-1R frozen: HI     n2      HII     ne  | exact n1+n2 cols: HI      n2       HII      ne | exact n2 only: HI   n2    HII')
for d in list(range(1,61,2)):
    k=d-1; xh=it.x[31,k]/it.x[1:32,k].sum()
    print(f'{d:4d} {it.x[133,k]:.3f} {xh:8.1e} {it.frac[2,k]:8.1e} | {y0[1,k]:9.2e} {y0[2,k]:9.2e} {y0[31,k]:9.2e} {y0[132,k]:9.2e} | {y1[1,k]:9.2e} {y1[2,k]:9.2e} {y1[31,k]:9.2e} {y1[132,k]:9.2e} | {y2[1,k]:9.2e} {y2[2,k]:9.2e} {y2[31,k]:9.2e}')
print('by group (fraction-significant only): frozen:', ck.group_report(y0,sig0))
print('   exact n1+n2:', ck.group_report(y1,sig0))
print('   exact n2 only:', ck.group_report(y2,sig0))
# exact vs frozen diagonals of n1 and n2 columns
lo1,di1,up1=it._pc
print('depth: n1 col diag frozen/exact (row 2); n2 col diag frozen/exact (row 3); n1 col entry in row n2 frozen/exact; row2 entries: n2 col frozen/exact')
for d in (5,10,15,18,20,25,30,33,35,37,40,45,50,60,80,100):
    k=d-1
    print(f'{d:4d}: {di0[k][1,1]:11.3e} {di1[k][1,1]:11.3e} | {di0[k][2,2]:11.3e} {di1[k][2,2]:11.3e} | {di0[k][2,1]:11.3e} {di1[k][2,1]:11.3e} | {di0[k][1,2]:11.3e} {di1[k][1,2]:11.3e}')
# conditioning of diagonal blocks restricted to H rows/cols + ne
for d in (10,18,25,33,35,40):
    k=d-1; idx=list(range(0,32))+[132]
    for nm,di in (('frozen',di0),('exact',di1)):
        A=di[k][np.ix_(idx,idx)]; r=1/np.abs(A).max(1); A2=A*r[:,None]; c=1/np.abs(A2).max(0); A2=A2*c[None,:]
        sv=np.linalg.svd(A2,compute_uv=False)
        print(f'  d{d} {nm}: equilibrated H block singular values min {sv[-1]:.2e} second {sv[-2]:.2e} max {sv[0]:.2e}')
np.savez(out/'it00_probe_analysis.npz',y0=y0,y1=y1,y2=y2)
