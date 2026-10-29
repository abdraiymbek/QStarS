import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck, cnk_trap
from qstars.nlte_newton import read_trace_record
out=ROOT/'results/nk167_trap/fixedT'; base=out/'it00_base'
x0,_=read_trace_record(base/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape; NT,ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=0.0; ck.WEIGHT_RAMP=2.0; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
z=np.load(out/'homotopy_frontpin.npz'); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
cols={n:(out/'cols_it00_base'/f'col_{n}_p',out/'cols_it00_base'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(out/'it00_base_tp',x0,cols,1e-4,True,hom,rj=base/'QSTARS_RJ.bin')
P=np.load(out/'it00_from_it00_base/probe/pairs.npz'); Z=P['z']; JZ=P['jz']
lo0,di0,up0=[m.copy() for m in it._pc]
dirs=[Z[:,:,i] for i in range(12)]; jzs=[JZ[:,:,i] for i in range(12)]
_,meta=cnk_trap.probe_directions([2,3],it.active,6)
it.install_exact_columns(dirs,jzs,meta,6); lo1,di1,up1=it._pc
names={0:'H-',1:'n1',2:'n2',3:'n3',4:'n4',31:'HII',32:'HeI1',101:'HeII1',132:'ne',133:'T'}
for d in (18,35,59,90):
    k=d-1
    print(f'=== depth {d}: T={it.x[133,k]*10:.2f} kK  n1={it.x[1,k]:.3e} n2={it.x[2,k]:.3e} nHII={it.x[31,k]:.3e} ne={it.x[132,k]:.3e}')
    print('  residuals: H- %.3e | row2 %.3e | n2 %.3e | n3 %.3e | sum rows 2..31 %.3e | HII row %.3e | ne row %.3e'%(it.R[0,k],it.R[1,k],it.R[2,k],it.R[3,k],it.R[1:31,k].sum(),it.R[31,k],it.R[132,k]))
    for r,nm in ((1,'row 2 (H I ground eq)'),(2,'row 3 (n2)'),(31,'row 32 (H II)')):
        fr=di0[k][r]; ex=di1[k][r]
        idx=np.argsort(-np.abs(fr))[:7]
        print(f'  {nm}: largest fixed-J entries: '+', '.join(f'{names.get(c,"c%d"%(c+1))}:{fr[c]:.3e}' for c in idx))
        idx=np.argsort(-np.abs(ex))[:7]
        print(f'  {nm}: with exact n1,n2 cols : '+', '.join(f'{names.get(c,"c%d"%(c+1))}:{ex[c]:.3e}' for c in idx))
    # exact n1 column: all rows
    c1=di1[k][:,1]; idx=np.argsort(-np.abs(c1))[:8]
    print('  exact n1 column, largest: '+', '.join(f'r{names.get(r,str(r+1))}:{c1[r]:.3e}' for r in idx))
    c0=di0[k][:,1]; idx=np.argsort(-np.abs(c0))[:8]
    print('  fixed-J n1 column, largest: '+', '.join(f'r{names.get(r,str(r+1))}:{c0[r]:.3e}' for r in idx))
    print('  neighbours (exact n1 col, row 2): lo[d+1] %.3e  up[d-1] %.3e ; diag %.3e'%(lo1[k+1][1,1],up1[k-1][1,1],di1[k][1,1]))
