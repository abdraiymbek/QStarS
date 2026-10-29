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
dirs=[Z[:,:,i] for i in range(Z.shape[2])]; jzs=[it.Jdir.matvec(d_)+drad[:,:,i] for i,d_ in enumerate(dirs)]
sig0=it.sig.copy(); W0=it.W.copy()
print('columns installed', it.install_exact_columns(dirs,jzs,[(3,0)]*len(dirs),4))
y=it.M.solve(it.R)
for thr in (1e-2,1e-3):
    for rounds in (1,8):
        it.W=W0.copy()
        n=it.extend_significance(thr,rounds)
        add=it.sig&~sig0
        grp={'H-':(0,1),'HI':(1,31),'HII':(31,32),'HeI':(32,101),'HeII':(101,131),'HeIII':(131,132)}
        print(f'--- threshold {thr:g}, rounds {rounds}: added {n}; by group: '+', '.join(f'{g} {int(add[a:b].sum())}' for g,(a,b) in grp.items()))
        for d0,d1 in ((1,25),(26,46),(47,80),(81,120),(121,167)):
            a_=add[:,d0-1:d1]
            if a_.any():
                yy=np.abs(y[:,d0-1:d1][a_]); fr=it.frac[:,d0-1:d1][a_]
                print(f'    depths {d0:3d}-{d1:3d}: {int(a_.sum()):5d} levels; frac min {fr.min():.1e} median {np.median(fr):.1e}; |M\'^-1R| median {np.median(yy):.2e} max {yy.max():.2e}')
        print(f'    max|M\'^-1R| over old sig {np.abs(y[sig0]).max():.3e}; over new sig {np.abs(y[it.sig]).max():.3e}; weighted norm old W {np.linalg.norm(W0*y):.3e} new W {np.linalg.norm(it.W*y):.3e}')
        # where are the largest flagged values
        yy=np.where(add,np.abs(y),0); idx=np.dstack(np.unravel_index(np.argsort(-yy.ravel())[:6],yy.shape))[0]
        print('    largest flagged:', [(int(r+1),int(c+1),float(f'{yy[r,c]:.3g}'),float(f'{it.frac[r,c]:.2g}')) for r,c in idx])
