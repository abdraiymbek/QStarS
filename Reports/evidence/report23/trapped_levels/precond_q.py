import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record, read_ba_blocks, BlockTridiag
out=ROOT/'results/nk167_pin/theta1.00'; base=out/'it05_trial_l1'
ali=ROOT/'results/two_photon_test/pin100_it05_linfac1'
x0,_=read_trace_record(out/'it00_base'/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape; NT,ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=0.0; ck.WEIGHT_RAMP=2.0; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
z=np.load(out/'homotopy_frontpin.npz'); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
cols={n:(out/f'cols_{base.name}'/f'col_{n}_p',out/f'cols_{base.name}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(base,x0,cols,1e-4,True,hom)
lo0,di0,up0=read_ba_blocks(base/'QSTARS_BA_TRACE.bin',NT,ND)
loa,dia,upa=read_ba_blocks(ali/'QSTARS_BA_TRACE.bin',NT,ND)
dC=dia-di0; lC=loa-lo0; uC=upa-up0          # native local-ALI increment (radiation response), all columns
dC[:,133,:]=0; lC[:,133,:]=0; uC[:,133,:]=0   # pinned T row: no RE row
P=np.load(ROOT/'results/two_photon_test/pin100_it05_n2probe/probe.npz'); Zp=P['z']; drad=P['drad'].copy(); deps=list(P['depths'])
drad[133]*=np.asarray(ck.pin_lam(hom))[...,None]
PZ=[Zp[:,:,i] for i in range(21)]; PJ=[it.Jdir.matvec(Zp[:,:,i])+drad[:,:,i] for i in range(21)]
itd=out/f'it05_from_{base.name}'
Zs,JZs=[],[]
for j in sorted(itd.glob('jv*')):
    if (j/'pairs.npz').exists():
        Q=np.load(j/'pairs.npz'); [(Zs.append(Q['z'][:,:,i]),JZs.append(Q['jz'][:,:,i])) for i in range(Q['z'].shape[2])]
lo,di,up=it._pc
def mk(dadd=None,ladd=None,uadd=None,cols_mask=None,rows_mask=None,n2exact=False):
    l,d,u=lo.copy(),di.copy(),up.copy()
    if dadd is not None:
        A=dadd.copy(); L=ladd.copy(); U=uadd.copy()
        if cols_mask is not None:
            A[:,:,~cols_mask]=0; L[:,:,~cols_mask]=0; U[:,:,~cols_mask]=0
        if rows_mask is not None:
            A[:,~rows_mask,:]=0; L[:,~rows_mask,:]=0; U[:,~rows_mask,:]=0
        d+=A; l+=L; u+=U
    if n2exact:
        for i,dd in enumerate(deps): d[dd-1][:,2]=PJ[i][:,dd-1]
    return BlockTridiag(l,d,u).factor()
allc=np.ones(NT,bool); pop=np.zeros(NT,bool); pop[:132]=True; Hc=np.zeros(NT,bool); Hc[0:32]=True
cands={'fixed-J (current M)':mk(),
       'fixed-J + exact n2 col (d26-46)':mk(n2exact=True),
       'ALI increment, all columns':mk(dC,lC,uC,allc),
       'ALI increment, population columns only':mk(dC,lC,uC,pop),
       'ALI increment, H columns only':mk(dC,lC,uC,Hc),
       'ALI increment, H columns, H rows':mk(dC,lC,uC,Hc,Hc),
       'ALI increment x0.5, all columns':mk(0.5*dC,0.5*lC,0.5*uC,allc)}
act=it.active; sig=it.sig
W0=it.W.copy(); W1=it.W.copy(); W1[2:31,:]=np.where(act[2:31,:],1.0,0.0)
def q(M,Z,JZ,W):
    return [float(np.linalg.norm(W*(M.solve(j)-zz))/np.linalg.norm(W*zz)) for zz,j in zip(Z,JZ)]
print(f'{len(Zs)} archived Krylov pairs, {len(PZ)} n2 unit probes; q = ||W(M^-1 J z - z)||/||W z|| (0 = perfect)')
for nm,M in cands.items():
    qa=np.array(q(M,Zs,JZs,W0)); qb=np.array(q(M,Zs,JZs,W1)); qp=np.array(q(M,PZ,PJ,W1))
    y=M.solve(it.R)
    print(f'{nm:42s} Krylov pairs, driver W: median {np.median(qa):8.3f} min {qa.min():7.3f} max {qa.max():9.3f} | H-excited W=1: median {np.median(qb):8.3f} max {qb.max():9.3f} | n2 probes: median {np.median(qp):8.3f} max {qp.max():8.3f} | max|M^-1R| sig {np.abs(y[sig]).max():9.3e}')
