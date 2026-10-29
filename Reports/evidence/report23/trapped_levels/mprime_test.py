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
P=np.load(ROOT/'results/two_photon_test/pin100_it05_n2probe/probe.npz'); Z=P['z']; drad=P['drad'].copy(); deps=P['depths']
drad[133]*=np.asarray(ck.pin_lam(hom))[...,None]
JZ=[it.Jdir.matvec(Z[:,:,i])+drad[:,:,i] for i in range(Z.shape[2])]
di2=it.M.di.copy()
for i,d in enumerate(deps): di2[d-1][:,2]=JZ[i][:,d-1]
M2=BlockTridiag(it.M.lo,di2,it.M.up).factor()
y=it.M.solve(it.R); y2=M2.solve(it.R)
print('depth   T(1e4K)  frac_n2   M^-1R: n2      M\'^-1R: n2   n3      n4      HI(row2)   HII      ne     | old: HI(row2) HII ne')
for d in range(26,47,1):
    k=d-1
    print(f'{d:4d} {it.x[133,k]:8.4f} {it.frac[2,k]:9.2e} {y[2,k]:11.3e} {y2[2,k]:11.3e} {y2[3,k]:8.2e} {y2[4,k]:8.2e} {y2[1,k]:10.3e} {y2[31,k]:9.2e} {y2[132,k]:9.2e} | {y[1,k]:9.2e} {y[31,k]:9.2e} {y[132,k]:9.2e}')
act=it.active; sig=(it.frac>1e-10)&act
print('max|M^-1R| sig old', np.abs(y[sig]).max(), ' new', np.abs(y2[sig]).max())
# coupling significance: column l at depth d -> largest row-scaled entry into significant rows other than itself
D=it.D
print('coupling significance c_l(d) = max_{i sig, i!=l} |M\'[i,l]| D_i  for H I levels 2..10 at d33/35/37/41 (exact n2 column)')
for d in (30,33,35,37,41,45,60,100):
    k=d-1; rowsig=np.where(sig[:,k])[0]; outl=[]
    for l in range(2,12):
        ii=[i for i in rowsig if i!=l]
        outl.append(float(np.max(np.abs(di2[k][ii,l])*D[ii,k])))
    print(f'  d{d:3d} frac n2 {it.frac[2,k]:.1e}: '+' '.join(f'{v:8.1e}' for v in outl))
# LS tests
itd=out/f'it05_from_{base.name}'
Zs,JZs=[],[]
for j in sorted(itd.glob('jv*')):
    if (j/'pairs.npz').exists():
        Q=np.load(j/'pairs.npz'); [(Zs.append(Q['z'][:,:,i]),JZs.append(Q['jz'][:,:,i])) for i in range(Q['z'].shape[2])]
print('archived pairs',len(Zs))
def solve(Msol,W,S,JS,label):
    Mb=Msol(it.R); Mi=[Msol(j) for j in JS]
    A=np.stack([(W*m).ravel() for m in Mi],1); b=(W*Mb).ravel()
    yv=np.linalg.lstsq(A,b,rcond=1e-12)[0]; s=sum(yi*zz for yi,zz in zip(yv,S)); r=(b-A@yv)
    res=Mb-sum(yi*m for yi,m in zip(yv,Mi))
    rn2=np.linalg.norm(res[2,25:46])/np.linalg.norm(Mb[2,25:46])
    print(f'{label:52s} |b| {np.linalg.norm(b):9.3e} LS res {np.linalg.norm(r)/np.linalg.norm(b):.3f}; n2-row resid -> {rn2:.3f}')
    print('      s_n2 d31/33/35/37/39/41 = '+'/'.join(f'{s[2,d-1]:.3g}' for d in (31,33,35,37,39,41))+f';  s_HII = '+'/'.join(f'{s[31,d-1]:.3g}' for d in (31,33,35,37,39,41))+'; s_ne = '+'/'.join(f'{s[132,d-1]:.3g}' for d in (31,33,35,37,39,41))+f'; max|s_sig| {np.abs(s[sig]).max():.3f}')
    return s
W0=it.W.copy(); W1=it.W.copy(); W1[2:31,:]=np.where(act[2:31,:],1.0,0.0)
PZ=[Z[:,:,i] for i in range(21)]
solve(it.M.solve,W0,Zs+PZ,JZs+JZ,'old M, driver weights, pairs+probes')
solve(it.M.solve,W1,Zs+PZ,JZs+JZ,'old M, H-excited weight 1, pairs+probes')
solve(M2.solve,W0,Zs+PZ,JZs+JZ,"M' (exact n2 col), driver weights, pairs+probes")
s=solve(M2.solve,W1,Zs+PZ,JZs+JZ,"M' (exact n2 col), H-excited weight 1, pairs+probes")
solve(M2.solve,W1,PZ,JZ,"M', H-excited weight 1, probes only")
# add the M'-preconditioned residual as a direction? needs a product: not available offline. Show its content instead.
b2=W1*y2
print('norm of W1 M\'^-1 R by row group: H-:%.3g HI ion row:%.3g n2:%.3g n3-30:%.3g HII:%.3g He:%.3g ne:%.3g T:%.3g'%(np.linalg.norm(b2[0]),np.linalg.norm(b2[1]),np.linalg.norm(b2[2]),np.linalg.norm(b2[3:31]),np.linalg.norm(b2[31]),np.linalg.norm(b2[32:132]),np.linalg.norm(b2[132]),np.linalg.norm(b2[133])))
