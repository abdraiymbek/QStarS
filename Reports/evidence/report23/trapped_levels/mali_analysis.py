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
xa,Ra=read_trace_record(ali/'QSTARS_RESID_TRACE.bin')
print('state identical', np.array_equal(xa,it.x), ' residual max rel diff', np.max(np.abs(Ra-it.R_native)/(np.abs(it.R_native)+1e-300)))
lo0,di0,up0=read_ba_blocks(base/'QSTARS_BA_TRACE.bin',NT,ND)
loa,dia,upa=read_ba_blocks(ali/'QSTARS_BA_TRACE.bin',NT,ND)
P=np.load(ROOT/'results/two_photon_test/pin100_it05_n2probe/probe.npz'); Z=P['z']; drad=P['drad'].copy(); deps=list(P['depths'])
drad[133]*=np.asarray(ck.pin_lam(hom))[...,None]
JZ=np.stack([it.Jdir.matvec(Z[:,:,i])+drad[:,:,i] for i in range(Z.shape[2])],2)
print('n2 column: depth, frozen diag, ALI diag, exact diag | n2-col l1 error of ALI local column vs exact (rows 1..132)')
for i,d in enumerate(deps):
    k=d-1; ex=JZ[:132,k,i]; al=dia[k][:132,2]; fr=di0[k][:132,2]
    print(f' d{d}: {fr[2]:11.3e} {al[2]:11.3e} {ex[2]:11.3e} | ALI err {np.abs(al-ex).sum()/np.abs(ex).sum():8.3f}  frozen err {np.abs(fr-ex).sum()/np.abs(ex).sum():9.1f}')
# all levels: ratio of ALI to frozen diagonal
ratio=np.zeros((132,ND))
for k in range(ND):
    f=np.diag(di0[k])[:132]; a_=np.diag(dia[k])[:132]
    ratio[:,k]=np.where(f!=0,a_/np.where(f!=0,f,1),1.0)
act=it.active[:132]
print('levels whose ALI diagonal is < 0.5 x the frozen diagonal somewhere (active rows): row, #depths, min ratio, depth of min, depth range')
for r in range(132):
    m=(ratio[r]<0.5)&act[r]
    if m.any():
        kk=np.where(m)[0]; j=kk[np.argmin(ratio[r][kk])]
        print(f'  row {r+1:4d}: {int(m.sum()):4d} depths, min ratio {ratio[r,j]:9.2e} at d{j+1}, range d{kk.min()+1}-d{kk.max()+1}, frac there {it.frac[r,j]:.1e}')
m=(ratio>2)&act
print('levels with ALI diag > 2x frozen:', int(m.sum()))
# preconditioner test with ALI blocks for population columns
lo,di,up=[m_.copy() for m_ in it._pc]
for k in range(ND):
    for c in list(range(0,132)):
        if c==31 or c>=101: continue      # keep the completed ion columns (POPION term) as in M
        di[k][:,c]=dia[k][:,c]; lo[k][:,c]=loa[k][:,c]; up[k][:,c]=upa[k][:,c]
Ma=BlockTridiag(lo,di,up).factor(); ya=Ma.solve(it.R)
lo2,di2,up2=[m_.copy() for m_ in it._pc]
for i,d in enumerate(deps): di2[d-1][:,2]=JZ[:,d-1,i]
M2=BlockTridiag(lo2,di2,up2).factor(); y2=M2.solve(it.R); y0=it.M.solve(it.R)
print('depth   M^-1R n2: frozen     exact-n2-col    ALI-pop-cols |  HII: frozen  exact-n2  ALI | HI(row2): frozen exact-n2 ALI')
for d in (28,30,31,32,33,34,35,36,37,38,39,40,41,43,45):
    k=d-1
    print(f'{d:4d} {y0[2,k]:11.3e} {y2[2,k]:11.3e} {ya[2,k]:11.3e} | {y0[31,k]:10.2e} {y2[31,k]:10.2e} {ya[31,k]:10.2e} | {y0[1,k]:10.2e} {y2[1,k]:10.2e} {ya[1,k]:10.2e}')
sig=it.sig
print('max|M^-1R| over fraction-significant: frozen %.3e exact-n2 %.3e ALI %.3e'%(np.abs(y0[sig]).max(),np.abs(y2[sig]).max(),np.abs(ya[sig]).max()))
print('ALI by group:', ck.group_report(ya, sig))
np.savez(ROOT/'results/two_photon_test/pin100_it05_linfac1_analysis.npz', ratio=ratio, ya=ya, y2=y2, y0=y0)
