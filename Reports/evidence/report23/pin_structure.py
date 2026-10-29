import sys, json, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record
out=Path(sys.argv[1]); base=Path(sys.argv[2]); k=int(sys.argv[3])
x0,_=read_trace_record(out/'it00_base'/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=1e-12; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
z=np.load(out/'homotopy_frontpin.npz'); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
cols={n:(out/f'cols_{base.name}'/f'col_{n}_p',out/f'cols_{base.name}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(base,x0,cols,1e-4,True,hom)
itd=out/f'it{k:02d}_from_{base.name}'
Z,JZ=[],[]
for j in sorted(itd.glob('jv*')):
    P=np.load(j/'pairs.npz'); [ (Z.append(P['z'][:,:,i]), JZ.append(P['jz'][:,:,i])) for i in range(P['z'].shape[2])]
act=it.active; W=it.W
Mb=it.M.solve(it.R); Mi=[it.M.solve(j) for j in JZ]
A=np.stack([(W*m).ravel() for m in Mi],1); b=(W*Mb).ravel()
y=np.linalg.lstsq(A,b,rcond=1e-12)[0]; r=(b-A@y).reshape(Mb.shape); b2=b.reshape(Mb.shape)
tot=np.sum(r**2); tot0=np.sum(b2**2)
print(f'LS residual {np.sqrt(tot/tot0):.4f}, {len(Z)} products')
ultra=(it.frac<=1e-20)&act; trace=(it.frac>1e-20)&(it.frac<=1e-10)&act; sig=(it.frac>1e-10)&act
rows=np.zeros_like(act); rows[132]=True; rT=np.zeros_like(act); rT[133]=True
for name,m in [('ultra',ultra&~rows&~rT),('trace',trace&~rows&~rT),('sig pops',sig&~rows&~rT),('n_e',rows&act),('T',rT&act)]:
    print(f'  {name:9s} final share {np.sum(r[m]**2)/tot:6.3f} initial {np.sum(b2[m]**2)/tot0:6.3f} reduced to {np.sqrt(np.sum(r[m]**2)/max(np.sum(b2[m]**2),1e-300)):.3f}')
for g,sl in {'H(1-32)':slice(0,32),'HeI(33-101)':slice(32,101),'HeII(102-131)':slice(101,131),'HeIII':slice(131,132)}.items():
    m=np.zeros_like(act); m[sl]=True; print(f'  {g:14s} final share {np.sum(r[m&act]**2)/tot:6.3f}')
prof=np.array([np.sum(r[:,d][act[:,d]]**2) for d in range(r.shape[1])])/tot
print('  top depths',[(int(d)+1,round(float(prof[d]),3)) for d in np.argsort(-prof)[:10]])
rr=(r**2)*act; i=np.unravel_index(np.argmax(rr),rr.shape); print('  largest single entry row',i[0]+1,'depth',i[1]+1,'frac %.2e'%it.frac[i])
