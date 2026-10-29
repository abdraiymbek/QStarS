import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record
out=Path(sys.argv[1]); base=Path(sys.argv[2]); uw=float(sys.argv[3])
x0,_=read_trace_record(out/'it00_base'/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape
ck.SIG_FRAC=1e-10; ck.TRACE_WEIGHT=1e-3; ck.ULTRA_FRAC=1e-20; ck.ULTRA_WEIGHT=uw; ck.MERIT_MODE='precond'; ck.EQ_FLOOR=0.0; ck.PC_T_FLOOR=1e-3
ck.SPLIT_REGIONS=ck.split_regions(ck.ND)
hp=out/'homotopy_frontpin.npz'
if hp.exists():
    z=np.load(hp); hom=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
else:
    import json; st=json.loads((out/'ptc_state.json').read_text()); hom=dict(mode='ptc',lam=1.0,delta=st['delta'])
cols={n:(out/f'cols_{base.name}'/f'col_{n}_p',out/f'cols_{base.name}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(base,x0,cols,1e-4,True,hom)
act=it.active; b=it.obj(it.R); tot=np.sum(b[act]**2)
print(f'{base.name}: merit {np.sqrt(tot):.4e}')
ultra=(it.frac<=1e-20)&act; trace=(it.frac>1e-20)&(it.frac<=1e-10)&act; sig=(it.frac>1e-10)&act
rows=np.zeros_like(act); rows[132]=True; rT=np.zeros_like(act); rT[133]=True
for name,m in [('ultra',ultra&~rows&~rT),('trace',trace&~rows&~rT),('sig pops',sig&~rows&~rT),('n_e',rows&act),('T',rT&act)]:
    print(f'  {name:9s} share of merit^2 {np.sum(b[m]**2)/tot:6.3f}  max|W M^-1 R| {np.max(np.abs(b[m])) if m.any() else 0:.3e}')
bb=(b**2)*act; prof=bb.sum(0)/tot; print('  top depths',[(int(d)+1,round(float(prof[d]),3)) for d in np.argsort(-prof)[:8]])
rowsh=bb.sum(1)/tot; print('  top rows',[(int(r)+1,round(float(rowsh[r]),3)) for r in np.argsort(-rowsh)[:8]])
