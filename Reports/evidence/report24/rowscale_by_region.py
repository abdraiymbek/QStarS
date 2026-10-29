import sys
from pathlib import Path
import numpy as np
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck
from qstars.nlte_newton import read_trace_record
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'rowscale', 0.0, 1e-3
O=ROOT/'results/nk167_trap/fixedT_b35'; ref='it01_trial_l0.5'; B=O/'branch_unconstrained_it2plus'
x0,_=read_trace_record(O/'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape
cols={n:(O/f'cols_{ref}'/f'col_{n}_p',O/f'cols_{ref}'/f'col_{n}_m') for n in ('ne','hii','T')}
it=ck.Iterate(O/(ref+'_tp'),x0,cols,1e-4,True,None,rj=O/ref/'QSTARS_RJ.bin'); it.extend_significance(0.0,rows=[3])
regs=[(0,22),(22,45),(45,75),(75,118),(118,150),(150,167)]
rows=list(range(0,133))   # populations and n_e (T row pinned)
print(f'{"state":28s}' + ''.join(f'{f"{a+1}-{b}":>10s}' for a,b in regs) + '     total')
for lab,d in (('start_base',O/'start_base'),('iterate 1 (it01_trial_l0.5)',O/'it01_trial_l0.5_tp'),('unconstr. it04',B/'it04_trial_l0.25_tp'),
              ('constrained 2',O/'it02_trial_l1_tp'),('constrained 3 (l1)',O/'it03_trial_l1'),('constrained 3 (l0.5)',O/'it03_trial_l0.5'),('constrained 4 (l1)',O/'it04_trial_l1'),('constrained 5 (l0.25)',O/'it05_trial_l0.25')):
    if not (d/'QSTARS_RESID_TRACE.bin').exists(): d=Path(str(d).replace('_tp',''))
    x,R=read_trace_record(d/'QSTARS_RESID_TRACE.bin'); w=it.obj(R)[:133]
    print(f'{lab:28s}' + ''.join(f'{np.linalg.norm(w[:,a:b]):10.3g}' for a,b in regs) + f'{np.linalg.norm(w):10.3g}')
print('--- deep (119-167) largest row-scaled residuals')
for lab,d in (('constrained 2',O/'it02_trial_l1_tp'),('constrained 3 (l1)',O/'it03_trial_l1')):
    x,R=read_trace_record(d/'QSTARS_RESID_TRACE.bin'); w=it.obj(R)
    sub=np.abs(w[:133,118:167]); idx=np.argsort(sub.ravel())[::-1][:8]
    print(lab, [(int(r)+1, int(c)+119, round(float(sub[r,c]),3)) for r,c in zip(*np.unravel_index(idx,sub.shape))])
s=np.load(O/'it02_from_it02_trial_l1'/'newton_step.npy')
sub=np.abs(s[:133,118:167]); idx=np.argsort(sub.ravel())[::-1][:8]
print('it02 step deep largest |s|', [(int(r)+1, int(c)+119, round(float(s[r,c+118]),3)) for r,c in zip(*np.unravel_index(idx,sub.shape))])
