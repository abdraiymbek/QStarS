"""Report 24: preconditioner variants measured against exact products at archived states."""
import sys, json, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
import cool_newton_krylov as ck, cnk_trap
from qstars.nlte_newton import read_trace_record, BlockTridiag
ck.SIG_FRAC,ck.TRACE_WEIGHT,ck.ULTRA_FRAC,ck.ULTRA_WEIGHT,ck.WEIGHT_RAMP=1e-10,1e-3,1e-20,0.0,2.0
ck.MERIT_MODE,ck.EQ_FLOOR,ck.PC_T_FLOOR='precond',0.0,1e-3
TAUS=[1e-6,1e-5,1e-4,1e-3,1e-2]
def states():
    F=ROOT/'results/nk167_trap/fixedT_b35'
    z=np.load(F/'homotopy_frontpin.npz'); homF=dict(Tfix=z['Tfix'],c=z['c'],w=z['w'],R0=None,mode='pin',lam=1.0)
    for it,b in ((1,'it01_trial_l0.5'),(2,'it02_trial_l1'),(3,'it03_trial_l0.25'),(4,'it04_trial_l0.25'),(5,'it05_trial_l1_d16'),(6,'it06_trial_l0.25')):
        yield f'fixedT it{it}',F,b,dict(homF),F/f'it{it:02d}_from_{b}'
    T=ROOT/'results/nk167_trap/free_b35'
    for it,b,dd in ((1,'it01_trial_l0.5','it01_from_it01_trial_l0.5'),(2,'it02_trial_l1','it02_from_it02_trial_l1'),(3,'it03_trial_l1','it03d0_from_it03_trial_l1')):
        yield f'free it{it}',T,b,dict(mode='ptc',lam=1.0,delta=0.06453),T/dd
rows=[]
x0=None
for name,O,b,hom,itd in states():
    pairs=[]
    for j in sorted(itd.glob('jv*')):
        if (j/'pairs.npz').exists():
            Q=np.load(j/'pairs.npz'); pairs+=[(Q['z'][:,:,i],Q['jz'][:,:,i]) for i in range(Q['z'].shape[2])]
    if not pairs or not (O/(b+'_ali')/'QSTARS_BA_TRACE.bin').exists(): print(name,'skipped'); continue
    if x0 is None:
        x0,_=read_trace_record(O/'it00_base'/'QSTARS_RESID_TRACE.bin'); ck.NT,ck.ND=x0.shape
    cols={n:(O/f'cols_{b}'/f'col_{n}_p',O/f'cols_{b}'/f'col_{n}_m') for n in ('ne','hii','T')}
    it=ck.Iterate(O/(b+'_tp'),x0,cols,1e-4,True,hom,rj=O/b/'QSTARS_RJ.bin'); it.extend_significance(0.0,rows=[3])
    lo,dif,up=[m.copy() for m in it._pc]
    it.install_ali_increment(O/(b+'_ali')); _,dia,_=it._pc; lo2,_,up2=it._pc
    W=it.W
    def evalM(M):
        q=np.array([np.linalg.norm(W*(M.solve(jz)-z))/np.linalg.norm(W*z) for z,jz in pairs])
        y=M.solve(it.R); return np.median(q),np.percentile(q,90),q.max(),np.linalg.norm(W*y*it.active)
    res={}
    res['fixedJ']=evalM(BlockTridiag(lo,dif,up).factor())
    res['ALI']=evalM(it.M)
    for t in TAUS:
        d2,n=cnk_trap.floor_blocks(dia,t,it.active); res[f'ALI floor {t:g} ({n})']=evalM(BlockTridiag(lo2,d2,up2).factor())
    print(f'== {name}: {len(pairs)} pairs')
    for k,v in res.items():
        print(f'   {k:26s} q median {v[0]:7.3f}  p90 {v[1]:8.3f}  max {v[2]:10.2f} | ||W M^-1 R|| {v[3]:10.3e}')
    rows.append((name,res))
json.dump([(n,{k:list(map(float,v)) for k,v in r.items()}) for n,r in rows],open(ROOT/'results/nk167_trap/floor_q.json','w'),indent=1)
