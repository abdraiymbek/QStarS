"""Exploratory: the third stop (DRYRUN 2, step 1624). dt continuation on the step's own mesh, sigma_min along it,
where Newton stalls, fold.classify (registered code) at that state, and the outer-cell table."""
import os,sys,json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'));sys.path.insert(0,HERE)
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import numpy as np
from qstars import quasistatic as qs,cavity,fold
from qstars.nonlinear import pattern_jacobian
import validate_ball_track as vbt
YEAR=365.25*86400;MSUN=1.98847e33
st,config=qs.load_checkpoint(os.path.join(HERE,'DRYRUN_continuation_counterfactual_ladder_from1598_400','checkpoint.npz'))
eos,opacity=vbt.physics(config);s=qs.settings(config);q_=s['quasistatic'];zones=config['numerics']['zones']
print('state',st.steps,st.bh_mass/MSUN,len(st.dm))
def problem(dt):
    pred=qs.predict(st,dt*YEAR,eos,config,'evolve',None,s)
    cut,ci=cavity.slice_base(st,pred['swept'],eos)
    pre,rg=qs.regrid_base(st,cut,eos,zones,q_['base_mass_fraction'])
    post,rm=qs.remesh_base(st,pre,eos,zones,q_['split'],q_['merge'],protect=rg['region_cells'])
    return qs.Problem(st,post,dt*YEAR,eos,opacity,config,'evolve',None,pred,s),rm
x=None;rows=[]
for dt in np.geomspace(1e-4,7.04,30):
    P,rm=problem(dt);n=P.n
    start=P.x0 if x is None else np.r_[x[:3*n+2],P.x0[3*n+2]]
    z,info=qs.newton(P.evaluate,start,P.pattern,q_['tolerance'],q_['max_iterations'],q_['polish_iterations'],q_['roundoff_residual'],q_['roundoff_step'])
    J=pattern_jacobian(P.evaluate,z,P.pattern).toarray();U,S,Vt=np.linalg.svd(J);v=Vt[-1]
    f=P.evaluate(z);k=int(np.argmax(abs(f)))
    th=z[n+1:2*n+1]-P.x0[n+1:2*n+1]
    rows.append(dict(dt_years=float(dt),events=rm['merges']+rm['splits'],success=bool(info['success']),iterations=info['iterations'],residual=info['residual'],
                     worst_row=k,sigma_min=float(S[-1]),sigma_2=float(S[-2]),vector_lnT_cell=int(np.argmax(abs(v[n+1:2*n+1]))),vector_lnT_weight=float(np.max(abs(v[n+1:2*n+1]))),
                     mapped_dlnT=float(np.max(abs(th))),cell=int(np.argmax(abs(th)))))
    print(json.dumps(rows[-1]),flush=True)
    if not info['success']:break
    x=z
json.dump(rows,open(os.path.join(HERE,'probe1624_dt_continuation.json'),'w'),indent=1)
out=fold.classify(st,eos,opacity,config)
print('fold.classify at step 1624:',out['label'],'|',out['reason'])
print('  arc points:',[(round(p['bh_mass']/MSUN,3),round(p['t_lambda'],3),p['det_sign'],f"{p['min_singular']:.3g}") for r in out['runs'] for p in r['points']])
json.dump(qs._json(dict(label=out['label'],reason=out['reason'],points=[dict(direction=r['direction'],bh_mass_msun=p['bh_mass']/MSUN,t_lambda=p['t_lambda'],det_sign=p['det_sign'],min_singular=p['min_singular']) for r in out['runs'] for p in r['points']])),
          open(os.path.join(HERE,'fold_classify_dryrun_step1624.json'),'w'),indent=1)
