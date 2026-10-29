"""First look at the DRYRUN's second blocker (step 1598): where Newton's residual stalls, Jacobian conditioning."""
import os,sys
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'patched','src'));sys.path.insert(0,os.path.join(HERE,'patched','scripts'))
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import numpy as np,json
from qstars import quasistatic as qs,cavity
from qstars.nonlinear import pattern_jacobian
import validate_ball_track as vbt
YEAR=365.25*86400;MSUN=1.98847e33
st,config=qs.load_checkpoint(os.path.join(HERE,'DRYRUN_continuation_counterfactual_400','checkpoint.npz'))
eos,opacity=vbt.physics(config)
s=qs.settings(config);q_=s['quasistatic'];zones=config['numerics']['zones']
print('state step',st.steps,'Mbh',st.bh_mass/MSUN,'zones',len(st.dm))
def blocks(P,f):
    n=P.n
    names=[('momentum',0,n),('energy',n,2*n),('closure',2*n,3*n-1),('piston',3*n-1,3*n),('cavity',3*n,3*n+1),('photosphere',3*n+1,3*n+2),('growth',3*n+2,3*n+3)]
    return {k:(float(np.max(abs(f[a:b]))),int(np.argmax(abs(f[a:b])))) for k,a,b in names}
for dt in (1e-3,1.,30.9):
    pred=qs.predict(st,dt*YEAR,eos,config,'evolve',None,s)
    cut,ci=cavity.slice_base(st,pred['swept'],eos)
    pre,rg=qs.regrid_base(st,cut,eos,zones,q_['base_mass_fraction'])
    post,rm=qs.remesh_base(st,pre,eos,zones,q_['split'],q_['merge'],protect=rg['region_cells'])
    P=qs.Problem(st,post,dt*YEAR,eos,opacity,config,'evolve',None,pred,s)
    x,info=P.solve()
    f=P.evaluate(x)
    print(f"dt={dt} yr: events splits={rm['splits']} merges={rm['merges']} n={P.n}; success={info['success']} it={info['iterations']} residual={info['residual']:.3g} msg={info['message']}")
    print('   residual per block at final iterate (max, index):',blocks(P,f))
    f0=P.evaluate(P.x0);print('   residual per block at predictor:',blocks(P,f0))
    for label,z in (('predictor',P.x0),('final iterate',x)):
        J=pattern_jacobian(P.evaluate,z,P.pattern).toarray()
        U,S,Vt=np.linalg.svd(J);v=Vt[-1];n=P.n
        peak={k:(int(np.argmax(abs(b))),round(float(np.max(abs(b))),3)) for k,b in (('y',v[:n+1]),('lnT',v[n+1:2*n+1]),('L',v[2*n+1:3*n]),('tail',v[3*n:]))}
        print(f'   J at {label}: sigma_min={S[-1]:.3g} next={S[-2]:.3g} cond={S[0]/S[-1]:.3g} smallest right vector peaks {peak}')
    th=x[P.n+1:2*P.n+1]-P.x0[P.n+1:2*P.n+1];k=int(np.argmax(abs(th)))
    print(f'   largest lnT move from predictor at the final iterate: cell {k} of {P.n} ({th[k]:+.3g}); mass below surface {(post.dm.sum()-np.cumsum(post.dm)[k]+.5*post.dm[k])/MSUN:.4g} Msun')
