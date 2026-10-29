"""dt -> 0 continuation of the step on the event mesh and on the pre-event mesh; Jacobian conditioning."""
from common import *
from locate_lib import internals
from qstars.nonlinear import pattern_jacobian
state,config,eos,opacity=load()
out={}
print('dt continuation (same accepted state, step 1255):')
print(f"{'dt/yr':>10} | {'merge: mapped dlnT':>18} {'cell':>4} {'it':>3} | {'no merge: mapped dlnT':>21} {'cell':>4} {'it':>3} | {'event jump (merge vs no-merge, cell 167 T)':>30}")
rows=[]
for dt in [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.,10.,30.,100.,172.1438701239381]:
    R={}
    for rm in (True,False):
        I=internals(state,dt,eos,opacity,config,remesh=rm)
        P,x=I['P'],I['x'];n=P.n
        th=x[n+1:2*n+1];th0=P.x0[n+1:2*n+1]
        k=int(np.argmax(abs(th-th0)))
        R[rm]=dict(m=float(np.max(abs(th-th0))),k=k,it=I['info']['iterations'],T167=float(P.post.T[167]*np.exp(th[167])),P=P,x=x)
    jump=np.log(R[True]['T167']/R[False]['T167'])
    rows.append(dict(dt=dt,merge=R[True]['m'],merge_cell=R[True]['k'],merge_it=R[True]['it'],nomerge=R[False]['m'],nomerge_cell=R[False]['k'],nomerge_it=R[False]['it'],event_dlnT167=jump))
    print(f"{dt:10.3g} | {R[True]['m']:18.5g} {R[True]['k']:4d} {R[True]['it']:3d} | {R[False]['m']:21.5g} {R[False]['k']:4d} {R[False]['it']:3d} | {jump:+.5f}")
    if dt in (1.,172.1438701239381):out[dt]=R
json.dump(rows,open(SCR/'continuity_rows.json','w'),indent=1)
print()
print('Jacobian (scaled residual rows, as Newton sees them) at the converged step solutions:')
for dt,R in out.items():
    for rm in (True,False):
        P,x=R[rm]['P'],R[rm]['x']
        J=pattern_jacobian(P.evaluate,x,P.pattern).toarray()
        sv=np.linalg.svd(J,compute_uv=False)
        U,S,Vt=np.linalg.svd(J)
        v=Vt[-1];n=P.n
        blocks=dict(y=np.abs(v[:n+1]),lnT=np.abs(v[n+1:2*n+1]),L=np.abs(v[2*n+1:3*n]),tail=np.abs(v[3*n:]))
        where={k:(int(np.argmax(b)),float(np.max(b))) for k,b in blocks.items()}
        sign,logdet=np.linalg.slogdet(J)
        print(f"dt={dt:g} merge={rm}: n={n} sigma_max={sv[0]:.3g} sigma_min={sv[-1]:.3g} next={sv[-2]:.3g},{sv[-3]:.3g} cond={sv[0]/sv[-1]:.3g} det_sign={sign:+.0f} smallest-vector peak per block {where}")
