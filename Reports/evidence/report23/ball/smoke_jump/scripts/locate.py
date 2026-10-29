from common import *
from qstars import cavity
MSUN=1.98847e33
state,config,eos,opacity=load()
s=qs.settings(config);q_=s['quasistatic'];zones=config['numerics']['zones']
def internals(old,dt_years,remesh=True):
    dt=dt_years*YEAR
    pred=qs.predict(old,dt,eos,config,'evolve',None,s)
    swept=pred['swept']
    for passes in range(1,4):
        cut,ci=cavity.slice_base(old,swept,eos)
        premesh,regrid=qs.regrid_base(old,cut,eos,zones,q_['base_mass_fraction'])
        if remesh:post,rm=qs.remesh_base(old,premesh,eos,zones,q_['split'],q_['merge'],protect=regrid['region_cells'])
        else:post,rm=premesh,dict(splits=0,merges=0)
        P=qs.Problem(old,post,dt,eos,opacity,config,'evolve',None,pred,s)
        x,info=P.solve()
        assert info['success'],info['message']
        delta=P.delta(x)
        if abs(delta)<=q_['pass_fraction']*post.dm[0]:break
        swept+=delta
    return dict(pred=pred,cut=cut,premesh=premesh,regrid=regrid,post=post,remesh=rm,P=P,x=x,info=info,passes=passes)
n_old=len(state.dm)
for dt in [1e-3,1.,172.1438701239381]:
    I=internals(state,dt)
    P,x,post=I['P'],I['x'],I['post'];n=P.n
    th=x[n+1:2*n+1];th0=P.x0[n+1:2*n+1]
    k=int(np.argmax(abs(th-th0)))
    mu=np.cumsum(post.dm)-.5*post.dm;M=post.dm.sum()
    print(f'dt={dt} n_old={n_old} n_post={n} premesh n={len(I["premesh"].dm)} merges={I["remesh"]["merges"]} splits={I["remesh"]["splits"]} region={I["regrid"]["region_cells"]} passes={I["passes"]}')
    print(f'  max |th-x0| at cell {k}: {th[k]-th0[k]:+.5f}  (th={th[k]:+.5f}, x0={th0[k]:+.5f}); mass below surface {(M-mu[k])/MSUN:.4g} Msun, from base {mu[k]/MSUN:.4g}; T={post.T[k]*np.exp(th[k]):.5g} rho={post.rho[k]:.4g}')
    d=th-th0
    top=np.argsort(-abs(d))[:6]
    print('  top cells',[(int(i),round(float(d[i]),5)) for i in top])
    # where did the merge happen
    pm=I['premesh']
    if len(pm.dm)!=n:
        a=np.cumsum(pm.dm);b=np.cumsum(post.dm)
        # first differing face
        j=0
        while j<n and abs(a[j]-b[j])<=1e-12*a[-1]:j+=1
        print(f'  first differing cell (post) {j}; premesh dm[{j}:{j+3}]/Msun={pm.dm[j:j+3]/MSUN}, post dm[{j}]={post.dm[j]/MSUN}; post cells {n}; mass below surface of merged cell {(M-mu[j])/MSUN:.4g} Msun; T_post={post.T[j]:.5g} premesh T={pm.T[j:j+2]}')
        print(f'  targets at premesh: ', qs.targets(state,pm,eos,zones)[j-1:j+3]/MSUN, 'dm',pm.dm[j-1:j+3]/MSUN)
