import pickle,numpy as np
b=pickle.load(open('v_bitwise_before.pkl','rb'));a=pickle.load(open('v_bitwise_after.pkl','rb'));r=pickle.load(open('v_repro.pkl','rb'))
for dt in (1.,172.1438701239381):
    B=b[('registered',dt)]
    for ev in ('registered','measured','counterfactual'):
        A=a[(ev,dt)];diffs=[k for k in ('T','r','dm','L','rc') if not np.array_equal(A[k],B[k])]+[k for k in ('ledger','aux','bh','cav','rho_ph') if A[k]!=B[k]]
        dk=[k for k in B['diag'] if not (B['diag'][k]==A['diag'].get(k) or (isinstance(B['diag'][k],float) and np.isnan(B['diag'][k])))]
        print(dt,ev,'state diffs',diffs,'diag diffs',dk,'extra diag',sorted(set(A['diag'])-set(B['diag'])),'metric diffs',[k for k in B['metrics'] if B['metrics'][k]!=A['metrics'][k]],'delta_actual',A['metrics']['delta_actual'])
    nm=r[(0.45,dt)]['metrics'];cf=a[('counterfactual',dt)]['metrics']
    print('  counterfactual metrics == merge-suppressed (0.45) step metrics:',all(cf[k]==nm[k] for k in cf))
