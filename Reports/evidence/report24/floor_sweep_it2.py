#!/usr/bin/env python3
"""Report 24: constrained iterate 2 of fixedT_b35. (1) The least squares with the preconditioner floor at 1e-2 ... 3e-1 on the
run's 19 products: linear residual, coefficient of the deep hydrogen piece (direction 15), its preconditioned action outside
its depths, the largest deep step component. (2) The same directions' exact action, row-scaled (no inverse), split into the
direct (fixed-J) part and the radiation part, inside and outside their depths. Output: floor_sweep_it2.txt."""
src=open(str(__import__('pathlib').Path(__file__).with_name('leakage_check.py'))).read().split('S1, J1 = pairs')[0]
g={'__file__':str(__import__('pathlib').Path(__file__).with_name('leakage_check.py'))}
exec(compile(src,'lc','exec'),g)
import numpy as np
O=g['O']
S,JS=g['pairs']([O/'it02_from_it02_trial_l1'])
for tau in (1e-2,3e-2,1e-1,3e-1):
    it=g['iterate'](O,'it02_trial_l1',tau)
    A=np.stack([it.obj(j).ravel() for j in JS],axis=1); b=it.obj(it.R).ravel()
    y=np.linalg.lstsq(A,b,rcond=1e-12)[0]; s=sum(yi*si for yi,si in zip(y,S))
    a15=A[:,15].reshape(S[0].shape); sup=np.abs(S[15]).sum(0)>0
    deep=np.max(np.abs(s[:133,118:167][it.sig[:133,118:167]]))
    print(f'tau {tau:5.0e}: lin res {np.linalg.norm(A@y-b)/np.linalg.norm(b):.3f} y15 {y[15]:8.2f} dir15 outside {np.linalg.norm(a15[:,~sup])/np.linalg.norm(a15):.1%} deep max|s| {deep:.3f} max|s| sig {np.max(np.abs(s[it.sig])):.2f}')
import cool_newton_krylov as ck
ck.MERIT_MODE='rowscale'
it=g['iterate'](O,'it02_trial_l1',1e-2)
for i in (15,16,17):
    jz=JS[i]; a=it.obj(jz); sup=np.abs(S[i]).sum(0)>0
    k=np.argmax(np.abs(a).max(0))
    print(f'dir {i}: row-scaled exact action outside support {np.linalg.norm(a[:,~sup])/np.linalg.norm(a):.1%}, largest at d{k+1} row {np.argmax(np.abs(a[:,k]))+1}; |Wz| {np.linalg.norm(it.W*S[i]):.3g}')
    jd=it.Jdir.matvec(S[i]); rad=jz-jd; ad=it.obj(jd); ar=it.obj(rad)
    print(f'        direct part outside {np.linalg.norm(ad[:,~sup]):.3g} inside {np.linalg.norm(ad[:,sup]):.3g}; radiation part outside {np.linalg.norm(ar[:,~sup]):.3g} inside {np.linalg.norm(ar[:,sup]):.3g}')
