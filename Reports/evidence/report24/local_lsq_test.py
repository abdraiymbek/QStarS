#!/usr/bin/env python3
"""Report 24: region-local least squares (each direction's action kept only within its depth support +- margin) against
the full least squares, constrained iterates 1 and 2 of fixedT_b35 (floor 1e-2). Prints the true linear residual, the
largest coefficient and step components. Output of the run: local_lsq_test.txt."""
src=open(str(__import__('pathlib').Path(__file__).with_name('leakage_check.py'))).read().split('S1, J1 = pairs')[0]
g={'__file__':str(__import__('pathlib').Path(__file__).with_name('leakage_check.py'))}
exec(compile(src,'lc','exec'),g)
import numpy as np
O=g['O']; x0=g['x0']
def test(b_, dirs, floor, label):
    it=g['iterate'](O,b_,floor); S,JS=g['pairs'](dirs)
    A=np.stack([it.obj(j).ravel() for j in JS],axis=1); b=it.obj(it.R).ravel(); nb=np.linalg.norm(b)
    ev=np.arange(117,166,2); W=it.W*it.active
    for m in (None, 0, 3, 8):
        if m is None:
            Am=A
        else:
            Am=A.copy()
            for i,z in enumerate(S):
                dep=np.nonzero(np.abs(z).sum(0))[0]; mask=np.zeros(x0.shape,bool); mask[:,max(0,dep[0]-m):dep[-1]+1+m]=True
                Am[:,i]*=mask.ravel()
        y=np.linalg.lstsq(Am,b,rcond=1e-12)[0]; s=sum(yi*si for yi,si in zip(y,S))
        deep=np.max(np.abs(s[:133,118:167][it.sig[:133,118:167]]))
        print(f'{label} margin {str(m):4s}: true lin res {np.linalg.norm(A@y-b)/nb:.3f}  max|y| {np.max(np.abs(y)):8.2f}  max|s| sig {np.max(np.abs(s[it.sig])):7.2f}  deep max|s| {deep:6.3f}  ||Ws|| {np.linalg.norm(W*s):7.2f}')
test('it02_trial_l1',[O/'it02_from_it02_trial_l1'],1e-2,'c-it2')
test('it01_trial_l0.5',[O/'it01f0_from_it01_trial_l0.5',O/'it01_from_it01_trial_l0.5'],1e-2,'c-it1')
