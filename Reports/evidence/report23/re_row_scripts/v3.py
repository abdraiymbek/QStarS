import sys, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record, read_ba_blocks
Q = '/Users/abraiymbek/Projects/QStarS/'
nk = Q + 'results/nk167_pin/theta1.00/it05_trial_l1'
x, R = read_trace_record(nk + '/QSTARS_RESID_TRACE.bin')
lo, di, up = read_ba_blocks(nk + '/QSTARS_BA_TRACE.bin', 134, 167)
h15 = 6.62607e-12; E = 3.288087 / np.arange(1, 31) ** 2
F = 4 * np.pi * 1e-10
net = {33: 4.3412e-10, 35: 1.1734e-09, 37: 2.6593e-09}; D = {33: 5.760e-10, 35: 1.411e-09, 37: 3.468e-09}
for d in (33, 35, 37):
    k = d - 1
    Sx_exc = h15 * ((E[0] - E[1:]) * R[2:31, k]).sum()
    ion = h15 * E[0] * R[1, k]
    rowl1 = [np.abs(lo[k, r]).sum() + np.abs(di[k, r]).sum() + np.abs(up[k, r]).sum() for r in range(134)]
    print('d%d RTc=%.4e  hSum(E1-En)Rn(n=2..30)=%.4e (n=2 alone %.4e)  hE1*R[row2]=%.3e (R[row2]=%.3e, /l1=%.1e)  -net+Sx+D (ion sign -) = %.4e ; (ion sign +) = %.4e ; no D: %.4e' % (
        d, F * R[133, k], Sx_exc, h15 * (E[0] - E[1]) * R[2, k], ion, R[1, k], R[1, k] / rowl1[1], -net[d] + Sx_exc - ion + D[d], -net[d] + Sx_exc + ion + D[d], -net[d] + Sx_exc - ion))
    print('    rows R/l1: H-(1) %.2e  ion(2) %.2e  n=2(3) %.2e  n=3(4) %.2e  n=4 %.2e  n=5 %.2e ; HII row32 R=%.3e ; n_e row R=%.3e ; n=2 row l1 %.3e ; di[n2,n2]=%.4e di[n2,n1]=%.4e' % (
        R[0, k] / rowl1[0], R[1, k] / rowl1[1], R[2, k] / rowl1[2], R[3, k] / rowl1[3], R[4, k] / rowl1[4], R[5, k] / rowl1[5], R[31, k], R[132, k], rowl1[2], di[k, 2, 2], di[k, 2, 1]))
    # local fixed-J Newton step for n=2 alone and block solve over H I excited levels (rows 3..31) at frozen everything else
    t = np.arange(2, 31)
    s = np.linalg.solve(di[k][np.ix_(t, t)], R[t, k])
    print('    frozen-J local Newton fractional step (x->x(1-s)) n=2: %.3e (scalar R/diag %.3e); n=3 %.3e' % (s[0], R[2, k] / di[k, 2, 2], s[1]))
