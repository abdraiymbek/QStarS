import sys, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record
Q = '/Users/abraiymbek/Projects/QStarS/results/nk167_pin/theta1.00/'
xp, Rp = read_trace_record(Q + 'cols_it05_trial_l1/col_T_p/QSTARS_RESID_TRACE.bin')
xm, Rm = read_trace_record(Q + 'cols_it05_trial_l1/col_T_m/QSTARS_RESID_TRACE.bin')
x, R = read_trace_record(Q + 'it05_trial_l1/QSTARS_RESID_TRACE.bin')
for d in (33, 35, 37):
    k = d - 1
    dl = np.log(xp[133, k] / xm[133, k])
    print('d%d  T+ %.6f T- %.6f  dlnT=%.3e  dR_T/dlnT (frozen-J FD) = %.4f ; dR_3/dlnT = %.3f ; R_T base %.4f ; n2 +/- %.5e %.5e' % (d, xp[133, k], xm[133, k], dl, (Rp[133, k] - Rm[133, k]) / dl, (Rp[2, k] - Rm[2, k]) / dl, R[133, k], xp[2, k], xm[2, k]))
