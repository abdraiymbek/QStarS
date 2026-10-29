import sys, glob, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record, read_ba_blocks, read_jv_rad
Q = '/Users/abraiymbek/Projects/QStarS/results/nk167_pin/theta1.00/'
x, R = read_trace_record(Q + 'it05_trial_l1/QSTARS_RESID_TRACE.bin')
lo, di, up = read_ba_blocks(Q + 'it05_trial_l1/QSTARS_BA_TRACE.bin', 134, 167)
for jv in sorted(glob.glob(Q + 'it05_from_it05_trial_l1/jv*')):
    P = np.load(jv + '/pairs.npz'); z = P['z']; jz = P['jz']
    s, drad = read_jv_rad(jv + '/QSTARS_JV_RAD.bin')
    print(jv.split('/')[-1], 'z shape', z.shape, 'dirs equal stored:', np.allclose(s, z))
    for i in range(z.shape[2]):
        for d in (33, 35, 37):
            k = d - 1
            # native fixed-J block action restricted to Ly-alpha partners (n1,n2 columns) on the n=2 row
            lya = di[k, 2, 1] * z[1, k, i] + di[k, 2, 2] * z[2, k, i]
            nat = lo[k, 2] @ z[:, k - 1, i] + di[k, 2] @ z[:, k, i] + up[k, 2] @ z[:, k + 1, i]
            natT = lo[k, 133] @ z[:, k - 1, i] + di[k, 133] @ z[:, k, i] + up[k, 133] @ z[:, k + 1, i]
            if abs(z[2, k, i]) + abs(z[1, k, i]) < 1e-12: continue
            print('   dir %d d%d: z_n1=%+.3e z_n2=%+.3e z_n3=%+.3e z_HII=%+.3e z_ne=%+.3e z_T=%+.3e | n=2 row: native fixed-J action %+.4e (Lya n1,n2 part %+.4e)  drad %+.4e  full Jz %+.4e  [R_3=%.1f] | RE row: native %+.4e drad %+.4e Jz %+.4e [R_T=%.3f]' % (
                i, d, z[1, k, i], z[2, k, i], z[3, k, i], z[31, k, i], z[132, k, i], z[133, k, i], nat, lya, drad[2, k, i], jz[2, k, i], R[2, k], natT, drad[133, k, i], jz[133, k, i], R[133, k]))
st = np.load(Q + 'it05_from_it05_trial_l1/newton_step.npy')
print('accepted Newton step it05: n2 at d33,35,37:', st[2, 32], st[2, 34], st[2, 36], ' n1:', st[1, 34], ' HII:', st[31, 34], ' ne:', st[132, 34], 'T:', st[133, 34])
