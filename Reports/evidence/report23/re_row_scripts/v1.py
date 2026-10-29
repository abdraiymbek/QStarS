import sys, struct, re, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record, read_ba_blocks
Q = '/Users/abraiymbek/Projects/QStarS/'
nk = Q + 'results/nk167_pin/theta1.00/it05_trial_l1'
gd = Q + 'results/gencool_167/pin100_it05'
x, R = read_trace_record(nk + '/QSTARS_RESID_TRACE.bin')
xg, Rg = read_trace_record(gd + '/QSTARS_RESID_TRACE.bin')
print('shape', x.shape, 'max rel pop diff nk vs gencool dir', np.max(np.abs(xg / x - 1)), 'RE diff d30-41', np.max(np.abs(Rg[133, 29:41] / R[133, 29:41] - 1)))
lo, di, up = read_ba_blocks(nk + '/QSTARS_BA_TRACE.bin', 134, 167)
def gencool_rows(path, label):
    L = open(path).read().split('\n'); out = {}
    i = 0; cols = None
    while i < len(L):
        s = L[i].strip()
        if s.startswith('Depth'):
            cols = [int(v) for v in L[i + 1].split()]; i += 2; continue
        if s.startswith(label):
            j = i + 1; rows = []
            while j < len(L) and re.match(r'^[-+0-9.Ee ]+$', L[j].strip() or 'x'):
                rows.append([float(v) for v in L[j].split()]); j += 1
            for c, dpt in enumerate(cols):
                out[dpt] = [r[c] for r in rows]
        i += 1
    return out
net = gencool_rows(gd + '/GENCOOL', 'Net Cooling Rate'); net2 = gencool_rows(nk + '/GENCOOL', 'Net Cooling Rate')
col = gencool_rows(gd + '/GENCOOL', 'HI Collisional Cooling'); bf = gencool_rows(gd + '/GENCOOL', 'HI Bound-Free Cooling')
hP = 6.62607e-27; h15 = hP * 1e15; E = 3.288087 / np.arange(1, 31) ** 2; nuA = E[0] - E[1]
F = 4 * np.pi * 1e-10
print(' d   T[K]   HI/H     n1        n2       n2/nH     R_T     c(l1)   R_T/c    4pie-10R_T  GENnet(new) GENnet(nk)  R_3      hnuR3/RTc  dRT/dln n1  dRT/dln n2  dRT/dlnT(native)  HIcol   HIbf(n=2)')
for d in range(30, 42):
    k = d - 1
    c = np.abs(lo[k, 133]).sum() + np.abs(di[k, 133]).sum() + np.abs(up[k, 133]).sum()
    nH = x[1:32, k].sum()
    print('%3d %7.0f %6.4f %9.3e %9.3e %8.1e %8.3f %9.1f %8.2e %10.3e %10.3e %10.3e %8.1f %7.3f %10.2f %10.2f %10.3e %10.3e %10.3e' % (
        d, x[133, k] * 1e4, x[1:31, k].sum() / nH, x[1, k], x[2, k], x[2, k] / nH, R[133, k], c, R[133, k] / c, F * R[133, k], net[d][0], net2[d][0],
        R[2, k], h15 * nuA * R[2, k] / (F * R[133, k]), di[k, 133, 1], di[k, 133, 2], di[k, 133, 133], col[d][0], bf[d][1]))
np.save('/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/synth/x_pin.npy', x)
