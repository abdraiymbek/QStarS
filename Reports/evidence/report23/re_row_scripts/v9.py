import sys, re, os, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record
Q = '/Users/abraiymbek/Projects/QStarS/results/'
gd = Q + 'gencool_167/pin100_it05'; nk = Q + 'nk167_pin/theta1.00/it05_trial_l1'
x, R = read_trace_record(nk + '/QSTARS_RESID_TRACE.bin'); nt, nd = x.shape
def blocks(path, label, nrow_max=40):
    L = open(path).read().split('\n'); out = {}; cols = None; i = 0
    while i < len(L):
        s = L[i].strip()
        if s.startswith('Depth'):
            cols = [int(v) for v in L[i + 1].split()]; i += 2; continue
        if s.startswith(label) and cols:
            j = i + 1; rr = []
            while j < len(L) and L[j].strip() and re.match(r'^[-+0-9.Ee ]+$', L[j].strip()) and len(rr) < nrow_max:
                rr.append([float(v) for v in L[j].split()]); j += 1
            for c, dpt in enumerate(cols): out[dpt] = [r[c] for r in rr]
        i += 1
    return out
ph = blocks(gd + '/HIPRRR', 'HI Photoionization Rates'); rc = blocks(gd + '/HIPRRR', 'HI Recombination Rates')
# two-photon up rate per atom (first transition block = H I 1s-2s)
L = open(gd + '/TWO_PHOT_SUM').read().split('\n')
i0 = next(i for i, s in enumerate(L) if s.strip().startswith('Up rate'))
vals = []
for s in L[i0 + 1:]:
    if s.strip() and not re.match(r'^[-+0-9.Ee ]+$', s.strip()): break
    vals += [float(v) for v in s.split()]
up = np.array(vals[:nd])
s_cum = np.frombuffer(open(nk + '/QSTARS_STEP.bin', 'rb').read(), '<f8').reshape((nt, nd), order='F').copy()
F = np.ones(nd)
print(' d   T[K]   n2       P2=photoion(n=2) [cm-3 s-1] per-atom[s-1]  phot(n>=3)  recomb(n>=2)  R_3     sumR(n>=2)  2g_up(OLD)  pool feed w/o 2g   F=1+feed/P2')
for d in range(22, 50):
    k = d - 1
    P2 = ph[d][1]; Pn = sum(ph[d][2:]); rec = sum(rc[d][1:])
    pool = R[2:31, k].sum(); nup = up[k] * x[1, k]; feed = pool - nup
    Fd = 1 + max(feed, 0) / P2
    if 26 <= d <= 46: F[k] = Fd
    print('%3d %7.0f %9.3e %10.3e %10.3e %10.3e %10.3e %8.1f %9.1f %9.1f %9.1f %10.1f' % (d, x[133, k] * 1e4, x[2, k], P2, P2 / x[2, k], Pn, rec, R[2, k], pool, nup, feed, Fd))
snew = s_cum.copy()
snew[2, :] = 1.0 - F * (1.0 - s_cum[2, :])
out = '/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/synth/n2boost_state'
os.makedirs(out, exist_ok=True)
open(out + '/QSTARS_STEP.bin', 'wb').write(np.asarray(snew, '<f8').tobytes(order='F'))
np.savetxt(out + '/F_by_depth.txt', np.c_[np.arange(1, nd + 1), F], fmt=['%d', '%.6e'])
print('wrote', out, 'size', os.path.getsize(out + '/QSTARS_STEP.bin'), 'max F', F.max(), 'n rows changed', int((snew != s_cum).sum()))
