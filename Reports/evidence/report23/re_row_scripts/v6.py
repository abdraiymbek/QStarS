import sys, re, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record
Q = '/Users/abraiymbek/Projects/QStarS/results/'
def rows(path, label):
    L = open(path).read().split('\n'); out = {}; cols = None; i = 0
    while i < len(L):
        s = L[i].strip()
        if s.startswith('Depth'):
            cols = [int(v) for v in L[i + 1].split()]; i += 2; continue
        if s.startswith(label):
            j = i + 1; rr = []
            while j < len(L) and L[j].strip() and re.match(r'^[-+0-9.Ee ]+$', L[j].strip()):
                rr.append([float(v) for v in L[j].split()]); j += 1
            for c, dpt in enumerate(cols): out[dpt] = [r[c] for r in rr]
        i += 1
    return out
F = 4 * np.pi * 1e-10
for p in ('nk167_hot/it00_base', 'nk167_hot/it01_trial_r0.0625', 'nk167_hot/it04_trial_r0.03125', 'nk167_hot/it17_trial_l1', 'nk167_pin/theta1.00/it05_trial_l1'):
    g = Q + p + '/GENCOOL'; x, R = read_trace_record(Q + p + '/QSTARS_RESID_TRACE.bin')
    net = rows(g, 'Net Cooling Rate'); col = rows(g, 'HI Collisional Cooling'); bf = rows(g, 'HI Bound-Free Cooling'); ffh = rows(g, 'HI (ion) Free-Free'); h0bf = rows(g, 'H0 Bound-Free'); h0ff = rows(g, 'H0 (ion) Free-Free'); adv = rows(g, 'Adiabatic cooling rate (V term)'); adt = rows(g, 'Adiabatic cooling rate (dTdR term)')
    print(p)
    print('   d    T[K]     ne        n2       GENnet      HIcol     HIbf(sum)   HIbf n=1   HIbf n=2   HIff      H-bf      H-ff     adV      adT    | RTc[erg]   hnuR3')
    for d in (31, 33, 35, 37, 39):
        k = d - 1
        print('  %3d %7.0f %9.3e %9.3e %10.3e %10.3e %10.3e %10.3e %10.3e %9.2e %9.2e %9.2e %8.1e %8.1e | %10.3e %10.3e' % (d, x[133, k] * 1e4, x[132, k], x[2, k], net[d][0], col[d][0], sum(bf[d]), bf[d][0], bf[d][1], ffh[d][0], h0bf[d][0], h0ff[d][0], adv.get(d,[float("nan")])[0], adt.get(d,[float("nan")])[0], F * R[133, k], 6.62607e-12 * 3.288087 * 0.75 * R[2, k]))
