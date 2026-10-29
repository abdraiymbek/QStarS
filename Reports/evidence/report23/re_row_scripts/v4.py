import sys, glob, os, re, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record
Q = '/Users/abraiymbek/Projects/QStarS/'
def net35(path, d=35):
    try:
        L = open(path).read().split('\n')
    except Exception:
        return float('nan'), float('nan')
    cols = None; out = float('nan'); col = float('nan')
    for i, s in enumerate(L):
        t = s.strip()
        if t.startswith('Depth'):
            cols = [int(v) for v in L[i + 1].split()]
        if cols and d in cols:
            if t.startswith('Net Cooling Rate'):
                out = float(L[i + 1].split()[cols.index(d)])
            if t.startswith('HI Collisional Cooling'):
                col = float(L[i + 1].split()[cols.index(d)])
    return out, col
dirs = []
for pat in ('results/nk167_hot/it*_trial_*', 'results/nk167_pin/theta1.00/it0*_trial_l1', 'results/nk167_pin/theta1.15/it*_trial_*', 'results/nk167_hot/it00_base', 'results/nk167_pin/theta1.00/it00_base'):
    dirs += sorted(glob.glob(Q + pat))
F = 4 * np.pi * 1e-10; hnuA = 6.62607e-12 * 3.288087 * 0.75
print('%-58s %8s %10s %10s %9s %9s %9s %11s %11s %s' % ('pass', 'T35[K]', 'n1', 'n2', 'HI/H', 'R_T', 'R_3', 'RTc[erg]', 'GENnet', 'GENCOOL mtime vs trace'))
for p in dirs:
    tr = p + '/QSTARS_RESID_TRACE.bin'
    if not os.path.exists(tr): continue
    x, R = read_trace_record(tr); k = 34
    g = p + '/GENCOOL'; n, c = net35(g)
    age = (os.path.getmtime(g) - os.path.getmtime(tr)) if os.path.exists(g) else float('nan')
    print('%-58s %8.0f %10.3e %10.4e %9.4f %9.3f %9.1f %11.3e %11.3e dt=%.0fs' % (p.replace(Q + 'results/', ''), x[133, k] * 1e4, x[1, k], x[2, k], x[1:31, k].sum() / x[1:32, k].sum(), R[133, k], R[2, k], F * R[133, k], n, age))
