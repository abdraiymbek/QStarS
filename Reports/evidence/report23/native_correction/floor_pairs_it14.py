"""Report 23: measured floating-point floor of the native correction from identical-state base passes on disk
(results/nk167_hot/it14_trial_l1 = build 28, 2 threads; results/bench22/matched/base_b28_t8, base_b30_t8, base_b30_t2), against the
random-sign model with the withdrawn EPS 3.2e-12, and the measured noise pattern (dR / row l1 at it14) rescaled to the row l1
of it24 and pushed through it24's native M. Output: floor_pairs_it14.txt. Read-only; seconds on one core."""
import sys, json, numpy as np
from pathlib import Path
ROOT = Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'scripts'))
import native_correction_check as ncc
from qstars.nlte_newton import BlockTridiag, read_ba_blocks, read_trace_record
A = ROOT/'results/nk167_hot/it14_trial_l1'
P = [ROOT/'results/bench22/matched'/n for n in ('base_b28_t8','base_b30_t8','base_b30_t2')]
def load(b):
    x, R = read_trace_record(b/'QSTARS_RESID_TRACE.bin'); return x, R
xa, Ra = load(A); nt, nd = xa.shape
lo, di, up = read_ba_blocks(A/'QSTARS_BA_TRACE.bin', nt, nd)
gross = (np.abs(lo)+np.abs(di)+np.abs(up)).sum(axis=2).T
act = Ra != 0
M = BlockTridiag(lo, di, up).factor()
rep = ncc.read_replace(A/'QSTARS_REPLACE_TRACE', nd)
frac, cls = ncc.classes(xa, Ra, di, rep, 1e-10)
Rs = {'it14': Ra}
for p in P:
    x, R = load(p); assert np.array_equal(x, xa), p
    print(p.name, 'active equal', np.array_equal(R != 0, act)); Rs[p.name] = R
def report(y, tag):
    ya = ncc.audit_normalized(y)
    d = {k: (float(np.abs(y[m]).max()) if m.any() else 0) for k, m in cls.items()}
    print(tag, {k: '%.2g' % v for k, v in d.items()}, 'n>1e-3 all:', int(np.sum(ya[act] > 1e-3)), 'audit-norm >1e-3:', int(np.sum((ya > 1e-3) & act)))
names = list(Rs)
dRs = {}
for i in range(len(names)):
    for j in range(i+1, len(names)):
        dR = Rs[names[i]] - Rs[names[j]]
        eps = float(np.max(np.abs(dR[act]) / gross[act]))
        dRs[(names[i], names[j])] = dR
        print(f'{names[i]} vs {names[j]}: eps_row_l1 {eps:.2e}; bench rowwise {np.max(np.max(np.abs(dR),1)[np.max(np.abs(Rs[names[i]]),1)>0]/np.max(np.abs(Rs[names[i]]),1)[np.max(np.abs(Rs[names[i]]),1)>0]):.2e}')
        report(M.solve(dR), '   direct |M^-1 dR| (it14 M)')
for eps in (3.2e-12,):
    fl = ncc.noise_floor(M, lo, di, up, act, eps, 8)
    report(fl, f'model eps {eps:g}')
# pattern transfer to it24
B = ROOT/'results/nk167_hot/it24_trial_l1'
xb, Rb = load(B); lo2, di2, up2 = read_ba_blocks(B/'QSTARS_BA_TRACE.bin', nt, nd)
g2 = (np.abs(lo2)+np.abs(di2)+np.abs(up2)).sum(axis=2).T
M2 = BlockTridiag(lo2, di2, up2).factor(); act2 = Rb != 0
frac2, cls2 = ncc.classes(xb, Rb, di2, ncc.read_replace(B/'QSTARS_REPLACE_TRACE', nd), 1e-10)
for key, dR in dRs.items():
    r = np.where(act & act2, dR / np.where(gross > 0, gross, 1) * g2, 0.0)
    y = np.abs(M2.solve(r)); ya = ncc.audit_normalized(y)
    print('it24 transfer', key, {k: '%.2g' % (float(y[m].max()) if m.any() else 0) for k, m in cls2.items()}, 'n>1e-3', int(np.sum((ya > 1e-3) & act2)))
