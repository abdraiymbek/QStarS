#!/usr/bin/env python3
"""Report 26: what exact front columns could do inside the least squares, from the stored exact products of one
free_r25 iterate (no new passes). The colored probes (one row on every 3rd depth of the window) are sums of the dense
single-depth columns, so the least squares over Krylov + colored probes is a lower bound on what the dense columns
give at that iterate. Prints the relative fit and the remaining residual share by depth region. Read-only.

    span_check.py ITERATE_DIR    (e.g. it03_from_it03_trial_l0.125 under results/nk167_trap/free_r25)"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'rowscale', 1e-15, 1e-3
O = ROOT / 'results/nk167_trap/free_r25'
x0, _ = read_trace_record(ROOT / 'results/nk167_trap/fixedT_b35/it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
REG = [(1, 9), (10, 18), (19, 32), (33, 47), (48, 54), (55, 61), (62, 100), (101, 167)]


def pairs(dirs):
    S, JS = [], []
    for d in dirs:
        p = np.load(d / 'pairs.npz')
        S += [p['z'][:, :, i] for i in range(p['z'].shape[2])]; JS += [p['jz'][:, :, i] for i in range(p['z'].shape[2])]
    return S, JS


def main():
    itd = O / sys.argv[1]
    base = itd.name.split('_from_')[1]
    cols = {n: (itd / f'col_{n}_p', itd / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    resdir = O / (base + '_tp') if (O / (base + '_tp')).exists() else O / base
    it = ck.Iterate(resdir, x0, cols, 1e-4, True, None, rj=O / base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    bo = it.obj(it.R)
    beta = np.linalg.norm(bo)
    kry = pairs(sorted(itd.glob('jv*')))
    fr = pairs(sorted(itd.glob('front*')))
    r2_0 = (bo ** 2)

    def fit(S, JS, label):
        A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
        y, *_ = np.linalg.lstsq(A, bo.ravel(), rcond=1e-12)
        r = (bo.ravel() - A @ y).reshape(bo.shape)
        s = sum(yi * si for yi, si in zip(y, S))
        print(f'{label}: {len(S)} directions, relative residual {np.linalg.norm(r) / beta:.4f}, max|s| {np.abs(s).max():.3f}')
        print('   region share of the base residual^2 -> after the fit: ' + ', '.join(
            f'd{a}-{b} {r2_0[:, a - 1:b].sum() / beta**2:.3f}->{(r[:, a - 1:b]**2).sum() / beta**2:.3f}' for a, b in REG))
    fit(*kry, 'Krylov')
    if fr[0]:
        fit(kry[0] + fr[0], kry[1] + fr[1], 'Krylov + colored front probes')
        fit(*fr, 'colored front probes alone')


if __name__ == '__main__':
    main()
