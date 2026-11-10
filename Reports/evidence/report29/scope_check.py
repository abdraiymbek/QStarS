#!/usr/bin/env python3
"""Report 29, step 1 (option 3, analysis only): at the best state, which levels fall below a population-fraction
threshold, what share of the residual (frozen it16 measure) and of the needed Newton corrections they carry, and how
close the remaining set is to the convergence test (|s| < 1e-3). Flags levels below the threshold that are NOT safe to
hand to a fixed-J cleanup because they are radiatively coupled: H I n=2 (Lya-trapped), and any level whose fixed-J
diagonal differs from the local-ALI diagonal by more than a factor 2 (the local radiation response matters). Read-only.
    scope_check.py ITERATE_DIR BASE_DIR"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
from qstars.nlte_newton import read_trace_record, read_ba_blocks  # noqa: E402
itdir, base = Path(sys.argv[1]), Path(sys.argv[2])
cn.configure()
w = np.load(ROOT / 'Reports/evidence/v1_exec_20261106/frozen_norm_it16.npz')['w']
x, R = read_trace_record(cn.residual_dir(base) / 'QSTARS_RESID_TRACE.bin')
s = np.load(itdir / 'newton_step.npy')
_, dfix, _ = read_ba_blocks(cn.residual_dir(base) / 'QSTARS_BA_TRACE.bin', 134, 167)
ali = base.with_name(base.name + '_ali')
dali = read_ba_blocks(ali / 'QSTARS_BA_TRACE.bin', 134, 167)[1] if (ali / 'QSTARS_BA_TRACE.bin').exists() else None
H = x[1:32].sum(0); He = x[32:132].sum(0)
frac = np.ones_like(x); frac[0:32] = x[0:32] / H; frac[32:132] = x[32:132] / He
act = R != 0
r2 = (w * R[:134]) ** 2; tot = r2.sum()
sig = act & (frac > 1e-10)
print(f'{base.name}: frozen {np.sqrt(tot):.4e}; significant variables (fraction > 1e-10) {int(sig.sum())}')
for th in (1e-8, 1e-7, 1e-6):
    low = sig & (frac < th)
    keep = sig & ~low
    diag_f = np.abs(np.einsum('kii->ik', dfix))
    coupled = np.zeros_like(low)
    if dali is not None:
        diag_a = np.abs(np.einsum('kii->ik', dali))
        ratio = np.where(diag_a > 0, diag_f / np.maximum(diag_a, 1e-300), 1.0)
        coupled = low & ((ratio > 2) | (ratio < 0.5))
    coupled[2] |= low[2]                                     # H I n=2 always counts as radiatively coupled
    v_keep, v_low = np.abs(s[keep]), np.abs(s[low])
    print(f' threshold {th:g}: {int(low.sum())} levels below it carry {r2[:134][low[:134]].sum() / tot:.1%} of the residual^2; '
          f'of them {int(coupled.sum())} radiatively coupled (fixed-J vs local-ALI diagonal off by > x2, or H I n=2)')
    print(f'   below: step median {np.median(v_low):.1e}, <1e-3 {np.mean(v_low < 1e-3):.0%} | remaining {int(keep.sum())}: '
          f'median {np.median(v_keep):.1e}, 90% {np.quantile(v_keep, 0.9):.1e}, max {v_keep.max():.2e}, <1e-3 {np.mean(v_keep < 1e-3):.0%}')
