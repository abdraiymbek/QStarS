#!/usr/bin/env python3
"""Report 27: progress in the convergence criterion itself: the unbounded Newton step (newton_step.npy of an iterate
directory) over significant variables (fraction > 1e-10, active), by depth region: max, median and the share of
variables already below the 1e-3 test. Read-only.   step_convergence.py ITERATE_DIR [ITERATE_DIR ...]"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_trace_record  # noqa: E402
REG = [(1, 9), (10, 18), (19, 32), (33, 47), (48, 54), (55, 61), (62, 100), (101, 167)]
for d in map(Path, sys.argv[1:]):
    base = d.parent / d.name.split('_from_', 1)[1]
    rd = base.with_name(base.name + '_tp') if (base.with_name(base.name + '_tp') / 'QSTARS_RESID_TRACE.bin').exists() else base
    x, R = read_trace_record(rd / 'QSTARS_RESID_TRACE.bin')
    s = np.load(d / 'newton_step.npy')
    H = x[1:32].sum(0); He = x[32:132].sum(0)
    frac = np.ones_like(x); frac[0:32] = x[0:32] / H; frac[32:132] = x[32:132] / He
    sig = (frac > 1e-10) & (R != 0)
    v = np.abs(s)[sig]
    print(f'{d.parent.name}/{d.name}: all significant: max {v.max():.2e}, median {np.median(v):.1e}, '
          f'below 1e-3 {np.mean(v < 1e-3):.1%} of {v.size}')
    print('   ' + '; '.join(f'd{a}-{b} med {np.median(np.abs(s[:, a-1:b])[sig[:, a-1:b]]):.0e} '
                           f'<1e-3 {np.mean(np.abs(s[:, a-1:b])[sig[:, a-1:b]] < 1e-3):.0%}' for a, b in REG))
