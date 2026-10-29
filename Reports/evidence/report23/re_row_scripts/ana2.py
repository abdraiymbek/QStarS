#!/usr/bin/env python3
"""History of n=2 / RE rows over a run's accepted iterates (read-only)."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_ba_blocks, read_trace_record
T1 = 4e-10 * np.pi; PC = 6.62607015e-12; dE21 = 3.288087 - 0.822021
dirs = [Path(p) for p in sys.argv[1:]]
ds = [30, 33, 35, 37, 39, 41]
for D in dirs:
    if not (D / 'QSTARS_RESID_TRACE.bin').exists():
        print('missing', D); continue
    x, R = read_trace_record(D / 'QSTARS_RESID_TRACE.bin')
    NT, ND = x.shape
    line = f'{D.parent.name}/{D.name:22s}'
    have_ba = (D / 'QSTARS_BA_TRACE.bin').exists()
    if have_ba:
        lo, di, up = read_ba_blocks(D / 'QSTARS_BA_TRACE.bin', NT, ND)
        l1 = lambda i, k: np.abs(lo[k][i]).sum() + np.abs(di[k][i]).sum() + np.abs(up[k][i]).sum()
    print(line)
    for d in ds:
        k = d - 1
        s = f'   d{d}: T={x[133,k]:.4f} n2={x[2,k]:.3e} n2/n1={x[2,k]/x[1,k]:.3e} R3={R[2,k]:10.3e} R_T={R[133,k]:10.3e} T1*R_T={T1*R[133,k]:10.3e} hv21*R3={PC*dE21*R[2,k]:10.3e} ratio={T1*R[133,k]/(PC*dE21*R[2,k]):7.3f}'
        if have_ba:
            s += f' R3/l1={R[2,k]/l1(2,k):9.2e} R_T/l1={R[133,k]/l1(133,k):9.2e} l1_3/(n2*A21)={l1(2,k)/(x[2,k]*4.696e8):6.3f}'
        print(s)
