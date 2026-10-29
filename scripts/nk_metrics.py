#!/usr/bin/env python3
"""Iterate-independent convergence metrics for a cool-closure Newton-Krylov path.

For each accepted base (it00_base, then the accepted trials in order) report:
  ||M0^-1 R||  with the FIXED x0 direct Jacobian M0 (comparable across iterates),
  max |M0^-1 R| over levels holding > 1e-6 of their element (significant levels),
  raw ||R||, and the same split into depth regions.
"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_trace_record, read_ba_blocks, BlockTridiag
NT, ND = 134, 84
REG = {'outer d1-9': (0, 9), 'cool d10-30': (9, 30), 'mid d31-60': (30, 60), 'inner d61-84': (60, 84)}


def significant(x):
    H = x[1:32].sum(0); He = x[32:132].sum(0)
    f = np.ones_like(x)
    f[0:32] = x[0:32] / H; f[32:132] = x[32:132] / He
    return f > 1e-6


def main(run, bases):
    run = Path(run)
    lo, di, up = read_ba_blocks(run / 'it00_base' / 'QSTARS_BA_TRACE.bin', NT, ND)
    M0 = BlockTridiag(lo, di, up).factor()
    for b in bases:
        x, R = read_trace_record(run / b / 'QSTARS_RESID_TRACE.bin')
        y = M0.solve(R); y[R == 0] = 0
        sig = significant(x) & (R != 0)
        reg = {k: float(np.max(np.abs(y[:, a:c][sig[:, a:c]]))) for k, (a, c) in REG.items()}
        print(f'{b:24s} ||R||={np.linalg.norm(R):.3e} ||M0^-1R||={np.linalg.norm(y):.3e} '
              f'max_sig|M0^-1R|={np.max(np.abs(y[sig])):.3e}  ' + ' '.join(f'{k}:{v:.2e}' for k, v in reg.items()))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
