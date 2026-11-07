#!/usr/bin/env python3
"""Report 26: where the row-scaled residual of the stage-2 run sits (by depth and by variable group), in the driver's
acceptance measure (--merit rowscale --accept-rowscale --eq-floor 1e-15, RE row included), with the row weights and
scales of a reference iterate of free_r25. Read-only.

    resid_map.py REF_ITERATE_DIR PASS [PASS ...]   (directories under results/nk167_trap/free_r25)"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'rowscale', 1e-15, 1e-3
O = ROOT / 'results/nk167_trap/free_r25'
x0, _ = read_trace_record(ROOT / 'results/nk167_trap/fixedT_b35/it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
GROUPS = {'H-': [0], 'H I n=1/bal': [1], 'H I n=2': [2], 'H I n>=3': list(range(3, 31)), 'H num': [31],
          'He I': list(range(32, 101)), 'He II': list(range(101, 131)), 'He num': [131], 'charge': [132], 'RE': [133]}


def main():
    ref = O / sys.argv[1]
    base = ref.name.split('_from_')[1]
    cols = {n: (ref / f'col_{n}_p', ref / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    resdir = O / (base + '_tp') if (O / (base + '_tp')).exists() else O / base
    it = ck.Iterate(resdir, x0, cols, 1e-4, True, None, rj=O / base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    w = (it.row_weights() * it.D * it.active)[:134]
    for p in sys.argv[2:]:
        d = O / (p + '_tp') if (O / (p + '_tp')).exists() else O / p
        x, R = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
        r2 = (w * R[:134]) ** 2
        tot = r2.sum()
        print(f'{p}: ||W D R|| = {np.sqrt(tot):.5f}; front window (H II < 0.5, +2) '
              f'{tuple(v + 1 for v in cnk_trap.front_window(x, 0.5, 2))}')
        print('  share by group: ' + ', '.join(f'{g} {r2[i].sum() / tot:.3f}' for g, i in GROUPS.items()))
        dep = r2.sum(0) / tot
        order = np.argsort(dep)[::-1]
        print('  largest depths (1-based: share): ' + ', '.join(f'd{j + 1}:{dep[j]:.3f}' for j in order[:20]))
        for a, b in [(1, 9), (10, 18), (19, 32), (33, 47), (48, 54), (55, 61), (62, 100), (101, 167)]:
            print(f'  d{a}-{b}: {dep[a - 1:b].sum():.3f}  (groups: ' + ', '.join(
                f'{g} {r2[i, a - 1:b].sum() / tot:.3f}' for g, i in GROUPS.items() if r2[i, a - 1:b].sum() / tot > 0.01) + ')')
        cells = np.argsort(r2.ravel())[::-1][:15]
        print('  largest cells (row 1-based, depth 1-based, share): ' + ', '.join(
            f'({c // ck.ND + 1},{c % ck.ND + 1}) {r2.ravel()[c] / tot:.3f}' for c in cells))


if __name__ == '__main__':
    main()
