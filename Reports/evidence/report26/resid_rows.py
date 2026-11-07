#!/usr/bin/env python3
"""Report 26: the rows (1-based) that hold the row-scaled residual in each depth region of one free_r25 pass, with the
population fraction of the matching variable. Read-only.   resid_rows.py REF_ITERATE_DIR PASS"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import resid_map as rm  # noqa: E402  (same settings and reference weights)
ck = rm.ck


def main():
    ref = rm.O / sys.argv[1]
    base = ref.name.split('_from_')[1]
    cols = {n: (ref / f'col_{n}_p', ref / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    resdir = rm.O / (base + '_tp') if (rm.O / (base + '_tp')).exists() else rm.O / base
    it = ck.Iterate(resdir, rm.x0, cols, 1e-4, True, None, rj=rm.O / base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    w = (it.row_weights() * it.D * it.active)[:134]
    d = rm.O / (sys.argv[2] + '_tp') if (rm.O / (sys.argv[2] + '_tp')).exists() else rm.O / sys.argv[2]
    x, R = rm.read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
    r2 = (w * R[:134]) ** 2; tot = r2.sum()
    for a, b in [(1, 9), (19, 47), (55, 61), (62, 100), (101, 167)]:
        rs = r2[:, a - 1:b].sum(1) / tot
        top = np.argsort(rs)[::-1][:8]
        print(f'd{a}-{b} ({rs.sum():.3f}): ' + ', '.join(
            f'row {i + 1} {rs[i]:.3f} (frac {np.median(it.frac[i, a - 1:b]):.0e}, |rel| {np.sqrt(r2[i, a-1:b].max()):.1e})'
            for i in top))
    print('T (K) at d1,19,33,47,55,61,100,130,167:', ', '.join(f'{x[133, j - 1]:.0f}' for j in (1, 19, 33, 47, 55, 61, 100, 130, 167)))


if __name__ == '__main__':
    main()
