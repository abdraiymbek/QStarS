#!/usr/bin/env python3
"""Transfer a converged 84-radius cool state onto the nested 167-radius grid.

The 167-radius restart bisects every 84-radius interval (shared radii, density
and velocity are identical; verified). Populations, n_e and T are copied at the
shared radii and interpolated linearly in (ln r, ln x) at the midpoints. The
result is written as the perturbation-hook step relative to the 167 restart's
own state (POPS_target = POPS_restart*(1-s)); nothing physical is fitted. The
84-radius equation-selection table is mapped to 167 depths (midpoint k uses
the selection of its outer shared neighbour).
"""
import argparse, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.cmfgen_backend import read_rvtj  # noqa: E402
from qstars.nlte_newton import read_trace_record, write_fortran  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('state84', type=Path, help='base pass dir holding the 84-radius state (RESID trace)')
    ap.add_argument('restart84', type=Path, help='84-radius restart (for radii)')
    ap.add_argument('base167', type=Path, help='unperturbed 167-radius base pass (RESID trace = restart state)')
    ap.add_argument('restart167', type=Path)
    ap.add_argument('replace84', type=Path)
    ap.add_argument('out', type=Path)
    a = ap.parse_args()
    x84, _ = read_trace_record(a.state84 / 'QSTARS_RESID_TRACE.bin')
    x167, _ = read_trace_record(a.base167 / 'QSTARS_RESID_TRACE.bin')
    r84 = read_rvtj(a.restart84 / 'RVTJ')['radius'][::-1]       # native order: outer first
    r167 = read_rvtj(a.restart167 / 'RVTJ')['radius'][::-1]
    if not np.array_equal(r167[::2], r84):
        raise SystemExit('167 grid does not nest the 84 grid')
    nt = x84.shape[0]
    tgt = np.empty_like(x167)
    tgt[:, ::2] = x84
    lr = np.log(r167)
    for k in range(1, len(r167), 2):
        j = k // 2
        w = (lr[k] - lr[k - 1]) / (lr[k + 1] - lr[k - 1])
        tgt[:, k] = np.exp((1 - w) * np.log(x84[:, j]) + w * np.log(x84[:, j + 1]))
    a.out.mkdir(parents=True, exist_ok=True)
    write_fortran(a.out / 'QSTARS_STEP_167.bin', 1.0 - tgt / x167)
    np.save(a.out / 'target_state_167.npy', tgt)
    rows = [l.split() for l in a.replace84.read_text().splitlines() if l.strip()]
    flags = {int(r[0]): r[1:] for r in rows}
    with open(a.out / 'replace_167.txt', 'w') as f:
        for k in range(1, len(r167) + 1):
            j = (k + 1) // 2
            f.write(f'{k:5d}' + ''.join(f'{int(v):2d}' for v in flags[j]) + '\n')
    print('target min/max', tgt.min(), tgt.max(), ' step range', (1 - tgt / x167).min(), (1 - tgt / x167).max())


if __name__ == '__main__':
    main()
