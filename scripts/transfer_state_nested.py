#!/usr/bin/env python3
"""Transfer a cool state (and its equation-selection table) onto a grid that bisects it (Report 21).

Generalizes scripts/transfer_state_167.py (84 -> 167) to any coarse/fine pair with fine[::2] = coarse
(167 -> 333). Populations, n_e and T are copied at the shared radii and interpolated linearly in
(ln r, ln x) at the midpoints; the result is the perturbation-hook step relative to the fine
restart's own state (POPS_target = POPS_restart*(1-s)). Nothing physical is fitted. The table
maps midpoint k to the selection of its outer shared neighbour, as for 84 -> 167.
Run with --table-only to write only the fine table (before a converged coarse state exists).
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
    ap.add_argument('restart_coarse', type=Path)
    ap.add_argument('restart_fine', type=Path)
    ap.add_argument('replace_coarse', type=Path)
    ap.add_argument('out', type=Path)
    ap.add_argument('--state', type=Path, help='coarse base pass holding the state (RESID trace)')
    ap.add_argument('--base-fine', type=Path, help='unperturbed fine base pass (RESID trace = fine restart state)')
    ap.add_argument('--table-only', action='store_true')
    a = ap.parse_args()
    rc = read_rvtj(a.restart_coarse / 'RVTJ')['radius'][::-1]       # native order: outer first
    rf = read_rvtj(a.restart_fine / 'RVTJ')['radius'][::-1]
    if not np.array_equal(rf[::2], rc):
        raise SystemExit('fine grid does not nest the coarse grid')
    nf = len(rf)
    a.out.mkdir(parents=True, exist_ok=True)
    rows = [l.split() for l in a.replace_coarse.read_text().splitlines() if l.strip()]
    flags = {int(r[0]): r[1:] for r in rows}
    with open(a.out / f'replace_{nf}.txt', 'w') as f:
        for k in range(1, nf + 1):
            f.write(f'{k:5d}' + ''.join(f'{int(v):2d}' for v in flags[(k + 1) // 2]) + '\n')
    if a.table_only:
        return
    for t in (a.state / 'QSTARS_RESID_TRACE.bin', a.base_fine / 'QSTARS_RESID_TRACE.bin'):
        # record 0 is what is read: a record inherited from a restart (results/nk333/superseded) is refused
        hdr = np.frombuffer(t.read_bytes()[:12], '<i4')
        if t.stat().st_size != 12 + 16 * int(hdr[1]) * int(hdr[2]):
            raise SystemExit(f'{t}: not a single-record trace')
    xc, _ = read_trace_record(a.state / 'QSTARS_RESID_TRACE.bin')
    xf, _ = read_trace_record(a.base_fine / 'QSTARS_RESID_TRACE.bin')
    tgt = np.empty_like(xf)
    tgt[:, ::2] = xc
    lr = np.log(rf)
    for k in range(1, nf, 2):
        j = k // 2
        w = (lr[k] - lr[k - 1]) / (lr[k + 1] - lr[k - 1])
        tgt[:, k] = np.exp((1 - w) * np.log(xc[:, j]) + w * np.log(xc[:, j + 1]))
    write_fortran(a.out / f'QSTARS_STEP_{nf}.bin', 1.0 - tgt / xf)
    np.save(a.out / f'target_state_{nf}.npy', tgt)
    (a.out / 'source_state.txt').write_text(f'{a.state.resolve()}\n')
    print('target min/max', tgt.min(), tgt.max(), ' step range', (1 - tgt / xf).min(), (1 - tgt / xf).max())


if __name__ == '__main__':
    main()
