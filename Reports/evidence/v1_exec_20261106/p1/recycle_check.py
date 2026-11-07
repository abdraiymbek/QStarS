#!/usr/bin/env python3
"""V1 plan P2 (reuse test, offline): the dense-window columns of an EARLIER base recycled at this base (radiation part
from their own base, direct part exact here: the driver's --recycle-dense), the conservation-consistent box step built
on them, and the same step scored with this base's EXACT dense columns for the same cells (unit-cell directions, so the
cell identity maps one span onto the other). Read-only.   recycle_check.py ITERATE_DIR OLD_ITERATE_DIR"""
import sys, copy
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached, O  # noqa: E402
from qstars.nlte_newton import read_jv_rad  # noqa: E402
import box_exact  # noqa: E402

d, old = sys.argv[1], sys.argv[2]
P = cached(d)                                   # exact: dense first, then Krylov
it = P.it
cell = lambda z: tuple(np.argwhere(z != 0)[0]) if np.count_nonzero(z) == 1 else None  # noqa: E731
exact = {cell(P.S[i]): i for i in range(P.ndense)}
Sr, Jr = [], []
for q in sorted((O / old).glob('dense[0-9][0-9]')):
    zs = np.load(q / 'pairs.npz')['z']
    _, drad = read_jv_rad(q / 'QSTARS_JV_RAD.bin')
    for i in range(zs.shape[2]):
        Sr.append(zs[:, :, i]); Jr.append(it.Jdir.matvec(zs[:, :, i]) + drad[:, :, i])
nk = len(P.S) - P.ndense
R = copy.copy(P)
R.S = Sr + P.S[P.ndense:]
R.Sm = np.stack([s.ravel() for s in R.S], axis=1)
R.A = np.concatenate([np.stack([it.obj(j).ravel() for j in Jr], axis=1), P.A[:, P.ndense:]], axis=1)
R.ndense = len(Sr)
E = box_exact.consistency_rows(R)
R.mask = P.mask.copy(); R.mask[[131, 132]] = it.x[[131, 132]] > 0
y, status = box_exact.solve(R, E)
# the same step in this base's exact span: recycled cell c -> exact column of cell c (cells absent here are dropped)
ye = np.zeros(P.A.shape[1]); missing = 0
for i, z in enumerate(Sr):
    j = exact.get(cell(z))
    if j is None:
        missing += 1
    else:
        ye[j] = y[i] * z[cell(z)] / P.S[j][cell(z)]
ye[P.ndense:] = y[R.ndense:]
s_r, s_e = R.step(y), P.step(ye)
print(f'{d} with the dense columns of {old}: box-consistent fit predicted by the recycled columns {R.lin(y):.4f} ({status}); '
      f'the same step with this base\'s exact columns {P.lin(ye):.4f} ({missing} recycled cells absent here, step mismatch '
      f'{np.abs(s_r - s_e).max():.2e})')
np.save(Path(__file__).with_name(f'y_REC_{d}_from_{old}.npy'), y)
np.save(Path(__file__).with_name(f'ye_REC_{d}_from_{old}.npy'), ye)
