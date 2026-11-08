#!/usr/bin/env python3
"""Report 28: predicted fit of the bounded step against the radius of the major components (low-population levels at
factor 10, T at 0.06 x radius), from stored products: how large a step the remaining linear gain needs. Read-only."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
itdir, base = Path(sys.argv[1]), Path(sys.argv[2])
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(cn.residual_dir(base), x0b, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
it.extend_significance(0.0, rows=[3])
b = it.obj(it.R).ravel(); beta = np.linalg.norm(b)
S, JS = [], []
for q in sorted(itdir.glob('jv[0-9][0-9]p[0-9]')):
    z = np.load(q / 'pairs.npz'); S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
A = np.stack([it.obj(j).ravel() for j in JS], 1)
Sb, Ab, M, kb = ck.span_basis(S, A, 1e-6)
c0 = np.linalg.lstsq(Ab, b, rcond=1e-12)[0]
s0 = sum(ci * si for ci, si in zip(c0, Sb))
_, _, _, _, m = cnk_trap.lowpop_bounds(it.frac, it.active, 0.25, 1e-6, 10.0)
maj = it.active & ~m
print(f'unbounded fit {np.linalg.norm(Ab @ c0 - b) / beta:.4f}; unbounded step over major components: max {np.abs(s0[maj]).max():.2f}, '
      f'99th pct {np.quantile(np.abs(s0[maj]), 0.99):.3f}, median {np.median(np.abs(s0[maj])):.1e}; where max: row {np.unravel_index(np.argmax(np.where(maj, np.abs(s0), 0)), s0.shape)}')
for r in (0.25, 1.0, 3.0, 10.0, 30.0):
    rad = np.full(it.x.shape, r); rad[133] = min(0.06 * r, 0.5)
    lo, hi, knee, cap, m = cnk_trap.lowpop_bounds(it.frac, it.active, rad, 1e-6, 10.0)
    s, _, c = ck.bounded_step((Sb, Ab, b), it.active, (lo, hi), it.active, 'box')
    print(f'  major radius {r:g}: predicted fit {np.linalg.norm(Ab @ c - b) / beta:.4f}')
