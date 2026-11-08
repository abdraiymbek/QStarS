#!/usr/bin/env python3
"""Report 28: predicted fit of the bounded step with factor-F bounds for low-population levels (cnk_trap.lowpop_bounds),
from one iterate's stored products (truncated span basis, box solver). Saves the steps for native tests. Read-only.
    lowpop_offline.py ITERATE_DIR BASE_DIR"""
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
print(f'{itdir.name}: basis {kb}; unbounded fit {np.linalg.norm(Ab @ c0 - b) / beta:.4f}')
rad = np.full(it.x.shape, 0.25); rad[133] = 0.015
steps = {}
for th, F in ((0.0, 1.0), (1e-8, 10.0), (1e-6, 10.0), (1e-6, 100.0), (1e-4, 10.0)):
    if th == 0:
        bnd, m = rad, np.zeros_like(it.active)
    else:
        lo, hi, knee, cap, m = cnk_trap.lowpop_bounds(it.frac, it.active, rad, th, F)
        bnd = (lo, hi)
    s, _, c = ck.bounded_step((Sb, Ab, b), it.active, bnd, it.active, 'box')
    tag = f'th{th:g}_F{F:g}'
    steps[tag] = s
    print(f'  {tag}: {int(m.sum())} low-population components; predicted fit {np.linalg.norm(Ab @ c - b) / beta:.4f}; '
          f'max|s| over them {np.abs(s[m]).max() if m.any() else 0:.2f}, elsewhere {np.abs(s[~m & it.active]).max():.3f}')
np.savez(Path(__file__).with_name(f'lowpop_steps_{itdir.name}.npz'), **steps)
