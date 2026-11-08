#!/usr/bin/env python3
"""Report 27: bounded-step fit vs the bound, from one iterate's stored products (truncated span basis). Variants: uniform
radius r (T r x 0.06), and radius r with near-trace variables (fraction < 1e-6 of the element) left to the component
limiter (bound 0.5, the exact range of the linear update map). Saves each step for native tests. Read-only."""
import sys, json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
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
steps = {}
import os
TS = [float(v) for v in os.environ.get('RT_T', '').split()] or [None]
for r in (0.25, 0.35, 0.5):
  for tr in TS:
    for trace_free in (False, True):
        rad = np.full(it.x.shape, r); rad[133] = r * 0.06 if tr is None else tr
        if trace_free:
            rad = np.where(it.frac < 1e-6, np.maximum(rad, 0.5), rad)
        s, _, c = ck.bounded_step((Sb, Ab, b), it.active, rad, it.active, 'box')
        tag = f'r{r:g}_T{rad[133, 0]:g}{"_tracefree" if trace_free else ""}'
        steps[tag] = s
        bind = np.abs(s) / rad > 0.999
        print(f'  {tag}: predicted fit {np.linalg.norm(Ab @ c - b) / beta:.4f}; components at the bound {int((bind & it.active).sum())}, '
              f'of them near-trace (<1e-6) {int((bind & it.active & (it.frac < 1e-6)).sum())}')
np.savez(Path(__file__).with_name(f'radius_steps_{itdir.name}.npz'), **steps)
