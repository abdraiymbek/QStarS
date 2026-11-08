#!/usr/bin/env python3
"""Report 27: which depth regions can the iterate's span reduce at all? For each region, the unbounded least squares of
that region's residual alone (other regions ignored), in the truncated orthonormal span basis, and the same with the
0.25 box; plus the global fit's effect per region. Read-only.   region_reach.py ITERATE_DIR BASE_DIR"""
import sys
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
b = it.obj(it.R); beta = np.linalg.norm(b)
S, JS = [], []
for q in sorted(itdir.glob('jv[0-9][0-9]p[0-9]')):
    z = np.load(q / 'pairs.npz'); S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
A = np.stack([it.obj(j).ravel() for j in JS], 1)
Sb, Ab, M, kb = ck.span_basis(S, A, 1e-6)
Ab3 = Ab.reshape(b.shape + (kb,))
cg = np.linalg.lstsq(Ab, b.ravel(), rcond=1e-12)[0]
rg = (b.ravel() - Ab @ cg).reshape(b.shape)
print(f'{itdir.name}: {len(S)} directions, basis {kb}; global unbounded fit {np.linalg.norm(rg) / beta:.4f}')
for a_, b_ in cn.REGIONS:
    sl = (slice(None), slice(a_ - 1, b_))
    br = b[sl].ravel(); nr = np.linalg.norm(br)
    Ar = Ab3[sl].reshape(-1, kb)
    cr = np.linalg.lstsq(Ar, br, rcond=1e-12)[0]
    own = np.linalg.norm(br - Ar @ cr) / nr
    print(f'  d{a_}-{b_}: share {nr**2 / beta**2:.3f}; region alone reducible to {own:.3f} of itself; '
          f'global fit leaves {np.linalg.norm(rg[sl]) / nr:.3f}')
