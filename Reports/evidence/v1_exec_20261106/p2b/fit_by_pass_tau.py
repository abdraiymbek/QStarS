#!/usr/bin/env python3
"""Council block: in-box fit by number of product passes for one v1_p2 iterate, in the truncated orthonormal span basis
(tau) with the bound the block uses; decides the --krylov depth for the next block (cut only if passes 3-8 add < 0.05).
Read-only.   fit_by_pass_tau.py ITERATE_DIR BASE_DIR [--bound 0.25 --bound-t 0.015 --tau 1e-6]"""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument('itdir', type=Path); ap.add_argument('base', type=Path)
ap.add_argument('--bound', type=float, default=0.25); ap.add_argument('--bound-t', type=float, default=0.015)
ap.add_argument('--tau', type=float, default=1e-6)
a = ap.parse_args()
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
cols = {n: (a.itdir / f'col_{n}_p', a.itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(cn.residual_dir(a.base), x0b, cols, 1e-4, True, None, rj=a.base / 'QSTARS_RJ.bin')
it.extend_significance(0.0, rows=[3])
b = it.obj(it.R).ravel(); beta = float(np.linalg.norm(b))
rad = np.full(it.x.shape, a.bound); rad[133] = a.bound_t
S, JS, rows = [], [], []
for k, pn in enumerate(sorted({p.name[:4] for p in a.itdir.glob('jv[0-9][0-9]p[0-9]')})):
    for q in sorted(a.itdir.glob(f'{pn}p[0-9]')):
        z = np.load(q / 'pairs.npz'); S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
    A = np.stack([it.obj(j).ravel() for j in JS], 1)
    Sb, Ab, M, kb = ck.span_basis(S, A, a.tau)
    cu = np.linalg.lstsq(Ab, b, rcond=1e-12)[0]
    _, _, cb = ck.bounded_step((Sb, Ab, b), it.active, rad, it.active, 'box')
    r = dict(passes=k + 1, directions=len(S), basis=kb, unbounded=float(np.linalg.norm(Ab @ cu - b) / beta),
             box=float(np.linalg.norm(Ab @ cb - b) / beta))
    rows.append(r)
    print(f'{k + 1} passes, {len(S)} directions ({kb} basis): unbounded {r["unbounded"]:.4f}, box {r["box"]:.4f}', flush=True)
Path(__file__).with_name(f'fit_by_pass_{a.itdir.name}.json').write_text(json.dumps(rows, indent=1))
