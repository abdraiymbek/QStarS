#!/usr/bin/env python3
"""V1 plan P2: where a native trial departs from the linear model of its step (rows, depths), for the v1_p2 products.
    p2_gap.py RUN ITER BASE PASSES TRIAL_DIR LAM [solver]. Read-only."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
run, name, base, npass, tdir, lam = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), int(sys.argv[4]), Path(sys.argv[5]), float(sys.argv[6])
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
itdir = run / name
cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(cn.residual_dir(base), x0b, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
it.extend_significance(0.0, rows=[3])
b = it.obj(it.R).ravel(); beta = np.linalg.norm(b)
passes = sorted({p.name[:4] for p in itdir.glob('jv[0-9][0-9]p[0-9]')})[:npass]
S, JS = [], []
for pn in passes:
    for q in sorted(itdir.glob(f'{pn}p[0-9]')):
        z = np.load(q / 'pairs.npz'); S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
A = np.stack([it.obj(j).ravel() for j in JS], 1)
rad = np.full(it.x.shape, 0.5); rad[133] = 0.03
s, _, y = ck.bounded_step((S, A, b), it.active, rad, it.active, sys.argv[7] if len(sys.argv) > 7 else 'qr')
xt, Rt = read_trace_record(tdir / 'QSTARS_RESID_TRACE.bin')
s_real = 1 - xt / it.x
pred = (b - lam * A @ y).reshape(it.x.shape) - it.obj(it.Jdir.matvec(np.where(it.active, s_real - lam * s, 0)))
nat = it.obj(Rt)
e2 = ((nat - pred)[:134]) ** 2; tot = e2.sum()
print(f'{tdir.name}: native {np.linalg.norm(nat[:134]) / beta:.4f}, predicted {np.linalg.norm(pred[:134]) / beta:.4f}, '
      f'||native - predicted|| / merit {np.sqrt(tot) / beta:.4f}')
print('  by group: ' + ', '.join(f'{g} {e2[r].sum() / tot:.3f}' for g, r in cn.GROUPS.items() if e2[r].sum() / tot > 0.005))
print('  by region: ' + ', '.join(f'd{a_}-{b_} {e2[:, a_ - 1:b_].sum() / tot:.3f}' for a_, b_ in cn.REGIONS))
c = np.argsort(e2.ravel())[::-1][:10]
print('  cells (row, depth, share, native, predicted, base): ' + '; '.join(
    f'({i // ck.ND + 1},{i % ck.ND + 1}) {e2.ravel()[i] / tot:.3f} {nat.ravel()[i]:+.3f} {pred.ravel()[i]:+.3f} {b[i]:+.3f}' for i in c))
j = c[0] % ck.ND
print(f'  step at depth {j + 1}: ' + ', '.join(f'{g} {np.abs(lam * s[r, j]).max():.3f}' for g, r in cn.GROUPS.items()))
