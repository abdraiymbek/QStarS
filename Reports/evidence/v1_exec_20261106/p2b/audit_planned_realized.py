#!/usr/bin/env python3
"""V1 plan, council block step 0: planned vs realized vs native for the stored trials of one v1_p2 iterate. Rebuilds
the driver's bounded step (span basis tau, box/QR best, |s| <= 0.5, T 0.03) from the iterate's stored products, checks
it against the logged linear residual, then for each trial lambda: the realized update minus lambda*s (by group), and
the native residual against the linear model of the REALIZED update (exact products for lambda*s, direct blocks for the
difference), by region and row group. Read-only.
    audit_planned_realized.py ITERATE_DIR BASE_DIR [--tau 1e-6] [--bound 0.5 --bound-t 0.03]"""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument('itdir', type=Path); ap.add_argument('base', type=Path)
ap.add_argument('--tau', type=float, default=1e-6); ap.add_argument('--bound', type=float, default=0.5)
ap.add_argument('--bound-t', type=float, default=0.03)
a = ap.parse_args()
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
cols = {n: (a.itdir / f'col_{n}_p', a.itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(cn.residual_dir(a.base), x0b, cols, 1e-4, True, None, rj=a.base / 'QSTARS_RJ.bin')
it.extend_significance(0.0, rows=[3])
b = it.obj(it.R).ravel(); beta = float(np.linalg.norm(b))
S, JS = [], []
for q in sorted(a.itdir.glob('jv[0-9][0-9]p[0-9]')):
    z = np.load(q / 'pairs.npz'); S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
A = np.stack([it.obj(j).ravel() for j in JS], 1)
Sb, Ab, M, kb = ck.span_basis(S, A, a.tau)
rad = np.full(it.x.shape, a.bound); rad[133] = a.bound_t
s, _, c = ck.bounded_step((Sb, Ab, b), it.active, rad, it.active, 'box')
Js = (Ab @ c).reshape(it.x.shape)                      # obj(J s)
out = dict(itdir=str(a.itdir), basis=kb, directions=len(S), planned_fit_l1=float(np.linalg.norm(Ab @ c - b) / beta))
print(f'{a.itdir.name}: {len(S)} directions, basis {kb}; rebuilt bounded step: linear fit {out["planned_fit_l1"]:.4f}')
k = int(a.itdir.name[2:4])
for tdir in sorted(a.itdir.parent.glob(f'it{k + 1:02d}_trial_l*')):
    if tdir.name.endswith(('_tp', '_ali')) or not (tdir / 'QSTARS_RESID_TRACE.bin').exists():
        continue
    lam = float(tdir.name.split('_l')[-1])
    xt, Rt = read_trace_record(tdir / 'QSTARS_RESID_TRACE.bin')
    sr = np.where(it.active, 1.0 - xt / it.x, 0.0)
    d = sr - lam * s
    plan = b.reshape(it.x.shape) - lam * Js
    real = plan - it.obj(it.Jdir.matvec(d))
    nat = it.obj(Rt)
    gap = (nat - real)[:134] ** 2; gt = gap.sum()
    row = dict(lam=lam, planned=float(np.linalg.norm(plan[:134]) / beta), realized=float(np.linalg.norm(real[:134]) / beta),
               native=float(np.linalg.norm(nat[:134]) / beta), gap=float(np.sqrt(gt) / beta),
               max_d=float(np.abs(d).max()), d_by_group={g: float(np.abs(d[r]).max()) for g, r in cn.GROUPS.items()},
               gap_regions={f'd{a_}-{b_}': float(gap[:, a_ - 1:b_].sum() / gt) for a_, b_ in cn.REGIONS},
               gap_groups={g: float(gap[r].sum() / gt) for g, r in cn.GROUPS.items()})
    out[f'l{lam:g}'] = row
    print(f'  lambda {lam:g}: planned {row["planned"]:.4f}  realized {row["realized"]:.4f}  NATIVE {row["native"]:.4f}  '
          f'|native - realized model| {row["gap"]:.4f}; max|realized - planned step| {row["max_d"]:.2e} '
          f'({max(row["d_by_group"], key=row["d_by_group"].get)})')
    print('     gap by region: ' + ', '.join(f'{k_} {v:.2f}' for k_, v in row['gap_regions'].items() if v > 0.02)
          + ' | by group: ' + ', '.join(f'{k_} {v:.2f}' for k_, v in row['gap_groups'].items() if v > 0.02))
Path(__file__).with_name(f'audit_{a.itdir.name}.json').write_text(json.dumps(out, indent=1))
