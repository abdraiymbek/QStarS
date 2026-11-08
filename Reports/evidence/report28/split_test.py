#!/usr/bin/env python3
"""Report 28: does splitting the preconditioned-residual pieces by population class give the least squares the freedom
it lacks? At it32 (the run's preconditioner: local ALI, no floor/PTC), the region x species-group pieces of W M^-1 R
(as --split) and the same pieces split into major / low-population (fraction < 1e-6, cnk_trap.lowpop_bounds mask)
parts, all made conservation-consistent; one exact product pass; fits of the residual by each set, alone and with the
it32 span (truncated basis), unbounded and in the 0.25 box with factor-10 low-population bounds."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'Reports/evidence/v1_exec_20261106/p1'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_jv_rad  # noqa: E402
from p1_native import EXE, REPLACE, X0DIR  # noqa: E402
RUN = ROOT / 'results/nk167_trap/v1_p2'
itdir, base = RUN / 'it32_from_it32_trial_l1', RUN / 'it32_trial_l1'
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
ck.SPLIT_REGIONS = ck.split_regions(167)
cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(cn.residual_dir(base), x0b, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
it.install_ali_increment(base.with_name(base.name + '_ali'))
it.extend_significance(0.0, rows=[3])
v = it.W * it.M.solve(it.R)
groups = ck.split_groups(True)
_, _, _, _, low = cnk_trap.lowpop_bounds(it.frac, it.active, 0.25, 1e-6, 10.0)
cons = lambda z: cnk_trap.project_linear(it.x, z) * it.active  # noqa: E731
lumped = [cons(p) for p in ck.split_pieces(v, ck.SPLIT_REGIONS, groups, 0.0)]
split = []
for p in ck.split_pieces(v, ck.SPLIT_REGIONS, groups, 0.0):
    for msk in (low, ~low):
        q = np.where(msk, p, 0.0)
        if np.any(q):
            split.append(cons(q))
lumped = [z for z in lumped if np.any(z)]
dirs = lumped + split
print(f'{len(lumped)} lumped pieces, {len(split)} split pieces')
out = RUN / 'report28_split_test'
R = ck.Runner(X0DIR, EXE, REPLACE, False, 4, RUN / 'report28_split_test.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
R.no_obs = True
R.run(out, it.s_cum, replay=it.rj, dirs=np.stack(dirs, 2), ba=False, store_j=False, two_phot_j=base / 'QSTARS_RJ.bin', threads=8)
_, drad = read_jv_rad(out / 'QSTARS_JV_RAD.bin')
J = [it.Jdir.matvec(z) + drad[:, :, i] for i, z in enumerate(dirs)]
b = it.obj(it.R).ravel(); beta = np.linalg.norm(b)
S, JS = [], []
for q in sorted(itdir.glob('jv[0-9][0-9]p[0-9]')):
    zz = np.load(q / 'pairs.npz'); S += [zz['z'][:, :, i] for i in range(zz['z'].shape[2])]; JS += [zz['jz'][:, :, i] for i in range(zz['z'].shape[2])]
rad = np.full(it.x.shape, 0.25); rad[133] = 0.015
lo, hi, _, _, _ = cnk_trap.lowpop_bounds(it.frac, it.active, rad, 1e-6, 10.0)


def fits(Sx, Jx, label):
    A = np.stack([it.obj(j).ravel() for j in Jx], 1)
    Sb, Ab, M, kb = ck.span_basis(Sx, A, 1e-6)
    c0 = np.linalg.lstsq(Ab, b, rcond=1e-12)[0]
    _, _, cb = ck.bounded_step((Sb, Ab, b), it.active, (lo, hi), it.active, 'box')
    print(f'{label}: {len(Sx)} directions ({kb} basis): unbounded {np.linalg.norm(Ab @ c0 - b) / beta:.4f}, '
          f'bounded {np.linalg.norm(Ab @ cb - b) / beta:.4f}', flush=True)


nl = len(lumped)
fits(dirs[:nl], J[:nl], 'lumped pieces alone')
fits(dirs[nl:], J[nl:], 'split pieces alone')
fits(S, JS, 'it32 span')
fits(S + dirs[:nl], JS + J[:nl], 'it32 span + lumped')
fits(S + dirs[nl:], JS + J[nl:], 'it32 span + split')
