#!/usr/bin/env python3
"""Report 27: can a preconditioner give directions that reduce the d55-61 and d62-100 residual, which the v1_p2 span no
longer reduces (it26: d62-100 alone reducible only to 0.863 of itself)? Directions: the region x species-group pieces of
W M^-1 R (as the driver's --split pieces), made conservation-consistent, for three preconditioners at it26:
  fix  fixed-J direct blocks (PC_T_FLOOR as the run);  ali  + local-ALI increment;  cur  + PTC 0.03 (pc) + SVD floor 0.01 (the run's).
One exact product pass (RUN/report27_pc_test); then, per region, the least-squares reduction of that region's own residual
by each set alone and by the it26 span plus each set. Passes only at nice 10."""
import sys, os
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
itdir, base = RUN / 'it26_from_it26_trial_l0.5', RUN / 'it26_trial_l0.5'
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}


def mk():
    it = ck.Iterate(cn.residual_dir(base), x0b, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    return it


its = {}
its['fix'] = mk()
it = mk(); it.install_ali_increment(base.with_name(base.name + '_ali')); its['ali'] = it
it = mk(); it.install_ali_increment(base.with_name(base.name + '_ali')); it.add_pop_ptc(0.03, scale='pc', jacobian=False)
it.floor_preconditioner(0.01); it.extend_significance(0.0, rows=[3]); its['cur'] = it
REG = {'d55-61': (54, 61), 'd62-100': (61, 100)}
groups = ck.split_groups(True)
dirs, meta = [], []
ref = its['fix']
for name, itp in its.items():
    v = itp.W * itp.M.solve(itp.R)
    for rn, (d0, d1) in REG.items():
        for g0, g1 in groups:
            z = np.zeros_like(v); z[g0:g1, d0:d1] = v[g0:g1, d0:d1]
            z = cnk_trap.project_linear(ref.x, z) * ref.active
            if np.any(z):
                dirs.append(z); meta.append((name, rn, f'{g0}-{g1}'))
print(f'{len(dirs)} directions')
out = RUN / 'report27_pc_test'
R = ck.Runner(X0DIR, EXE, REPLACE, False, 3, RUN / 'report27_pc_test.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
R.no_obs = True
R.run(out, ref.s_cum, replay=ref.rj, dirs=np.stack(dirs, 2), ba=False, store_j=False, two_phot_j=base / 'QSTARS_RJ.bin', threads=3)
_, drad = read_jv_rad(out / 'QSTARS_JV_RAD.bin')
J = [ref.obj(ref.Jdir.matvec(z) + drad[:, :, i]) for i, z in enumerate(dirs)]
b = ref.obj(ref.R)
S, JS = [], []
for q in sorted(itdir.glob('jv[0-9][0-9]p[0-9]')):
    zz = np.load(q / 'pairs.npz'); S += [zz['z'][:, :, i] for i in range(zz['z'].shape[2])]; JS += [zz['jz'][:, :, i] for i in range(zz['z'].shape[2])]
Asp = np.stack([ref.obj(j).ravel() for j in JS], 1)
_, Ab, _, _ = ck.span_basis(S, Asp, 1e-6)
for rn, (d0, d1) in REG.items():
    sl = (slice(None), slice(d0, d1))
    br = b[sl].ravel(); nr = np.linalg.norm(br)
    rows = []
    for name in its:
        idx = [i for i, m in enumerate(meta) if m[0] == name]
        Ar = np.stack([J[i][sl].ravel() for i in idx], 1)
        c = np.linalg.lstsq(Ar, br, rcond=1e-12)[0]
        Ac = np.concatenate([Ab.reshape(b.shape + (-1,))[sl].reshape(len(br), -1), Ar], 1)
        c2 = np.linalg.lstsq(Ac, br, rcond=1e-12)[0]
        rows.append(f'{name}: alone {np.linalg.norm(br - Ar @ c) / nr:.3f}, with span {np.linalg.norm(br - Ac @ c2) / nr:.3f}')
    print(f'{rn} (residual reducible to, fraction of itself): ' + '; '.join(rows))
