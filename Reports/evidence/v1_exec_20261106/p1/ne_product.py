#!/usr/bin/env python3
"""V1 plan P1: exact product of the projection's correction to the BOXP step (P S y - S y; its n_e part alone; the rest)
at the it16 base, to test whether the radiation response of the n_e change is the first-order term the model missed.
One native product pass (3 directions) in results/nk167_trap/v1_p1/<base>/ne_product."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached, ck, read_trace_record, ROOT  # noqa: E402
from p1_native import EXE, REPLACE, X0DIR  # noqa: E402
from box_proj import projected_problem  # noqa: E402
from qstars.nlte_newton import read_jv_rad  # noqa: E402

d = sys.argv[1]
P = cached(d); it = P.it; Q = projected_problem(P)
y = np.load(Path(__file__).with_name(f'y_BOXP_{d}.npy'))
corr = Q.step(y) - P.step(y)
ne = np.zeros_like(corr); ne[132] = corr[132]
rest = corr - ne
zs = [corr, ne, rest]
out = ROOT / 'results/nk167_trap/v1_p1' / P.base.name / 'ne_product'
R = ck.Runner(X0DIR, EXE, REPLACE, False, 8, out.parent / 'ne_product.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
R.no_obs = True
tpj = P.base / 'QSTARS_RJ.bin'
R.run(out, it.s_cum, replay=it.rj, dirs=np.stack(zs, axis=2), ba=False, store_j=False, two_phot_j=tpj, threads=8)
_, drad = read_jv_rad(out / 'QSTARS_JV_RAD.bin')
np.savez(out / 'pairs.npz', z=np.stack(zs, 2), jz=np.stack([it.Jdir.matvec(z) + drad[:, :, i] for i, z in enumerate(zs)], 2))
for i, n in enumerate(['correction', 'n_e part', 'rest']):
    dr = it.obj(drad[:, :, i])
    print(f'{n}: ||obj(Jdir z)|| {np.linalg.norm(it.obj(it.Jdir.matvec(zs[i]))[:134]):.4e}, ||obj(radiation part)|| '
          f'{np.linalg.norm(dr[:134]):.4e} (row 2: {np.linalg.norm(dr[1]):.4e}) relative to the base merit {P.beta:.4e}')
# the native residual change of the lambda trials against the model with the exact correction product
Jc = it.obj(it.Jdir.matvec(corr) + drad[:, :, 0])
JSy = (P.A @ y).reshape(P.shape)
for lam in (0.25, 0.5, 1.0):
    xn, Rn = read_trace_record(ROOT / f'results/nk167_trap/v1_p1/{P.base.name}/boxp_l{lam:g}/QSTARS_RESID_TRACE.bin')
    pred_exact = P.b.reshape(P.shape) - lam * (JSy + Jc)
    pred_dir = (Q.b - lam * Q.A @ y).reshape(P.shape)
    nat = it.obj(Rn)
    print(f'lambda {lam:g}: native {np.linalg.norm(nat[:134]) / P.beta:.4f}, model with the exact correction product '
          f'{np.linalg.norm(pred_exact[:134]) / P.beta:.4f}, model with its direct part only {np.linalg.norm(pred_dir[:134]) / P.beta:.4f}; '
          f'native - exact-model {np.linalg.norm((nat - pred_exact)[:134]) / P.beta:.4f}')
