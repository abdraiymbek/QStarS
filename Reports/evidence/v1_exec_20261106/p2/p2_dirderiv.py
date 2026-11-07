#!/usr/bin/env python3
"""V1 plan P1 item 3 / P2: shrinking-step directional derivative on the actual update path of the v1_p2 5-pass bounded
step s at it16: (a) one exact product pass of s itself and of its pieces (populations without n_e/T, n_e, T), against
the sum of the stored products sum y_i J z_i; (b) native trials at small lambda: (R(x_lam) - R) / lam against J s_real.
Passes in RUN/dirderiv/. Run while nothing else needs exact timing."""
import sys, json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'p1'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record, read_jv_rad  # noqa: E402
from p1_native import realize, EXE, REPLACE, X0DIR  # noqa: E402

run, name, base = ROOT / 'results/nk167_trap/v1_p2', 'it16_from_it16_trial_l0.5', ROOT / 'results/nk167_trap/free_r25/it16_trial_l0.5'
x0b = cn.configure(); ck.STRUCT_ACTIVE = True
itdir = run / name
cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(cn.residual_dir(base), x0b, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
it.extend_significance(0.0, rows=[3])
b = it.obj(it.R).ravel(); beta = np.linalg.norm(b)
passes = sorted({p.name[:4] for p in itdir.glob('jv[0-9][0-9]p[0-9]')})[:5]
S, JS = [], []
for pn in passes:
    for q in sorted(itdir.glob(f'{pn}p[0-9]')):
        z = np.load(q / 'pairs.npz'); S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
A = np.stack([it.obj(j).ravel() for j in JS], 1)
rad = np.full(it.x.shape, 0.5); rad[133] = 0.03
s, _, y = ck.bounded_step((S, A, b), it.active, rad, it.active, 'qr')
Jsum = sum(yi * ji for yi, ji in zip(y, JS))
pieces = {'s': s}
p_ = s.copy(); p_[[132, 133]] = 0.0; pieces['pops'] = p_
p_ = np.zeros_like(s); p_[132] = s[132]; pieces['ne'] = p_
p_ = np.zeros_like(s); p_[133] = s[133]; pieces['T'] = p_
out = run / 'dirderiv'; out.mkdir(exist_ok=True)
R = ck.Runner(X0DIR, EXE, REPLACE, False, 3, out / 'dirderiv.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
R.no_obs = True
tpj = base / 'QSTARS_RJ.bin'
x0, _ = read_trace_record(ROOT / 'results/nk167_trap/free_r25/it00_base/QSTARS_RESID_TRACE.bin')
ct = cnk_trap.conservation_totals(x0)
lams = (0.02, 0.005)
trial_states = {lam: realize(it, lam * s, x0, ct)[0] for lam in lams}


def job(k):
    if k == 'prod':
        R.run(out / 'product', it.s_cum, replay=it.rj, dirs=np.stack(list(pieces.values()), 2), ba=False, store_j=False,
              two_phot_j=tpj, threads=3)
        return k
    R.run(out / f'trial_l{k:g}', 1.0 - trial_states[k] / x0, threads=2, warm=tpj, two_phot_j=tpj)
    return k


with ThreadPoolExecutor(3) as ex:
    list(ex.map(job, ['prod', *lams]))
_, drad = read_jv_rad(out / 'product' / 'QSTARS_JV_RAD.bin')
Jp = {k: it.Jdir.matvec(v) + drad[:, :, i] for i, (k, v) in enumerate(pieces.items())}
nrm = lambda v: float(np.linalg.norm(it.obj(v)[:134]) / beta)  # noqa: E731
res = dict(product_of_s_vs_sum_of_products=nrm(Jp['s'] - Jsum), product_of_s=nrm(Jp['s']), sum_of_products=nrm(Jsum),
           pieces={k: nrm(v) for k, v in Jp.items()}, pieces_sum_vs_s=nrm(Jp['pops'] + Jp['ne'] + Jp['T'] - Jp['s']))
print(json.dumps(res, indent=1))
print(f'max|y| {np.abs(y).max():.3e}, ||y * |z_i| || {np.linalg.norm(y * np.array([np.linalg.norm(z) for z in S])):.3e} against ||s|| {np.linalg.norm(s):.3e}')
for lam in lams:
    xt, Rt = read_trace_record(out / f'trial_l{lam:g}' / 'QSTARS_RESID_TRACE.bin')
    sr = np.where(it.active, 1.0 - xt / it.x, 0.0)
    fd = (Rt - it.R) / lam                          # model: R(x(1 - lam s)) = R - lam J s
    corr = it.Jdir.matvec(sr / lam - s)             # realized minus planned, direct blocks
    for lab, Jm in (('fresh product of s', Jp['s']), ('sum of stored products', Jsum)):
        e = it.obj(fd + Jm + corr)[:134]
        print(f'lambda {lam:g}: ||FD + J s|| / merit with the {lab}: {np.linalg.norm(e) / beta:.4f} '
              f'(||FD|| / merit {np.linalg.norm(it.obj(fd)[:134]) / beta:.4f})')
