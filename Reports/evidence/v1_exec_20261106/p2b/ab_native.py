#!/usr/bin/env python3
"""Council block, matched A/B: a bounded step rebuilt from the SAME stored products of one v1_p2 iterate (span basis
tau, box/QR best, |s| <= BOUND, T BOUND_T), realized by the driver's trial pipeline (p1_native.realize: linear map with
the 0.693 limit, He III trace balance, conservation projection) and run natively as a driver trial pass (warm ALI and
two-photon J from the base's stored J). Prints the planned linear fit, the native merit in the iterate's local
acceptance measure (relative to the base) and the frozen it16 measure. Passes: ITERATE_DIR/../ab_<tag>.

    ab_native.py ITERATE_DIR BASE_DIR TAG --bound 0.5 --bound-t 0.03 --lams 0.5 [--tau 1e-6]"""
import argparse, json, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'p1'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
from p1_native import realize, EXE, REPLACE, X0DIR, FROZEN  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument('itdir', type=Path); ap.add_argument('base', type=Path); ap.add_argument('tag')
ap.add_argument('--bound', type=float, required=True); ap.add_argument('--bound-t', type=float, required=True)
ap.add_argument('--lams', type=float, nargs='+', required=True); ap.add_argument('--tau', type=float, default=1e-6)
ap.add_argument('--threads', type=int, default=3)
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
x0, _ = read_trace_record(ROOT / 'results/nk167_trap/free_r25/it00_base/QSTARS_RESID_TRACE.bin')
ct = cnk_trap.conservation_totals(x0)
R = ck.Runner(X0DIR, EXE, REPLACE, False, a.threads, a.itdir.parent / 'ab_native.log',
              [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
R.no_obs = True
tpj = a.base / 'QSTARS_RJ.bin'
z = np.load(FROZEN)
w = it.row_weights() * it.D * it.active
jobs = []
for lam in a.lams:
    xt, pst, _ = realize(it, lam * s, x0, ct)
    jobs.append((lam, a.itdir.parent / f'ab_{a.tag}_l{lam:g}', 1.0 - xt / x0,
                 float(np.linalg.norm(lam * (Ab @ c) - b) / beta), pst))


def run(j):
    t0 = time.time(); R.run(j[1], j[2], threads=a.threads, warm=tpj, two_phot_j=tpj)
    return read_trace_record(j[1] / 'QSTARS_RESID_TRACE.bin')[1], time.time() - t0


with ThreadPoolExecutor(len(jobs)) as ex:
    res = list(ex.map(run, jobs))
rows = []
for j, (Rt, sec) in zip(jobs, res):
    r = dict(tag=a.tag, bound=a.bound, bound_t=a.bound_t, lam=j[0], basis=kb, predicted=j[3],
             native=float(np.linalg.norm((w * Rt)[:134]) / beta), native_frozen=float(np.linalg.norm(z['w'] * Rt[:134])),
             projection=j[4], seconds=sec)
    rows.append(r)
    print(f'{a.tag} lambda {j[0]:g}: predicted {r["predicted"]:.4f}, NATIVE {r["native"]:.4f} (frozen {r["native_frozen"]:.6e}) '
          f'[{sec:.0f} s]', flush=True)
f = Path(__file__).with_name('ab_native.json')
prev = json.loads(f.read_text()) if f.exists() else []
f.write_text(json.dumps(prev + rows, indent=1))
