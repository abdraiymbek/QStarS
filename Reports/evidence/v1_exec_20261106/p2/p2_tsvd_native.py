#!/usr/bin/env python3
"""V1 plan P2: native trials (driver trial pipeline) of one candidate step saved by p2_tsvd.py.
    p2_tsvd_native.py STEPS.npy INDEX TAG [--lams 1 0.5 0.25]"""
import sys, json, time, argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p2_tsvd import load, RUN, BASE  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'p1'))
from p1_native import realize, EXE, REPLACE, X0DIR, FROZEN  # noqa: E402
import cool_newton_krylov as ck, cnk_trap  # noqa: E402,E401
from qstars.nlte_newton import read_trace_record  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument('steps'); ap.add_argument('index', type=int); ap.add_argument('tag')
ap.add_argument('--lams', type=float, nargs='+', default=[1.0, 0.5, 0.25]); a = ap.parse_args()
it, *_ = load(1)
s = np.load(a.steps)[:, :, a.index]
b = it.obj(it.R); beta = float(np.linalg.norm(b[:134])); w = it.row_weights() * it.D * it.active
x0, _ = read_trace_record(Path(X0DIR).parents[0] / 'nk167_trap/free_r25/it00_base/QSTARS_RESID_TRACE.bin')
ct = cnk_trap.conservation_totals(x0)
R = ck.Runner(X0DIR, EXE, REPLACE, False, 3, RUN / 'offline_trials.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
R.no_obs = True
tpj = BASE / 'QSTARS_RJ.bin'; z = np.load(FROZEN)
jobs = []
for lam in a.lams:
    xt, pst, _ = realize(it, lam * s, x0, ct)
    jobs.append((lam, RUN / 'offline_trials' / f'{a.tag}_l{lam:g}', 1.0 - xt / x0, pst))
def run(j):
    t0 = time.time(); R.run(j[1], j[2], threads=3, warm=tpj, two_phot_j=tpj)
    return read_trace_record(j[1] / 'QSTARS_RESID_TRACE.bin')[1], time.time() - t0
with ThreadPoolExecutor(3) as ex:
    res = list(ex.map(run, jobs))
rows = []
for j, (Rt, sec) in zip(jobs, res):
    r = dict(tag=a.tag, lam=j[0], native=float(np.linalg.norm((w * Rt)[:134]) / beta), native_frozen=float(np.linalg.norm(z['w'] * Rt[:134])),
             projection=j[3], seconds=sec)
    rows.append(r); print(f'{a.tag} lambda {j[0]:g}: NATIVE {r["native"]:.4f} (frozen {r["native_frozen"]:.6e}), projection '
                          + ', '.join(f'{k} {v:.1e}' for k, v in j[3].items()) + f' [{sec:.0f} s]', flush=True)
f = Path(__file__).with_name('tsvd_native.json')
prev = json.loads(f.read_text()) if f.exists() else []
f.write_text(json.dumps(prev + rows, indent=1))
