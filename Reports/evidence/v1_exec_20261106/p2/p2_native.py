#!/usr/bin/env python3
"""V1 plan P2: native trials of the bounded step built from the v1_p2 run's products (all passes or the first N), with the
driver's --struct-active measure, through the driver's trial pipeline (p1_native.realize). The step: the better of the
LM-QR and box solutions (cool_newton_krylov.bounded_step), or the one named. Passes in RUN_DIR/offline_trials/<tag>.

    p2_native.py RUN_DIR ITERATE_NAME BASE_PATH [--passes N] [--solver best|qr|box] [--lams 1 0.5 0.25] [--dry-run]"""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path); ap.add_argument('name'); ap.add_argument('base', type=Path)
    ap.add_argument('--passes', type=int, default=0); ap.add_argument('--solver', default='best')
    ap.add_argument('--lams', type=float, nargs='+', default=[1.0, 0.5, 0.25]); ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--workers', type=int, default=3); ap.add_argument('--threads', type=int, default=3)
    a = ap.parse_args()
    x0b = cn.configure()
    ck.STRUCT_ACTIVE = True
    itdir = a.run / a.name
    cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(cn.residual_dir(a.base), x0b, cols, 1e-4, True, None, rj=a.base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    b = it.obj(it.R).ravel(); beta = float(np.linalg.norm(b))
    rad = np.full(it.x.shape, 0.5); rad[133] = 0.03
    passes = sorted({p.name[:4] for p in itdir.glob('jv[0-9][0-9]p[0-9]') if (p / 'pairs.npz').exists()})
    if a.passes:
        passes = passes[:a.passes]
    S, JS = [], []
    for pn in passes:
        for q in sorted(itdir.glob(f'{pn}p[0-9]')):
            z = np.load(q / 'pairs.npz')
            S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
    A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
    cands = {}
    for sv in (['qr', 'box'] if a.solver == 'best' else [a.solver]):
        s_, mu_, y_ = ck.bounded_step((S, A, b), it.active, rad, it.active, sv)
        cands[sv] = (s_, y_, float(np.linalg.norm(A @ y_ - b) / beta))
    sv = min(cands, key=lambda k: cands[k][2])
    s, y, fit = cands[sv]
    print(f'{len(passes)} passes, {len(S)} directions: ' + ', '.join(f'{k} {v[2]:.4f}' for k, v in cands.items())
          + f' -> {sv}', flush=True)
    x0, _ = read_trace_record(ROOT / 'results/nk167_trap/free_r25/it00_base/QSTARS_RESID_TRACE.bin')
    cons_tot = cnk_trap.conservation_totals(x0)
    w = it.row_weights() * it.D * it.active
    Js = (A @ y).reshape(it.x.shape)
    rows, jobs = [], []
    out = a.run / 'offline_trials'
    for lam in a.lams:
        xt, pst, _ = realize(it, lam * s, x0, cons_tot)
        d = np.where(it.active, (1.0 - xt / it.x) - lam * s, 0.0)
        r_pl = b.reshape(it.x.shape) - lam * Js
        r_re = r_pl - it.obj(it.Jdir.matvec(d))
        tag = f'p{len(passes)}_{sv}_l{lam:g}'
        rows.append(dict(tag=tag, lam=lam, solver=sv, passes=len(passes), planned=float(np.linalg.norm(r_pl[:134]) / beta),
                         realized=float(np.linalg.norm(r_re[:134]) / beta), projection=pst))
        jobs.append((tag, out / tag, 1.0 - xt / x0))
        print(f'{tag}: planned {rows[-1]["planned"]:.4f}, realized {rows[-1]["realized"]:.4f}, projection '
              + ', '.join(f'{k} {v:.1e}' for k, v in pst.items()), flush=True)
    if a.dry_run:
        return 0
    R = ck.Runner(X0DIR, EXE, REPLACE, False, a.threads, out / 'offline_trials.log',
                  [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
    R.no_obs = True
    out.mkdir(parents=True, exist_ok=True)
    tpj = a.base / 'QSTARS_RJ.bin'
    z = np.load(FROZEN)

    def run(j):
        tag, dd, st = j
        t0 = time.time()
        try:
            R.run(dd, st, threads=a.threads, warm=tpj, two_phot_j=tpj)
        except RuntimeError as e:
            return None, str(e), time.time() - t0
        return read_trace_record(dd / 'QSTARS_RESID_TRACE.bin')[1], None, time.time() - t0
    with ThreadPoolExecutor(a.workers) as ex:
        res = list(ex.map(run, jobs))
    for row, (Rt, err, sec) in zip(rows, res):
        row['seconds'] = sec
        if Rt is None:
            row['error'] = err; print(f'{row["tag"]}: native pass failed ({err})'); continue
        row['native'] = float(np.linalg.norm((w * Rt)[:134]) / beta)
        row['native_frozen'] = float(np.linalg.norm(z['w'] * Rt[:134]))
        print(f'{row["tag"]}: planned {row["planned"]:.4f} realized {row["realized"]:.4f} NATIVE {row["native"]:.4f} '
              f'(frozen {row["native_frozen"]:.6e}) [{sec:.0f} s]', flush=True)
    f = Path(__file__).with_name(f'native_{a.run.name}.json')
    prev = json.loads(f.read_text()) if f.exists() else dict(run=str(a.run), base_merit=beta, rows=[])
    prev['rows'] += rows
    f.write_text(json.dumps(prev, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
