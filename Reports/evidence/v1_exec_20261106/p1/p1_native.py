#!/usr/bin/env python3
"""V1 plan P1, item 3: planned step, realized update and native response of candidate steps at one archived base.

For each candidate y (a .npy coefficient vector over the iterate's span, or LM_NE / LM_QR from steps_<it>.npz) and each
lambda: the trial state is made exactly as the driver makes it (cool_newton_krylov.py trial(): linear update map with
the 0.693 component limit, He III trace balance against the realized update, conservation projection onto x0's
totals), and the native pass is run as a driver trial pass (warm ALI from the base's J, the base's stored J as the
two-photon J). Reported, in the base's local acceptance measure (row-scaled, its weights; merit_rel = /base merit):
  planned    ||w (R - lam J s)||                     linear model of the planned step (exact products)
  realized   ||w (R - lam J s - Jdir d)||            d = realized - planned step (its action by the direct blocks)
  native     ||w R(x_trial)||                        the trial pass
and the frozen reference measure of the native residual. Passes go to results/nk167_trap/v1_p1/<base>/<tag>.

    p1_native.py ITERATE_DIR TAG=SOURCE [TAG=SOURCE ...] [--lams 1 0.5] [--workers 4] [--dry-run]"""
import argparse, json, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached, ck, cnk_trap, cn, read_trace_record, ROOT  # noqa: E402

RUN = ROOT / 'results/nk167_trap/free_r25'
EXE = ROOT / 'results/nk167_trap/exe_b36/cmfgen_dev.exe'
REPLACE = ROOT / 'results/nk167_adv/transfer/replace_167.txt'
X0DIR = ROOT / 'results/cmfgen_ali_lambda167'
FROZEN = ROOT / 'Reports/evidence/v1_exec_20261106/frozen_norm_it16.npz'
CAP = 0.693


def realize(it, s, x0, cons_tot):
    """The driver's trial state for the planned step s (trial(), --update-map linear, --trace-ion-balance, --conserve)."""
    xt = it.x * cnk_trap.update_factor(s, 'linear', CAP)
    s_real = 1.0 - xt / it.x
    s_real[131] = 0.0
    Jr = it.Jdir.matvec(s_real)
    mask_, sc_ = cnk_trap.trace_ion_balance(it.R, Jr, it.Jdir.di, it.frac, 1.0)
    f_ = cnk_trap.update_factor(sc_, 'linear', CAP)
    xt[131, mask_] = it.x[131, mask_] * f_[mask_]
    xt, pst = cnk_trap.project_conservation(xt, cons_tot)
    return xt, pst, int(mask_.sum())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('itdir'); ap.add_argument('cands', nargs='+')
    ap.add_argument('--lams', type=float, nargs='+', default=[1.0, 0.5]); ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--threads', type=int, default=2); ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    P = cached(a.itdir)
    it = P.it
    x0, _ = read_trace_record(RUN / 'it00_base' / 'QSTARS_RESID_TRACE.bin')
    cons_tot = cnk_trap.conservation_totals(x0)
    w = P.w
    nr = 134
    phi = float(np.linalg.norm((w * it.R)[:nr]))
    here = Path(__file__).resolve().parent
    zst = np.load(here / f'steps_{a.itdir}.npz')
    jobs, rows = [], []
    outdir = ROOT / 'results/nk167_trap/v1_p1' / P.base.name
    Q = None
    for c in a.cands:
        tag, src = c.split('=', 1)
        y = zst[src] if src in zst.files else np.load(here / src)
        PP = P
        if 'BOXP' in src:                                # coefficients over the projected span (box_proj.py)
            if Q is None:
                from box_proj import projected_problem
                Q = projected_problem(P)
            PP = Q
        s = PP.step(y)
        Js = (PP.A @ y)                                  # obj(J s), flattened (rowscale: w * J s)
        for lam in a.lams:
            xt, pst, ntr = realize(it, lam * s, x0, cons_tot)
            s_real = 1.0 - xt / it.x
            d = np.where(it.active, s_real - lam * s, 0.0)
            r_pl = (P.b - lam * Js).reshape(P.shape)
            r_re = r_pl - it.obj(it.Jdir.matvec(d))
            lab = f'{tag}_l{lam:g}'
            row = dict(tag=lab, lam=lam, fit_planned=float(np.linalg.norm(r_pl[:nr]) / phi),
                       fit_realized=float(np.linalg.norm(r_re[:nr]) / phi), max_d=float(np.max(np.abs(d))),
                       max_s=float(np.max(np.abs(lam * s[PP.mask]))), projection=pst, trace_ion_depths=ntr)
            rows.append(row)
            jobs.append((lab, outdir / lab, 1.0 - xt / x0))
            print(f'{lab}: planned {row["fit_planned"]:.4f}, realized {row["fit_realized"]:.4f}, max|d| {row["max_d"]:.3f}, '
                  f'projection {", ".join(f"{k} {v:.1e}" for k, v in pst.items())}', flush=True)
    if a.dry_run:
        return 0
    R = ck.Runner(X0DIR, EXE, REPLACE, False, a.threads, outdir / 'p1_native.log',
                  [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
    R.no_obs = True
    outdir.mkdir(parents=True, exist_ok=True)
    tpj = P.base / 'QSTARS_RJ.bin'
    z = np.load(FROZEN)

    def run(j):
        lab, d, st = j
        t0 = time.time()
        try:
            R.run(d, st, threads=a.threads, warm=tpj, two_phot_j=tpj)
        except RuntimeError as e:
            return lab, None, str(e), time.time() - t0
        xr, Rt = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
        return lab, Rt, None, time.time() - t0
    with ThreadPoolExecutor(a.workers) as ex:
        res = list(ex.map(run, jobs))
    for row, (lab, Rt, err, sec) in zip(rows, res):
        row['seconds'] = sec
        if Rt is None:
            row['native'] = None; row['error'] = err
            print(f'{lab}: native pass failed ({err})')
            continue
        row['native'] = float(np.linalg.norm((w * Rt)[:nr]) / phi)
        row['native_frozen'] = float(np.linalg.norm(z['w'] * Rt[:nr]))
        print(f'{lab}: planned {row["fit_planned"]:.4f}  realized {row["fit_realized"]:.4f}  NATIVE {row["native"]:.4f} '
              f'(frozen {row["native_frozen"]:.6e})  [{sec:.0f} s]', flush=True)
    out = a.out or here / f'native_{a.itdir}.json'
    prev = json.loads(out.read_text()) if out.exists() else dict(iterate=a.itdir, base_merit=phi, rows=[])
    prev['rows'] += rows
    out.write_text(json.dumps(prev, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
