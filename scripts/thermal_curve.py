#!/usr/bin/env python3
"""Local thermal-equilibrium curve at one depth of a cool-closure state (diagnostic).

Question: does the radiative-equilibrium (RE) equation at depth d have a root in T when the
local populations are in statistical equilibrium with the CURRENT radiation field?
For each trial temperature T = f * T_d (all other depths and the stored converged J fixed),
the populations and n_e at depth d are Newton-solved with frozen-J replay passes (the native
fixed-J block at d is the exact local Jacobian for frozen J apart from the level-dissolution
terms, which only slow the local iteration). The RE residual is then reported relative to the
RE row's gross scale (row l1-norm of the block): net heating minus cooling over gross.
Frozen J is the optically thin limit of the local response; it is a diagnostic of the local
thermal balance, not a solution of the coupled problem.
"""
import argparse, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_ba_blocks, read_trace_record  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('base', type=Path, help='base pass (state + stored J)')
    ap.add_argument('depth', type=int, help='native depth index (1 = outermost)')
    ap.add_argument('out', type=Path)
    ap.add_argument('--factors', type=float, nargs='+', default=[0.5, 0.65, 0.8, 1.0, 1.25, 1.6, 2.0, 2.6, 3.2])
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--x0-trace', type=Path, default=ROOT / 'results/nk_freeT/it00_base/QSTARS_RESID_TRACE.bin')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--maxit', type=int, default=8)
    ap.add_argument('--tol', type=float, default=1e-6)
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'),
                    help='extra VADAT settings for the diagnostic passes (e.g. INC_AD T)')
    a = ap.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    base = a.base.resolve()
    x0, _ = read_trace_record(a.x0_trace)
    nt, nd = x0.shape
    ck.NT, ck.ND = nt, nd
    d = a.depth - 1
    xb, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    R = ck.Runner(a.x0, a.exe, a.replace, False, 1, out / 'driver.log', [tuple(v) for v in a.vadat])

    def point(f):
        x = xb.copy()
        x[nt - 1, d] *= f
        hist = []
        for it in range(a.maxit):
            p = R.run(out / f'f{f:g}_i{it}', 1.0 - x / x0, replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False, threads=1)
            xp, Rp = read_trace_record(p / 'QSTARS_RESID_TRACE.bin')
            _, di, _ = read_ba_blocks(p / 'QSTARS_BA_TRACE.bin', nt, nd)
            blk = di[d]
            idx = [i for i in range(nt - 1) if np.any(blk[i])]
            s = np.linalg.solve(blk[np.ix_(idx, idx)], Rp[idx, d])
            re_scale = float(np.abs(blk[nt - 1]).sum())
            rec = dict(it=it, T=float(xp[nt - 1, d] * 1e4), max_local_step=float(np.max(np.abs(s))),
                       RE_residual=float(Rp[nt - 1, d]), RE_scale=re_scale,
                       RE_relative=float(Rp[nt - 1, d] / re_scale) if re_scale > 0 else None,
                       ne=float(xp[nt - 2, d]), neutral_H=float(xp[1:31, d].sum() / xp[1:32, d].sum()))
            hist.append(rec)
            if rec['max_local_step'] < a.tol:
                break
            # linear update where it is safe, capped to a factor 1e4 per local iteration either way
            fs = np.where(1.0 - s > 0.1, 1.0 - s, np.exp(-np.clip(s, -np.log(1e4), np.log(1e4))))
            x[idx, d] = x[idx, d] * np.clip(fs, 1e-4, 1e4)
        return f, hist
    def safe(f):
        try:
            return point(f)
        except RuntimeError as e:          # record, do not hide: a rejected/failed local solve
            return f, [dict(it=-1, failed=str(e)[:300], T=float('nan'), max_local_step=float('inf'),
                            RE_relative=float('nan'), ne=float('nan'), neutral_H=float('nan'))]
    with ThreadPoolExecutor(len(a.factors)) as ex:
        res = list(ex.map(safe, a.factors))
    summary = []
    for f, hist in sorted(res):
        last = hist[-1]
        summary.append(dict(factor=f, **last, iterations=len(hist), converged=last['max_local_step'] < a.tol))
        print(f"f={f:5.2f} T={last['T']:8.1f} K  RE/gross={last['RE_relative']:+.3e}  ne={last['ne']:.3e}  "
              f"neutral H={last['neutral_H']:.4f}  local step {last['max_local_step']:.1e} ({len(hist)} its)")
    (out / 'thermal_curve.json').write_text(json.dumps(dict(base=str(base), depth=a.depth, vadat=a.vadat, points=summary,
                                                            histories={str(f): h for f, h in res}), indent=1))


if __name__ == '__main__':
    main()
