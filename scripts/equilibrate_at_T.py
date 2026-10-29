#!/usr/bin/env python3
"""Prepare a starting state on a prescribed temperature profile (initial guess only).

T is set to a given profile at selected depths; the populations and n_e at ALL depths are
then Newton-iterated in statistical equilibrium with the stored converged J frozen (replay
passes; the native fixed-J blocks are block-diagonal, so one pass updates every depth). T is
held fixed during this preparation. The result is only a starting point for the T-released
Newton-Krylov solve, which releases T and recomputes J; nothing here is part of a solution.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_ba_blocks, read_trace_record  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('base', type=Path, help='base pass (state + stored J)')
    ap.add_argument('out', type=Path)
    ap.add_argument('--t-floor', type=float, required=True, help='T floor in K applied at --depths (initial guess)')
    ap.add_argument('--depths', type=int, nargs=2, required=True, metavar=('FIRST', 'LAST'), help='native depths, inclusive')
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--x0-trace', type=Path, default=ROOT / 'results/nk_freeT/it00_base/QSTARS_RESID_TRACE.bin')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--maxit', type=int, default=8)
    ap.add_argument('--tol', type=float, default=1e-4)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'),
                    help='VADAT settings applied to every pass (as in the Newton-Krylov run)')
    a = ap.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    base = a.base.resolve()
    x0, _ = read_trace_record(a.x0_trace)
    nt, nd = x0.shape
    ck.NT, ck.ND = nt, nd
    x, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    d0, d1 = a.depths[0] - 1, a.depths[1]
    x = x.copy()
    x[nt - 1, d0:d1] = np.maximum(x[nt - 1, d0:d1], a.t_floor / 1e4)
    R = ck.Runner(a.x0, a.exe, a.replace, False, a.threads, out / 'driver.log', [tuple(v) for v in a.vadat])
    hist = []
    for it in range(a.maxit):
        p = R.run(out / f'i{it}', 1.0 - x / x0, replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False)
        xp, Rp = read_trace_record(p / 'QSTARS_RESID_TRACE.bin')
        _, di, _ = read_ba_blocks(p / 'QSTARS_BA_TRACE.bin', nt, nd)
        smax = 0.0
        for d in range(nd):
            blk = di[d]
            idx = [i for i in range(nt - 1) if np.any(blk[i])]
            s = np.linalg.solve(blk[np.ix_(idx, idx)], Rp[idx, d])
            smax = max(smax, float(np.max(np.abs(s))))
            fs = np.where(1.0 - s > 0.1, 1.0 - s, np.exp(-np.clip(s, -np.log(1e4), np.log(1e4))))
            x[idx, d] = xp[idx, d] * np.clip(fs, 1e-4, 1e4)
        hist.append(dict(it=it, max_local_step=smax))
        R.log(f'equilibrate sweep {it}: max local frozen-J step {smax:.3e}')
        if smax < a.tol:
            break
    np.save(out / 'state.npy', x)
    from qstars.nlte_newton import write_fortran
    write_fortran(out / 'QSTARS_STEP.bin', 1.0 - x / x0)
    (out / 'equilibrate.json').write_text(json.dumps(dict(base=str(base), t_floor=a.t_floor, depths=a.depths,
                                                          history=hist), indent=1))


if __name__ == '__main__':
    main()
