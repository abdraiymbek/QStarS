#!/usr/bin/env python3
"""Measure candidate preconditioners against EXACT Jacobian products at an NK state.

The fixed-J direct Jacobian M (completed as in Report 20 §2.2) omits the radiation
response. In optically thick continua the response nearly cancels the direct
rate derivative (detailed balance with J ~ S), so M overstates those rows by up to
~1e4 and the left-preconditioned Krylov solve stagnates. The native local-ALI
blocks (QS_LIN_FAC > 0) add the diagonal approximate-lambda response. Report 18
showed they are NOT an accurate Jacobian; here they are tested only as a
preconditioner, against the exact products (z, J z) already computed at the state:

    q(M) = || W (M^-1 J z - z) || / || W z ||      (0 = perfect preconditioner)

M_f = M_complete + (BA_f - BA_0): the native ALI increment added to the completed
fixed-J Jacobian. Passes run at the base state with the ALI warm-started from the
base J (converged field reproduced to ~1e-9), BA dumped, stopped after the solve.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import BlockTridiag, read_ba_blocks, read_trace_record  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path)
    ap.add_argument('base', help='accepted base pass name inside run')
    ap.add_argument('k', type=int)
    ap.add_argument('--facs', type=float, nargs='*', default=[0.5, 1.0])
    ap.add_argument('--moment', action='store_true',
                    help='also test the native moment-equation linearization (QS_ALI_LIN=F, Report 18 "native moment")')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--free-t', action='store_true')
    a = ap.parse_args()
    run = a.run.resolve(); base = run / a.base
    itdir = sorted(run.glob(f'it{a.k:02d}_from_{a.base}'))[0]
    x0, _ = read_trace_record(run / 'it00_base' / 'QSTARS_RESID_TRACE.bin')
    ck.NT, ck.ND = x0.shape
    ck.SIG_FRAC, ck.TRACE_WEIGHT = 1e-10, 1e-3
    names = ['ne', 'hii'] + (['T'] if a.free_t else [])
    it = ck.Iterate(base, x0, {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in names}, 1e-4, a.free_t)
    xk, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    s_k = 1.0 - xk / x0
    out = itdir / 'precond_probe'; out.mkdir(exist_ok=True)
    ap_out = out / ('probe_' + '_'.join([f'fac{f:g}' for f in a.facs] + (['moment'] if a.moment else [])) + '.json')
    from concurrent.futures import ThreadPoolExecutor

    variants = [(f'fac{f:g}', [('QS_LIN_FAC', f'{f:g}')]) for f in a.facs]
    if a.moment:
        variants.append(('moment', [('QS_ALI_LIN', 'F')]))

    def one(v):
        name, sets = v
        R = ck.Runner(ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3', a.exe, ROOT / 'results/jv/replace_fac1_x0.txt',
                      not a.free_t, a.threads, out / 'probe.log', sets)
        return name, R.run(out / name, s_k, ba=True, store_j=False, warm=base / 'QSTARS_RJ.bin')
    with ThreadPoolExecutor(len(variants)) as ex:
        passes = list(ex.map(one, variants))
    lo0, di0, up0 = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', ck.NT, ck.ND)
    pairs = []
    for d in sorted(itdir.glob('jv*')):
        if (d / 'pairs.npz').exists():
            f = np.load(d / 'pairs.npz')
            pairs += [(f['z'][:, :, i], f['jz'][:, :, i]) for i in range(f['z'].shape[2])]
    precs = {'fixedJ_complete': it.M}
    for name, d in passes:
        _, R_f = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
        lo, di, up = read_ba_blocks(d / 'QSTARS_BA_TRACE.bin', ck.NT, ck.ND)
        print(f'{name}: residual vs base {np.linalg.norm(R_f - it.R) / np.linalg.norm(it.R):.2e}; '
              f'radiation increment ||dBA||/||BA0|| diag {np.linalg.norm(di - di0) / np.linalg.norm(di0):.2e} '
              f'off-diag {np.linalg.norm(lo) + np.linalg.norm(up):.2e}')
        Mf = BlockTridiag(it.Jdir.lo + (lo - lo0), it.Jdir.di + (di - di0), it.Jdir.up + (up - up0)).factor()
        precs[name] = Mf
    res = {}
    W = it.W * it.active
    for name, M in precs.items():
        q = [float(np.linalg.norm(W * (M.solve(jz) - z)) / np.linalg.norm(W * z)) for z, jz in pairs]
        b = (W * M.solve(it.R)).ravel()
        A = np.stack([(W * M.solve(jz)).ravel() for _, jz in pairs], 1)
        y, *_ = np.linalg.lstsq(A, b, rcond=1e-12)
        res[name] = dict(q=q, lsq_rel_residual_same_directions=float(np.linalg.norm(b - A @ y) / np.linalg.norm(b)),
                         merit=float(np.linalg.norm(b)))
        print(f'{name:18s} q per product {[round(v, 3) for v in q]}  lsq residual (same {len(pairs)} directions) '
              f'{res[name]["lsq_rel_residual_same_directions"]:.3f}  ||W M^-1 R||={res[name]["merit"]:.3e}')
    ap_out.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
