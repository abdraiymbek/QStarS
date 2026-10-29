#!/usr/bin/env python3
"""Multi-preconditioned Krylov probe at an NK state.

Directions W P^-1 R for the alternative native linearizations P (moment, local ALI 0.5
and 1; each = completed fixed-J Jacobian + its native radiation increment, from the
precond_probe passes) get EXACT products in one batched tangent pass. The least squares
||W M^-1 (R - J S y)|| (M = completed fixed-J, the driver's objective) is then evaluated
with and without them, next to the directions the driver already used.
Only exact products enter; a poor P merely contributes a useless direction.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import BlockTridiag, read_ba_blocks, read_jv_rad, read_trace_record  # noqa: E402


def precs_at(it, itdir, names):
    lo0, di0, up0 = read_ba_blocks(it.dir / 'QSTARS_BA_TRACE.bin', ck.NT, ck.ND)
    out = {}
    for n in names:
        d = itdir / 'precond_probe' / n
        lo, di, up = read_ba_blocks(d / 'QSTARS_BA_TRACE.bin', ck.NT, ck.ND)
        out[n] = BlockTridiag(it.Jdir.lo + (lo - lo0), it.Jdir.di + (di - di0), it.Jdir.up + (up - up0)).factor()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path)
    ap.add_argument('base')
    ap.add_argument('k', type=int)
    ap.add_argument('--precs', nargs='+', default=['moment', 'fac1', 'fac0.5'])
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--threads', type=int, default=4)
    a = ap.parse_args()
    run = a.run.resolve(); base = run / a.base
    itdir = sorted(run.glob(f'it{a.k:02d}_from_{a.base}'))[0]
    x0, _ = read_trace_record(run / 'it00_base' / 'QSTARS_RESID_TRACE.bin')
    ck.NT, ck.ND = x0.shape
    ck.SIG_FRAC, ck.TRACE_WEIGHT = 1e-10, 1e-3
    it = ck.Iterate(base, x0, {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii')}, 1e-4, False)
    P = precs_at(it, itdir, a.precs)
    dirs = [it.W * P[n].solve(it.R) for n in a.precs]
    (itdir / 'mpk_probe').mkdir(exist_ok=True)
    R = ck.Runner(ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3', a.exe, ROOT / 'results/jv/replace_fac1_x0.txt',
                  True, a.threads, itdir / 'mpk_probe.log')
    d = R.run(itdir / 'mpk_probe' / 'jv', it.s_cum, replay=it.rj, dirs=np.stack(dirs, 2), ba=False, store_j=False)
    _, drad = read_jv_rad(d / 'QSTARS_JV_RAD.bin')
    new = [(z, it.Jdir.matvec(z) + drad[:, :, i]) for i, z in enumerate(dirs)]
    np.savez(d / 'pairs.npz', z=np.stack([p[0] for p in new], 2), jz=np.stack([p[1] for p in new], 2))
    old = []
    for j in sorted(itdir.glob('jv*')):
        f = np.load(j / 'pairs.npz')
        old += [(f['z'][:, :, i], f['jz'][:, :, i]) for i in range(f['z'].shape[2])]
    wm = lambda r: it.W * it.M.solve(r)
    b = wm(it.R).ravel(); nb = np.linalg.norm(b)

    def fit(pairs):
        A = np.stack([wm(jz).ravel() for _, jz in pairs], 1)
        y, *_ = np.linalg.lstsq(A, b, rcond=1e-12)
        return float(np.linalg.norm(b - A @ y) / nb)
    res = {'driver_directions_all': fit(old)}
    first = old[:20]
    res['driver_first_pass'] = fit(first)
    for i, n in enumerate(a.precs):
        res[f'{n}_alone'] = fit([new[i]])
        res[f'first_pass_plus_{n}'] = fit(first + [new[i]])
    res['first_pass_plus_all_new'] = fit(first + new)
    res['all_plus_all_new'] = fit(old + new)
    for k, v in res.items():
        print(f'{k:32s} {v:.3f}')
    (itdir / 'mpk_probe' / 'result.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
