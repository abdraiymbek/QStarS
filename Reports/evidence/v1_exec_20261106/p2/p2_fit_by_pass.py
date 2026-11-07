#!/usr/bin/env python3
"""V1 plan P2: fit of the least squares by the number of sequential global product passes at one base (the v1_p2 run),
unbounded, LM-QR bounded and box bounded (|s| <= 0.5, T 0.03), from the run's stored exact products (pairs.npz of
jvNNpK), with the run's --struct-active measure. Read-only.

    p2_fit_by_pass.py RUN_DIR ITERATE_NAME BASE_PATH     (e.g. results/nk167_trap/v1_p2 it16_from_it16_trial_l0.5
                                                          results/nk167_trap/free_r25/it16_trial_l0.5)"""
import json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402


def main():
    run, name, base = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    x0 = cn.configure()
    ck.STRUCT_ACTIVE = True
    itdir = run / name
    cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(cn.residual_dir(base), x0, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    b = it.obj(it.R).ravel(); beta = np.linalg.norm(b)
    rad = np.full(it.x.shape, 0.5); rad[133] = 0.03
    passes = sorted({p.name[:4] for p in itdir.glob('jv[0-9][0-9]p[0-9]') if (p / 'pairs.npz').exists()})
    S, JS, rows = [], [], []
    for k, pn in enumerate(passes):
        for q in sorted(itdir.glob(f'{pn}p[0-9]')):
            z = np.load(q / 'pairs.npz')
            S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]
            JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
        A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
        t0 = time.time()
        y0 = np.linalg.lstsq(A, b, rcond=1e-12)[0]
        _, mu, yq = ck.bounded_step((S, A, b), it.active, rad, it.active, 'qr')
        _, _, yb = ck.bounded_step((S, A, b), it.active, rad, it.active, 'box')
        r = dict(passes=k + 1, directions=len(S), unbounded=float(np.linalg.norm(A @ y0 - b) / beta),
                 lm_qr=float(np.linalg.norm(A @ yq - b) / beta), box=float(np.linalg.norm(A @ yb - b) / beta),
                 seconds=time.time() - t0)
        rows.append(r)
        print(f'{k + 1:2d} passes, {len(S):4d} directions: unbounded {r["unbounded"]:.4f}, LM-QR {r["lm_qr"]:.4f}, '
              f'box {r["box"]:.4f}', flush=True)
    out = Path(__file__).with_name(f'fit_by_pass_{run.name}_{name}.json')
    out.write_text(json.dumps(dict(run=str(run), iterate=name, base=str(base), base_merit=float(beta), rows=rows), indent=1))


if __name__ == '__main__':
    sys.exit(main())
