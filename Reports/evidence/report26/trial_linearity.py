#!/usr/bin/env python3
"""Report 26: why the trials of a free_r25 iterate miss the linear prediction. Rebuilds the iterate's least squares
(stored exact products, row-scaled objective) and its |s| <= 0.5 (T 0.03) bounded step, then for every trial state
compares, in the acceptance measure:
  actual    ||w R(x_trial)||                       (the trial pass)
  planned   ||w (R - lam J s_b)||                  (exact products: the linear model of the planned step)
  realized  ||w (R - lam J s_b - Jdir d)||         d = realized - planned fractional step (projection, He III balance,
                                                   limiter), its action by the direct blocks only
and the same by depth region. Read-only.   trial_linearity.py ITERATE_DIR"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import span_check as sc  # noqa: E402  (same settings)
ck = sc.ck
REG = sc.REG


def main():
    itd = sc.O / sys.argv[1]
    base = itd.name.split('_from_')[1]
    k = int(itd.name[2:4])
    cols = {n: (itd / f'col_{n}_p', itd / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    resdir = sc.O / (base + '_tp') if (sc.O / (base + '_tp')).exists() else sc.O / base
    it = ck.Iterate(resdir, sc.x0, cols, 1e-4, True, None, rj=sc.O / base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    w = (it.row_weights() * it.D * it.active)[:134]
    bo = it.obj(it.R)
    S, JS = sc.pairs(sorted(itd.glob('jv*')))
    A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
    rad = np.full(it.x.shape, 0.5); rad[133] = 0.03
    sb, mu, y = ck.lm_step_bound((S, A, bo.ravel()), it.active.astype(float), rad, it.active.copy())
    Jsb = sum(yi * ji for yi, ji in zip(y, JS))
    m0 = np.linalg.norm(w * it.R[:134])
    print(f'{itd.name}: base merit {m0:.5f}; bounded step mu {mu:.2e}, max|s| {np.abs(sb).max():.3f}, '
          f'linear fit {np.linalg.norm(w * (it.R - Jsb)[:134]) / m0:.4f}')
    for lab in ('l1', 'l0.5', 'l0.25', 'l0.125'):
        tdir = sc.O / f'it{k + 1:02d}_trial_{lab}'
        if not (tdir / 'QSTARS_RESID_TRACE.bin').exists():
            continue
        lam = float(lab[1:])
        xt, Rt = sc.read_trace_record(tdir / 'QSTARS_RESID_TRACE.bin')
        s_real = 1.0 - xt / it.x
        dlt = np.where(it.active, s_real - lam * sb, 0.0)
        r_pl = it.R - lam * Jsb
        r_re = r_pl - it.Jdir.matvec(dlt)
        ra, rp, rr = (w * v[:134] for v in (Rt, r_pl, r_re))
        print(f' {lab}: actual {np.linalg.norm(ra) / m0:.4f}  planned {np.linalg.norm(rp) / m0:.4f}  '
              f'realized(direct) {np.linalg.norm(rr) / m0:.4f}  max|d| {np.abs(dlt).max():.3f} '
              f'(rows {sorted(set((np.argsort(np.abs(dlt).ravel())[::-1][:5] // ck.ND + 1).tolist()))})')
        print('    by region, actual / planned / realized (share of the base merit^2): ' + ', '.join(
            f'd{a}-{b} {np.sum(ra[:, a - 1:b]**2) / m0**2:.3f}/{np.sum(rp[:, a - 1:b]**2) / m0**2:.3f}/'
            f'{np.sum(rr[:, a - 1:b]**2) / m0**2:.3f}' for a, b in REG))


if __name__ == '__main__':
    main()
