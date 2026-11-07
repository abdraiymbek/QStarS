#!/usr/bin/env python3
"""Report 26: the --geodesic computation of the driver, offline on one free_r25 iterate (stored products and trials):
second-order term Q from the trials with lambda <= 0.5, how well the iterate's subspace represents J t = Q, and the
predicted residual of the curved path lam s + lam^2 t if the third-order remainder of the fit is kept. With
--dense DIR (an iterate directory with dense* passes at the SAME base) the dense columns join the subspace. Read-only.

    geodesic_check.py ITERATE_DIR"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import span_check as sc  # noqa: E402
import cnk_trap  # noqa: E402
ck = sc.ck


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
    use, Es, Rts = [], [], {}
    for lab in ('l1', 'l0.5', 'l0.25', 'l0.125'):
        tdir = sc.O / f'it{k + 1:02d}_trial_{lab}'
        lam = float(lab[1:])
        xt, Rt = sc.read_trace_record(tdir / 'QSTARS_RESID_TRACE.bin')
        Rts[lam] = Rt
        if lam <= 0.5:
            dl = np.where(it.active, (1.0 - xt / it.x) - lam * sb, 0.0)
            Es.append(Rt - (it.R - lam * Jsb - it.Jdir.matvec(dl))); use.append(lam)
    Q = cnk_trap.quadratic_term(use, Es)
    # cubic coefficient of the same fit, for the prediction
    lams = np.asarray(use); e = [Ei / l ** 2 for Ei, l in zip(Es, lams)]
    cC = np.linalg.pinv(np.stack([np.ones_like(lams), lams], 1))[1]
    C = sum(c * ei for c, ei in zip(cC, e))
    q = it.obj(Q).ravel()
    yt = np.linalg.lstsq(A, q, rcond=1e-12)[0]
    Jt = sum(yi * ji for yi, ji in zip(yt, JS))
    t = sum(yi * si for yi, si in zip(yt, S))
    wv = it.W * it.active
    print(f'{itd.name}: |obj(Q)|/merit {np.linalg.norm(q) / m0:.3f}, fit of J t = Q {np.linalg.norm(q - A @ yt) / np.linalg.norm(q):.3f}, '
          f'||W t||/||W s|| {np.linalg.norm(wv * t) / np.linalg.norm(wv * sb):.3f}, max|t| {np.abs(t).max():.3f}')
    for lam in (1.0, 0.5, 0.25):
        quad = it.R - lam * Jsb + lam ** 2 * Q + lam ** 3 * C          # the fitted model of the plain trial
        geo = it.R - lam * Jsb + lam ** 2 * (Q - Jt) + lam ** 3 * C     # same, with the correction (cross terms left out)
        print(f'  lam {lam}: actual {np.linalg.norm(w * Rts[lam][:134]) / m0:.4f}, fitted model {np.linalg.norm(w * quad[:134]) / m0:.4f}, '
              f'geodesic estimate {np.linalg.norm(w * geo[:134]) / m0:.4f}, linear {np.linalg.norm(w * (it.R - lam * Jsb)[:134]) / m0:.4f}')


if __name__ == '__main__':
    main()
