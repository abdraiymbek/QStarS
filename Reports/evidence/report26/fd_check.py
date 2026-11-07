#!/usr/bin/env python3
"""Report 26: finite-difference check of one Jacobian column at a free_r25 base (one variable at one depth), against
the driver's direct blocks (Jdir) and the exact product (Jdir + native radiation part). Three passes, in
results/nk167_trap/free_r25/fdcheck_<base>/<tag>:
  frozen  : frozen-J replay at x (1 - eps e)           -> the direct column the driver's Jdir should hold
  full    : full pass (converged J) at x (1 - eps e)   -> the true column (what a trial sees)
  product : product pass of direction e at x           -> the native radiation part of the column
Prints, for the equations most changed, the frozen FD, Jdir, the full FD and Jdir + radiation part, in units of each
equation's gross scale (driver row scale D).

    fd_check.py ITERATE_DIR ROW DEPTH [EPS]     (1-based ROW and DEPTH; ITERATE_DIR under free_r25, its base the state)"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import span_check as sc  # noqa: E402
ck = sc.ck
ROOT = sc.ROOT


def main():
    itd = sc.O / sys.argv[1]; row, dep = int(sys.argv[2]), int(sys.argv[3])
    eps = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-3
    base = sc.O / itd.name.split('_from_')[1]
    cols = {n: (itd / f'col_{n}_p', itd / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(base.with_name(base.name + '_tp'), sc.x0, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
    R = ck.Runner(ROOT / 'results/cmfgen_ali_lambda167', ROOT / 'results/nk167_trap/exe_b36/cmfgen_dev.exe',
                  ROOT / 'results/nk167_adv/transfer/replace_167.txt', False, 2, sc.O / f'fdcheck_{base.name}' / 'log.txt',
                  [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
    out = sc.O / f'fdcheck_{base.name}'; out.mkdir(exist_ok=True)
    rj = base / 'QSTARS_RJ.bin'
    e = np.zeros(it.x.shape); e[row - 1, dep - 1] = 1.0
    xk = it.x
    s_k = 1.0 - xk / sc.x0
    s_e = 1.0 - xk * (1.0 - eps * e) / sc.x0
    tag = f'r{row}_d{dep}_e{eps:g}'
    jobs = {
        'frozen': lambda: R.run(out / f'{tag}_frozen', s_e, replay=rj, ba=False, store_j=False, two_phot_j=rj),
        'full': lambda: R.run(out / f'{tag}_full', s_e, warm=rj, two_phot_j=rj),
        'product': lambda: R.run(out / f'r{row}_d{dep}_product', s_k, replay=rj, dirs=e[:, :, None], ba=False,
                                 store_j=False, two_phot_j=rj),
    }
    with ThreadPoolExecutor(3) as ex:
        res = dict(zip(jobs, ex.map(lambda f: f(), jobs.values())))
    _, Rf = sc.read_trace_record(res['frozen'] / 'QSTARS_RESID_TRACE.bin')
    _, Ru = sc.read_trace_record(res['full'] / 'QSTARS_RESID_TRACE.bin')
    _, drad = ck.read_jv_rad(res['product'] / 'QSTARS_JV_RAD.bin')
    D = it.D
    fd_frozen = -(Rf - it.R) / eps * D
    fd_full = -(Ru - it.R) / eps * D
    jdir = it.Jdir.matvec(e) * D
    prod = jdir + drad[:, :, 0] * D
    big = np.argsort(np.abs(fd_full).ravel())[::-1][:12]
    big = list(dict.fromkeys(list(big) + list(np.argsort(np.abs(prod).ravel())[::-1][:6])))
    print(f'{base.name}: column of row {row} at depth {dep} (eps {eps:g}); entries x row scale D (equation, depth)')
    print(f'{"eq":>4s} {"d":>4s} {"frozen FD":>11s} {"Jdir":>11s} {"full FD":>11s} {"Jdir+rad":>11s}')
    for c in big:
        i, j = divmod(int(c), ck.ND)
        print(f'{i + 1:4d} {j + 1:4d} {fd_frozen[i, j]:+11.4e} {jdir[i, j]:+11.4e} {fd_full[i, j]:+11.4e} {prod[i, j]:+11.4e}')
    for nm, v in (('frozen FD - Jdir', fd_frozen - jdir), ('full FD - (Jdir + rad)', fd_full - prod)):
        print(f'{nm}: max {np.abs(v).max():.3e} at (eq {np.argmax(np.abs(v)) // ck.ND + 1}, d {np.argmax(np.abs(v)) % ck.ND + 1}); '
              f'norm {np.linalg.norm(v):.3e} vs column norm {np.linalg.norm(fd_full):.3e}')


if __name__ == '__main__':
    main()
