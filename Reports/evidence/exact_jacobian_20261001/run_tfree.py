"""T-released (radiative-equilibrium row ACTIVE) validation of the complete Jacobian action at x0.

Earlier "FIX_T=F" passes (results/jv/oracle_t, b21baseT) inherited FIXED_T=.TRUE. from the
restart's BAMATPNT, so the radiative-equilibrium row was the identity with zero residual
(T column validated, T row not). Here every pass goes through cool_newton_krylov.Runner with
fix_t=False, which removes BAMATPNT/BAMAT (release_t).

base       : x0, BA dump, stored J, products for the 4 directions of results/jv/dirs_t.npy
col_*_p/m  : frozen-J replay passes, n_e / H II / T columns (h = 1e-4)
<dir>_a/_b : exact-residual passes at x0(1 - eps s) (full formal-ALI, no warm start)
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402

EXE = ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe'
X0 = ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3'
REP = ROOT / 'results/jv/replace_fac1_x0.txt'
NAMES = ['T', 'ne', 'rand', 'mix']
EPS = {'T': (1e-4, 1e-5), 'ne': (1e-3, 1e-4), 'rand': (1e-3, 1e-4), 'mix': (1e-4, 1e-5)}


def main():
    R = ck.Runner(X0, EXE, REP, False, 2, HERE / 'driver.log')
    d = np.load(ROOT / 'results/jv/dirs_t.npy')
    nt, nd = d.shape[:2]
    base = R.run(HERE / 'base', np.zeros((nt, nd)), dirs=d, ba=True, store_j=True, threads=4)
    R.log('base done')
    jobs = []
    for n, row in (('ne', 133), ('hii', 32), ('T', 134)):
        for sg in (+1, -1):
            e = np.zeros((nt, nd)); e[row - 1, :] = sg * 1e-4
            jobs.append(dict(dest=HERE / f'col_{n}_{"p" if sg > 0 else "m"}', s=1.0 - (1.0 + e),
                             replay=base / 'QSTARS_RJ.bin'))
    for k, n in enumerate(NAMES):
        for lab, eps in zip('ab', EPS[n]):
            jobs.append(dict(dest=HERE / f'{n}_{lab}', s=eps * d[:, :, k], replay=None))

    def one(j):
        R.run(j['dest'], j['s'], replay=j['replay'], ba=False, store_j=False, threads=2)
        R.log(f'{j["dest"].name} done')
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(one, jobs))
    R.log('all done')


if __name__ == '__main__':
    main()
