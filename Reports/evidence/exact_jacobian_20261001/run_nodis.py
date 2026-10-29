"""Hypothesis test for the RE-row error of the complete action (Report 20 §4.3): repeat the
T-released validation with level dissolution OFF (VADAT DO_DIS=F). If the native opacity
derivatives' missing dissolution dependence causes the RE-row error, the error and the need
for the POPION completion should disappear. Same x0, same directions (n_e, random) as
results/jv/tfree; one oracle step each. This is a diagnostic of the Jacobian only; the
production model keeps dissolution."""
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


def main():
    R = ck.Runner(X0, EXE, REP, False, 1, HERE / 'driver.log', [('DO_DIS', 'F')])
    d = np.load(ROOT / 'results/jv/dirs_t.npy')
    nt, nd = d.shape[:2]
    base = R.run(HERE / 'base', np.zeros((nt, nd)), dirs=d, ba=True, store_j=True, threads=2)
    R.log('base done')
    jobs = []
    for n, row in (('ne', 133), ('hii', 32), ('T', 134)):
        for sg in (+1, -1):
            e = np.zeros((nt, nd)); e[row - 1, :] = sg * 1e-4
            jobs.append((HERE / f'col_{n}_{"p" if sg > 0 else "m"}', -e, base / 'QSTARS_RJ.bin'))
    for k, n, eps in ((1, 'ne', 1e-3), (2, 'rand', 1e-3), (0, 'T', 1e-4)):
        jobs.append((HERE / f'{n}_a', eps * d[:, :, k], None))

    def one(j):
        R.run(j[0], j[1], replay=j[2], ba=False, store_j=False, threads=1)
        R.log(f'{j[0].name} done')
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(one, jobs))
    R.log('all done')


if __name__ == '__main__':
    main()
