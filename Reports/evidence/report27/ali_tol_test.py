#!/usr/bin/env python3
"""Report 27: is the irreducible d62-100 residual the native ALI's convergence error? The accepted state it26_trial_l0.5
is recomputed exactly as its trial pass was (same step, warm J, two-photon J) but with QS_ALI_TOL 1e-13 (and QS_ALI_MAX
20000) instead of 1e-10; the residuals are compared in the frozen measure by region and row group. Pass at nice 10."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'Reports/evidence/v1_exec_20261106/p1'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record, read_fortran  # noqa: E402
from p1_native import EXE, REPLACE, X0DIR, FROZEN  # noqa: E402
RUN = ROOT / 'results/nk167_trap/v1_p2'
st = RUN / 'it26_trial_l0.5'
tol = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-13
env = json.loads((st / 'run_manifest.json').read_text())['extra_env']
x0b = cn.configure()
x, R0 = read_trace_record(st / 'QSTARS_RESID_TRACE.bin')
step = read_fortran(st / 'QSTARS_STEP.bin', x.shape)
R = ck.Runner(X0DIR, EXE, REPLACE, False, 4, RUN / 'report27_ali_tol.log',
              [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F'), ('QS_ALI_TOL', f'{tol:.3E}'), ('QS_ALI_MAX', '20000')])
R.no_obs = True
dest = RUN / f'report27_ali_tol{tol:.0e}'
R.run(dest, step, threads=4, warm=env['QSTARS_ALI_WARM'], two_phot_j=env['QSTARS_TWO_PHOT_J'])
x1, R1 = read_trace_record(dest / 'QSTARS_RESID_TRACE.bin')
assert np.array_equal(x1, x), 'state differs'
w = np.load(FROZEN)['w']
a2, b2, d2 = (w * R0[:134]) ** 2, (w * R1[:134]) ** 2, (w * (R1 - R0)[:134]) ** 2
print(f'frozen measure: tol 1e-10 {np.sqrt(a2.sum()):.6e}, tol {tol:g} {np.sqrt(b2.sum()):.6e}, of the difference {np.sqrt(d2.sum()):.3e}')
for a_, b_ in cn.REGIONS:
    s = slice(a_ - 1, b_)
    print(f'  d{a_}-{b_}: {np.sqrt(a2[:, s].sum()):.3e} -> {np.sqrt(b2[:, s].sum()):.3e} (difference {np.sqrt(d2[:, s].sum()):.3e})')
for g, r in cn.GROUPS.items():
    if a2[r].sum() > 0.01 * a2.sum():
        print(f'  {g}: {np.sqrt(a2[r].sum()):.3e} -> {np.sqrt(b2[r].sum()):.3e}')
print(open(dest / 'QSTARS_ALI_STATS').read().splitlines()[-1])
