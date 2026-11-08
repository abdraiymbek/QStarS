#!/usr/bin/env python3
"""Report 27: one local relaxation sweep of the low-population levels with the LOCAL-ALI blocks (QS_LIN_FAC=1: the
local radiation response included, so trapped levels are not solved at fixed J; the fixed-J version diverged x1500).
Levels with fraction <= THRESH are solved with CMFGEN's structural pairing (native_correction_check.structural_step) on
the ALI blocks of the base; the update is limited to a factor 2 (|log f| <= 0.693, as the solver's component limit),
then projected onto x0's H/He totals and the charge sum; one native pass (warm and two-photon J of the base, as a
driver trial); frozen and local merit before/after.   ali_relax_test.py BASE_DIR THRESH TAG"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'Reports/evidence/v1_exec_20261106/p1'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
import native_correction_check as ncc  # noqa: E402
from qstars.nlte_newton import read_trace_record, read_ba_blocks  # noqa: E402
from p1_native import EXE, REPLACE, X0DIR, FROZEN  # noqa: E402
base, th, tag = Path(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
D0 = int(sys.argv[4]) if len(sys.argv) > 4 else 1          # first relaxed depth (1-based); fronts above it untouched
x0b = cn.configure()
tp, ali = base.with_name(base.name + '_tp'), base.with_name(base.name + '_ali')
x, R = read_trace_record(tp / 'QSTARS_RESID_TRACE.bin')
xa, _ = read_trace_record(ali / 'QSTARS_RESID_TRACE.bin')
assert np.array_equal(xa, x)
lo, di, up = read_ba_blocks(ali / 'QSTARS_BA_TRACE.bin', 134, 167)
rep = ncc.read_replace(REPLACE, 167)
s, solved, held, bw = ncc.structural_step(x, R, lo, di, up, rep, th)
s[:, :D0 - 1] = 0.0; solved[:, :D0 - 1] = False
print(f'{int(solved.sum())} variables relaxed (fraction <= {th:g}), {int(held.sum())} held; backward error {bw:.1e}; '
      f'max|s| {np.abs(s[solved]).max():.3e}, median {np.median(np.abs(s[solved])):.1e}')
f = cnk_trap.update_factor(s, 'linear', 0.693)
xn = x * f
x0, _ = read_trace_record(ROOT / 'results/nk167_trap/free_r25/it00_base/QSTARS_RESID_TRACE.bin')
xn, pst = cnk_trap.project_conservation(xn, cnk_trap.conservation_totals(x0))
print('projection', pst, 'limited', int((np.abs(np.log(np.maximum(f, 1e-300))) >= 0.692).sum()))
Rn_ = ck.Runner(X0DIR, EXE, REPLACE, False, 4, base.parent / 'report27_relax.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
Rn_.no_obs = True
d = base.parent / f'report27_relax_{tag}'
Rn_.run(d, 1.0 - xn / x0, threads=4, warm=base / 'QSTARS_RJ.bin', two_phot_j=base / 'QSTARS_RJ.bin')
for p in (tp, d):
    sc = cn.score(FROZEN, p, x0b)
    print(f'{p.name}: frozen {sc["frozen"]:.6e}, local {sc["local"]:.6e}; regions ' + ', '.join(f'{k} {v:.3f}' for k, v in sc['regions'].items()))
