#!/usr/bin/env python3
"""Report 25: progress of the constrained stage-1 solve (fixedT_b35) in the two measures of Report 24, by depth region,
plus conservation against x0's totals. Read-only.

    progress.py PASS [PASS ...]      (pass directory names under fixedT_b35; a <name>_tp replay is used when present)

fixed : ||W M^-1 R|| with the preconditioner and weights of it01_trial_l0.5 (local-ALI increment, no floor), as
        scripts/fixed_metric.py and Reports/evidence/report24/collect_page_data.py;
rows  : ||W D R|| over populations and n_e, D = 1 / row l1-norm of it01's direct blocks (rowscale_by_region.py);
rows* : the same with --eq-floor 1e-15, as the run since iterate 8 (equations whose gross scale is below 1e-15 of the largest at their depth, e.g.
        the He II / He III balance where He III is 1e-11..1e-43 of He, are not scaled up);
cons  : largest |factor - 1| of the H, He and charge projections (cnk_trap.project_conservation)."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'precond', 0.0, 1e-3
O = ROOT / 'results/nk167_trap/fixedT_b35'; ref = 'it01_trial_l0.5'
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
cols = {n: (O / f'cols_{ref}' / f'col_{n}_p', O / f'cols_{ref}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(O / (ref + '_tp'), x0, cols, 1e-4, True, None, rj=O / ref / 'QSTARS_RJ.bin')
D, Wr = it.D.copy(), it.W.copy()
rs_ = 1.0 / np.where(D > 0, D, np.inf)
Df = np.where(rs_ > 0, 1.0 / np.maximum(rs_, 1e-15 * rs_.max(axis=0, keepdims=True)), 0.0)   # --eq-floor 1e-15 (the run)
it.install_ali_increment(O / (ref + '_ali')); it.extend_significance(0.0, rows=[3])
Wr = it.W                                     # rowscale_by_region.py: same weights (rows=[3]), no ALI needed for D
regs = [(0, 22), (22, 45), (45, 75), (75, 118), (118, 150), (150, 167)]
tot = cnk_trap.conservation_totals(x0)


def find(b):
    for d in (O / (b + '_tp'), O / b, O / 'superseded/it05_trials_pre_front' / b):
        if (d / 'QSTARS_RESID_TRACE.bin').exists():
            return d
    return None


def main():
    hdr = ''.join(f'{f"{a + 1}-{b}":>8s}' for a, b in regs)
    print(f'{"pass":22s}  fixed{hdr}  | rows {hdr} | rows* | cons')
    for b in sys.argv[1:]:
        d = find(b)
        if d is None:
            print(f'{b:22s}  (no residual)'); continue
        x, R = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
        w = it.W * it.M.solve(R) * it.active
        r = (Wr * D * R)[:133]
        _, st = cnk_trap.project_conservation(x, tot)
        print(f'{b:22s} {np.linalg.norm(w):6.1f}' + ''.join(f'{np.linalg.norm(w[:, a:b_]):8.1f}' for a, b_ in regs)
              + f'  | {np.linalg.norm(r):5.2f}' + ''.join(f'{np.linalg.norm(r[:, a:b_]):8.3f}' for a, b_ in regs)
              + f' | {np.linalg.norm((Wr * Df * R)[:133]):5.3f} | {max(st.values()):.0e}')


if __name__ == '__main__':
    main()
