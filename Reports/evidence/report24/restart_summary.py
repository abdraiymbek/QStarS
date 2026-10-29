#!/usr/bin/env python3
"""Report 24: the constrained restart of fixedT_b35 (from iterate 1, 2026-10-22 01:20) in the fixed measure.
Reads history.json and every it*_trial_* directory of the run; for each trial: native outcome, the driver's natural merit,
the fixed measure ||W M_ref^-1 R|| (preconditioner of it01_trial_l0.5 with its local-ALI increment, as fixed_metric.py),
its depth-region split, and conservation against x0. Writes restart.json (used by build_page.py) and prints a table."""
import html
import json
import re
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
E = Path(__file__).resolve().parent
O = ROOT / 'results/nk167_trap/fixedT_b35'; ref = 'it01_trial_l0.5'
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
cols = {n: (O / f'cols_{ref}' / f'col_{n}_p', O / f'cols_{ref}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(O / (ref + '_tp'), x0, cols, 1e-4, True, None, rj=O / ref / 'QSTARS_RJ.bin')
it.install_ali_increment(O / (ref + '_ali')); it.extend_significance(0.0, rows=[3])
tot = cnk_trap.conservation_totals(x0)
regs = [(0, 22), (22, 45), (45, 75), (75, 118), (118, 150), (150, 167)]


def measure(d):
    x, R = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
    w = it.W * it.M.solve(R) * it.active
    fh = x[cnk_trap.H_ROWS].sum(0) / tot[0]; fhe = x[cnk_trap.HE_ROWS].sum(0) / tot[1]
    ch = (x[cnk_trap.HII_ROW] + x[cnk_trap.HEII_ROWS].sum(0) + 2 * x[cnk_trap.HEIII_ROW] - x[cnk_trap.HMI_ROW]) / x[cnk_trap.NE_ROW]
    worst = float(max(np.max(np.abs(fh - 1)), np.max(np.abs(fhe - 1)), np.max(np.abs(ch - 1))))
    return dict(total=float(np.linalg.norm(w)), by_region=[float(np.linalg.norm(w[:, a:b])) for a, b in regs],
                conservation_worst=worst, he_over_he0=fhe.tolist(), charge_over_ne=ch.tolist())


hist = json.loads((O / 'history.json').read_text())
log = (O / 'stdout.log').read_text().splitlines()
rows, traj = [], [(1, 223.99)]
for r in hist:
    k = r['it']
    for lab, merit in r.get('trials', {}).items():
        d = O / f'it{k + 1:02d}_trial_{lab}'
        row = dict(iterate=k, label=lab, natural=None if merit == float('inf') else merit, accepted=(lab == r.get('accepted')))
        if (d / 'QSTARS_RESID_TRACE.bin').exists():
            row.update(measure(d)); row['native'] = 'completed'
        else:
            nw = d / 'QSTARS_ALI_NEWTON'
            alb = [float(m) for m in re.findall(r'max_scat_over_chi\s+([0-9.E+-]+)', nw.read_text())] if nw.exists() else []
            ok = re.findall(r' ok ([TF]) ', nw.read_text()) if nw.exists() else []
            if alb and 'F' in ok:
                row['native'] = f'failed: maser, CHI_SCAT/CHI {max(alb):.2f}, no solution'
            else:
                msg = [ln for ln in log if f'it{k:02d} trial {lab}: rejected' in ln]
                row['native'] = 'failed' + (': ' + re.sub(r'.*no residual \(', '', msg[-1]).rstrip(')')[:60] if msg else '')
        rows.append(row)
    acc = r.get('accepted')
    if acc and r.get('merit_new', np.inf) < r.get('merit', 0):
        m = [q for q in rows if q['iterate'] == k and q['label'] == acc]
        if m and 'total' in m[0]:
            traj.append((k + 1, m[0]['total']))
print(f'{"iter":>4s} {"trial":>8s} {"native":>10s} {"natural":>10s} {"fixed":>8s} {"H fr":>7s} {"deep":>7s} {"cons":>8s}')
for q in rows:
    deep = np.linalg.norm(q['by_region'][4:]) if 'by_region' in q else float('nan')
    print(f"{q['iterate']:4d} {q['label']:>8s} {q['native'][:10]:>10s} {q['natural'] if q['natural'] else float('inf'):10.4g} "
          f"{q.get('total', float('nan')):8.2f} {q['by_region'][1] if 'by_region' in q else float('nan'):7.2f} {deep:7.2f} "
          f"{q.get('conservation_worst', float('nan')):8.1e}{'  accepted' if q['accepted'] else ''}")
trs = ''.join(
    f"<tr><td>{q['iterate']}</td><td>λ = {html.escape(q['label'][1:].replace('_d16', ' / 16'))}</td><td>{html.escape(q['native'])}</td>"
    f"<td>{'—' if 'total' not in q else f'{q['total']:.1f}'}</td>"
    f"<td>{'—' if 'by_region' not in q else f'{q['by_region'][1]:.1f}'}</td>"
    f"<td>{'—' if 'by_region' not in q else f'{float(np.linalg.norm(q['by_region'][4:])):.1f}'}</td>"
    f"<td>{'—' if 'conservation_worst' not in q else f'{q['conservation_worst']:.0e}'}</td>"
    f"<td>{'accepted' if q['accepted'] else ''}</td></tr>" for q in rows)
table = ('<div class="tablewrap"><table class="trials"><thead><tr><th>Iterate</th><th>Trial</th><th>Native pass</th>'
         '<th>Fixed measure</th><th>H fronts</th><th>Deep (119–167)</th><th>Worst conservation error</th><th></th></tr></thead>'
         f'<tbody>{trs}</tbody></table></div>')
old = [(1, 223.99), (2, 152.61), (3, 130.29), (4, 118.20)]
it04p = measure(O / 'branch_unconstrained_it2plus' / 'it05_trial_l0')['total']   # projected iterate-4 base (00:29 restart)
out = dict(rows=rows, trajectory_old=old, trajectory_new=traj, table_html=table, it04_projected=it04p,
           yticks=[50, 100, 200, 400])
(E / 'restart.json').write_text(json.dumps(out) + '\n')
print('trajectory (constrained):', traj)
