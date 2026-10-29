#!/usr/bin/env python3
"""Report 24: numbers behind the summary page (Reports/evidence/report24/page_data.json), all read from run outputs.
fixed measure by depth region (preconditioner of fixedT_b35/it01_trial_l0.5 with its local-ALI increment), conservation
of each pass against x0's totals, the albedo profile of the failing trial frequency, the SVD-floor q statistics and the
Newton-fallback spectral radii."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'precond', 0.0, 1e-3
E = ROOT / 'Reports/evidence/report24'
O = ROOT / 'results/nk167_trap/fixedT_b35'; ref = 'it01_trial_l0.5'
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
cols = {n: (O / f'cols_{ref}' / f'col_{n}_p', O / f'cols_{ref}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(O / (ref + '_tp'), x0, cols, 1e-4, True, None, rj=O / ref / 'QSTARS_RJ.bin')
it.install_ali_increment(O / (ref + '_ali')); it.extend_significance(0.0, rows=[3])
regs = [(0, 22, '1–22'), (22, 45, '23–45'), (45, 75, '46–75'), (75, 118, '76–118'), (118, 150, '119–150'), (150, 167, '151–167')]
tot = cnk_trap.conservation_totals(x0)
passes = [a for a in sys.argv[1:]] or ['start_base', 'it01_trial_l0.5', 'it02_trial_l1', 'it03_trial_l0.25', 'it04_trial_l0.25']
data = dict(regions=[r[2] for r in regs], passes=[])
for b in passes:
    d = O / (b + '_tp') if (O / (b + '_tp') / 'QSTARS_RESID_TRACE.bin').exists() else O / b
    for br in ('branch_unconstrained_it2plus', 'branch_b35_noconserve'):     # moved at the restarts of 2026-10-22
        if not (d / 'QSTARS_RESID_TRACE.bin').exists():
            d = O / br / (b + '_tp') if (O / br / (b + '_tp') / 'QSTARS_RESID_TRACE.bin').exists() else O / br / b
    x, R = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
    w = it.W * it.M.solve(R) * it.active
    fh = x[cnk_trap.H_ROWS].sum(0) / tot[0]; fhe = x[cnk_trap.HE_ROWS].sum(0) / tot[1]
    ch = (x[cnk_trap.HII_ROW] + x[cnk_trap.HEII_ROWS].sum(0) + 2 * x[cnk_trap.HEIII_ROW] - x[cnk_trap.HMI_ROW]) / x[cnk_trap.NE_ROW]
    data['passes'].append(dict(name=b, dir=str(d.relative_to(ROOT)), total=float(np.linalg.norm(w)),
                               by_region=[float(np.linalg.norm(w[:, r0:r1])) for r0, r1, _ in regs],
                               by_group={g: float(np.linalg.norm(w[r0:r1])) for g, r0, r1 in
                                         (('H I', 1, 31), ('H II', 31, 32), ('He', 32, 132), ('n_e', 132, 133))},
                               he_over_he0=fhe.tolist(), h_over_h0=fh.tolist(), charge_over_ne=ch.tolist()))
alb = np.loadtxt(E / 'albedo_profile_it04_trial_l0.5_f35609.txt')
data['albedo'] = dict(depth=alb[:, 0].astype(int).tolist(), scat_over_chi=alb[:, 4].tolist())
data['floor_q'] = json.loads((E / 'floor_q.json').read_text())
nd = E / 'newton_dump_it04_trial_l0.5.json'
data['newton_dump'] = json.loads(nd.read_text()) if nd.exists() else None
(E / 'page_data.json').write_text(json.dumps(data) + '\n')
for p in data['passes']:
    print(f"{p['name']:24s} total {p['total']:8.2f} regions " + ' '.join(f'{v:7.2f}' for v in p['by_region']) +
          f"  He/He0 [{min(p['he_over_he0']):.3f},{max(p['he_over_he0']):.3f}] charge [{min(p['charge_over_ne']):.3f},{max(p['charge_over_ne']):.3f}]")
