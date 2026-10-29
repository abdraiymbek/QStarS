#!/usr/bin/env python3
"""Report 24: update Reports/evidence/completion_gates.json (indent=2). Re-run when the closure status changes.
usage: update_completion_gates.py CLOSURE_STATUS_TEXT"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
p = ROOT / 'Reports/evidence/completion_gates.json'
g = json.loads(p.read_text())
closure = sys.argv[1]
if 'previous_release_blocker_report23' not in g:
    g['previous_release_blocker_report23'] = g.get('release_blocker')
g['assessment_date_utc'] = '2026-10-22'
g['latest_report'] = 'Reports/Report_24.md'
g['automated_tests_passed'] = 341
g['automated_tests_failed'] = 0
g['automated_tests_skipped'] = 37
g['automated_tests_scope'] = ('2026-10-22, full suite once (nice 15, shared machine): compiled MESA/FreeEOS library tests and '
                              'native-compiler tests skip unless their environment variables are configured.')
g['release_blocker'] = (
    'Cool closure at 167 radii, NOT CONVERGED. Report 24: the Report 23 solves violated the exact number and charge '
    'conservation rows (He up to 8x, H 2.3x, charge 34x) because the per-component step limiter clips populations '
    'separately; the least squares used distant directions through the nearly singular local-ALI front blocks (coefficient '
    '172) and moved T in the fixed-T stage (coefficient 358); full trial steps created supercritical infrared masers '
    '(Lambda spectral radius 1.063, 1.058). Fixes: --conserve (trial projection onto the constraints), --pc-svd-floor 1e-2, '
    '--freeze-t, build 36 (Newton solve of the scattering fixed point where the local ALI stalls). Status: ' + closure)
g['report24_numerical_checks'] = {
    'builds': {'36': json.loads((ROOT / 'backends/cmfgen/ali_newton_patch_manifest.json').read_text())['reproducibility']['executable_sha256']},
    'conservation_violation_it04_fixedT_b35': {'he_total_over_x0': [0.128, 1.236], 'h_total_over_x0': [0.880, 2.350],
                                               'charge_sum_over_ne': [0.059, 33.9]},
    'fixed_measure_by_region_start_vs_it04': {'h_fronts_23_45': [350.0, 77.7], 'deep_119_167': [0.22, 72.4]},
    'maser_lambda_spectral_radius': [1.063, 1.058],
    'precond_floor_1e-2_worst_q': {'raw_local_ali_max': 236.8, 'floored_max': 9.76},
    'lsq_far_field_coefficient_it01': {'coefficient': 172.0, 'action_outside_support': 0.967},
    'build36_neutrality_one_thread': 'residual and BA traces byte-identical to build 35 (converged pass, fallback unused)',
    'build36_neutrality_two_threads_rel_resid_diff': 6.8e-13,
    'evidence': ['Reports/Report_24.md', 'Reports/evidence/report24/', 'results/nk167_trap/fixedT_b35/README.txt',
                 'results/nk167_trap/b36_tests/']}
oc = g.get('outer_campaign_completion', {})
oc.update(outcome='NUMERICAL BLOCKER (constrained re-solve of stage 1 in progress)', released_cool_closure_converged=False)
g['outer_campaign_completion'] = oc
p.write_text(json.dumps(g, indent=2) + '\n')
print('updated', p)
