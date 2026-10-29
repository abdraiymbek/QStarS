#!/usr/bin/env python3
"""Report 23: update Reports/evidence/completion_gates.json (indent=2, as at HEAD). Re-run when the closure status changes.
usage: update_completion_gates.py CLOSURE_STATUS_TEXT"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
p = ROOT / 'Reports/evidence/completion_gates.json'
g = json.loads(p.read_text())
closure = sys.argv[1]
if 'previous_release_blocker_report22' not in g:
    g['previous_release_blocker_report22'] = g.get('release_blocker')
g['assessment_date_utc'] = '2026-10-14'
g['latest_report'] = 'Reports/Report_23.md'
g['automated_tests_passed'] = 334
g['automated_tests_failed'] = 0
g['automated_tests_skipped'] = 37
g['automated_tests_scope'] = ('2026-10-14, full suite once (single thread): compiled MESA/FreeEOS library tests and native-compiler '
                              'tests skip unless their environment variables are configured.')
g['release_blocker'] = (
    'T-released cool closure at 167 radii, NOT CONVERGED. Diagnosed as numerical (Report 23 section 2): the Newton-Krylov '
    'measure used fixed-J Jacobian blocks that hide the error of every optically thick transition (H I n=1-9, He I and He II '
    'resonance levels), and every single-iteration pass ran the OLD_DEFAULT two-photon approximation (about 75 % of the RE row '
    'at the H I front). Earlier "converged" cool states are converged only in that measure. Fixes: build 32 (two-photon USE_RAD '
    'from stored J), local-ALI preconditioner for every column (build 33 makes it affordable), trapped-level driver options. '
    'Status: ' + closure)
g['report23_numerical_checks'] = {
    'builds': {'32': '290d1490f87e29e3fa269fffaa1ed803c392065feb0f890c6403399579a0fe3e',
               '33': '6d72dbadb10e3e785751b9c3686e2db681a38ba46136fe7a95583cf2ca7977d6',
               '34': 'd3e227737611c3bc473eab3510668921951f83a8140cb37272f254d1179829d2',
               '35': 'c69836d31fd0425b93ea0cf529e988cdd356ad9a6ea50d0c6949492355a3a6dd'},
    'build34_ff_gaunt_dT': {'tangent_over_fd_deta_dlnT_167_d80_90_before': [1.97, 1.96, 1.54, 1.52, 1.31], 'after': 1.0,
                            're_row_T_column_error_d80_90': {'build33': 3.5e-2, 'build34': 1.19e-5, 'oracle': 1.24e-5},
                            'rows_above_3x_oracle_84_radii': {'build33': 111, 'build34': 2}},
    'build35_inner_boundary_dDTDR_dT_error': {'build34': 6.97e-4, 'build35': 7.3e-6},
    'build32_feature_off_vs_build31': 'byte-identical residual, BA, products (84 and 167 radii, one thread)',
    'build32_two_photon_up_rate_d35': {'old_default': 186.64, 'use_rad': 3.081, 'band_exact_replica': 3.08097},
    'build32_energy_identity_without_D_d30_41': [2e-5, 9e-4],
    'build32_tangent_population_directions_max_error': 5e-6,
    'build33_local_ali_pass': {'instructions_ratio': 17.8, 'cycles_ratio': 20.2, 'ba_identical_entries': '8995956 of 8995956',
                               'lin_fac0_vs_build32': 'byte-identical'},
    'n2_exact_over_fixedJ_diagonal_d26_46': [3.05e-4, 6.5e-4],
    'ali_n2_column_vs_exact_local_l1_error': [0.001, 0.051],
    'preconditioner_q_median': {'fixed_J': 0.88, 'fixed_J_plus_exact_n2': 0.88, 'local_ali_all_columns': 0.40},
    'max_sig_Minv_R_pin_state': {'fixed_J': 7.7e-3, 'exact_n2_column': 1.63},
    'departure_coefficient_b1_front_ptc_state_max': 6.9e10,
    'fixed_T_trapped_solver_natural_measure': [139.8, 97.3, 45.7],
    't_released_build34_fixed_measure': [224.0, 130.7, 87.3],
    'ball2011_stage7': {'amendments_registered_in_order': [3, 4, 2, 1], 'row0_converged_unique': ['C1', 'C0', 'C2', 'C5', 'C6', 'C7'],
                        'smoke_gate': 'failed at M_BH = 46.0 Msun after 1255 accepted steps (numerical blocker)',
                        'smoke_per_step_energy_residual_max': 1.184e-14, 'smoke_per_step_mass_residual_max': 3.562e-16,
                        'q2_row0_max_density_inversion': {'N200': 9.07e-3, 'N400': 5.36e-3, 'N800': 2.75e-3, 'limit': 1e-3},
                        'amendment5': 'rejected (unbounded tolerance relaxation)'},
    'evidence': ['Reports/Report_23.md', 'Reports/evidence/report23/trapped_levels/README.md',
                 'Reports/evidence/report23/ball/stage7_commands.log', 'results/two_photon_test/b32/README.txt',
                 'results/fast_ali_b33/README.txt']}
oc = g.get('outer_campaign_completion', {})
oc.update(outcome='NUMERICAL BLOCKER (trapped-transition measure; closure re-solve in progress)', released_cool_closure_converged=False,
          published_full_evolution_replicas=0)
g['outer_campaign_completion'] = oc
g['published_full_evolution_replication'] = False
g['published_full_evolution_replicas'] = 0
ir = g.get('intermediate_reports_without_gate_update', [])
g['intermediate_reports_without_gate_update'] = ir
p.write_text(json.dumps(g, indent=2) + '\n')
print('updated', p)
