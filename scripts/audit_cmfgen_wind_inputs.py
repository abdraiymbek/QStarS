#!/usr/bin/env python3
"""Read-only Gate-4 input-chain audit: export g_rad versus native HYDRO g_RAD.

Plan: benchmarks/cmfgen_selected_wind_plan.json (registered before this audit).
Only same-pass exports qualify: HYDRO and QSTARS_SPECTRAL.bin written within
120 s with matching radius grids (USE_FIXED_J passes leave HYDRO stale).
Archived states are not converged released closures: this checks the input
chain, never a wind. No native code is run; nothing is written but the JSON.
"""
import argparse, json, time
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import sha256
from qstars.cmfgen_wind import load_wind_inputs, grad_input_audit, frozen_sources
from qstars.constants import G, MSUN

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ['cmfgen_ali_lambda84_long', 'cmfgen_ali_lambda333', 'cmfgen_alilin_fac0', 'cmfgen_cool_T9_bootstrap',
           'cmfgen_cool_T6_refined', 'cmfgen_mesh_84', 'cmfgen_mesh_167', 'cmfgen_mesh_333',
           'cmfgen_recovery_T9_fixedJ']  # last: stale-HYDRO negative control


def audit(d, threshold):
    spec, hyd = d/'QSTARS_SPECTRAL.bin', d/'HYDRO'
    row = dict(run=str(d.relative_to(ROOT)))
    if not spec.exists() or not hyd.exists(): return dict(row, status='missing_files')
    dt = abs(spec.stat().st_mtime-hyd.stat().st_mtime); row['mtime_difference_s'] = dt
    if dt > 120: return dict(row, status='not_same_pass', reason='HYDRO and spectral export from different passes')
    q = load_wind_inputs(d)
    if q['hydro_sync_radius'] > 2e-6: return dict(row, status='not_same_pass', reason='radius grids differ')
    mass = float(q['vadat']['MASS'])*MSUN
    a = grad_input_audit(q, mass=mass); arr = a.pop('arrays')
    s = frozen_sources(q, mass, 'radius'); x = s.nodes; a2, da2, gam, _ = s(x)
    n = (2*a2-da2-G*mass*(1-gam)/np.exp(x))*np.exp(x)/(G*mass)
    return dict(row, status='audited', files=q['files'], vadat={k: q['vadat'].get(k) for k in ('MASS', 'LSTAR', 'RSTAR', 'MDOT', 'VINF')},
                hydrogen_mass_fraction=q['hydrogen_mass_fraction'],
                gencool_channel_closure_max=float(np.max(q['cooling']['channel_closure_relative'])),
                **a, passes_L_corrected_threshold=a['L_corrected_relative_max'] <= threshold,
                illustrative_not_evidence=dict(
                    note='frozen radius-mode numerator on this (non-converged/prescribed-flow) state; not a hydro step',
                    scaled_numerator_at_base=float(n[0]), sign_changes=int(np.sum(np.sign(n[1:]) != np.sign(n[:-1]))),
                    gamma_max=float(gam.max()), radius_of_gamma_max=float(np.exp(x[np.argmax(gam)]))),
                per_depth=dict(radius=arr['radius'].tolist(), raw=arr['raw'].tolist(), L_corrected=arr['L_corrected'].tolist(),
                               momentum_relative=arr['momentum_relative'].tolist(),
                               native_percent_error=[None if not np.isfinite(v) else float(v) for v in arr['native_percent_error']]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('runs', nargs='*', default=DEFAULT)
    ap.add_argument('--out', default=str(ROOT/'Reports/evidence/report23/grad_input_audit.json'))
    a = ap.parse_args()
    plan = ROOT/'benchmarks/cmfgen_selected_wind_plan.json'
    threshold = json.loads(plan.read_text())['acceptance_thresholds']['spectral_vs_native_grad_relative_L_corrected_except_ND']
    t = time.time(); rows = []
    for name in a.runs:
        d = Path(name) if Path(name).is_absolute() else ROOT/'results'/name
        try: rows.append(audit(d, threshold))
        except (ValueError, OSError, KeyError) as exc: rows.append(dict(run=name, status='unreadable', reason=str(exc)))
    out = dict(audit='Gate 4 input chain: spectral g_rad vs native HYDRO g_RAD (Report 23)', created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               plan=str(plan.relative_to(ROOT)), plan_sha256=sha256(plan), threshold_L_corrected=threshold,
               module_sha256=sha256(ROOT/'src/qstars/cmfgen_wind.py'), script_sha256=sha256(Path(__file__)),
               definition='native g_RAD = kappa_F L(ND)/(4 pi r^2 c) (hydro_terms_v5.f:214); compared with g_spec after scaling by L_comoving(r)/L_comoving(r_ND); depth ND excluded and reported separately',
               evidence_scope='archived component/prescribed-flow states; none is a converged T-released closure; no hydro step was run',
               wall_seconds=None, runs=rows)
    out['wall_seconds'] = time.time()-t
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=1)+'\n')
    for r in rows:
        print(r['run'], r['status'], *(f"{k}={r[k]:.3e}" for k in ('raw_relative_max', 'L_corrected_relative_max', 'L_and_density_corrected_relative_max', 'momentum_relative_max') if k in r))


if __name__ == '__main__':
    main()
