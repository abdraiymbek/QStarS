#!/usr/bin/env python3
"""Complete native verification iteration at a Newton-Krylov state, then audits.

The state is reached from the archived restart x0 through the read-only
perturbation hook (the accepted base pass's QSTARS_STEP.bin). One full native
iteration runs with the spectral export and WITHOUT stopping after the solve,
so RVTJ, GENCOOL, OUTGEN and QSTARS_SPECTRAL.bin are written. The native
correction of that iteration is the fixed-J frozen-J step at the state (reported
separately from the exact Newton step). Then scripts/audit_cmfgen_cool.py runs.

V1 plan P3 (Report 26 review): the executable, restart and equation table are required (the defaults were builds 23
and the 84-radius tables). With --two-phot-j STORED_J (the solve's --two-phot-j: the state's own QSTARS_RJ.bin) the pass
uses TWO_METH=USE_RAD with that conjugate-frequency J, as every pass of the solve did, stores its own J, and writes
two_photon_selfconsistency.json: the relative difference between the J it was given and the J it produced (a
converged state reproduces its own stored J).
"""
import argparse, json, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'backends' / 'cmfgen'))
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))


def two_photon_selfconsistency(given, produced):
    """Relative difference of the J a verification pass produced against the stored J it was given (same grid)."""
    import numpy as np
    from qstars.nlte_newton import read_rj
    ml0, nu0, j0 = read_rj(given)
    ml1, nu1, j1 = read_rj(produced)
    if not (np.array_equal(ml0, ml1) and np.array_equal(nu0, nu1)):
        return dict(same_grid=False)
    rel = np.abs(j1 - j0) / np.maximum(np.abs(j0), 1e-300)
    return dict(same_grid=True, given=str(given), produced=str(produced), max_rel=float(rel.max()),
                median_rel=float(np.median(rel)), p99_rel=float(np.quantile(rel, 0.99)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('state', type=Path, help='accepted base pass directory (has QSTARS_STEP.bin)')
    ap.add_argument('dest', type=Path)
    ap.add_argument('--x0', type=Path, required=True)
    ap.add_argument('--exe', type=Path, required=True)
    ap.add_argument('--replace', type=Path, required=True)
    ap.add_argument('--two-phot-j', type=Path, default=None,
                    help='stored conjugate-frequency J of the state (its QSTARS_RJ.bin): TWO_METH=USE_RAD as in the solve')
    ap.add_argument('--free-t', action='store_true')
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'),
                    help='extra VADAT settings (e.g. the doubled-frequency grid), as used by the NK run')
    a = ap.parse_args()
    from formal_ali_trial import prepare
    from run_component import run
    dest = a.dest.resolve()
    sets = [('QS_ALI_LIN', 'T'), ('FIX_T', 'F' if a.free_t else 'T'), ('QS_LIN_FAC', '0')] + [tuple(v) for v in a.vadat]
    tpj = a.two_phot_j.resolve() if a.two_phot_j else None
    if tpj is not None:
        cap = a.exe.resolve().with_suffix('.capabilities.json')
        if not (cap.is_file() and json.loads(cap.read_text()).get('qstars_two_phot_stored_j') is True):
            raise SystemExit(f'--two-phot-j needs an executable recorded with qstars_two_phot_stored_j: {a.exe}')
        if any(k == 'TWO_METH' for k, _ in sets):
            raise SystemExit('--two-phot-j sets TWO_METH=USE_RAD itself')
        sets.append(('TWO_METH', 'USE_RAD'))
    prepare(a.x0, dest, 1, 1e-10, 2000, sets)          # touches QSTARS_EXPORT
    shutil.copy2(a.state / 'QSTARS_STEP.bin', dest / 'QSTARS_STEP.bin')
    if a.free_t:
        from cool_newton_krylov import release_t
        release_t(dest)
    env = {'QSTARS_PERTURB_FILE': str(dest / 'QSTARS_STEP.bin'), 'QSTARS_PERTURB_EPS': '1.0',
           'QSTARS_REPLACE_FORCE': str(a.replace.resolve())}
    if tpj is not None:
        env.update(QSTARS_TWO_PHOT_J=str(tpj), QSTARS_RJ_STORE='1')
    r = run(dest, a.exe.resolve(), 40000, a.threads, extra_env=env)
    print(json.dumps({k: r[k] for k in ('status', 'error', 'seconds')}))
    if tpj is not None and (dest / 'QSTARS_RJ.bin').exists():
        (dest / 'two_photon_selfconsistency.json').write_text(json.dumps(two_photon_selfconsistency(tpj, dest / 'QSTARS_RJ.bin'), indent=1))
    out = dest.with_name(dest.name + '_audit')
    p = subprocess.run([sys.executable, str(ROOT / 'scripts/audit_cmfgen_cool.py'), str(dest), str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    (dest / 'audit_stdout.txt').write_text(p.stdout + p.stderr)
    print(p.stdout[-4000:], p.stderr[-2000:])


if __name__ == '__main__':
    main()
