#!/usr/bin/env python3
"""Complete native verification iteration at a Newton-Krylov state, then audits.

The state is reached from the archived restart x0 through the read-only
perturbation hook (the accepted base pass's QSTARS_STEP.bin). One full native
iteration runs with the spectral export and WITHOUT stopping after the solve,
so RVTJ, GENCOOL, OUTGEN and QSTARS_SPECTRAL.bin are written. The native
correction of that iteration is the fixed-J frozen-J step at the state (reported
separately from the exact Newton step). Then scripts/audit_cmfgen_cool.py runs.
"""
import argparse, json, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'backends' / 'cmfgen'))
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('state', type=Path, help='accepted base pass directory (has QSTARS_STEP.bin)')
    ap.add_argument('dest', type=Path)
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--free-t', action='store_true')
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'),
                    help='extra VADAT settings (e.g. the doubled-frequency grid), as used by the NK run')
    a = ap.parse_args()
    from formal_ali_trial import prepare
    from run_component import run
    dest = a.dest.resolve()
    sets = [('QS_ALI_LIN', 'T'), ('FIX_T', 'F' if a.free_t else 'T'), ('QS_LIN_FAC', '0')] + [tuple(v) for v in a.vadat]
    prepare(a.x0, dest, 1, 1e-10, 2000, sets)          # touches QSTARS_EXPORT
    shutil.copy2(a.state / 'QSTARS_STEP.bin', dest / 'QSTARS_STEP.bin')
    if a.free_t:
        from cool_newton_krylov import release_t
        release_t(dest)
    env = {'QSTARS_PERTURB_FILE': str(dest / 'QSTARS_STEP.bin'), 'QSTARS_PERTURB_EPS': '1.0',
           'QSTARS_REPLACE_FORCE': str(a.replace.resolve())}
    r = run(dest, a.exe.resolve(), 40000, a.threads, extra_env=env)
    print(json.dumps({k: r[k] for k in ('status', 'error', 'seconds')}))
    out = dest.with_name(dest.name + '_audit')
    p = subprocess.run([sys.executable, str(ROOT / 'scripts/audit_cmfgen_cool.py'), str(dest), str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    (dest / 'audit_stdout.txt').write_text(p.stdout + p.stderr)
    print(p.stdout[-4000:], p.stderr[-2000:])


if __name__ == '__main__':
    main()
