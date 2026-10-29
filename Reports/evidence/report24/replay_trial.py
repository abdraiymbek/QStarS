#!/usr/bin/env python3
"""Report 24: re-run one line-search trial pass of cool_newton_krylov.py (same state, warm J, two-photon J, VADAT and
replacement list as the driver) with another executable, and report how the pass ended.

    replay_trial.py TRIAL OUT --exe EXE [--x0 X0] [--replace FILE] [--threads N] [--env K=V ...]
TRIAL is a trial directory written by the driver (its QSTARS_STEP.bin and run_manifest.json define the pass).
"""
import argparse, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_fortran  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('trial', type=Path); ap.add_argument('out', type=Path); ap.add_argument('--exe', type=Path, required=True)
ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda167')
ap.add_argument('--threads', type=int, default=2); ap.add_argument('--env', action='append', default=[])
ap.add_argument('--vadat', nargs=2, action='append', default=[['INC_AD', 'T'], ['INC_ADV', 'T'], ['STOP_IF_MW', 'F']])
ap.add_argument('--path-map', nargs=2, action='append', default=[], metavar=('OLD', 'NEW'),
                help='replace the prefix OLD by NEW in a recorded input path that no longer exists (moved run directories)')
a = ap.parse_args()
trial = a.trial.resolve(); out = a.out.resolve(); out.parent.mkdir(parents=True, exist_ok=True)
env = json.loads((trial / 'run_manifest.json').read_text())['extra_env']
for k_, v_ in list(env.items()):
    for old, new in a.path_map:
        if isinstance(v_, str) and v_.startswith(old) and not Path(v_).exists():
            env[k_] = new + v_[len(old):]
nt = int(next(l.split()[0] for l in (a.x0 / 'MODEL').read_text().splitlines() if 'Total number of variables' in l))
nd = int(next(l.split()[0] for l in (a.x0 / 'MODEL').read_text().splitlines() if 'Number of depth points' in l))
ck.NT, ck.ND = nt, nd
R = ck.Runner(a.x0, a.exe, Path(env['QSTARS_REPLACE_FORCE']), False, a.threads, out.parent / (out.name + '.log'),
              [tuple(v) for v in a.vadat])
R.no_obs = 'QSTARS_NO_OBS_FORMAL' in env
extra = dict(kv.split('=', 1) for kv in a.env) or None
st = read_fortran(trial / 'QSTARS_STEP.bin', (nt, nd))
res = dict(trial=str(trial), exe=str(a.exe.resolve()), env=extra)
try:
    R.run(out, st, threads=a.threads, warm=Path(env['QSTARS_ALI_WARM']) if 'QSTARS_ALI_WARM' in env else None,
          two_phot_j=Path(env['QSTARS_TWO_PHOT_J']) if 'QSTARS_TWO_PHOT_J' in env else None, extra_env=extra)
    res['completed'] = True
except RuntimeError as e:
    res['completed'] = False; res['error'] = str(e)[:500]
for f in ('QSTARS_ALI_STATS', 'QSTARS_ALI_NEWTON'):
    p = out / f
    res[f] = p.read_text().splitlines() if p.exists() else None
(out.parent / (out.name + '.json')).write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps({k: (v if k != 'QSTARS_ALI_NEWTON' or v is None else v[-12:]) for k, v in res.items()}, indent=1))
