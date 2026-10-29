#!/usr/bin/env python3
"""Stage 7 (Report 23): register Ball amendment n from its draft, citing preserved failing runs.

    register_amendment.py N RUN_DIR [RUN_DIR ...]

Copies benchmarks/ball2011_track_plan_amendment_N.draft.json to ..._N.json (the draft stays), sets status to
'registered <UTC>', fills preserved_failing_run with the cited directories (relative to the repository root) and the
result field (relax_result or smoke_gate_result) with a summary read from each run's relax.json / smoke.json, changes
nothing else, then writes ..._N.registration.json (sha256 of the file, plan sha256, UTC time) and verifies the
amendment with the driver's own check."""
import datetime, hashlib, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'scripts'))
import validate_ball_track as v  # noqa: E402

n = int(sys.argv[1]); runs = [Path(r) for r in sys.argv[2:]]
draft = ROOT / f'benchmarks/ball2011_track_plan_amendment_{n}.draft.json'
out = ROOT / f'benchmarks/ball2011_track_plan_amendment_{n}.json'
if out.exists():
    raise SystemExit(f'{out} exists: already registered')
a = json.loads(draft.read_text())
now = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def summary(d):
    d = (ROOT / d) if not d.is_absolute() else d
    if (d / 'relax.json').exists():
        r = json.loads((d / 'relax.json').read_text())
        seeds = [dict(seed_radius_rsun=s.get('seed_radius_rsun'), status=s.get('status'), reason=s.get('reason'),
                      accepted=s.get('accepted_steps'), rejections=s.get('rejections')) for s in r.get('seeds', [])]
        return dict(run=str(d.relative_to(ROOT)), kind='relax', zones=r.get('zones'), status=r.get('status'),
                    classification=r.get('classification'),
                    seed_spread=max(r['uniqueness']['max_abs_dln'].values()) if isinstance(r.get('uniqueness'), dict) else None,
                    seeds=seeds)
    if (d / 'smoke.json').exists():
        r = json.loads((d / 'smoke.json').read_text())
        keep = ('status', 'passed', 'accepted_steps', 'median_dt_years', 'rejected_fraction', 'rejection_causes',
                'termination', 'reason', 'final_bh_mass_msun')
        return dict(run=str(d.relative_to(ROOT)), kind='smoke', **{k: r[k] for k in keep if k in r})
    raise SystemExit(f'{d}: no relax.json or smoke.json')


sums = [summary(r) for r in runs]
a['status'] = f'registered {now}'
a['preserved_failing_run'] = [s['run'] for s in sums]
key = 'smoke_gate_result' if 'smoke_gate_result' in a else 'relax_result'
a[key] = sums
out.write_text(json.dumps(a, indent=1) + '\n')
rec = dict(amendment=str(out.relative_to(ROOT)), amendment_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
           plan_sha256='a10b4166222342acee53da9622885845ff7744914fb38316d8ea4f48fce8cc0e', registered_utc=now)
out.with_name(out.stem + '.registration.json').write_text(json.dumps(rec, indent=2) + '\n')
print(json.dumps(rec, indent=1))
print('driver check:', v.verify_amendment(str(out), v.verify_registration()['plan_sha256'])[0])
