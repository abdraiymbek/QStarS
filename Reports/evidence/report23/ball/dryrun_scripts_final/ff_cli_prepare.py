"""DRYRUN: scratch copies of the four amendment drafts with registration records, for exercising the CLI flow only.
They are NOT registrations: they live in scratch, say so in their status, and cite the 2026-10-14 relax failures for every
amendment (a real Amendment 4, 2 or 1 must cite its own preserved failing run)."""
import json,hashlib,sys
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS');OUT=Path(sys.argv[1])
plan_sha=json.loads((ROOT/'benchmarks/ball2011_track_plan.registration.json').read_text())['plan_sha256']
cited=json.loads((ROOT/'benchmarks/ball2011_track_plan_amendment_3.draft.json').read_text())['preserved_failing_run']
for n in (3,4,2,1):
    d=json.loads((ROOT/f'benchmarks/ball2011_track_plan_amendment_{n}.draft.json').read_text())
    d['status']='registered (DRYRUN scratch copy for the CLI flow; not a registration)'
    d['preserved_failing_run']=cited
    p=OUT/f'ball2011_track_plan_amendment_{n}.json';p.write_text(json.dumps(d,indent=2)+'\n')
    rec=dict(amendment=str(p),amendment_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),plan_sha256=plan_sha,
             registered_utc='DRYRUN',note='scratch copy, not a registration')
    p.with_name(p.stem+'.registration.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(p,rec['amendment_sha256'][:12])
