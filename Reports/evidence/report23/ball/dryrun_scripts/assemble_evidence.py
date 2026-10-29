"""Assemble the DRYRUN evidence (no row-0 observable values, no table comparison) into
Reports/evidence/report23/ball/ and fill dry_run_evidence of the Amendment 3 draft."""
import json,os,shutil,glob,math
ROOT='/Users/abraiymbek/Projects/QStarS'
OUT=ROOT+'/Reports/evidence/report23/ball';os.makedirs(OUT+'/dryrun_scripts',exist_ok=True)
driver=[json.loads(l) for l in open('DRYRUN_driver/DRYRUN_summary.jsonl')]
for r in driver:r['directory']='scratch:'+r['directory'];r.pop('row0',None)
def summary(d):
    p=d+'/DRYRUN_summary.json'
    if not os.path.exists(p):return None
    s=json.load(open(p));s.pop('row0',None);s['directory']='scratch:'+d;return s
library=[s for s in (summary(d) for d in sorted(glob.glob('dryrun_relax_*'))) if s]
def spread(a,b):
    x,y=(json.load(open(d+'/DRYRUN_summary.json')) for d in (a,b))
    if 'row0' not in x or 'row0' not in y:return None
    d={k:abs(math.log(x['row0'][k]/y['row0'][k])) for k in x['row0']}
    return dict(seeds=[a,b],zones_final=[x.get('zones_final'),y.get('zones_final')],max_abs_dln=max(d.values()),per_observable=d)
spreads=[]
for run,n in (('C1',200),('C0',400),('C2',800)):
    for tag in ('cyc0','cyc40','cyc40_A1'):
        a=f'dryrun_relax_{run}_n{n}_6000rsun_grey_{tag}';b=f'dryrun_relax_{run}_n{n}_3464rsun_grey_{tag}'
        if os.path.exists(a+'/DRYRUN_summary.json') and os.path.exists(b+'/DRYRUN_summary.json'):
            s=spread(a,b)
            if s:spreads.append(dict(run=run,zones=n,settings=tag,**s))
smokes=[s for s in (summary(d) for d in sorted(glob.glob('dryrun_smoke_*'))) if s]
evidence=dict(
    label='DRYRUN evidence for Stage 7 unblocking (Gate 6, Ball et al. 2011). Scratch runs made by calling library and driver functions with a scratch output root; nothing was written under benchmarks/results. No row-0 observable value is recorded here and none was compared with Table 1 or thesis Table 3.1; only seed-to-seed differences, step counts, timestep statistics and wall times.',
    made_utc='2026-10-14',machine_note='10 cores shared with 4 CMFGEN jobs (load average 9-18); single thread, nice 15',
    scratch_root='/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/ball_unblock (purgeable)',
    driver_path=dict(note='scripts/validate_ball_track.py relax_seeds and smoke called directly with amendment settings taken from the drafts (A3: seed surface grey + appended radii; A2: mesh cycles; A1: fixed-mass base + mapped metric), final source of this workflow',runs=driver),
    library_path=dict(note='quasistatic.seed + quasistatic.relax called directly, one seed per process (made while the source was being edited; numerics of seed and relaxation unchanged since)',runs=library,seed_spreads=spreads),
    smoke_like_probes=dict(note='quasistatic.run from a DRYRUN row-0 model with the plan track options, stopped after at most 300 accepted steps',runs=smokes))
json.dump(evidence,open(OUT+'/stage7_unblock_dryruns.json','w'),indent=1)
for f in ('dry.py','dry_relax.py','dry_smoke.py','dry_driver.py','spread.py','exp_seed.py','probe_ph.py','probe_ph2.py','probe_seed.py','probe_screens.py','probe_8485.py','toy_f1.py','toy_f1b.py','toy_f1c.py','batch_relax.sh','batch_driver.sh','assemble_evidence.py'):
    if os.path.exists(f):shutil.copy(f,OUT+'/dryrun_scripts/'+f)
print('evidence written',len(driver),len(library),len(spreads),len(smokes))
