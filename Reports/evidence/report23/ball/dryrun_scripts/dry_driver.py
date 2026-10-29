"""DRYRUN through the driver's own functions (relax_seeds, smoke) with a scratch output root.
Nothing is written under benchmarks/results; row-0 observables are never printed or compared with a table.
usage: dry_driver.py RUN [A3] [A2] [A1] [smoke] [F1]"""
from dry import *
import sys
run=sys.argv[1];flags=sys.argv[2:]
draft=lambda n,name=None:json.load(open(name or f'{ROOT}/benchmarks/ball2011_track_plan_amendment_{n}.draft.json'))['settings']
amendments=[]
if 'A1' in flags:amendments.append((dict(path='DRYRUN:amendment_1.draft'),draft(1)))
if 'A2' in flags:amendments.append((dict(path='DRYRUN:amendment_2.draft'),draft(2)))
if 'A3' in flags:amendments.append((dict(path='DRYRUN:amendment_3.draft'),draft(3,os.environ.get('A3_DRAFT'))))
parent='DRYRUN_driver';os.makedirs(parent,exist_ok=True)
registration=vbt.verify_registration()
config=vbt.build_config(PLAN,run,LIB,MESA_DIR,amendments);config['name']='DRYRUN_'+config['name']
eos,opacity=vbt.physics(config)
tag='+'.join(f for f in ('A3','A2','A1') if f in flags) or 'registered'
summary_path=f'{parent}/DRYRUN_summary.jsonl'
def note(record):
    with open(summary_path,'a') as f:f.write(json.dumps(quasistatic._json(record))+'\n')
    print('DRYRUN',json.dumps(quasistatic._json(record)),flush=True)
if 'smoke' not in flags:
    path=vbt.run_directory(parent,'DRYRUN_relax_'+tag,run,config,registration,[a for a in amendments])
    t=time.perf_counter()
    out=vbt.relax_seeds(PLAN,config,eos,opacity,path,config['replication']['seed_radii_rsun'],config['numerics']['zones'],log=lambda *a:None)
    wall=time.perf_counter()-t
    from qstars import io;io.finish_provenance(path)
    seeds=[dict(radius=e['seed_radius_rsun'],status=e['status'],reason=(e.get('reason') or '')[:90] if e['status']!='converged' else None,
                accepted=e.get('accepted_steps'),rejected=e.get('rejections'),cycles=len(e.get('mesh_cycles') or []),zones=e.get('row0_zones'),
                enforced=(e.get('row0_screens') or {}).get('enforced_screen_errors'),report_only=(e.get('row0_screens') or {}).get('report_only')) for e in out['seeds']]
    steps=sum((e.get('accepted_steps') or 0) for e in out['seeds'])
    note(dict(kind='relax',run=run,zones=config['numerics']['zones'],amendments=tag,status=out['status'],classification=out.get('classification'),
              seeds=seeds,wall_seconds=wall,seconds_per_accepted=wall/max(steps,1),
              uniqueness_max_abs_dln=max(out['uniqueness']['max_abs_dln'].values()) if 'uniqueness' in out else None,
              uniqueness_worst=max(out['uniqueness']['max_abs_dln'],key=out['uniqueness']['max_abs_dln'].get) if 'uniqueness' in out else None,
              directory=str(path)))
else:
    # the row-0 model of the latest DRYRUN relax made with the same amendments
    relaxes=sorted(d for d in os.listdir(parent) if d.endswith(f"DRYRUN_relax_{tag}_{run}_n{config['numerics']['zones']}"))
    row0=None
    for d in relaxes[::-1]:
        r=json.load(open(f'{parent}/{d}/relax.json'))
        if r['status']=='converged_unique':row0=f"{parent}/{d}/{r['row0_checkpoint']}";break
    fallbacks=['F1'] if 'F1' in flags else []
    # build_config inside vbt.smoke: same amendments, DRYRUN name is cosmetic (not in the signature)
    t=time.perf_counter()
    path,out=vbt.smoke(None,PLAN,registration,amendments,LIB,MESA_DIR,parent,row0,fallbacks=fallbacks)
    note(dict(kind='smoke',run=run,amendments=tag,fallbacks=fallbacks,row0=row0,wall_seconds=time.perf_counter()-t,
              **{k:out[k] for k in ('status','reason','accepted_steps','attempts','rejections','median_dt_years','rejected_fraction','passed','history_written','start_row_screen_errors')},
              directory=str(path)))
