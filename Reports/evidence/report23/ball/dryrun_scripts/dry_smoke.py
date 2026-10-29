"""DRYRUN smoke-like probe (scratch output only; at most a few hundred accepted steps).
usage: dry_smoke.py ROW0_DIR TAG MAX_STEPS [ptc] """
from dry import *
import sys,csv,statistics
row0_dir,tag,max_steps=sys.argv[1],sys.argv[2],int(sys.argv[3])
ptc='ptc' in sys.argv[4:]
state,config=quasistatic.load_checkpoint(row0_dir if row0_dir.endswith('.npz') else row0_dir+'/DRYRUN_row0.npz')
config['quasistatic']['ptc']=ptc
eos,opacity=vbt.physics(config)
out=f'dryrun_smoke_{tag}';os.makedirs(out,exist_ok=True)
for f in os.listdir(out):os.remove(os.path.join(out,f))
options=vbt.track_options(PLAN,config)
t=time.perf_counter()
_,report=quasistatic.run(state,eos,opacity,config,out,end_bh_mass=options['landings'][0],stop_after_steps=max_steps,
                         **{**options,'landings':options['landings'][:1],'max_wall_seconds':1500.})
wall=time.perf_counter()-t
rows=list(csv.DictReader(open(out+'/history.csv')))
dts=[float(r['step_dt']) for r in rows if r['event'] in ('step','landing')]
ev=[json.loads(x) for x in open(out+'/events.jsonl')]
rej=[e for e in ev if e['event']=='rejected'];acc=[e for e in ev if e['event']=='accepted']
kinds={}
for e in rej:
    key=e['message'].split(':')[0][:70] if 'delta_actual' not in e['message'] else 'delta_actual > 2 delta_target ('+e['message'].split('(')[-1]
    kinds[key]=kinds.get(key,0)+1
summary=dict(DRYRUN=True,tag=tag,ptc=ptc,status=report['status'],reason=report['reason'],accepted=len(dts),attempts=report['attempts'],
             rejections=report['rejections'],rejected_fraction=report['rejections']/max(report['attempts'],1),
             median_dt_years=statistics.median(dts) if dts else None,min_dt_years=min(dts,default=None),max_dt_years=max(dts,default=None),
             last_dt_years=dts[-1] if dts else None,final_bh_mass_msun=report['final_bh_mass_msun'],wall_seconds=wall,
             seconds_per_accepted=wall/max(len(dts),1),seconds_per_attempt=wall/max(report['attempts'],1),rejection_kinds=kinds,
             newton_iterations_median=statistics.median([e['newton_iterations'] for e in acc]) if acc else None,
             ptc_entered_accepted=sum(1 for e in acc if e.get('ptc_entered')),
             ptc_in_rejections=sum(1 for e in rej if 'F1 pseudo-transient' in e['message']),
             limiting_metric={},start_row_screen_errors=report.get('start_row_screen_errors'))
for e in acc:
    m=e['metrics'];k=max((v,k) for k,v in m.items() if k!='delta_actual')[1];summary['limiting_metric'][k]=summary['limiting_metric'].get(k,0)+1
json.dump(quasistatic._json(summary),open(out+'/DRYRUN_summary.json','w'),indent=1)
print('DRYRUN smoke summary',json.dumps(quasistatic._json(summary),indent=1))
