"""DRYRUN relaxation of row 0 with library functions (scratch output only).
usage: dry_relax.py RUN RADIUS SURFACE [mesh_cycles] [tag]"""
from dry import *
import sys
run,radius,surface=sys.argv[1],float(sys.argv[2]),sys.argv[3]
cycles=int(sys.argv[4]) if len(sys.argv)>4 else 0
tag=sys.argv[5] if len(sys.argv)>5 else ''
amend=dict(replication=dict(seed_surface=surface))
if cycles:amend['relaxation']=dict(mesh_cycles=cycles)
if 'A1' in tag:amend['quasistatic']=dict(base='fixed_mass',base_mass_fraction=.05);amend['stepcontrol']=dict(metric='mapped')
config,eos,opacity=setup(run,amend)
N=config['numerics']['zones']
out=f'dryrun_relax_{run}_n{N}_{radius:g}rsun_{surface}_cyc{cycles}{tag}'
os.makedirs(out,exist_ok=True)
log=open(out+'/DRYRUN_relax.jsonl','w')
t0=time.perf_counter()
try:
    s,info=quasistatic.seed(config,eos,opacity,N,5*MSUN,1e4*MSUN,radius*RSUN)
except Exception as e:
    print('DRYRUN',out,'seed_failed',type(e).__name__,e);json.dump(dict(DRYRUN=True,status='seed_failed',reason=f'{type(e).__name__}: {e}'),open(out+'/DRYRUN_summary.json','w'));sys.exit(0)
print('DRYRUN',out,'seed built %.1fs'%(time.perf_counter()-t0),{k:v for k,v in info.items() if 'surface' in k or k=='law_luminosity'},flush=True)
ledger.initialize(s,eos,cavity=True)
clock=[time.perf_counter()]
def note(kind,data):
    now=time.perf_counter();data=dict(data,wall=now-clock[0]);clock[0]=now
    log.write(json.dumps(data)+'\n');log.flush()
    if kind=='accepted':
        print('DRYRUN acc cyc%d dt=%.3g it=%d N=%d chg=%.2e bal=%.2e G=%.4f %.1fs'%(data['cycle'],data['dt_years'],data['newton_iterations'],data['zones'],data['change'] if data['change'] is not None else float('inf'),data['luminosity_balance'],data.get('photosphere_eddington_factor') or -1,data['wall']),flush=True)
    elif kind=='rejected':print('DRYRUN rej dt=%.3g %s'%(data['dt_years'],data['message'][:150]),flush=True)
    else:print('DRYRUN cycle',data,flush=True)
t1=time.perf_counter()
state,record=quasistatic.relax(s,eos,opacity,config,config['relaxation'],log=note)
wall=time.perf_counter()-t1
summary=dict(DRYRUN=True,run=run,zones=N,seed_radius_rsun=radius,surface=surface,mesh_cycles=cycles,status=record['status'],reason=record['reason'],
             accepted=len(record['steps']),rejected=len(record['rejections']),cycles=len(record['cycles']),wall_seconds=wall,
             seconds_per_accepted=wall/max(len(record['steps']),1))
quasistatic.save_checkpoint(out+'/DRYRUN_relaxed.npz',state,config)
if record['status']=='converged':
    row0=quasistatic.zero_time(state,eos)
    quasistatic.save_checkpoint(out+'/DRYRUN_row0.npz',row0,config)
    summary['row0']=vbt.row0_observables(row0,eos,opacity,config)
    summary['zones_final']=len(row0.dm)
json.dump(quasistatic._json(summary),open(out+'/DRYRUN_summary.json','w'),indent=1)
print('DRYRUN summary',{k:v for k,v in summary.items() if k!='row0'},flush=True)
