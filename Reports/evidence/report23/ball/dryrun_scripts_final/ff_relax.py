"""DRYRUN: the driver's relax_seeds with a scratch output root. usage: ff_relax.py RUN FLAGS (e.g. C1 A3,A4,A2)."""
from ff_common import *
run,flags=sys.argv[1],(sys.argv[2] if len(sys.argv)>2 else '')
registration=vbt.verify_registration();os.makedirs(OUT,exist_ok=True);vbt.import_ball_modules()
amendments=amend(flags)
config=vbt.build_config(PLAN,run,LIB,MESA_DIR,amendments);config['name']='DRYRUN_'+config['name']
eos,opacity=vbt.physics(config)
label=flags.replace(',','+') or 'registered'
path=vbt.run_directory(os.path.join(OUT,f'DRYRUN_relax_{label}'),'relax',run,config,registration,amendments)  # DRYRUN parent; kind relax (row0_record)
t=time.time()
try:out=vbt.relax_seeds(PLAN,config,eos,opacity,path,config['replication']['seed_radii_rsun'],config['numerics']['zones'],log=lambda *a:None)
finally:vbt.finish(path)
print(f'DRYRUN relax {run} [{label}] photosphere_bracket={config["quasistatic"]["photosphere_bracket"]} wall {time.time()-t:.1f} s')
print(summary(out),flush=True)
