"""DRYRUN: what Amendment 3's settings give when photosphere.solve is the one of the preserved (registered) source."""
from rp_common import *
import importlib.util
spec=importlib.util.spec_from_file_location('old_photosphere',HERE+'/snap/qstars/photosphere.py')
old=importlib.util.module_from_spec(spec);old.__package__='qstars';spec.loader.exec_module(old)
photosphere.solve=old.solve
registration=vbt.verify_registration();parent=HERE+'/DRYRUN_oldsolve';os.makedirs(parent,exist_ok=True)
for run in sys.argv[1:]:
    A3=amend(['A3']);config=vbt.build_config(PLAN,run,LIB,MESA_DIR,A3);config['name']='DRYRUN_'+config['name']
    eos,opacity=vbt.physics(config)
    path=vbt.run_directory(parent,'DRYRUN_relax_A3_oldsolve',run,config,registration,A3)
    out=vbt.relax_seeds(PLAN,config,eos,opacity,path,config['replication']['seed_radii_rsun'],config['numerics']['zones'],log=lambda *a:None)
    print('DRYRUN A3 settings + registered photosphere.solve,',run,'->',out['status'],out.get('classification'))
    for e in out['seeds']:print('   seed',e['seed_radius_rsun'],e['status'],(e.get('reason') or '')[:95],'acc',e.get('accepted_steps'),'rej',e.get('rejections'))
