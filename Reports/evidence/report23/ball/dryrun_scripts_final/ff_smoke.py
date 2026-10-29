"""DRYRUN: the driver's smoke() (registered smoke gate) from a dry-run row 0, scratch output. usage: ff_smoke.py ROW0 FLAGS [F1]"""
from ff_common import *
row0,flags=sys.argv[1],sys.argv[2];fallbacks=['F1'] if 'F1' in sys.argv[3:] else []
registration=vbt.verify_registration();os.makedirs(OUT,exist_ok=True)
vbt.import_ball_modules()
t=time.time()
path,out=vbt.smoke(None,PLAN,registration,amend(flags),LIB,MESA_DIR,OUT,row0,fallbacks=fallbacks)
print(f'DRYRUN smoke [{flags}] fallbacks {fallbacks} wall {time.time()-t:.1f} s -> {path.name}')
print('   ',{k:out[k] for k in ('status','reason','accepted_steps','attempts','rejections','median_dt_years','passed','rejection_causes')})
closure=json.load(open(path/'code_closure.json'));print('    code_closure loaded_outside',closure['loaded_outside_closure'],'sha',closure['sha256'][:16])
if 'fold' in sys.argv[3:]:
    from qstars import fold
    state,config=quasistatic.load_checkpoint(path/'checkpoint.npz');eos,opacity=vbt.physics(config)
    t=time.time();f=fold.classify(state,eos,opacity,config)
    print(f'DRYRUN fold.classify on the smoke end state: {f["label"]} | {f["reason"][:120]} | {time.time()-t:.1f} s')
