"""fold.classify (registered code, unpatched) on the last accepted state of the smoke run (step 1255). Scratch only."""
from common import *
import time
from qstars import fold
state,config,eos,opacity=load()
t=time.time()
def log(direction,run):print('direction',direction,'stop:',run['stop'],'points',len(run['points']),'fold',run['fold'] is not None,flush=True)
out=fold.classify(state,eos,opacity,config,log=log)
out['wall_seconds']=time.time()-t
slim={k:v for k,v in out.items() if k!='runs'}
slim['runs']=[dict(direction=r['direction'],limit=r['limit'],stop=r['stop'],points=len(r['points']),
                   bh_mass_range_msun=[min(p['bh_mass'] for p in r['points'])/1.98847e33,max(p['bh_mass'] for p in r['points'])/1.98847e33] if r['points'] else None,
                   sigma_min_range=[min(p.get('sigma_min',np.inf) for p in r['points']),max(p.get('sigma_min',0) for p in r['points'])] if r['points'] else None,
                   fold=r['fold']) for r in out['runs']]
json.dump(qs._json(slim),open(SCR/'fold_classify_step1255.json','w'),indent=1)
print(json.dumps(qs._json({k:slim[k] for k in ('label','reason','fold_confirmed','start_residual','wall_seconds')}),indent=1))
