from common import *
from qstars import fold
state,config,eos,opacity=load()
out=fold.classify(state,eos,opacity,config)
pts=[dict(direction=r['direction'],**{k:p[k] for k in p if k in ('lam','bh_mass','t_lambda','iterations','sigma_min','det_sign','min_pivot','min_singular')}) for r in out['runs'] for p in r['points']]
for p in pts:p['bh_mass_msun']=p.pop('bh_mass')/1.98847e33
for p in pts:print(p)
full=json.load(open(SCR/'fold_classify_step1255.json'));full['points']=pts
json.dump(qs._json(full),open(SCR/'fold_classify_step1255.json','w'),indent=1)
