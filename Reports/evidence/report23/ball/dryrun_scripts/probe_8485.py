"""DRYRUN: why the 8485 Rsun grey seed cannot take a relaxation step (no table comparison)."""
from dry import *
from qstars.mesh import q_targets,split_merge
import traceback
config,eos,opacity=setup('C1',dict(replication=dict(seed_surface='grey')))
for R in (8485.,7000.,6000.):
    s,info=quasistatic.seed(config,eos,opacity,200,5*MSUN,1e4*MSUN,R*RSUN)
    print('DRYRUN R',R,'Teff',info['surface_effective_temperature'],'T_last',s.T[-3:],'polytrope T_last',info['surface_polytrope_temperature'],'rho_last',s.rho[-3:])
    ledger.initialize(s,eos,cavity=True)
    try:
        goal=quasistatic.targets(s,s,eos,200)
        ratio=s.dm/goal
        print('   dm/target last cells',ratio[-5:],'max',ratio.max(),'argmax',ratio.argmax(),'min',ratio.min())
        st,inf=split_merge(s,eos,goal,2.,.5)
        print('   split_merge ok',inf['splits'],inf['merges'])
    except Exception as e:
        print('   split_merge error:',type(e).__name__,e)
        tb=traceback.extract_tb(e.__traceback__);print('   at',[(f.name,f.lineno) for f in tb][-4:])
