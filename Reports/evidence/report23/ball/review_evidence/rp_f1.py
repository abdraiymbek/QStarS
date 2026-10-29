"""DRYRUN: fallback F1 on the MESA Ball step (one evolve step from the dry-run row 0; nothing accepted or written)."""
from rp_common import *
import glob,copy
from qstars.numerics import StepRejected
cu=sorted(glob.glob(HERE+'/DRYRUN_out/*DRYRUN_relax_A3+A2_C1_n200'))[-1]
state,stored=quasistatic.load_checkpoint(cu+'/row0_6000rsun.npz');eos,opacity=vbt.physics(stored)
from qstars import stepcontrol
off=copy.deepcopy(stored);on=copy.deepcopy(stored);on['quasistatic']['ptc']=True
assert off['quasistatic']['ptc'] is False
fields=('r','T','dm','v','cell_volume','face_luminosity','radius_compensation')
for dt in (0.125,1.0):
    a,da=quasistatic.bondi_step(state,dt*YEAR,eos,opacity,off,mode='evolve')
    b,db=quasistatic.bondi_step(state,dt*YEAR,eos,opacity,on,mode='evolve')
    same=all(np.array_equal(getattr(a,f),getattr(b,f)) for f in fields) and (a.bh_mass,a.cavity_mass,a.photosphere_density)==(b.bh_mass,b.cavity_mass,b.photosphere_density) and json.dumps(a.ledger,sort_keys=True)==json.dumps(b.ledger,sort_keys=True)
    print('DRYRUN dt',dt,'yr: Newton converges, ptc off vs on bitwise equal:',same,'| residual',da['nonlinear_residual'],'| extra keys with ptc on:',sorted(set(db)-set(da)),'| differing diag keys:',[k for k in da if da[k]!=db[k] and not (da[k]!=da[k])])
# force the first Newton solve to report failure so that F1 is entered on the real equations
newton=quasistatic.newton;calls=[]
def first_fails(fun,x,*a,**k):
    calls.append(1)
    if len(calls)>1:return newton(fun,x,*a,**k)
    _,info=newton(fun,x,*a,**k)
    return x,dict(info,success=False,message='forced failure (review)')
for dt in (0.125,1.0):
    ref,dr=quasistatic.bondi_step(state,dt*YEAR,eos,opacity,off,mode='evolve')
    quasistatic.newton=first_fails;calls.clear()
    try:
        new,d=quasistatic.bondi_step(state,dt*YEAR,eos,opacity,on,mode='evolve')
        print('DRYRUN dt',dt,'forced Newton failure -> F1 entered',d.get('ptc_entered'),'handover',d.get('ptc_handover'),'ptc iterations',d.get('ptc_iterations'),'final residual',d['nonlinear_residual'],
              '| vs plain Newton solution max|dln r|',float(np.max(abs(np.log(new.r/ref.r)))),'max|dln T|',float(np.max(abs(np.log(new.T/ref.T)))),'| newton_iterations reported',d['newton_iterations'],'(plain:',dr['newton_iterations'],')',
              '| step mass/energy residual change',d['mass_residual_change'],d['energy_residual_change'])
    except StepRejected as e:print('DRYRUN dt',dt,'forced Newton failure -> rejected:',str(e)[:260])
    finally:quasistatic.newton=newton
    calls.clear();quasistatic.newton=first_fails
    try:quasistatic.bondi_step(state,dt*YEAR,eos,opacity,off,mode='evolve');print('   ptc off: accepted?!')
    except StepRejected as e:print('   ptc off, same forced failure -> rejected as registered:',str(e)[:100])
    finally:quasistatic.newton=newton
