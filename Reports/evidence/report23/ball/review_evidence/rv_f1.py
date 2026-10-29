"""DRYRUN reviewer: fallback F1 (quasistatic.ptc) on MESA steps where the plain Newton solve fails, and bitwise identity where it converges."""
from rv import *
from qstars.numerics import StepRejected
import copy
for path,dts in ((sys.argv[1],[50.,30.,20.]),(sys.argv[2],[1.,50.,400.,1000.])):
    state,config=quasistatic.load_checkpoint(path)
    eos,opacity=vbt.physics(config)
    on=copy.deepcopy(config);on.setdefault('quasistatic',{})['ptc']=True
    print('DRYRUN',path,'base',quasistatic.settings(config)['quasistatic']['base'])
    for dt in dts:
        res={}
        for name,c in (('newton',config),('F1',on)):
            t=time.perf_counter()
            try:
                new,d=quasistatic.bondi_step(state,dt*YEAR,eos,opacity,c,mode='evolve')
                res[name]=(new,d);msg='accepted it=%d res=%.2e passes=%d dE=%.2e dM=%.2e ptc_entered=%s ptc_it=%s'%(d['newton_iterations'],d['nonlinear_residual'],d['passes'],d['energy_residual_change'],d['mass_residual_change'],d.get('ptc_entered'),d.get('ptc_iterations'))
            except StepRejected as e:msg='REJECTED: '+str(e)[:330]
            print('   dt=%g %s: %s  [%.1fs]'%(dt,name,msg,time.perf_counter()-t),flush=True)
        if len(res)==2:
            a,b=res['newton'][0],res['F1'][0]
            print('      bitwise equal states:',all(np.array_equal(getattr(a,f),getattr(b,f)) for f in ('r','T','dm','cell_volume','face_luminosity','radius_compensation')) and a.bh_mass==b.bh_mass and a.cavity_mass==b.cavity_mass and json.dumps(a.ledger,sort_keys=True)==json.dumps(b.ledger,sort_keys=True))
        elif 'F1' in res and 'newton' not in res:
            # F1 rescued a step: compare with the same interval taken in two Newton half-steps? report observables drift only
            new,d=res['F1'];print('      F1 RESCUED a step Newton failed: delta_fraction %.3g, map_dlnr %.3g, map_dlnT %.3g, Gamma_ph %.4f'%(d['delta_fraction'],d['mapped_max_dlnr'],d['mapped_max_dlnT'],d['photosphere_eddington_factor']))
