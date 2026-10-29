"""DRYRUN: reproduce the registered seed blocker (no relaxation output under benchmarks/results)."""
from dry import *
from qstars.transport import luminosities
import copy
config,eos,opacity=setup('C1')
M=1e4*MSUN
for R in (2000,6000,1000,12000):
    t=time.perf_counter()
    try:
        s,info=quasistatic.seed(config,eos,opacity,200,5*MSUN,M,R*RSUN)
    except Exception as e:
        print('DRYRUN seed',R,'FAILED',type(e).__name__,str(e)[:100]);continue
    bc=copy.deepcopy(config);bc['transport']['convection']='off'
    L,_=luminosities(s,eos,opacity,bc,0.)
    kap=opacity(s.rho,s.T)
    print('DRYRUN seed',R,'built in %.1fs'%(time.perf_counter()-t),'law L/Lsun %.3e'%(info['law_luminosity']/LSUN),'T_last',s.T[-3:],'rho_last',s.rho[-3:],
          'L_N grey/Lsun %.3e'%(L[-1]/LSUN),'ratio',L[-1]/info['law_luminosity'],'kap_last',kap[-3:],'dm_last/Msun',s.dm[-3:]/MSUN,'pext',info['external_pressure'])
    try:
        rho=photosphere.solve(float(L[-1]),R*RSUN,M,eos,opacity)
        print('   photosphere at L_N exists',rho)
    except Exception as e:print('   photosphere at grey L_N:',type(e).__name__,e)
    ledger.initialize(s,eos,cavity=True)
    try:
        new,d=quasistatic.bondi_step(s,.01*YEAR,eos,opacity,config,mode='relax')
        print('   first relax step accepted')
    except Exception as e:print('   first relax step:',type(e).__name__,str(e)[:200])
