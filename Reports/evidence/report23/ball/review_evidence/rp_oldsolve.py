"""DRYRUN: does Amendment 3 depend on the unconditional photosphere.solve change?
Grey seeds built and stepped once with (a) the current solve, (b) the solve of the preserved source snapshot."""
from rp_common import *
import importlib.util
from qstars.numerics import StepRejected
spec=importlib.util.spec_from_file_location('old_photosphere',HERE+'/snap/qstars/photosphere.py',submodule_search_locations=None)
import qstars
old=importlib.util.module_from_spec(spec);old.__package__='qstars';spec.loader.exec_module(old)
new_solve=photosphere.solve
calls=[]
def traced(which):
    f=old.solve if which=='old' else new_solve
    def g(*a,**k):
        try:
            out=f(*a,**k);calls.append((which,'ok'));return out
        except Exception as e:
            calls.append((which,type(e).__name__+': '+str(e)[:80]));raise
    return g
for N,run in ((200,'C1'),(400,'C0'),(800,'C2')):
    config=vbt.build_config(PLAN,run,LIB,MESA_DIR,amend(['A3']));config['name']='DRYRUN_'+config['name']
    eos,opacity=vbt.physics(config)
    for R in config['replication']['seed_radii_rsun']:
        res={}
        for which in ('new','old'):
            photosphere.solve=traced(which);calls.clear()
            try:
                s,info=quasistatic.seed(config,eos,opacity,N,5*MSUN,1e4*MSUN,R*RSUN)
            except Exception as e:
                res[which]=('seed_failed',type(e).__name__,str(e)[:70]);continue
            ledger.initialize(s,eos,cavity=True)
            try:
                new,d=quasistatic.bondi_step(s,.01*YEAR,eos,opacity,config,mode='relax')
                res[which]=('step_ok',s.T.tobytes(),new.T.tobytes(),new.r.tobytes(),list(calls))
            except StepRejected as e:
                res[which]=('step_rejected',str(e)[:90],list(calls))
        same=res['new'][0]==res['old'][0] and (res['new'][0]!='step_ok' or res['new'][1:4]==res['old'][1:4])
        print('DRYRUN N',N,'R',R,'new:',res['new'][0],res['new'][-1] if res['new'][0]!='step_ok' else res['new'][4],'| old:',res['old'][0],res['old'][-1] if res['old'][0]!='step_ok' else res['old'][4],'| same outcome and bitwise:',same,flush=True)
photosphere.solve=new_solve
