"""No-merge DRYRUN final state (1667, M_BH 55.0): refine the front cells with a tiny-dt step, then compare Newton
iterations of steps at larger dt on the coarse vs refined mesh (unpatched code, in-memory remesh override only)."""
from v_common import *
from qstars.mesh import split_merge
ck=OUT/'DRYRUN_registered_controller_merge0_600/checkpoint.npz'
state,config=qs.load_checkpoint(ck)
eos,opacity=vbt.physics(config)
config["quasistatic"]["tolerance"]=float(sys.argv[1]) if len(sys.argv)>1 else config["quasistatic"].get("tolerance",1e-12)
print("DIAGNOSTIC tolerance",config["quasistatic"]["tolerance"])
orig=qs.remesh_base;mode=dict(spec=None)
def custom(old,cut,eos_,zones,split=2.,merge=.5,max_passes=40,protect=0):
    if mode['spec'] is None:return orig(old,cut,eos_,zones,split,merge,max_passes,protect)
    n=len(cut.dm);x=np.ones(n)
    for f,cells in mode['spec']:
        for c in cells:x[c]={2:2.5,4:4.5,8:9.}[f]
    st,info=split_merge(cut,eos_,cut.dm/x,split,merge)
    return st,dict({k:info[k] for k in ('splits','merges','energy_defect','entropy_change','merge_mixing_entropy','potential_reallocation','max_event_defect')},passes=1)
qs.remesh_base=custom
for factor in (4,8):
    mode['spec']=[(factor,range(154,164))]
    try:fine,d,m=attempt(state,1e-4,eos,opacity,config)
    except Exception as e:print('refine x',factor,'FAIL',str(e)[:150]);continue
    print(f'refine x{factor} at dt=1e-4 yr: zones {d["zones"]} it={d["newton_iterations"]} delta={m["delta_actual"]:.4g}',flush=True)
    mode['spec']=None
    for dt in (.5,1.,2.,5.,10.,30.):
        for name,st in (('coarse',state),('fine',fine)):
            try:
                new,d,m=attempt(st,dt,eos,opacity,config)
                print(f'  dt={dt:5.3g} {name:6s} ok it={d["newton_iterations"]:2d} delta={m["delta_actual"]:.4g} m/s={d["merges"]}/{d["splits"]}',flush=True)
            except Exception as e:print(f'  dt={dt:5.3g} {name:6s} FAIL {str(e)[:110]}',flush=True)
qs.remesh_base=orig
