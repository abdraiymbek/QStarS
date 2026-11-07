"""Goal 3 DIAGNOSTIC DRYRUN (not a registered run, not a registration change): the registered smoke track (C1, N = 200,
the stored configuration of the preserved smoke run, amendments 3, 4, 2, 1; registered controller) from the registered
row 0, with ONE in-memory change: the P2 split/merge targets are reduced where the opacity of the accepted state
changes steeply between neighbouring cells (front refinement). No file under src/ or benchmarks/ is touched; the
change is a monkeypatch of qstars.quasistatic.targets in this process only.

Front refinement (accepted state `old`, cells k): g_k = max |d ln kappa| over the two faces of cell k, smoothed by a
running max over +-W cells (so cells ahead of a moving front are refined before it arrives); the STARS-Q target of
cell k is multiplied by clip(C / g_k, 1/FMAX, 1), then mapped onto the sliced mesh exactly as targets() maps ln dm*.
The registered split (2) and merge (0.5) thresholds act unchanged on these targets, so every split and merge is still
an event measured by the registered controller.

usage: front_refine_dryrun.py NAME C W FMAX [TOL] [MAX_WALL_S] [START]
  C = 0 switches the refinement off (control: must reproduce the registered smoke run).
  TOL overrides quasistatic.tolerance (diagnostic only; default: the registered 1e-12).
  START: 'row0' (default) or a checkpoint path to continue from.
"""
from bm_common import *
import shutil
name=sys.argv[1];CF=float(sys.argv[2]);W=int(sys.argv[3]);FMAX=float(sys.argv[4])
TOL=float(sys.argv[5]) if len(sys.argv)>5 and sys.argv[5]!='-' else None
WALL=float(sys.argv[6]) if len(sys.argv)>6 else 2400.
START=sys.argv[7] if len(sys.argv)>7 else 'row0'
out=BIG/name
if out.exists():shutil.rmtree(out)
out.mkdir(parents=True)
_,config=qs.load_checkpoint(RUN/'checkpoint.npz')          # the registered smoke configuration, as run
ZONES=int(os.environ.get('BM_ZONES','200'))      # 400: the same track at the C0 resolution (row 0 of relax C0)
config['numerics']['zones']=ZONES
if ZONES!=200:config['name']=f'ball2011_C0_n{ZONES}';config['replication'].update(run='C0',plan_text=f'N = {ZONES}')
registered=copy.deepcopy(config);config=copy.deepcopy(config)
if TOL is not None:config['quasistatic']['tolerance']=TOL
if START=='row0':
    state,stored=qs.load_checkpoint(ROW0[ZONES])
    assert state.steps==0 and state.time==0 and vbt.signature(stored)==vbt.signature(registered),'row 0 does not match the smoke config'
else:
    state,stored=qs.load_checkpoint(START)
eos,opacity=vbt.physics(config)
orig_targets=qs.targets
stats=dict(calls=0,max_reduction=1.,reduced_cells_max=0)
def front_factor(old):
    kap=np.log(opacity(old.rho,old.T));gf=abs(np.diff(kap))
    g=np.maximum(np.r_[0.,gf],np.r_[gf,0.])
    if W>0:
        pad=np.r_[np.zeros(W),g,np.zeros(W)]
        g=np.max(np.lib.stride_tricks.sliding_window_view(pad,2*W+1),axis=1)
    return np.clip(CF/np.maximum(g,1e-300),1./FMAX,1.)
def targets_front(old,cut,eos_,zones):
    lnt=np.log(qs.q_targets(old,eos_,zones))+np.log(front_factor(old))
    centre_old=np.cumsum(old.dm)-.5*old.dm;centre=np.cumsum(cut.dm)-.5*cut.dm
    total_old,total=old.dm.sum(),cut.dm.sum()
    inner=np.interp(centre,centre_old,lnt)
    outer=np.interp((total-centre)[::-1],(total_old-centre_old)[::-1],lnt[::-1])[::-1]
    return np.exp(np.where(centre<=.5*total,inner,outer))
if CF>0:qs.targets=targets_front
plan=PLAN
options=vbt.track_options(plan,config)
options['max_wall_seconds']=WALL
meta=dict(name=name,front_C=CF,window=W,fmax=FMAX,newton_tolerance=config['quasistatic']['tolerance'],start=str(START),
          row0=str(ROW0.get(ZONES)),zones=ZONES,config_from=str(RUN/'checkpoint.npz'),max_wall_seconds=WALL,
          note='DIAGNOSTIC DRYRUN: not a registered run; registered controller; in-memory front refinement only')
(out/'DRYRUN_meta.json').write_text(json.dumps(meta,indent=1))
t0=time.time()
_,report=qs.run(state,eos,opacity,config,out,end_bh_mass=options['landings'][0],**{**options,'landings':options['landings'][:1]})
report['wall_seconds']=time.time()-t0
(out/'DRYRUN_report.json').write_text(json.dumps(qs._json(report),indent=1))
print(json.dumps({k:report.get(k) for k in ('status','reason','accepted_steps_this_run','attempts','rejections','final_step',
                                            'final_bh_mass_msun','final_time_years','final_dt_years','wall_seconds')},indent=1))
