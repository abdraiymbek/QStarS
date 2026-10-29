"""DRYRUN: smoke-like probe, quasistatic.run from a dry-run row 0 with stop_after_steps (scratch). usage: ff_probe.py ROW0 STEPS"""
from ff_common import *
row0,steps=sys.argv[1],int(sys.argv[2])
state,config=quasistatic.load_checkpoint(row0);eos,opacity=vbt.physics(config)
path=os.path.join(OUT,'DRYRUN_probe_'+os.path.basename(os.path.dirname(row0))+f'_{steps}');os.makedirs(path,exist_ok=True)
options=vbt.track_options(PLAN,config)
t=time.time()
state,report=quasistatic.run(state,eos,opacity,config,path,stop_after_steps=steps,**{**options,'landings':options['landings'][:1]},
                             end_bh_mass=options['landings'][0])
wall=time.time()-t
rows=vbt.read_history(path);acc=[r for r in rows if r['event'] in ('step','landing')]
dts=[float(r['step_dt']) for r in acc];res=[float(r['step_nonlinear_residual']) for r in acc]
print(f'DRYRUN probe N={len(state.dm)} {report["status"]} accepted {len(acc)} attempts {report["attempts"]} rejections {report["rejections"]} '
      f'wall {wall:.1f} s = {wall/max(len(acc),1):.3f} s/accepted; median dt {np.median(dts):.4g} yr; M_BH {state.bh_mass/MSUN:.4g}; '
      f'residual median {np.median(res):.2e} max {np.max(res):.2e}')
kn=[float(r['screen_radiation_knudsen_deep']) for r in rows if r.get('screen_radiation_knudsen_deep') not in ('',None)]
inv=[float(r['max_density_inversion']) for r in rows]
print(f'    deep Knudsen max {max(kn):.4f} (limit 0.1); row-0 max_density_inversion {inv[0]:.3e}; causes {vbt.rejection_causes(path)}')
