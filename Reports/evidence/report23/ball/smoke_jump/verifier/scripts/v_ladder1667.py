"""Diagnostic: refined (x4 front cells) vs coarse state at no-merge step 1667, each laddered through small dt to
re-equilibrate, then single steps at larger dt. Diagnostic tolerance (argv[1]) because refined meshes stall at a
~1e-11..1e-10 roundoff floor above the registered 1e-12."""
exec(open('v_refine1667.py').read().split("for factor in (4,8):")[0])
ladder=(1e-4,1e-3,1e-2,.1,.3)
mode['spec']=[(4,range(154,164))]
fine,d,m=attempt(state,ladder[0],eos,opacity,config);print('refined zones',d['zones'],'jump',m['delta_actual'])
mode['spec']=None
coarse,d,m=attempt(state,ladder[0],eos,opacity,config)
for dt in ladder[1:]:
    fine,df,mf=attempt(fine,dt,eos,opacity,config);coarse,dc,mc=attempt(coarse,dt,eos,opacity,config)
    print(f'ladder dt={dt:g}: fine it={df["newton_iterations"]} delta={mf["delta_actual"]:.3g} | coarse it={dc["newton_iterations"]} delta={mc["delta_actual"]:.3g}',flush=True)
for dt in (.5,1.,2.,5.,10.):
    for name,st in (('coarse',coarse),('fine',fine)):
        try:
            new,d,m=attempt(st,dt,eos,opacity,config)
            print(f'  dt={dt:5.3g} {name:6s} ok it={d["newton_iterations"]:2d} delta={m["delta_actual"]:.4g}',flush=True)
        except Exception as e:print(f'  dt={dt:5.3g} {name:6s} FAIL {str(e)[:110]}',flush=True)
