"""DRYRUN: does fold.classify run on a MESA state at all, and what does it cost? (row-0 state, N = 200; no table comparison)"""
from dry import *
from qstars import fold
import sys
state,config=quasistatic.load_checkpoint(sys.argv[1])
eos,opacity=vbt.physics(config)
t=time.perf_counter()
try:
    out=fold.classify(state,eos,opacity,config)
    print('DRYRUN fold.classify: label',out['label'],'|',out['reason'][:300],'| start_residual %.3g'%out['start_residual'],'| points',[len(r['points']) for r in out['runs']],'| stops',[r['stop'] for r in out['runs']],'| wall %.1fs'%(time.perf_counter()-t))
except Exception as e:
    import traceback;traceback.print_exc()
    print('DRYRUN fold.classify raised',type(e).__name__,e,'wall %.1fs'%(time.perf_counter()-t))
