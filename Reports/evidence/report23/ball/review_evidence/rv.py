"""DRYRUN reviewer helpers (scratch only; never writes under benchmarks/results)."""
import os,sys,json,time
for k in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
ROOT='/Users/abraiymbek/Projects/QStarS'
sys.path.insert(0,ROOT+'/scripts');sys.path.insert(0,ROOT+'/src')
import numpy as np
import validate_ball_track as vbt
from qstars import quasistatic,photosphere,ledger,cavity
from qstars.constants import MSUN,RSUN,LSUN,YEAR,G,C,PI,A,SIGMA
LIB=ROOT+'/backends/mesa/build_ball_abi7/libqstars_mesa.dylib'
MESA_DIR=json.load(open(ROOT+'/backends/mesa/build_ball_abi7/build.json'))['mesa_dir']
PLAN=json.loads(vbt.PLAN.read_text())
def draft(n):return json.load(open(f'{ROOT}/benchmarks/ball2011_track_plan_amendment_{n}.draft.json'))['settings']
def amendments(flags):
    out=[]
    for n in (1,2,3):
        if f'A{n}' in flags:out.append((dict(path=f'DRYRUN:amendment_{n}.draft'),draft(n)))
    return out
def setup(run='C1',flags=(),extra=None):
    config=vbt.build_config(PLAN,run,LIB,MESA_DIR,amendments(flags))
    config['name']='DRYRUN_'+config['name']
    for section,values in (extra or {}).items():config[section].update(values)
    eos,opacity=vbt.physics(config)
    return config,eos,opacity
