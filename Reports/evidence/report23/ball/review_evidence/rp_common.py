"""DRYRUN helpers of the plan-compliance reviewer (scratch only; never writes under benchmarks/results)."""
import os,sys,json,time
for k in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
ROOT='/Users/abraiymbek/Projects/QStarS'
SRC=os.environ.get('RP_SRC',ROOT+'/src');SCRIPTS=os.environ.get('RP_SCRIPTS',ROOT+'/scripts')
sys.path.insert(0,SCRIPTS);sys.path.insert(0,SRC)
import numpy as np
import validate_ball_track as vbt
from qstars import quasistatic,photosphere,ledger,cavity
from qstars.constants import MSUN,RSUN,LSUN,YEAR,G,C,PI,A,SIGMA
LIB=ROOT+'/backends/mesa/build_ball_abi7/libqstars_mesa.dylib'
MESA_DIR=json.load(open(ROOT+'/backends/mesa/build_ball_abi7/build.json'))['mesa_dir']
PLAN=json.loads(open(ROOT+'/benchmarks/ball2011_track_plan.json').read())
HERE=os.path.dirname(os.path.abspath(__file__))
def draft(n):return json.load(open(f'{ROOT}/benchmarks/ball2011_track_plan_amendment_{n}.draft.json'))['settings']
def amend(flags):
    out=[]
    for n in (3,2,1):
        if f'A{n}' in flags:out.append((dict(path=f'DRYRUN:amendment_{n}.draft'),draft(n)))
    return out
