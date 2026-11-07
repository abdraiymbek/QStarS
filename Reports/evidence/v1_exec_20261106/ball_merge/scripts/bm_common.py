"""Shared setup for the 46 Msun merge-event diagnostics (V1 plan P5 bullet 2).

Diagnostic only: reads the preserved registered smoke run, never writes under benchmarks/.
Large outputs go to results/ball_merge_20261107/, small ones next to the scripts (../data).
Adapted from Reports/evidence/report23/ball/smoke_jump/verifier/scripts/v_common.py.
"""
import os
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ.setdefault(_n,'1')
import sys,copy,json,math,time
from pathlib import Path
import numpy as np
ROOT=Path('/Users/abraiymbek/Projects/QStarS')
SRC=Path(os.environ.get('BM_SRC',str(ROOT)))
sys.path.insert(0,str(SRC/'src'));sys.path.insert(0,str(SRC/'scripts'))
RUNS=ROOT/'benchmarks/results/2026-10-14-ball-track'
RUN=RUNS/'20261014T051036.115106Z_smoke_C1_n200'
ROW0={200:RUNS/'20261014T045753.248504Z_relax_C1_n200/row0_6000rsun.npz',
      400:RUNS/'20261014T045958.492497Z_relax_C0_n400/row0_6000rsun.npz',
      800:RUNS/'20261014T050203.014930Z_relax_C2_n800/row0_6000rsun.npz'}
DATA=Path(__file__).resolve().parent.parent/'data';DATA.mkdir(exist_ok=True)
BIG=ROOT/'results/ball_merge_20261107';BIG.mkdir(parents=True,exist_ok=True)
YEAR=365.25*86400;MSUN=1.98847e33
from qstars import quasistatic as qs, stepcontrol as sc
import qstars
assert str(Path(qstars.__file__).resolve()).startswith(str(SRC.resolve())),qstars.__file__
import validate_ball_track as vbt
PLAN=json.loads((ROOT/'benchmarks/ball2011_track_plan.json').read_text())
FIRST_REJECTION_DT=172.1438701239381

def load(path=RUN/'checkpoint.npz'):
    state,config=qs.load_checkpoint(path)
    eos,opacity=vbt.physics(config)
    return state,config,eos,opacity

def attempt(state,dt_years,eos,opacity,config,mode='evolve'):
    """One bondi_step plus the registered controller metrics (stepcontrol.metrics)."""
    S=sc.settings(config);s=qs.settings(config);bm=s['cavity']['boundary_state']
    before=sc.observables(state,eos,bm)
    new,d=qs.bondi_step(state,dt_years*YEAR,eos,opacity,config,mode=mode)
    after=sc.observables(new,eos,bm)
    return new,d,sc.metrics(before,after,d,S)

def mbs(st):
    """Mass below the surface at cell centres (g)."""
    return st.dm.sum()-(np.cumsum(st.dm)-.5*st.dm)

def lnT_at(st,mass_below):
    """ln T of st, linear in mass below the surface between cell centres, at the given masses below the surface."""
    x=mbs(st)[::-1];y=np.log(st.T)[::-1]
    return np.interp(mass_below,x,y)

def with_merge(config,merge):
    c=copy.deepcopy(config);c['quasistatic']['merge']=merge;return c

def jdump(obj,path):
    Path(path).write_text(json.dumps(obj,indent=1,default=lambda o:o.item() if hasattr(o,'item') else str(o)))
