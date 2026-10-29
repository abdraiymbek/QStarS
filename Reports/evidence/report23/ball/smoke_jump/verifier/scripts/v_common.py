"""Verifier's own setup (scratch only; never writes under benchmarks/results)."""
import os
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import sys,copy,json,math,time
from pathlib import Path
import numpy as np
SRC=os.environ.get('VSRC','/Users/abraiymbek/Projects/QStarS')
ROOT=Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0,str(Path(SRC)/'src'));sys.path.insert(0,str(Path(SRC)/'scripts'))
RUN=ROOT/'benchmarks/results/2026-10-14-ball-track/20261014T051036.115106Z_smoke_C1_n200'
OUT=Path(__file__).resolve().parent
YEAR=365.25*86400;MSUN=1.98847e33
from qstars import quasistatic as qs, stepcontrol as sc
import qstars
assert str(Path(qstars.__file__).resolve()).startswith(str(Path(SRC).resolve())),qstars.__file__
import validate_ball_track as vbt

def load():
    state,config=qs.load_checkpoint(RUN/'checkpoint.npz')
    eos,opacity=vbt.physics(config)
    return state,config,eos,opacity

def attempt(state,dt_years,eos,opacity,config,mode='evolve'):
    S=sc.settings(config);s=qs.settings(config);bm=s['cavity']['boundary_state']
    before=sc.observables(state,eos,bm)
    new,d=qs.bondi_step(state,dt_years*YEAR,eos,opacity,config,mode=mode)
    after=sc.observables(new,eos,bm)
    return new,d,sc.metrics(before,after,d,S)

def below_surface_map(old,new):
    """ln T of `old` at the cell centres of `new`, in mass below the surface (as the outer-half predictor)."""
    mo=np.cumsum(old.dm)-.5*old.dm;mn=np.cumsum(new.dm)-.5*new.dm
    To,Tn=old.dm.sum(),new.dm.sum()
    return np.interp((Tn-mn)[::-1],(To-mo)[::-1],np.log(old.T)[::-1])[::-1]
