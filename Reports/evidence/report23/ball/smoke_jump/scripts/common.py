"""Shared setup for the smoke-jump diagnosis (scratch only; never writes under benchmarks/results)."""
import os
for _n in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_n]='1'
import sys,copy,json,math
from pathlib import Path
import numpy as np
ROOT=Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
RUN=ROOT/'benchmarks/results/2026-10-14-ball-track/20261014T051036.115106Z_smoke_C1_n200'
SCR=Path(__file__).resolve().parent
YEAR=365.25*86400
from qstars import quasistatic as qs, stepcontrol as sc
import validate_ball_track as vbt

def load():
    state,config=qs.load_checkpoint(RUN/'checkpoint.npz')
    eos,opacity=vbt.physics(config)
    return state,config,eos,opacity

def step(state,dt_years,eos,opacity,config,mode='evolve'):
    """One driver-equivalent attempt: bondi_step + controller metrics (as quasistatic.run)."""
    S=sc.settings(config);s=qs.settings(config);bm=s['cavity']['boundary_state']
    before=sc.observables(state,eos,bm)
    new,d=qs.bondi_step(state,dt_years*YEAR,eos,opacity,config,mode=mode)
    after=sc.observables(new,eos,bm)
    m=sc.metrics(before,after,d,S)
    return new,d,m
