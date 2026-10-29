"""DRYRUN helpers of the final fix round (scratch only; never writes under benchmarks/results)."""
import os,sys,json,time
for k in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
ROOT='/Users/abraiymbek/Projects/QStarS'
sys.path.insert(0,ROOT+'/scripts');sys.path.insert(0,ROOT+'/src')
import numpy as np
import validate_ball_track as vbt
from qstars import quasistatic,photosphere,ledger
from qstars.constants import MSUN,RSUN,LSUN,YEAR
LIB=ROOT+'/backends/mesa/build_ball_abi7/libqstars_mesa.dylib'
MESA_DIR=json.load(open(ROOT+'/backends/mesa/build_ball_abi7/build.json'))['mesa_dir']
PLAN=json.loads(open(ROOT+'/benchmarks/ball2011_track_plan.json').read())
HERE=os.path.dirname(os.path.abspath(__file__))
OUT=HERE+'/DRYRUN_out'
A4={'quasistatic':{'photosphere_bracket':'scan'}}
def settings(n):
    if n==4:return A4
    return json.load(open(f'{ROOT}/benchmarks/ball2011_track_plan_amendment_{n}.draft.json'))['settings']
def amend(flags):
    """[(record, settings)] in the order given, e.g. 'A3,A4,A2' (DRYRUN records: no registration)."""
    out=[]
    for f in [x for x in flags.split(',') if x]:
        n=int(f[1:]);out.append((dict(path=f'DRYRUN:amendment_{n}',sha256=f'dryrun-amendment-{n}'),settings(n)))
    return out
def merged(amendments):
    """Merge amendment settings section by section (as build_config does for one section at a time)."""
    return amendments
def summary(out):
    lines=[f"status {out['status']} {out.get('classification','')}"]
    for e in out['seeds']:
        lines.append(f"   seed {e['seed_radius_rsun']:g}: {e['status']} acc {e.get('accepted_steps')} rej {e.get('rejections')} "
                     f"cycles {len(e['mesh_cycles']) if isinstance(e.get('mesh_cycles'),list) else e.get('mesh_cycles')} "
                     f"{(e.get('reason') or '')[:110]}")
    if 'uniqueness' in out:
        u=out['uniqueness'];k=max(u['max_abs_dln'],key=u['max_abs_dln'].get)
        lines.append(f"   uniqueness seeds {u['seeds']} max|dln| {u['max_abs_dln'][k]!r} ({k})")
    return '\n'.join(lines)
