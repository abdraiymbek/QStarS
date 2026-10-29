"""Offline reconstruction of the cool_newton_krylov Iterate for results/nk167_hot iterates (no native runs)."""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa
from qstars.nlte_newton import read_trace_record, read_jv_rad  # noqa

RUN = ROOT / 'results/nk167_hot'
HIST = json.loads((RUN / 'history.json').read_text())


def setup_globals(ultra_frac=1e-20, pc_floor=1e-3):
    x0, _ = read_trace_record(RUN / 'it00_base' / 'QSTARS_RESID_TRACE.bin')
    ck.NT, ck.ND = x0.shape
    ck.SIG_FRAC = 1e-10; ck.TRACE_WEIGHT = 1e-3
    ck.ULTRA_FRAC = ultra_frac; ck.ULTRA_WEIGHT = 1e-8
    ck.MERIT_MODE = 'precond'; ck.EQ_FLOOR = 0.0; ck.PC_T_FLOOR = pc_floor
    ck.SPLIT_REGIONS = ck.split_regions(ck.ND)
    return x0


def load_iterate(k, x0, delta=None):
    rec = HIST[k]
    base = Path(rec['base'])
    if not base.is_absolute():
        base = ROOT / base
    base = RUN / base.name
    if delta is None:
        delta = rec.get('ptc_delta')
    hom = dict(mode='ptc', lam=1.0, delta=delta)
    cdir = RUN / f'cols_{base.name}'
    cols = {n: (cdir / f'col_{n}_p', cdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(base, x0, cols, 1e-4, True, hom)
    return it, rec


def load_pairs(k, base_name):
    itd = RUN / f'it{k:02d}_from_{base_name}'
    Z, JZ, pas = [], [], []
    for p, j in enumerate(sorted(itd.glob('jv*'))):
        P = np.load(j / 'pairs.npz')
        for i in range(P['z'].shape[2]):
            Z.append(P['z'][:, :, i]); JZ.append(P['jz'][:, :, i]); pas.append(p)
    return Z, JZ, np.array(pas), itd
