"""V1 plan P1: offline reconstruction of a free_r25 iterate's least squares from its stored exact products (helpers).

load(ITERATE_DIR) rebuilds, as the driver did at that iterate: the Iterate (rowscale measure, H I n=2 significant), the
span S (dense-window columns first, then the Krylov passes jv00p*, jv01p*, ... in order) with the exact products JS
(pairs.npz: Jdir z + native radiation part), A = obj(JS), b = obj(R), the bound radii (0.5, T 0.03) and mask. Read-only."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402

O = ROOT / 'results/nk167_trap/free_r25'
X0 = cn.configure()


def pairs(dirs):
    S, JS = [], []
    for d in dirs:
        p = np.load(d / 'pairs.npz')
        S += [p['z'][:, :, i] for i in range(p['z'].shape[2])]
        JS += [p['jz'][:, :, i] for i in range(p['z'].shape[2])]
    return S, JS


class Problem:
    def __init__(self, itdir, dense=True):
        self.itdir = O / itdir if not Path(itdir).is_absolute() else Path(itdir)
        self.k = int(self.itdir.name[2:4])
        self.base = O / self.itdir.name.split('_from_', 1)[1]
        self.it = cn.iterate_at(self.itdir, X0)
        it = self.it
        dd = sorted(self.itdir.glob('dense[0-9][0-9]')) if dense else []
        jv = sorted(self.itdir.glob('jv[0-9][0-9]p[0-9]'))
        self.ndense = sum(np.load(d / 'pairs.npz')['z'].shape[2] for d in dd)
        S, JS = pairs(dd + jv)
        self.S, self.JS = S, JS
        self.shape = it.x.shape
        self.Sm = np.stack([s.ravel() for s in S], axis=1)
        self.A = np.stack([it.obj(j).ravel() for j in JS], axis=1)
        self.b = it.obj(it.R).ravel()
        self.beta = float(np.linalg.norm(self.b))
        self.w = (it.row_weights() * it.D * it.active)
        self.rad = np.full(self.shape, 0.5); self.rad[133] = 0.03
        self.mask = it.active.copy()

    def lin(self, y):
        return float(np.linalg.norm(self.A @ y - self.b) / self.beta)

    def step(self, y):
        return (self.Sm @ y).reshape(self.shape)

    def bound_ratio(self, y):
        s = self.step(y)
        return float(np.max(np.abs(s[self.mask]) / self.rad[self.mask]))


def cached(itdir, dense=True):
    """Problem with A, b, Sm, w, mask, rad cached in the scratch directory (the loads take ~2 min)."""
    import os, pickle
    cdir = Path(os.environ.get('P1_CACHE', '/private/tmp/claude-501/p1_cache')); cdir.mkdir(parents=True, exist_ok=True)
    f = cdir / f'{itdir}{"" if dense else "_nodense"}.pkl'
    if f.exists():
        return pickle.loads(f.read_bytes())
    P = Problem(itdir, dense)
    f.write_bytes(pickle.dumps(P, protocol=5))
    return P
