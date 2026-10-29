#!/usr/bin/env python3
"""Recompute the H I 1s-2s two-photon upward rate from the pass's stored J (QSTARS_RJ.bin):
 OLD_DEFAULT (two_phot_rate.f:68-76: companion photon occupation exp(-h nu_B/kT_e))  vs
 USE_RAD     (two_phot_rate.f:77-84: companion photon occupation J(nu_B)/(2h nu_B^3/c^2)),
and the RE-vs-SE energy mismatch of the OLD_DEFAULT opacity (two_phot_opac_v3.f:119-131). Read-only."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_trace_record
D = Path(sys.argv[1])
RJ = Path(sys.argv[2]) if len(sys.argv) > 2 else D / 'QSTARS_RJ.bin'
x, R = read_trace_record(D / 'QSTARS_RESID_TRACE.bin')
NT, ND = x.shape
raw = np.fromfile(RJ, dtype=np.dtype([('ml', '<i4'), ('fl', '<f8'), ('rj', '<f8', (ND,))]), offset=12)
nu = raw['fl']; J = raw['rj']                     # nu decreasing, 1e15 Hz
assert np.all(np.diff(nu) < 0), 'frequency grid not monotonic'
fqw = np.empty_like(nu)                           # trapezoid weights in Hz
fqw[1:-1] = 0.5 * (nu[:-2] - nu[2:]); fqw[0] = 0.5 * (nu[0] - nu[1]); fqw[-1] = 0.5 * (nu[-2] - nu[-1])
fqw *= 1e15
TWOHCSQ = 0.0147452575; HDKT = 4.7994145; PC = 6.62607015e-12; T1 = 4e-10 * np.pi
nu0 = 3.288087 - 0.822021; A = 8.225 / 4.0; GU, GL = 8.0, 2.0
m = nu < nu0
v = nu[m]; w = fqw[m]; Jm = J[m]
y = v / nu0; U = y * (1 - y); FU = 4 * U
AY = 24.56 * A * (U * (1 - FU ** 0.8) + 0.88 * U ** 1.53 * FU ** 0.8)
print('check: 0.5*int AY dy =', 0.5 * np.sum(AY * w * 1e-15 / nu0), ' (adjusted A =', A, ')')
alpha = TWOHCSQ * v ** 3
vB = nu0 - v
# J at the companion frequency (log-log interpolation on the stored grid)
lognu = np.log(nu[::-1]);
ds = [28, 30, 33, 34, 35, 36, 37, 39, 41, 43, 45]
upc = 0.5e-15 * AY / nu0 * GU / GL
print(' d   T      n1        UP_old[calc]  UP_old[file]  UP_use_rad   old/use_rad  N_up_old   N_up_rad    <nu>/nu0(old)  RE-SE mismatch(old) [erg/cm3/s]   R3      (R3-N_up_old+N_up_rad)')
for d in ds:
    k = d - 1
    T = x[133, k]; n1 = x[1, k]
    Yj = Jm[:, k] / alpha
    JB = np.exp(np.interp(np.log(vB), lognu, np.log(J[::-1, k])))
    YB = JB / (TWOHCSQ * vB ** 3)
    old_int = upc * np.exp(-HDKT * vB / T) * Yj * w
    up_old = old_int.sum()
    up_rad = (upc * Yj * YB * w).sum()
    ymean = (y * old_int).sum() / up_old
    mism = PC * nu0 * n1 * 2.0 * ((y - 0.5) * old_int).sum()
    print(f'{d:3d} {T:6.4f} {n1:9.3e} {up_old:12.4e} {"":12s} {up_rad:12.4e} {up_old/up_rad:10.3e} {n1*up_old:10.3e} {n1*up_rad:10.3e} {ymean:10.3f}       {mism:12.4e}                 {R[2,k]:9.3e} {R[2,k]-n1*up_old+n1*up_rad:10.3e}')
