#!/usr/bin/env python3
"""Closed energy budget of the RE row at chosen depths (read-only).
 T1*R_T = -GENCOOL_net(rad+coll) - adiabatic_thermal + h*sum_j e_j R_j + M_2gamma
 M_2gamma = RE-minus-SE energy of the OLD_DEFAULT two-photon term, recomputed from the stored J.
Usage: ana6.py <pass_dir> [d_first d_last]"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_ba_blocks, read_trace_record
import gencool_breakdown as gb
D = Path(sys.argv[1])
d0, d1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (28, 45)
x, R = read_trace_record(D / 'QSTARS_RESID_TRACE.bin')
NT, ND = x.shape
T1 = 4e-10 * np.pi; PC = 6.62607015e-12
TWOHCSQ = 0.0147452575; HDKT = 4.7994145


def edges(path):
    out = []; lines = Path(path).read_text().splitlines(); n = None
    for i, s in enumerate(lines):
        if 'Number of energy levels' in s:
            n = int(s.split()[0])
        if 'Number of transitions' in s:
            j = i + 1
            while len(out) < n:
                t = lines[j].split(); j += 1
                if len(t) >= 5:
                    out.append(float(t[3]))
            break
    return np.array(out)


E = np.zeros(NT)
E[0] = edges(D / 'H0_F_OSCDAT')[0]; E[1:31] = edges(D / 'HI_F_OSCDAT')[:30]
ad = []
for s in (D / 'ADIABAT_CHK').read_text().splitlines():
    t = s.split()
    if len(t) == 9:
        try: ad.append([float(v) for v in t])
        except ValueError: pass
ad = np.array(ad)
g = gb.parse(D / 'GENCOOL')
G = lambda n, d: g.get(n, {}).get(d, 0.0)
have_ba = (D / 'QSTARS_BA_TRACE.bin').exists()
if have_ba:
    lo, di, up = read_ba_blocks(D / 'QSTARS_BA_TRACE.bin', NT, ND)
    l1 = lambda i, k: np.abs(lo[k][i]).sum() + np.abs(di[k][i]).sum() + np.abs(up[k][i]).sum()
# two-photon from stored J
raw = np.fromfile(D / 'QSTARS_RJ.bin', dtype=np.dtype([('ml', '<i4'), ('fl', '<f8'), ('rj', '<f8', (ND,))]), offset=12)
nu = raw['fl']; J = raw['rj']
fqw = np.empty_like(nu); fqw[1:-1] = 0.5 * (nu[:-2] - nu[2:]); fqw[0] = 0.5 * (nu[0] - nu[1]); fqw[-1] = 0.5 * (nu[-2] - nu[-1]); fqw *= 1e15
nu0 = E[1] - E[2]; A = 8.225 / 4.0; GU, GL = 8.0, 2.0
m = nu < nu0; v = nu[m]; w = fqw[m]; Jm = J[m]
y = v / nu0; U = y * (1 - y); FU = 4 * U
AY = 24.56 * A * (U * (1 - FU ** 0.8) + 0.88 * U ** 1.53 * FU ** 0.8)
alpha = TWOHCSQ * v ** 3; vB = nu0 - v; lognu = np.log(nu[::-1]); upc = 0.5e-15 * AY / nu0 * GU / GL
hv = PC * nu0
print(f'state {D}')
print(' d  T[1e4K]  R_T/l1    R3/l1  |  T1*R_T    = -GCnet_rc   -adi_th    +h*sum(e_j R_j)  +M_2g      (closure resid/T1R_T) |  R3      N2g_old  N2g_exact  C_coll<=   R3-N2g_old |  2g share of T1R_T   SE-limit T1*R_T(old)=-GC+M2g   SE-limit(exact 2g)=-GC   [raw/l1]')
for d in range(d0, d1 + 1):
    k = d - 1
    T = x[133, k]; n1 = x[1, k]
    net = G('Net Cooling Rate [ergs/cm**3/sec]', d) - G('Adiabatic cooling rate (V term) [ergs/cm**3/sec]', d) - G('Adiabatic cooling rate (dTdR term) [ergs/cm**3/sec]', d)
    adi = T1 * ad[k, 5:8].sum()
    dE = PC * (((E[1] - E[2:31]) * R[2:31, k]).sum() - E[1] * R[1, k] - E[0] * R[0, k])
    Yj = Jm[:, k] / alpha
    old = upc * np.exp(-HDKT * vB / T) * Yj * w
    JB = np.exp(np.interp(np.log(vB), lognu, np.log(J[::-1, k])))
    new = upc * Yj * (JB / (TWOHCSQ * vB ** 3)) * w
    M = hv * n1 * 2.0 * ((y - 0.5) * old).sum()
    Nold = n1 * old.sum(); Nnew = n1 * new.sum()
    lhs = T1 * R[133, k]
    pred = -net - adi + dE + M
    l1T = l1(133, k) if have_ba else np.nan; l13 = l1(2, k) if have_ba else np.nan
    se_old = -net - adi + M; se_new = -net - adi
    print(f'{d:3d} {T:7.4f} {R[133,k]/l1T:9.2e} {R[2,k]/l13:9.2e} | {lhs:10.3e} = {-net:10.3e} {-adi:10.2e} {dE:11.3e} {M:11.3e}   ({(lhs-pred)/lhs:9.1e}) | {R[2,k]:9.3e} {Nold:9.3e} {Nnew:9.3e} '
          f'{G("HI Collisional Cooling",d)/hv:9.3e} {R[2,k]-Nold:9.3e} |  {(hv*Nold+M)/lhs:6.3f}       {se_old:10.3e} [{se_old/T1/l1T:9.2e}]   {se_new:10.3e} [{se_new/T1/l1T:9.2e}]')
