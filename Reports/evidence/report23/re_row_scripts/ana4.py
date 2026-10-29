#!/usr/bin/env python3
"""Two-photon up/down transitions (TWO_PHOT_SUM) vs the n=2 and RE residuals (read-only)."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_trace_record
import gencool_breakdown as gb
D = Path(sys.argv[1])
x, R = read_trace_record(D / 'QSTARS_RESID_TRACE.bin')
NT, ND = x.shape
T1 = 4e-10 * np.pi; PC = 6.62607015e-12; dE21 = 3.288087 - 0.822021


def blocks(path):
    lines = Path(path).read_text().splitlines()
    out = {}; key = None; trans = None
    for s in lines:
        t = s.split()
        if len(t) >= 3 and t[1] in ('THI', 'THeI', 'THe2', 'FHI', 'FHeI', 'FHe2'):
            trans = t[1]; key = None; continue
        st = s.strip()
        if st in ('Down rate', 'Up rate', 'Total # of down transitions', 'Total # of upward transitions'):
            key = (trans, st); out[key] = []; continue
        if key is not None:
            try:
                out[key] += [float(v) for v in t]
            except ValueError:
                pass
    return {k: np.array(v) for k, v in out.items()}


b = blocks(D / 'TWO_PHOT_SUM')
for k, v in b.items():
    print(k, len(v))
up = b[('THI', 'Total # of upward transitions')]; dn = b[('THI', 'Total # of down transitions')]
upr = b[('THI', 'Up rate')]; dnr = b[('THI', 'Down rate')]
g = gb.parse(D / 'GENCOOL')
print(' d   T[1e4]  n1         n2        up_rate/atom  dn_rate    N_up[cm-3s-1] N_down     R3         N_up/R3   hv21*N_up   T1*R_T     hv21*N_up/(T1 R_T)  colH/hv21  LTEup=A*exp(-hv0/kT)  up/LTE')
for d in list(range(20, 60)) + [70, 90, 120, 160]:
    k = d - 1
    col = g['HI Collisional Cooling'][d]
    lte = 8.225 * np.exp(-11.8348 / x[133, k])
    print(f'{d:3d} {x[133,k]:7.4f} {x[1,k]:10.3e} {x[2,k]:10.3e} {upr[k]:11.3e} {dnr[k]:10.3e} {up[k]:11.3e} {dn[k]:10.3e} {R[2,k]:10.3e} {up[k]/R[2,k]:8.3f} '
          f'{PC*dE21*up[k]:11.3e} {T1*R[133,k]:11.3e} {PC*dE21*up[k]/(T1*R[133,k]):8.3f}        {col/(PC*dE21):10.3e} {lte:10.3e} {upr[k]/lte:10.3e}')
