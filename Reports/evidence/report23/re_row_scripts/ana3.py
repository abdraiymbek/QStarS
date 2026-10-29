#!/usr/bin/env python3
"""Energy identity over all depths + top energy-weighted rows (read-only)."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_trace_record
import gencool_breakdown as gb
sys.path.insert(0, str(Path(__file__).parent))

D = Path(sys.argv[1])
x, R = read_trace_record(D / 'QSTARS_RESID_TRACE.bin')
NT, ND = x.shape
T1 = 4e-10 * np.pi; PC = 6.62607015e-12


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
E[32:101] = edges(D / 'HeI_F_OSCDAT')[:69]; E[101:131] = edges(D / 'He2_F_OSCDAT')[:30]
ad = []
for s in (D / 'ADIABAT_CHK').read_text().splitlines():
    t = s.split()
    if len(t) == 9:
        try: ad.append([float(v) for v in t])
        except ValueError: pass
ad = np.array(ad)
g = gb.parse(D / 'GENCOOL')
net = np.array([g['Net Cooling Rate [ergs/cm**3/sec]'][d] for d in range(1, ND + 1)])
adv = np.array([g.get('Adiabatic cooling rate (V term) [ergs/cm**3/sec]', {}).get(d, 0.0) for d in range(1, ND + 1)])
adt = np.array([g.get('Adiabatic cooling rate (dTdR term) [ergs/cm**3/sec]', {}).get(d, 0.0) for d in range(1, ND + 1)])
colH = np.array([g['HI Collisional Cooling'][d] for d in range(1, ND + 1)])
net_rc = net - adv - adt
exc = PC * ((E[1] - E[2:31])[:, None] * R[2:31]).sum(0)
ion = PC * (-E[1] * R[1]); hm = PC * (-E[0] * R[0])
# He: e_int: HeI level j: E1-Ej ; HeII ground and above: E1(HeI) + (E1(HeII)-Ej); HeIII: E1(HeI)+E1(HeII)
# rows: 33..101 HeI levels (idx 32..100); 102 (idx 101) HeII ground replaced by ionization eq of He2 (= -R_HeIII); idx 131 HeIII row = number conservation
# R_HeIIground(full) = -(sum HeI levels + sum HeII excited + R_HeIII),  R_HeIII = -R[101]
RHe3 = -R[101]
RHe2g = -(R[32:101].sum(0) + R[102:131].sum(0) + RHe3)
he = PC * (((E[32] - E[32:101])[:, None] * R[32:101]).sum(0) + E[32] * RHe2g
           + ((E[32] + E[101] - E[102:131])[:, None] * R[102:131]).sum(0) + (E[32] + E[101]) * RHe3)
dE = exc + ion + hm + he
adi = T1 * ad[:, 5:8].sum(1)
pred = -net_rc - adi + dE
lhs = T1 * R[133]
print(' d   T[1e4]   n_e      T1*R_T      dEint(H)    dEint(He)   -GCnet_rc   adi_th     pred       gap=lhs-pred  gap/lhs  gap/hv21R3  gap/colH   lhs/(hv21 R3)')
for d in range(1, ND + 1):
    k = d - 1
    if d % 1 == 0 and (d <= 60 or d % 5 == 0):
        hv = PC * (E[1] - E[2]) * R[2, k]
        print(f'{d:3d} {x[133,k]:7.4f} {x[132,k]:9.2e} {lhs[k]:11.3e} {exc[k]+ion[k]+hm[k]:11.3e} {he[k]:11.3e} {-net_rc[k]:11.3e} {adi[k]:10.2e} {pred[k]:11.3e} {lhs[k]-pred[k]:11.3e} '
              f'{(lhs[k]-pred[k])/lhs[k]:8.3f} {(lhs[k]-pred[k])/hv:9.3f} {(lhs[k]-pred[k])/colH[k]:9.3f} {lhs[k]/hv:9.3f}')
for d in (30, 35, 41):
    k = d - 1
    eint = np.zeros(NT); eint[2:31] = E[1] - E[2:31]; eint[1] = -E[1]; eint[0] = -E[0]; eint[32:101] = E[32] - E[32:101]
    c = PC * eint * R[:, k]
    o = np.argsort(-np.abs(c))[:8]
    print(f'd{d}: top energy-weighted rows (row index 1-based, erg/cm3/s):', [(int(i) + 1, float(f'{c[i]:.3e}')) for i in o])
