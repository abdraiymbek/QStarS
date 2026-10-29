#!/usr/bin/env python3
"""Read-only numbers for the RE-row puzzle. Usage: ana1.py <pass_dir> [gencool_path]"""
import sys, re, json
from pathlib import Path
import numpy as np
ROOT = Path('/Users/abraiymbek/Projects/QStarS')
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_ba_blocks, read_trace_record
import gencool_breakdown as gb

np.set_printoptions(linewidth=250)
D = Path(sys.argv[1])
GC = Path(sys.argv[2]) if len(sys.argv) > 2 else D / 'GENCOOL'
x, R = read_trace_record(D / 'QSTARS_RESID_TRACE.bin')
NT, ND = x.shape
lo, di, up = read_ba_blocks(D / 'QSTARS_BA_TRACE.bin', NT, ND)
rs = (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T      # (NT, ND) raw native row l1
T1 = 1e-10 * 4 * np.pi
h = 6.62607015e-27
PC = 1e15 * h


def edges(path):
    out = []
    lines = Path(path).read_text().splitlines()
    n = None
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
E[0] = edges(D / 'H0_F_OSCDAT')[0]
E[1:31] = edges(D / 'HI_F_OSCDAT')[:30]
E[32:101] = edges(D / 'HeI_F_OSCDAT')[:69]
E[101:131] = edges(D / 'He2_F_OSCDAT')[:30]
print('E(H-)=', E[0], 'E(HI 1..4)=', E[1:5], 'E(HeI 1..3)=', E[32:35], 'E(He2 1)=', E[101])

# adiabatic terms
ad = []
for s in (D / 'ADIABAT_CHK').read_text().splitlines():
    t = s.split()
    if len(t) == 9:
        try:
            ad.append([float(v) for v in t])
        except ValueError:
            pass
ad = np.array(ad)
assert ad.shape[0] == ND, ad.shape
g = gb.parse(GC)
def G(name, d):
    return g[name].get(d, float('nan'))

ds = list(range(28, 46))
xH = x[1:32].sum(0); xHtot = x[0:32].sum(0)
print(f'state {D}\nGENCOOL {GC}')
print('\n--- A. state at the front ---')
print(' d      r[1e10cm]   T[1e4K] T(GENCOOL)   n_e        N_H        HI/H       n1          n2         n3        H-        HII')
for d in ds:
    k = d - 1
    print(f'{d:3d} {ad[k,0]:12.5e} {x[133,k]:9.5f} {G("Temperature [1.0E+4K]",d):9.5f} {x[132,k]:10.3e} {xHtot[k]:10.3e} {x[1:31,k].sum()/xHtot[k]:9.4f} '
          f'{x[1,k]:10.3e} {x[2,k]:10.3e} {x[3,k]:10.3e} {x[0,k]:9.2e} {x[31,k]:10.3e}')

print('\n--- B. RE row (row 134) ---')
print(' d   R_T[raw]     l1[raw BA]   R_T/l1     4pi1e-10*R_T   GENCOOL_net   -(GC_net)/T1    ratio(T1 R_T / -GCnet)  adi(A+B+C)*T1  D*T1       WORK*T1')
for d in ds:
    k = d - 1
    net = G('Net Cooling Rate [ergs/cm**3/sec]', d)
    print(f'{d:3d} {R[133,k]:12.4e} {rs[133,k]:12.4e} {R[133,k]/rs[133,k]:10.3e} {T1*R[133,k]:13.4e} {net:13.4e} {-net/T1:13.4e} '
          f'{T1*R[133,k]/(-net):12.3e}      {T1*ad[k,5:8].sum():11.3e} {T1*ad[k,8]:11.3e} {T1*ad[k,5:9].sum():11.3e}')

print('\n--- C. hydrogen rate rows: R/row-l1 (rows: 1=H-, 2=HI ground [replaced by ionization eq], 3=n2, 4=n3, ...; 32=number cons.) ---')
print(' d  ' + ''.join(f'{"row"+str(i+1):>11s}' for i in [0, 1, 2, 3, 4, 5, 6, 10, 20, 30, 31, 132]))
for d in ds:
    k = d - 1
    print(f'{d:3d} ' + ''.join(f'{(R[i,k]/rs[i,k] if rs[i,k]>0 else np.nan):11.2e}' for i in [0, 1, 2, 3, 4, 5, 6, 10, 20, 30, 31, 132]))
print(' raw residuals [cm^-3 s^-1]')
for d in ds:
    k = d - 1
    print(f'{d:3d} ' + ''.join(f'{R[i,k]:11.2e}' for i in [0, 1, 2, 3, 4, 5, 6, 10, 20, 30, 31, 132]))
print(' row l1 scales')
for d in ds:
    k = d - 1
    print(f'{d:3d} ' + ''.join(f'{rs[i,k]:11.2e}' for i in [0, 1, 2, 3, 4, 5, 6, 10, 20, 30, 31, 132]))
print(' max over HI rows 3..31 of |R|/l1, and which')
for d in ds:
    k = d - 1
    q = np.abs(R[2:31, k]) / rs[2:31, k]
    print(f'{d:3d} max {q.max():.3e} at row {2+int(q.argmax())+1}; l1-weighted mean {np.abs(R[2:31,k]).sum()/rs[2:31,k].sum():.3e}')

print('\n--- D. energy bookkeeping (erg/cm3/s) ---')
print('  dEint = PC*[ sum_{n>=2}(E1-En) R_n  - E1*row2 - E(H-)*R_Hminus ]  (H);  He: PC*sum (E1-Ej) R_j levels only (ion rows replaced unknown)')
print(' d   T1*R_T      -GCnet(all)  adi_th*T1   PC*sumH_exc   PC*(-E1*row2)  PC*Hminus   dEint_H     closure: T1*R_T -( -GCnet_rc - adi_th + dEint_H)   He_levels')
for d in ds:
    k = d - 1
    net = G('Net Cooling Rate [ergs/cm**3/sec]', d)
    adv = G('Adiabatic cooling rate (V term) [ergs/cm**3/sec]', d); adt = G('Adiabatic cooling rate (dTdR term) [ergs/cm**3/sec]', d)
    net_rc = net - (0 if np.isnan(adv) else adv) - (0 if np.isnan(adt) else adt)
    exc = PC * ((E[1] - E[2:31]) * R[2:31, k]).sum()
    ion = PC * (-E[1] * R[1, k])
    hm = PC * (-E[0] * R[0, k])
    dE = exc + ion + hm
    adi = T1 * ad[k, 5:8].sum()
    he = PC * ((E[32] - E[33:101]) * R[33:101, k]).sum()
    pred = -net_rc - adi + dE
    print(f'{d:3d} {T1*R[133,k]:11.3e} {-net:11.3e} {adi:11.3e} {exc:12.3e} {ion:12.3e} {hm:11.3e} {dE:11.3e}   pred {pred:11.3e}  diff {T1*R[133,k]-pred:11.3e}   {he:10.2e}')
