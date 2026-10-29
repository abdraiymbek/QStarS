import sys, numpy as np
from pathlib import Path
ROOT=Path('/Users/abraiymbek/Projects/QStarS'); sys.path.insert(0,str(ROOT/'src'))
from qstars.nlte_newton import read_trace_record
kB=1.380649e-16; h=6.62607015e-27; me=9.1093837e-28; eV=1.602176634e-12
def saha_n1(ne,npr,T):   # LTE H I ground population for given n_e, n_p, T (g1=2, g_p=1)
    lam=(h*h/(2*np.pi*me*kB*T))**1.5
    return ne*npr*lam*(2/2.0)*np.exp(13.598434*eV/(kB*T))
def n2_boltz(n1,T): return n1*4*np.exp(-10.19884*eV/(kB*T))
states={'x0 (cool start)':ROOT/'results/nk167_trap/fixedT/it00_base','nk167_hot it34':ROOT/'results/nk167_hot/it34_trial_l1','pin100 it05':ROOT/'results/nk167_pin/theta1.00/it05_trial_l1'}
for nm,p in states.items():
    x,R=read_trace_record(p/'QSTARS_RESID_TRACE.bin'); T=x[133]*1e4; ne=x[132]; npr=x[31]; n1=x[1]; n2=x[2]
    b1=n1/saha_n1(ne,npr,T); b2=n2/n2_boltz(saha_n1(ne,npr,T),T)
    print(f'== {nm}: depth, T[kK], n_e, xHII, b1 = n1/n1*(Saha), b2 = n2/n2*(LTE rel. to continuum), n2/n1 over Boltzmann')
    for d in (1,5,10,15,18,20,23,26,29,31,33,35,37,39,41,43,45,47,49,52,55,59,65,72,80,90,100,115,130,150,167):
        k=d-1; print(f'  {d:4d} {T[k]/1e3:7.2f} {ne[k]:9.2e} {npr[k]/x[1:32,k].sum():8.1e} {b1[k]:10.3e} {b2[k]:10.3e} {n2[k]/n2_boltz(n1[k],T[k]):10.3e}')
