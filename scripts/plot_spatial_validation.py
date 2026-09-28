#!/usr/bin/env python3
"""Plot retained spatial evidence only; does not run a simulation."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from qstars.io import load_checkpoint
p=argparse.ArgumentParser();p.add_argument('campaign',type=Path);a=p.parse_args()
x=json.loads(a.campaign.read_text());rows=x['rows']
fig,ax=plt.subplots(2,2,figsize=(11,8.2),layout='constrained')
ns=np.array([r['zones'] for r in rows]);err=np.array([r['initial_density_mass_log_L2'] for r in rows])
ax[0,0].loglog(ns,err,'o-',label='Independent continuum comparison')
ax[0,0].loglog(ns,err[-1]*(ns[-1]/ns)**2,'--',color='gray',label=r'$N^{-2}$')
ax[0,0].set(xlabel='Zones',ylabel=r'Mass-weighted RMS $\Delta\ln\rho$',title='Initial hydrostatic structure')
ax[0,0].legend(fontsize=8)
fields={'density':'Density','temperature':'Temperature','pressure':'Pressure','specific_internal_energy':'Specific energy','radius':'Radius'}
for field,label in fields.items():
 usable=[r for r in rows if 'final_comparison_to_finest' in r]
 ax[0,1].loglog([r['zones'] for r in usable],[r['final_comparison_to_finest']['fields'][field]['interpolated_log_Linf'] for r in usable],'o-',label=label)
ax[0,1].axhline(.01,color='black',ls=':',label='1% target')
ax[0,1].set(xlabel='Zones',ylabel='Maximum log difference to finest',title='Ten-year internal profiles')
ax[0,1].legend(fontsize=8,ncol=2)
states=[]
for row in rows:
 s,c=load_checkpoint(Path(row['path'])/'checkpoint.npz');states.append(s)
 q=(np.cumsum(s.dm)-s.dm/2)/s.dm.sum()
 ax[1,0].loglog(1-q,s.rho,label=str(row['zones'])+' zones')
ax[1,0].set(xlim=(1e-6,.01),xlabel='Exterior envelope mass fraction, 1−q',ylabel=r'Density [g cm$^{-3}$]',title='Outer layer retains the largest error')
ax[1,0].legend(fontsize=8)
ref=states[-1];qr=(np.cumsum(ref.dm)-ref.dm/2)/ref.dm.sum()
for row,s in zip(rows[:-1],states[:-1]):
 q=(np.cumsum(s.dm)-s.dm/2)/s.dm.sum();lo=max(q[0],qr[0]);hi=min(q[-1],qr[-1])
 grid=1-np.geomspace(max(1-hi,1e-8),min(1-lo,.1),2000)
 d=PchipInterpolator(q,np.log(s.rho))(grid)-PchipInterpolator(qr,np.log(ref.rho))(grid)
 ax[1,1].semilogx(1-grid,100*np.expm1(d),label=str(row['zones'])+' vs '+str(rows[-1]['zones']))
ax[1,1].axhline(0,color='gray',lw=.7);ax[1,1].axhline(1,color='black',ls=':',lw=.7);ax[1,1].axhline(-1,color='black',ls=':',lw=.7)
ax[1,1].set(xlabel='Exterior envelope mass fraction, 1−q',ylabel='Density difference [%]',title='Common cell-center domain; coverage in JSON')
ax[1,1].legend(fontsize=8)
for aa in ax.flat:aa.grid(alpha=.2)
fig.suptitle('QStarS: bounded hot-envelope convergence evidence',fontsize=15)
fig.savefig(a.campaign.parent/'spatial_validation.png',dpi=180)
print(a.campaign.parent/'spatial_validation.png')
