#!/usr/bin/env python3
"""Scientific component figures; inputs are an existing immutable campaign."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('campaign',type=Path);a=p.parse_args();root=a.campaign
colors=['black','tab:red','tab:green','tab:blue']
fig,axes=plt.subplots(3,1,figsize=(8,9),sharex=True,layout='constrained')
data=list(csv.DictReader(open('benchmarks/nakauchi2017_figure2.csv')))
for i,rho in enumerate([9.65e-9,9.71e-9,9.85e-9,10.5e-9]):
 path=root/('paper_'+str(rho)+'.npz')
 if not path.exists():continue
 f=np.load(path)
 for ax,key in zip(axes,['velocity','density','temperature']):
  ax.loglog(f['radius']/6.957e10,f[key],color=colors[i],label=f'{rho/1e-9:g}e-9 g/cm³')
  d=[x for x in data if x['field']==key and float(x['sonic_density'])==rho]
  ax.loglog([float(x['radius_rsun']) for x in d],[float(x['value']) for x in d],color=colors[i],ls='--',lw=1,alpha=.7)
  ax.set_ylabel({'velocity':'Velocity [cm/s]','density':'Density [g/cm³]','temperature':'Temperature [K]'}[key]);ax.grid(alpha=.15)
axes[0].legend(fontsize=8,ncol=2);axes[-1].set_xlabel('Radius [solar radii]')
fig.suptitle('Nakauchi et al. Fig. 2: solid = QStarS FreeEOS/MESA; dashed = published\nDifferent microphysics; escaping red branch is not reproduced',fontsize=12)
fig.savefig(root/'published_wind_profiles.png',dpi=160);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
for row,case in enumerate(['recombining','diffusion_exit']):
 f=np.load(root/(case+'.npz'));x=(f['radius']/1.152e15-1)*100
 ax=axes[row,0];ax.plot(x,f['mach'],label='Gas Mach');ax.axhline(1,color='grey',ls='--');ax.set_ylabel('Gas Mach number')
 ax2=ax.twinx();ax2.plot(x,f['temperature']/1e4,color='tab:red',ls=':',label='T');ax2.set_ylabel('T / 10000 K',color='tab:red')
 ax.set_title('Recombining / second sonic approach' if row==0 else 'Diffusion validity endpoint')
 ax=axes[row,1];ax.plot(x,f['hydrogen_ionized'],label='H ionized fraction');ax.plot(x,f['eddington_factor'],label='Radiative Eddington ratio');ax.legend(fontsize=8);ax.set_title('Ionization and radiative force')
 ax=axes[row,2]
 for key,label in [('radiative_luminosity','Diffusion'),('radiation_enthalpy_luminosity','Radiation enthalpy'),('gas_internal_luminosity','Gas internal'),('gas_work_luminosity','Gas pressure work'),('kinetic_luminosity','Kinetic')]:
  ax.semilogy(x,f[key],label=label)
 ax.legend(fontsize=7);ax.set_title('Positive energy channels [erg/s]')
 for ax in axes[row]:ax.set_xlabel('100 × (r / r sonic − 1)');ax.grid(alpha=.15)
fig.suptitle('Spherical transonic controls: 10000 solar masses, real recombination, no MLT\nPrescribed sonic states; neither curve is a certified stellar mass-loss boundary',fontsize=12)
fig.savefig(root/'recombination_components.png',dpi=160);plt.close(fig)
print(root/'published_wind_profiles.png');print(root/'recombination_components.png')
