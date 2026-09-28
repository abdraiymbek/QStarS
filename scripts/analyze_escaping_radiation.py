#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('main',type=Path);p.add_argument('extent',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
main=json.loads((a.main/'campaign.json').read_text());extent=json.loads((a.extent/'campaign.json').read_text())
checks=json.loads((a.main/'checks.json').read_text())+json.loads((a.extent/'checks.json').read_text())
metrics=dict(main_directory=str(a.main),extent_directory=str(a.extent),main_campaign=main,extent_campaign={k:v for k,v in extent.items() if k!='cases'},
 checks=len(checks),passing_checks=sum(x['passed'] for x in checks),failing_checks=[x for x in checks if not x['passed']],
 published_table=json.loads((a.main/'published_table_comparison.json').read_text()),
 fine_wind=json.loads((a.main/'wind_n801.json').read_text()),extended_wind=json.loads((a.extent/'outer_100000.json').read_text()))
(a.output/'analysis.json').write_text(json.dumps(metrics,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,3,figsize=(15,8.4),layout='constrained')
colors=['#ad4b38','#b68223','#29888b','#235a93']
for n,color in zip([201,401,801,1601],colors):
 d=np.load(a.main/f'grey_n{n}.npz');r=d['r']
 ax[0,0].loglog(r,r*r*d['J'],color=color,lw=1.3,label=f'{n} log nodes')
 ax[0,1].semilogx(r,d['f'],color=color,lw=1.3)
 ax[0,2].semilogx(r,100*(r*r*d['angular_H']-1),color=color,lw=1.2)
t=np.genfromtxt(a.main/'hummer1971_table1.csv',delimiter=',',names=True)
ax[0,0].scatter(t['r'],t['r']**2*t['J'],facecolors='none',edgecolors='black',s=30,zorder=5,label='Published Table I')
ax[0,1].scatter(t['r'],t['f'],facecolors='none',edgecolors='black',s=30,zorder=5)
ax[0,0].set(xlabel='Reference radius r',ylabel='r² J / H₀',title='Published grey radiation profile');ax[0,0].legend(fontsize=8)
ax[0,1].axhline(1/3,color='0.6',ls='--');ax[0,1].set(xlabel='Reference radius r',ylabel='Eddington factor K/J',title='Angular closure through escape')
ax[0,2].axhline(0,color='0.6',lw=.7);ax[0,2].set(xlabel='Reference radius r',ylabel='Angular luminosity error (%)',title='Independent flux: no renormalization')
for name,label,color in [('wind_n801','Small drag','#235a93'),('wind_drag_negative_control','Larger drag','#ad4b38')]:
 d=np.load(a.main/f'{name}.npz');ax[1,0].semilogx(d['r'],100*(np.sqrt(d['w']/d['reference_w'])-1),label=label,color=color)
ax[1,0].axhline(-1,color='0.4',ls='--',label='1% comparison limit');ax[1,0].set(xlabel='Wind radius r/R',ylabel='Velocity vs analytic limit (%)',title='Asymptotic wind test without fitting');ax[1,0].legend(fontsize=8)
for path,label,color in [(a.main/'wind_n801.npz','Outer 100 R','#ad4b38'),(a.main/'wind_outer300.npz','Outer 300 R','#b68223'),(a.extent/'outer_100000.npz','Outer 100,000 R','#235a93')]:
 d=np.load(path);mask=d['r']<=1000;ax[1,1].loglog(d['r'][mask],(d['r']**2*d['J'])[mask],color=color,label=label)
ax[1,1].set(xlabel='Wind radius r/R',ylabel='r² J / (R² Hbase)',title='Outer boundary changes radiation structure');ax[1,1].legend(fontsize=8)
rs=[100,300,10000,30000,100000];photos=[]
for r in rs:
 file=a.main/('wind_n801.json' if r==100 else 'wind_outer300.json') if r<1000 else a.extent/f'outer_{r}.json'
 q=json.loads(file.read_text());photos.append(q.get('photosphere',q.get('photospheric_radius_ratio')))
ax[1,2].semilogx(rs,photos,'o-',color='#235a93');ax[1,2].axhline(photos[-1],color='0.6',ls='--')
ax[1,2].set(xlabel='Vacuum boundary r/R',ylabel='Optical-depth photosphere r/R',title='Photosphere requires extent convergence')
fig.suptitle('QStarS Report 11 — grey component validation\nNo material-temperature, ionization or envelope-launch validation',fontsize=15)
fig.savefig(a.output/'escaping_validation.png',dpi=170)
print(a.output)
