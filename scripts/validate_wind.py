#!/usr/bin/env python3
"""Published wind profiles, numerical convergence, and independent ODE integrals."""
import csv,copy,time
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from qstars import io
from qstars.wind_reference import solve_wind

out=io.create_directory('results','fiacconi_wind_validation')
io.provenance(out,{'family':'fiacconi2016_wind','source':'Figure 4, equations 12-22'})
data=list(csv.DictReader(open('benchmarks/fiacconi2016_figure4.csv')))
rows=[];profiles={}
for alpha in (.4,7.64):
    for mach in (.75,.8,.85):
        start=time.perf_counter();f,s=solve_wind(alpha,5.,100.,mach,samples=4001)
        f['density_ratio']=f['radius_ratio']**-2*np.sqrt(f['w'][0]/f['w'])
        profiles[(alpha,mach)]=f
        np.savez_compressed(out/f'alpha_{alpha}_mach_{mach}.npz',**f)
        s['runtime_seconds']=time.perf_counter()-start;rows.append(s)
    base=profiles[(alpha,.8)]
    for label,kwargs in [('tolerance',{'rtol':1e-10}),('outer_boundary',{'delta':1e-7}),('independent_integrator',{'method':'DOP853'})]:
        f,s=solve_wind(alpha,5.,100.,.8,samples=4001,**kwargs)
        grid=np.geomspace(1,1e4,501)
        s['comparison']=label;s['profile_log_Linf']={}
        for field in ('velocity_ratio','u','luminosity','advected_luminosity','mach'):
            a=PchipInterpolator(np.log(base['radius_ratio']),np.log(base[field]))(np.log(grid))
            b=PchipInterpolator(np.log(f['radius_ratio']),np.log(f[field]))(np.log(grid))
            s['profile_log_Linf'][field]=float(max(abs(a-b)))
        rows.append(s)
comparisons=[]
for d in data:
    alpha=float(d['alpha']);field=d['field'];r=float(d['radius_ratio']);value=float(d['value']);unc=float(d['log_uncertainty'])
    predictions=[]
    for mach in (.75,.8,.85):
        f=profiles[(alpha,mach)]
        predictions.append(float(PchipInterpolator(np.log(f['radius_ratio']),np.log(f[field]))(np.log(r))))
    distance=max(0,min(predictions)-np.log(value),np.log(value)-max(predictions))
    comparisons.append(dict(**d,log_distance_outside_prescribed_mach_band=float(distance),passed=bool(distance<=unc)))
summary=dict(scope='Stationary FR16 radiation-dominated wind component, not coupled evolution',
 cases=rows,figure4= dict(points=len(comparisons),passed=sum(x['passed'] for x in comparisons),comparisons=comparisons),
 v1_validated=False)
io.json_write(out/'validation.json',summary)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for row,alpha in enumerate((.4,7.64)):
    for col,fields in enumerate((('velocity_ratio','sound_speed_ratio','mach'),('optical_depth_proxy','u','density_ratio'),('kinetic_luminosity','advected_luminosity','luminosity'))):
        ax=axes[row,col]
        for k,field in enumerate(fields):
            color=f'C{k}';f=profiles[(alpha,.8)];r=f['radius_ratio']
            values=np.array([profiles[(alpha,m)][field] for m in (.75,.8,.85)])
            ax.fill_between(r,values.min(0),values.max(0),color=color,alpha=.15)
            ax.loglog(r,f[field],color=color,label=field.replace('_',' '))
            selected=[d for d in comparisons if float(d['alpha'])==alpha and d['field']==field]
            ax.scatter([float(d['radius_ratio']) for d in selected],[float(d['value']) for d in selected],s=12,marker='x',color=color)
        ax.set(xlim=(1,1e4),xlabel='r / R base',title=f'alpha = {alpha}, E = 5, beta = 100')
        ax.grid(alpha=.2);ax.legend(fontsize=7)
        if col==0:ax.set_ylim(.01,20)
        elif col==1:ax.set_ylim(1e-3,1e4)
        else:ax.set_ylim(.01,20)
fig.suptitle('Fiacconi–Rossi Figure 4: lines = solution, x = digitized figure, bands = stated base-Mach range')
fig.savefig(out/'wind_validation.png',dpi=160)
print(out);print(summary['figure4']['passed'],'/',len(comparisons),'digitized points within uncertainty')
