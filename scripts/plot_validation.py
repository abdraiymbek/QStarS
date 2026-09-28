"""Plot saved campaign data; no solver calls."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(sys.argv[1]);data=json.loads(p.read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
for ax,axis,key,xlabel,ykey,ylabel in [
 (axes[0],'mesh','zones','Zones','internal','Internal energy difference / finest'),
 (axes[1],'timestep','dt_years','Maximum timestep [yr]','luminosity_erg_s','Luminosity difference / finest')]:
 rows=[r for r in data['convergence'] if r['axis']==axis and r.get('comparable',True)]
 points=[(r[key],r.get('relative_difference_to_finest',{}).get(ykey,0)) for r in rows]
 points=[(x,y) for x,y in points if y>0]
 if points:ax.loglog(*np.array(points).T,'o-',color='#146a79')
 ax.set(xlabel=xlabel,ylabel=ylabel)
rows=[r for r in data['convergence'] if r['axis']=='solver']
axes[2].loglog([r['solver_tolerance'] for r in rows],[abs(r['energy_residual_relative']) for r in rows],'o-',color='#784bb2')
axes[2].set(xlabel='Nonlinear tolerance',ylabel='Absolute energy residual / initial U')
fig.suptitle('Measured convergence · zero-pressure hot-envelope test',fontsize=15)
fig.savefig(p.parent/'convergence.png',dpi=160);plt.close(fig)
# Published scalars are explicitly structure-level comparisons.
fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
rows=[x for x in data['checks'] if x['name'].startswith('Ball 2011 n=')]
for i,row in enumerate(rows):
 axes[0].plot(i,row['actual'],'o',color='#146a79',label='Integrated ODE' if i==0 else None)
 axes[0].errorbar(i+.08,row['reference'],yerr=row['tolerance'],fmt='x',color='#df8832',label='Published rounding interval' if i==0 else None)
axes[0].set(xticks=range(3),xticklabels=['n = 2','n = 3','n = 4'],ylabel='Maximum BH / total mass',title='Ball loaded polytropes')
axes[0].legend(fontsize=8)
curve=json.loads((p.parent/'saturated_sequence.json').read_text())
axes[1].plot([x['Ki'] for x in curve],[x['bh_fraction'] for x in curve],color='#146a79')
axes[1].plot(.4655453,.617828,'o',color='#df8832',label='Printed cutoff interpretation')
axes[1].set(xlabel='Dimensionless inner parameter Ki',ylabel='BH / total mass',title='Saturated convection + adiabatic envelope')
axes[1].legend(fontsize=8)
fig.suptitle('Published structure equations · not stellar evolutionary tracks',fontsize=14)
fig.savefig(p.parent/'structure_benchmarks.png',dpi=160);plt.close(fig)
