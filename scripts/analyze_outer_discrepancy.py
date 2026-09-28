#!/usr/bin/env python3
"""Read-only campaign analysis; original measurements and labels are retained."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('discrepancy',type=Path);p.add_argument('transport',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
raw=a.discrepancy/'cases.json';cases=json.loads(raw.read_text())
for c in cases:
    if 'endpoint_rsun' not in c and 'photo_rsun' in c:c['endpoint_rsun']=c['photo_rsun']
    if 'endpoint_rsun' in c:
        c.pop('photo_rsun',None)
        event=c.get('terminations',{}).get('outward',c.get('termination'))
        if event=='paper_photosphere':c['photo_rsun']=c['endpoint_rsun']
        c['endpoint_kind']=event
# An initial field was named photo_rsun even on stalled cases. Correct only a
# derived copy. Do not silently rewrite the raw campaign or its frozen driver.
(a.output/'normalized_cases.json').write_text(json.dumps(dict(source=str(raw),source_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),correction='photo_rsun restricted to actual paper_photosphere events; endpoint_rsun covers all endpoints',cases=cases),indent=2)+'\n')
data=list(csv.DictReader(open('benchmarks/nakauchi2017_figure2.csv')));by={c['label']:c for c in cases}
comparison=[]
for d in data:
    f=np.load(a.transport/f"paper_{float(d['sonic_density'])}.npz")
    x=np.log(float(d['radius_rsun'])*6.957e10);dx=float(d['log_r_uncertainty']);y=np.log(float(d['value']));dy=float(d['log_y_uncertainty'])
    axis=np.log(f['radius']);lo=max(x-dx,axis[0]);hi=min(x+dx,axis[-1]);row=dict(d)
    if lo>hi:row['status']='outside_model_coverage'
    else:
        pred=PchipInterpolator(axis,np.log(f[d['field']]))(np.linspace(lo,hi,9));distance=max(0.,min(pred)-y,y-max(pred))
        row.update(status='compared',log_distance_outside_radius_band=float(distance),passed=bool(distance<=dy))
    comparison.append(row)
(a.output/'published_profile_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
reddata=[d for d in data if float(d['sonic_density'])==9.71e-9 and d['field']=='temperature']
end=max(reddata,key=lambda d:float(d['radius_rsun']));rend=float(end['radius_rsun']);dr=float(end['log_r_uncertainty'])
comp=json.loads((a.discrepancy/'comparisons.json').read_text())
num=[c for c in comp if c['label'] in ['tight','radial_step','radau','sonic_offset','critical_derivative']]
rounds=[c for c in cases if c['label'].startswith('round_')]
summary=dict(published_vector_endpoint_rsun=rend,published_radius_stroke_band_rsun=[rend*np.exp(-dr),rend*np.exp(dr)],
    published_vector_endpoint_temperature_K=float(end['value']),
    red_model_endpoint_rsun=by['paper_gn93']['endpoint_rsun'],
    red_model_max_profile_log_refinement=max(max(c['max_profile_log_errors'].values()) for c in num),
    red_model_max_endpoint_relative_refinement=max(abs(c['endpoint_relative_change']) for c in num),
    density_rounding_photo_rsun=[by[f'round_axis2_{s:+}']['photo_rsun'] for s in [-1,1]],
    conditional_rounding_samples=len(rounds),
    conditional_rounding_photospheres=sum(c['endpoint_kind']=='paper_photosphere' for c in rounds),
    conditional_rounding_stalled=sum(c['endpoint_kind']=='approaching_second_sonic_point' for c in rounds),
    conditional_rounding_escape_criterion_passed=sum(c['endpoint_kind']=='paper_photosphere' and c['endpoint_velocity_over_escape']>1 for c in rounds),
    conditional_rounding_is_not_a_bound=True,
    published_profile_counts=dict(compared=sum(c['status']=='compared' for c in comparison),passed=sum(c.get('passed',False) for c in comparison),uncovered=sum(c['status']=='outside_model_coverage' for c in comparison)),
    outcome='OUTER FAIL',reproduction=False,causal_explanation_complete=False)
(a.output/'analysis.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for ax,key,ylabel in [(axes[0,0],'velocity','Velocity [km/s]'),(axes[0,1],'temperature','Temperature [K]')]:
    ds=sorted([d for d in data if float(d['sonic_density'])==9.71e-9 and d['field']==key],key=lambda d:float(d['radius_rsun']))
    scale=1e5 if key=='velocity' else 1.
    ax.loglog([float(d['radius_rsun']) for d in ds],[float(d['value'])/scale for d in ds],color='tab:red',lw=2,label='Published red branch')
    for label,title,color in [('default','Default MESA','tab:blue'),('paper_gn93','OPAL / AF94','tab:orange'),('mass_number_diagnostic','Mass-number experiment','tab:green')]:
        f=np.load(a.discrepancy/(label+'.npz'));ax.loglog(f['radius']/6.957e10,f[key]/scale,label=title,color=color,lw=1.4)
    ax.set_xlabel('Radius [solar radii]');ax.set_ylabel(ylabel);ax.grid(alpha=.2)
axes[0,0].legend(fontsize=8);axes[0,1].set_ylim(5000,1e6)
ax=axes[1,0];ec=[c for c in cases if c['label'].startswith('euler_')]
ax.semilogx([c['epsilon'] for c in ec],[c['endpoint_rsun'] for c in ec],'o-',label='Physical-variable Euler')
ax.axhline(by['paper_gn93']['endpoint_rsun'],color='tab:orange',label='Converged DOP853 / Radau')
ax.axhspan(*summary['published_radius_stroke_band_rsun'],color='tab:red',alpha=.2,label='Published stroke band')
ax.set_xlabel('Euler fractional step');ax.set_ylabel('Comparison photosphere [solar radii]');ax.legend(fontsize=8);ax.grid(alpha=.2)
ax=axes[1,1]
rr=[by[f'round_axis{k}_{s:+}'] for k in range(4) for s in [-1,1]]
for j,c in enumerate(rr):
    photo=c['endpoint_kind']=='paper_photosphere'
    ax.scatter(j,c['endpoint_rsun'],marker='o' if photo else 'x',color='tab:blue' if photo else 'black',s=55)
ax.scatter([],[],marker='o',color='tab:blue',label='Comparison photosphere');ax.scatter([],[],marker='x',color='black',label='Stalled cutoff')
ax.axhspan(*summary['published_radius_stroke_band_rsun'],color='tab:red',alpha=.2)
ax.axhline(by['paper_gn93']['endpoint_rsun'],color='tab:orange',ls=':')
ax.set_xticks(range(8),['M−','M+','r−','r+','ρ−','ρ+','T−','T+']);ax.set_ylabel('Endpoint [solar radii]')
ax.set_xlabel('Conditional last-digit perturbation; not fitted inputs');ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
fig.suptitle('Unresolved stationary benchmark: controlled changes do not reproduce the escaping branch\nDiffusion photospheres are comparison endpoints, not validated radiation boundaries',fontsize=12)
fig.savefig(a.output/'outer_discrepancy.png',dpi=170);plt.close(fig)
print(json.dumps(summary,indent=2))
