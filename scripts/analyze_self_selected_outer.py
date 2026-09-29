#!/usr/bin/env python3
"""Summarize frozen Report 12 components; no physical wind promotion."""
import json,argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qstars.thermal_transfer import opacity_means
p=argparse.ArgumentParser();p.add_argument('--main',type=Path,required=True);p.add_argument('--refinement',type=Path,required=True)
p.add_argument('--out',type=Path,default=Path('Reports/evidence/self_selected_20260928/analysis'));a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
read=lambda p:json.loads(p.read_text());load=lambda p:dict(np.load(p))
m,r=a.main,a.refinement;main=read(m/'campaign.json');ref=read(r/'campaign.json');checks=read(m/'checks.json')+read(r/'checks.json')
thermal=load(r/'cluster1601.npz');rad=thermal['radius'];ka=np.array([.2,2.])[:,None]*np.ones((2,len(rad)))
means=opacity_means(ka,ka*.2,thermal['gas_temperature'],[0,20000,np.inf],thermal['J'],thermal['H'])
adiabatic=load(r/'dense4001.npz');ad=read(r/'dense4001.json');critical=read(m/'native_critical_states.json')
recovery=read(r/'rejected_trial_recovery.json')
subset=lambda prefix:max(x['value'] for x in checks if x['name'].startswith(prefix))
report=dict(outcome='FAIL',main_campaign=main,refinement_campaign=ref,
 checks=len(checks),passing=sum(x['passed'] for x in checks),failed=[x for x in checks if not x['passed']],
 fine_thermal=dict(energy_defect=read(r/'cluster1601.json')['energy_integral_relative_defect'],
                  max_profile_difference_801_1601=subset('cluster801_vs_cluster1601'),
                  exterior_difference_60_180=subset('cluster_outer60_vs_cluster_outer180'),
                  exterior180_energy_defect=read(r/'cluster_outer180.json')['energy_integral_relative_defect'],
                  Tgas_outer=float(thermal['gas_temperature'][-1]),Trad_outer=float(thermal['radiation_temperature'][-1]),
                  max_flux_to_rosseland=float(max(means['flux_total']/means['rosseland_total']))),
 adiabatic=dict(mdot=ad['mdot'],sonic_radius=ad['sonic_radius'],base_velocity=ad['base_velocity'],terminal_speed=ad['terminal_speed'],
                profile_difference_2001_4001=subset('dense2001_vs_dense4001'),exterior_profile_difference=subset('dense_exterior_vs_dense4001')),
 native_critical_speed_ratio_range=[min(x['adiabatic_to_isothermal_speed_ratio'] for x in critical),max(x['adiabatic_to_isothermal_speed_ratio'] for x in critical)],
 recovered_mass_flux_error=recovery['mdot_relative_error'],rejected_trials=recovery['rejected_trials'],
 max_native_slice_mass_flux_error=max(read(f)['mdot_relative_error'] for f in m.glob('diffusion_match_*.json')))
(a.out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
fig,ax=plt.subplots(2,3,figsize=(14,8.4),layout='constrained')
ax[0,0].loglog(rad,thermal['gas_temperature'],label='Material temperature')
ax[0,0].loglog(rad,thermal['radiation_temperature'],label='Radiation energy temperature')
ax[0,0].set(title='Frozen two-band thermal equilibrium',xlabel='r / core radius',ylabel='Temperature [K]');ax[0,0].legend(fontsize=9)
for key,label in [('planck_absorption','Planck absorption'),('energy_absorption','Energy absorption'),('rosseland_total','Rosseland total'),('flux_total','Flux total')]:ax[0,1].semilogx(rad,means[key],label=label)
ax[0,1].set(title='Distinct opacity means',xlabel='r / core radius',ylabel='Specified-opacity units');ax[0,1].legend(fontsize=9)
for names,path,label in [(['thermal_'+str(n) for n in [101,201,401,801]],m,'Uniform log radius'),(['cluster'+str(n) for n in [401,801,1601]],r,'Core-clustered radius')]:
 n=[int(x.split('_')[-1]) if '_' in x else int(x[7:]) for x in names]
 errors=[read(path/(name+'.json'))['energy_integral_relative_defect'] for name in names]
 ax[0,2].loglog(n,errors,'o-',label=label)
ax[0,2].axhline(.0005,color='black',ls='--',label='Energy target')
ax[0,2].set(title='Independent radiation energy defect',xlabel='Radial points',ylabel='Relative defect');ax[0,2].legend(fontsize=9)
x=adiabatic['radius']/adiabatic['radius'][0]
ax[1,0].loglog(x,adiabatic['velocity']/ad['sonic_speed'],label='v / sonic speed')
ax[1,0].loglog(x,adiabatic['density']/adiabatic['density'][0],label='Density / base')
t=adiabatic['pressure']/adiabatic['density'];ax[1,0].loglog(x,t/t[0],label='Temperature / base')
ax[1,0].axvline(ad['sonic_radius']/adiabatic['radius'][0],color='grey',ls=':')
ax[1,0].set(title='Adiabatic eigenvalue control (no radiation)',xlabel='r / base radius',ylabel='Normalized profiles');ax[1,0].legend(fontsize=9)
values=[max(abs(np.array(i['residual']))) for i in recovery['iterations']]
ax[1,1].semilogy(np.arange(1,len(values)+1),values,'o-',ms=3)
ax[1,1].axhline(1e-8,color='black',ls='--')
ax[1,1].set(title=f'Diffusion root: {recovery["rejected_trials"]} rejected trials recovered',xlabel='Accepted Newton step',ylabel='Largest logarithmic residual')
ns=[16,32,64,128];error=[read(m/f'frequency_{n}.json')['temperature_relative_error'] for n in ns]
ax[1,2].loglog(ns,error,'o-');ax[1,2].axhline(.001,color='black',ls='--',label='Thermal target')
ax[1,2].set(title='Independent thin-opacity limit: κν ∝ ν',xlabel='Frequency bands',ylabel='Temperature error');ax[1,2].legend(fontsize=9)
ax[1,2].set_xticks(ns,labels=[str(n) for n in ns]);ax[1,2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for panel in ax.ravel():panel.grid(alpha=.2)
fig.suptitle('QStarS Report 12 — component validation; physical outer closure FAIL',fontsize=16)
fig.savefig(a.out/'self_selected_validation.png',dpi=160);plt.close(fig)
print(json.dumps({k:v for k,v in report.items() if k not in ['main_campaign','refinement_campaign','failed']},indent=2))
