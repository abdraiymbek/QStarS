#!/usr/bin/env python3
"""Plot saved reference profiles; never confuse a converged branch with validity."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args()
report=json.loads((a.directory/'validation.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for row,color in zip([r for r in report['cases'] if r['case'].endswith('_baseline')],['#2471a3','#c65c20']):
    q=np.load(a.directory/(row['case']+'.npz'));label=f"{row['inputs']['teff']:.0f} K, g={row['inputs']['gravity']:g}"
    axes[0].loglog(q['tau'],q['density'],color=color,label=label)
    thick=q['tau']>=2/3
    axes[1].loglog(q['tau'][thick],q['convective_velocity'][thick]/q['sound_speed'][thick],color=color)
    depth=(q['inward_depth'][thick]-q['inward_depth'][300])/row['radius_for_geometry_screen_cm']
    axes[2].semilogx(q['tau'][thick],100*depth,color=color)
axes[0].set_ylabel('Density (g cm$^{-3}$)');axes[0].legend(fontsize=8)
axes[1].axhline(1,color='black',ls='--',lw=1);axes[1].set_ylim(1e-5,10)
axes[1].set_ylabel('Convective Mach number');axes[1].text(1,1.4,'Subsonic MLT limit',fontsize=8)
axes[2].axhline(1,color='black',ls='--',lw=1);axes[2].set_ylabel('Depth below photosphere / radius (%)')
for ax in axes:ax.set_xlabel('Optical depth');ax.grid(alpha=.2)
fig.suptitle('Converged outer-layer candidates fail physical validity — not stellar evolution',fontsize=11)
fig.savefig(a.directory/'atmosphere_validity.png',dpi=180);plt.close(fig)
print(a.directory/'atmosphere_validity.png')
