#!/usr/bin/env python3
"""Compare completed thermal iterates, including material and atomic profiles."""
import argparse,json
from pathlib import Path
import numpy as np

def compare(root):
    rows=[]
    for left,right in [('hhe_3','thermal_tight'),('thermal_rel','thermal_equilibrium')]:
        if not (root/right/'profiles.npz').exists():continue
        a=np.load(root/left/'profiles.npz');b=np.load(root/right/'profiles.npz')
        np.testing.assert_allclose(a['radius'],b['radius'],rtol=1e-12,atol=0)
        relative={}
        for key in ['material_temperature','material_electron_density','gas_internal_energy',
                    'E_comoving','F_comoving','P_comoving','radiative_acceleration','kappa_flux']:
            relative[key]=float(np.max(abs(a[key]-b[key])/np.maximum(abs(b[key]),1e-100)))
        populations={}
        for ion in ['HI','HII','HeI','HeII','HeIII']:
            norm=b['material_hydrogen_density' if ion in ['HI','HII'] else 'material_helium_density']
            populations[ion]=float(np.max(abs(a['population_'+ion]-b['population_'+ion])/norm))
        rows.append(dict(left=left,right=right,relative_profile_differences=relative,
            population_absolute_change_per_nucleus=populations,
            profile_limit=.01,profile_status='PASS' if max(relative.values())<=.01 else 'FAIL',
            population_limit=.01,population_status='PASS' if max(populations.values())<=.01 else 'FAIL',
            scope='same prescribed flow and transfer equations; iteration tolerance only. Per-nucleus population norm does not bound relative errors in trace levels.'))
    result=dict(comparisons=rows,production_certified=False)
    (root/'thermal_profiles.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);compare(**vars(p.parse_args()))
