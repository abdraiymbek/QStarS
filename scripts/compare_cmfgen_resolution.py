#!/usr/bin/env python3
"""Compare common-radius internal radiation fields; retain initialization defects."""
import argparse,json
from pathlib import Path
import numpy as np

def compare(root):
    comparisons=[]
    pairs=[('transport_rel','transport_frequency','frequency'),('transport_rel','mesh_84','restart_initialization'),
           ('mesh_84','mesh_167','mesh'),('mesh_167','mesh_333','mesh'),('mesh_333','mesh_665_fixed','mesh')]
    for left,right,kind in pairs:
        if not (root/left/'profiles.npz').exists() or not (root/right/'profiles.npz').exists():continue
        a=np.load(root/left/'profiles.npz');b=np.load(root/right/'profiles.npz');values={}
        for key in ['temperature','density','E_comoving','F_comoving','P_comoving','radiative_acceleration','kappa_flux']:
            target=np.interp(a['radius'],b['radius'],b[key])
            values[key]=float(np.max(abs(a[key]-target)/np.maximum(abs(target),1e-100)))
        vals=[values[k] for k in ['E_comoving','F_comoving','P_comoving','radiative_acceleration','kappa_flux']]
        comparisons.append(dict(left=left,right=right,kind=kind,values=values,radiation_profile_limit=.01,
            radiation_profile_status='PASS' if max(vals)<=.01 else 'FAIL',material_state_limit=.001,
            material_state_status='PASS' if max(values['temperature'],values['density'])<=.001 else 'FAIL',
            scope='incoming T/rho and radiation operator only; final Lambda-corrected populations are not equilibria'))
    # Tiny incoming-state differences can cause large net-source cancellation errors.
    # Do not infer a thermal/eigenvalue convergence certificate from these comparisons.
    result=dict(comparisons=comparisons,thermal_convergence_certified=False,wind_certified=False)
    (root/'resolution.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);compare(**vars(p.parse_args()))
