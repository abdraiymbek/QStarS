#!/usr/bin/env python3
"""Explicit cool H/He controls; prescribed flow, never a wind boundary."""
import argparse,json
from pathlib import Path
import numpy as np
from prepare_hhe import make_hhe
from prepare_reference import keyword
from qstars.cmfgen_backend import sha256


def prepare(source,destination,atomic,case='T9',iterations=40,factor=1,initialization='auto'):
    if factor not in [1,2,4] or initialization not in ['auto','full']:
        raise ValueError('Unsupported mesh or initialization control')
    plan_path=Path(__file__).resolve().parents[2]/'benchmarks/cmfgen_cool_validation_plan.json'
    plan=json.loads(plan_path.read_text());c=plan['cases'][case];flow=plan['prescribed_control_flow']
    make_hhe(source,destination,atomic,iterations)
    # Table radii are in solar radii; use the IAU nominal conversion, explicitly.
    R=c['radius_Rsun']*6.957; vmax=flow['terminal_velocity_km_s'];v0=flow['base_velocity_km_s']
    n=flow['initial_radii'];x=1+np.r_[0,np.geomspace(1e-5,flow['outer_radius_ratio']-1,n-1)]
    if factor>1:
        x=np.exp(np.r_[np.concatenate([np.linspace(a,b,factor,endpoint=False) for a,b in zip(np.log(x[:-1]),np.log(x[1:]))]),np.log(x[-1])])
    n=len(x)
    for key,val in [('ND',n),('NP',n+flow['core_rays'])]:keyword(destination/'MODEL_SPEC',key,str(val))
    v=v0+(vmax-v0)*(1-1/x);sig=(vmax-v0)/x/v-1
    (destination/'RVSIG_COL').write_text(f'{n} !Number of depth points\n'+''.join(
        f'{r:.14e} {vv:.14e} {ss:.14e}\n' for r,vv,ss in zip((R*x)[::-1],v[::-1],sig[::-1])))
    for k,val in {'RSTAR':R,'RMAX':x[-1],'VINF':vmax,'BETA':1,'MDOT':c['mdot_Msun_yr'],
        'LSTAR':c['luminosity_Lsun'],'TEFF':c['inner_flux_temperature_K']/1e4,
        'MASS':1e6,'LOGG':np.log10(6.67430e-8*1e6*1.98847e33/(R*1e10)**2),
        'LIN_INT':'F','IT_ON_T':'T','DC_METH':'LTE','T_EXC':0.0,'RD_IN_R_GRID':'F',
        'INC_TWO':'T','INC_RAY':'T','USE_J_REL':'T','INCL_REL':'T','INCL_ADV_TRANS':'T',
        'AUTO_ADD':'F','DIE_H0':'F,F','VTURB':50,'EPS_TERM':.0001}.items():keyword(destination/'VADAT',k,str(val))
    # H- is named H0 by the current upstream code (manual and rd_model_spec).
    keyword(destination/'MODEL_SPEC','H0_ISF','1,1,1')
    hmi=atomic.resolve()/'HYD/MI/13apr15'
    extra=[]
    for name,target in {'H0_F_OSCDAT':'hmi_osc','H0_COL_DATA':'col_guess',
                        'PHOTH0_A':'hmi_phot','H0_FF':'hmi_ff'}.items():
        p=hmi/target
        if not p.is_file():raise FileNotFoundError(p)
        (destination/name).symlink_to(p);extra.append(dict(link=name,path=str(p),sha256=sha256(p)))
    manifest=destination/'atomic_manifest.json'
    manifest.write_text(json.dumps(json.loads(manifest.read_text())+extra,indent=2)+'\n')
    # No old hydro table or grey correction is used to claim a cool equilibrium.
    # Retain copied initial guesses, but disable importing an O-star grey correction.
    for name in ['GREY_SCL_FAC_IN','GREY_SCL_FACIN']:
        p=destination/name
        if p.exists():p.rename(destination/(name+'.unused_hot_guess'))
    (destination/'QSTARS_EXPORT').touch()
    if initialization=='full':
        for k in ['DO_LAM_IT','DO_LAM_AUTO','DO_T_AUTO']:keyword(destination/'IN_ITS',k,'F')
        (destination/'IT_SPECIFIER').write_text('FULL 1\n')
    definition=json.loads((destination/'component_definition.json').read_text())
    definition.update(input_flow='explicit cool beta-one control; see cool_control.json and RVSIG_COL',
        atomic_levels={'Hminus':1,'HI':30,'HeI':69,'HeII':30},
        two_photon=True,rayleigh=True,coherent_electron_scattering=True,
        thermal_solution_required=True,published_model_reproduction=False)
    (destination/'component_definition.json').write_text(json.dumps(definition,indent=2)+'\n')
    (destination/'cool_control.json').write_text(json.dumps(dict(
        plan=plan,plan_sha256=sha256(plan_path),case=case,factor=factor,initialization=initialization,
        radius_conversion_cm=6.957e10,mass_role='unused prescribed-flow bookkeeping, not a hydrostatic envelope',
        initial_population_estimate='native LTE with native grey temperature initialization',
        hminus_collision_warning='official col_guess: guessed collisional ionization, not independently validated',
        production=False,published_model_reproduction=False),indent=2)+'\n')
    return destination


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('destination',type=Path)
    p.add_argument('--atomic',type=Path,required=True);p.add_argument('--case',choices=['T9','T6'],default='T9')
    p.add_argument('--factor',type=int,choices=[1,2,4],default=1);p.add_argument('--initialization',choices=['auto','full'],default='auto')
    p.add_argument('--iterations',type=int,default=40);print(prepare(**vars(p.parse_args())))
