#!/usr/bin/env python3
"""Cool-state coverage, extended material derivatives and independent ledgers."""
import argparse,json
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.integrate import cumulative_trapezoid
from validate_cmfgen_component import audit
from qstars.cmfgen_backend import read_rvtj,read_hhe_populations,read_atomic_levels,read_spectral,sha256
from qstars.nlte_material import hhe_resolved_material

def validate(directory,output):
    result=audit(directory,output,hhe=True)
    plan_path=Path('benchmarks/cmfgen_cool_validation_plan.json');plan=json.loads(plan_path.read_text());th=plan['thresholds']
    s=read_rvtj(directory/'RVTJ');p=read_hhe_populations(directory/'RVTJ')
    levels={ion:read_atomic_levels(directory/name,p[ion].shape[0]) for ion,name in [
        ('HI','HI_F_OSCDAT'),('HeI','HeI_F_OSCDAT'),('HeII','He2_F_OSCDAT'),('Hminus','H0_F_OSCDAT')] if ion in p}
    m=hhe_resolved_material(s,p,levels);q=m.values();rad=read_spectral(directory/'QSTARS_SPECTRAL.bin').moments(coherent_scattering=True)
    neutral=p['HI'].sum(0)/s['hydrogen_density'];T=s['temperature'];r=s['radius'];rho=s['density']
    mask=(T>=th['cool_temperature_min_K'])&(T<=th['cool_temperature_max_K'])&(neutral>=th['recombination_neutral_fraction_min'])&(neutral<=th['recombination_neutral_fraction_max'])
    eps=th['finite_difference_step_ln'];dy=np.zeros_like(m.number_per_gram)
    ih=m.labels.index('HI:0');ip=m.labels.index('HII:0')
    d=.1*np.minimum(m.number_per_gram[ih],m.number_per_gram[ip]);dy[ih]=-d;dy[ip]=d
    plus=replace(m,density=m.density*np.exp(eps*.2),temperature=T*np.exp(-eps*.3),number_per_gram=m.number_per_gram+eps*dy).values()
    minus=replace(m,density=m.density*np.exp(-eps*.2),temperature=T*np.exp(eps*.3),number_per_gram=m.number_per_gram-eps*dy).values()
    direction=m.differential(.2,-.3,dy)
    errors={k:float(np.max(abs((plus[k]-minus[k])/(2*eps)-direction[k])/np.maximum(abs(q[k]),1e-100))) for k in ['pressure','internal_energy']}
    # Radial material first law as an integral, not a rescaled Bernoulli identity.
    mdot=float(np.mean(s['mass_flux']));L=rad['luminosity_lab'];Lscale=max(abs(L[0]),abs(L[-1]))
    work=mdot*cumulative_trapezoid(q['pressure']/rho,np.log(rho),initial=0)
    heat=cumulative_trapezoid(4*np.pi*r*r*rad['heating_comoving'],r,initial=0)
    firstlaw=mdot*(q['internal_energy']-q['internal_energy'][0])-work-heat
    atomic_flux=mdot*(q['atomic_energy']-q['atomic_energy'][0])
    radexchange=cumulative_trapezoid(4*np.pi*r*r*rad['energy_to_matter_lab'],r,initial=0)
    raderror=(L-L[0]+radexchange)/Lscale
    initial=dict(np.load(output/'profiles.npz'))
    initial.update(cool_recombination_mask=mask,neutral_hydrogen_fraction=neutral,
       hminus_fraction=p.get('Hminus',np.zeros((1,len(r)))).sum(0)/s['hydrogen_density'],
       material_firstlaw_cumulative_fraction=firstlaw/Lscale,atomic_advected_energy_fraction=atomic_flux/Lscale,
       radiation_cumulative_fraction=raderror,gas_pressure_dlnT_frozen=q['pressure_dlnT'],
       gas_energy_dlnT_frozen=q['energy_dlnT'])
    np.savez_compressed(output/'cool_profiles.npz',**initial)
    run=json.loads((directory/'run_result.json').read_text())
    cap=Path(run['executable']).with_suffix('.capabilities.json')
    capability=json.loads(cap.read_text()) if cap.exists() else {}
    ray_fixed=capability.get('executable_sha256')==run['executable_sha256'] and capability.get('hminus_rayleigh_fixed') is True
    checks=[dict(name='cool_partial_recombination_coverage',value=int(mask.sum()),minimum=th['minimum_recombining_depths'],status='PASS' if mask.sum()>=th['minimum_recombining_depths'] else 'FAIL'),
       dict(name='native_Hminus_Rayleigh_repair',value=ray_fixed,status='PASS' if ray_fixed else 'FAIL'),
       dict(name='extended_material_directional_derivatives',value=max(errors.values()),limit=th['material_derivative_relative'],status='PASS' if max(errors.values())<th['material_derivative_relative'] else 'FAIL'),
       dict(name='cumulative_radiation_energy',value=float(max(abs(raderror))),limit=th['cumulative_radiation_energy_relative'],status='PASS' if max(abs(raderror))<th['cumulative_radiation_energy_relative'] else 'FAIL')]
    model=(directory/'MODEL').read_text()
    synchronized=all(c['status']=='PASS' for c in result['checks'] if c['name'] in ['last_nonlinear_correction','spectral_material_temperature'])
    summary=dict(component=result,additional_checks=checks,derivative_normalization='absolute directional error divided by P or u; no relaxed-response claim',
        derivative_errors=errors,cool_depths=int(mask.sum()),hminus_fraction_max=float(max(initial['hminus_fraction'])),
        material_firstlaw_cumulative_fraction_max=float(max(abs(firstlaw))/Lscale),
        material_radiation_synchronized_to_component_tolerance=synchronized,
        material_firstlaw_scope='prescribed flow diagnostic only; final spectrum precedes the final material update. Do not interpret mixed-iteration gas/heat defects physically when nonlinear convergence fails. Radiative-equilibrium controls omit gas/population advection.',
        atomic_advected_energy_fraction_max=float(max(abs(atomic_flux))/Lscale),
        radiation_cumulative_fraction_max=float(max(abs(raderror))),
        equilibrium_response_validated=False,global_wind_eigenvalue_solved=False,
        published_cool_reference_reproduced=False,production_certified=False,
        plan_sha256=sha256(plan_path),driver_sha256=sha256(__file__))
    (output/'cool_validation.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (output/'cool_driver.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k:v for k,v in summary.items() if k not in ['component']},indent=2))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('output',type=Path);validate(**vars(p.parse_args()))
