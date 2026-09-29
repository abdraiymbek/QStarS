#!/usr/bin/env python3
"""Independent audit of completed native CMFGEN data; failures never reclassified."""
import argparse,json,re
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from qstars import io
from qstars.constants import C
from qstars.cmfgen_backend import read_rvtj,read_spectral,read_cooling,read_hhe_populations,read_atomic_levels,hhe_material,sha256

def audit(directory,output,reference=None,hhe=False):
    p=directory.resolve();output.mkdir(parents=True,exist_ok=False)
    plan=json.loads(Path('benchmarks/cmfgen_validation_plan.json').read_text());t=plan['thresholds'];checks=[];metrics={}
    io.provenance(output,dict(plan=plan,native_directory=str(p),reference=str(reference),hhe=hhe))
    (output/'driver.py').write_bytes(Path(__file__).read_bytes())
    def check(name,value,limit,scope):
        value=float(value);checks.append(dict(name=name,value=value if np.isfinite(value) else None,limit=limit,
            status='PASS' if np.isfinite(value) and value<=limit else 'FAIL',scope=scope))
    s=read_rvtj(p/'RVTJ');spec=read_spectral(p/'QSTARS_SPECTRAL.bin')
    model=(p/'MODEL').read_text()
    coh=re.search(r'(?m)^\s*([TF])\s+\[COH_ES\]',model)
    if not coh or coh[1]!='T':raise ValueError('Energy diagnostic currently requires native coherent electron scattering')
    q=spec.moments(coherent_scattering=True)
    if not np.allclose(q['radius'],s['radius'],rtol=1e-7,atol=0):raise ValueError('Material and spectral grids differ')
    frac=lambda a,b:float(np.max(np.abs(a-b)/np.maximum(np.abs(b),1e-100)))
    for name in ['J','H','K']:
        check('spectral_'+name,frac(q[name],s[name]),t['spectral_moments_relative'],'native printed vs independently integrated spectral moment')
    check('spectral_flux_mean',frac(q['kappa_flux'],s['kappa_flux']),t['spectral_flux_mean_relative'],'independent extinction-force quadrature')
    # The final spectrum precedes any final population/temperature update.
    check('spectral_material_temperature',frac(q['temperature'],s['temperature']),t['iteration_max_fractional_correction'],'state synchronization')
    check('mass_flux_uniformity',float(np.ptp(s['mass_flux'])/abs(np.mean(s['mass_flux']))),1e-6,'prescribed-flow input consistency only')
    cool=read_cooling(p/'GENCOOL',len(s['radius']))
    check('electron_balance',np.max(abs(cool['net_percent']))/100,t['electron_balance_relative'],'electron heating/cooling; not total material energy')
    L=q['luminosity_lab'];exchange=4*np.pi*trapezoid(q['radius']**2*q['energy_to_matter_lab'],q['radius'])
    check('global_radiation_energy',abs(L[-1]-L[0]+exchange)/max(abs(L[0]),abs(L[-1])),t['global_radiation_energy_relative'],
          'independent trapezoid volume four-force vs lab luminosity; mesh error and transfer approximations included')
    corrections=re.findall(r'Maximm changes as returned by SOLVEBA_V13 is\s+([0-9.Ee+\-]+)',(p/'OUTGEN').read_text())
    check('last_nonlinear_correction',float(corrections[-1])/100 if corrections else float('nan'),t['iteration_max_fractional_correction'],'all native variables; full convergence needed')
    metrics.update(radii=len(s['radius']),frequencies=len(spec.frequency),temperature_min=float(s['temperature'].min()),temperature_max=float(s['temperature'].max()),
       luminosity_lab_inner=float(L[0]),luminosity_lab_outer=float(L[-1]),volume_energy_to_matter=float(exchange),
       prescribed_mdot=float(np.mean(s['mass_flux'])),heating_max_abs=float(max(abs(q['heating_comoving']))),
       warnings=(p/'WARNINGS').read_text() if (p/'WARNINGS').exists() else '')
    arrays={**{k:v for k,v in q.items() if isinstance(v,np.ndarray)},**{('material_'+k):v for k,v in s.items()},**cool}
    if hhe:
        pops=read_hhe_populations(p/'RVTJ')
        levels={k:read_atomic_levels(p/n,pops[k].shape[0]) for k,n in [('HI','HI_F_OSCDAT'),('HeI','HeI_F_OSCDAT'),('HeII','He2_F_OSCDAT')]}
        if 'Hminus' in pops:levels['Hminus']=read_atomic_levels(p/'H0_F_OSCDAT',pops['Hminus'].shape[0])
        mat=hhe_material(s,pops,levels)
        for k,v in mat['conservation_relative'].items():check('population_'+k,v,t['population_charge_relative'] if k=='charge' else t['population_nuclei_relative'],'atomic conservation on native output')
        metrics['energy_gauge']=mat['energy_gauge'];metrics['material_derivatives_available']=False
        metrics['neutral_hydrogen_fraction_range']=[float(x) for x in [min(pops['HI'].sum(axis=0)/s['hydrogen_density']),max(pops['HI'].sum(axis=0)/s['hydrogen_density'])]]
        arrays.update({k:v for k,v in mat.items() if isinstance(v,np.ndarray)})
        arrays.update({('population_'+k):v for k,v in pops.items()})
    if reference:
        ref=read_rvtj(reference/'RVTJ')
        # No fitted radius, luminosity or temperature normalization.
        for name,key in [('temperature','reference_temperature_relative'),('electron_density','reference_electron_density_relative'),('H','reference_flux_moment_relative'),('kappa_flux','reference_force_relative')]:
            target=np.interp(s['radius'],ref['radius'],ref[name])
            check('reference_'+name,frac(s[name],target),t[key],'supplied reference across release versions; not quasi-star validation')
            arrays['reference_'+name]=target
        check('reference_acceleration',frac(s['kappa_flux']*4*np.pi*s['H']/C,ref['kappa_flux']*4*np.pi*ref['H']/C),t['reference_force_relative'],'force profile compatibility with supplied reference')
        check('reference_grid',frac(s['radius'],ref['radius']),1e-7,'unchanged supplied radial grid')
    np.savez_compressed(output/'profiles.npz',**arrays)
    result=dict(scope='component only; prescribed flow, no wind eigenvalue',backend='CMFGEN 18-Jun-2025',checks=checks,metrics=metrics,
       all_requested_component_gates_pass=False,production_certified=False,wind_eigenvalue_solved=False,
       unperformed_in_this_run=['radial/frequency convergence','export invariance','cool recombination regime','full coupled energy/population response'],
       this_run_checks_pass=all(c['status']=='PASS' for c in checks),
       passed=sum(c['status']=='PASS' for c in checks),failed=sum(c['status']=='FAIL' for c in checks),
       input_hashes={n:sha256(p/n) for n in ['RVTJ','MODEL','GENCOOL','OUTGEN','QSTARS_SPECTRAL.bin']},
       plan_sha256=sha256('benchmarks/cmfgen_validation_plan.json'))
    result['source_snapshot_verified']=io.finish_provenance(output)
    if not result['source_snapshot_verified']:raise RuntimeError('Analysis package changed while running')
    (output/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
    return result
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('directory',type=Path);a.add_argument('output',type=Path)
    a.add_argument('--reference',type=Path);a.add_argument('--hhe',action='store_true');audit(**vars(a.parse_args()))
