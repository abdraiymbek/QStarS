#!/usr/bin/env python3
"""Quantify an outer-boundary candidate, including explicit validity failures."""
import argparse,dataclasses,json,time,resource
from pathlib import Path
import numpy as np
from qstars import io
from qstars.constants import G,PI,C,MSUN
from qstars.eos import Composition
from qstars.freeeos_reference import FreeEOSReference
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.atmosphere_reference import integrate_atmosphere
p=argparse.ArgumentParser();p.add_argument('--freeeos-library',required=True);p.add_argument('--mesa-library',required=True);p.add_argument('--mesa-dir',required=True);a=p.parse_args()
out=io.create_directory('results','atmosphere_reference');io.provenance(out,vars(a))
m=MesaMicrophysics(Composition(X=.7,Y=.3),a.mesa_library,a.mesa_dir)
f=FreeEOSReference(a.freeeos_library,radiation_constant=m.radiation_constant)
rows=[];profiles={};start=time.perf_counter()
for teff,g in ((8000.,3.),(6000.,1.)):
 for variant,kwargs in [('baseline',{}),('tolerance',{'rtol':1e-9}),('outer_tau',{'tau_outer':1e-7}),
                        ('outer_pressure',{'gas_pressure_outer':1e-7}),('finite_pressure',{'gas_pressure_outer':1.}),
                        ('ML1',{'mlt_model':'ML1'}),
                        ('constant_opacity',{'constant_opacity':.34})]:
  name=f'T{int(teff)}_g{g}_{variant}';t0=time.perf_counter()
  inputs=dict(teff=teff,gravity=g,gas_pressure_outer=1e-6,tau_outer=1e-6,rtol=1e-7,points=301);inputs.update(kwargs)
  mass=1e4*MSUN;radius=np.sqrt(G*mass/g)
  row=dict(case=name,inputs=inputs,enclosed_mass_msun=1e4,radius_for_geometry_screen_cm=radius)
  try:
   q=integrate_atmosphere(f,m,**inputs);profiles[name]=q
   np.savez_compressed(out/f'{name}.npz',**dataclasses.asdict(q))
   area=4*PI*radius**2
   mach=q.convective_velocity/q.sound_speed
   hydro=(q.pressure-q.pressure[0])-g*q.column_mass
   row.update(status='integrated_reference',photosphere={key:float(getattr(q,key)[300]) for key in ('gas_pressure','density','opacity','pressure','temperature')},
              max_convective_mach=float(max(mach)),max_total_flux_eddington_factor=float(max(q.total_flux_eddington_factor)),
              max_radiative_eddington_factor=float(max(q.radiative_eddington_factor)),
              tau_at_peak_convective_mach=float(q.tau[np.argmax(mach)]),
              blob_optical_depth_at_peak_mach=float(2*q.opacity[np.argmax(mach)]*q.pressure[np.argmax(mach)]/g),
              total_depth_over_radius=float(q.inward_depth[-1]/radius),
              photosphere_to_base_depth_over_radius=float((q.inward_depth[-1]-q.inward_depth[300])/radius),
              column_mass_g_cm2=float(q.column_mass[-1]),layer_mass_fraction=float(area*q.column_mass[-1]/mass),
              column_energy_erg_cm2=float(q.column_energy[-1]),
              hydrostatic_integral_relative_residual=float(max(abs(hydro))/max(q.pressure)),
              supersonic_mlt=bool(max(mach)>=1),thin_geometry_one_percent=bool(q.inward_depth[-1]/radius<=.01))
   if variant=='constant_opacity':
    tau=q.tau[:301];F=f.radiation_constant*C*teff**4/4
    expected=inputs['gas_pressure_outer']+(tau-inputs['tau_outer'])*(g/.34-F/C)
    row['analytic_grey_gas_pressure_relative_max']=float(max(abs(q.gas_pressure[:301]-expected))/max(expected))
   if variant in ('tolerance','outer_tau','outer_pressure'):
    b=profiles[f'T{int(teff)}_g{g}_baseline'];common=np.geomspace(max(q.tau[0],b.tau[0]),min(q.tau[-1],b.tau[-1]),3001)
    row['max_log_profile_difference_vs_baseline']={key:float(max(abs(np.interp(np.log(common),np.log(q.tau),np.log(getattr(q,key)))-np.interp(np.log(common),np.log(b.tau),np.log(getattr(b,key)))))) for key in ('gas_pressure','temperature','density','opacity')}
    row['photosphere_relative_difference_vs_baseline']={key:float(getattr(q,key)[300]/getattr(b,key)[300]-1) for key in ('gas_pressure','temperature','density','opacity')}
  except (ValueError,RuntimeError) as error:row.update(status='failed_reference',error=str(error))
  row['runtime_seconds']=time.perf_counter()-t0;rows.append(row);print(json.dumps(row),flush=True)
report=dict(scope='Uncoupled candidate atmosphere; validity failures prevent promotion to production.',
            cases=rows,runtime_seconds=time.perf_counter()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            mesa_build=m.manifest,freeeos_build=f.manifest,production_validated=False)
report['source_snapshot_verified']=io.finish_provenance(out);io.json_write(out/'validation.json',report);print(out)
