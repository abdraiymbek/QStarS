from qstars import io
from qstars.spherical_wind import *
from qstars.constants import MSUN,RSUN,YEAR
from qstars.freeeos_reference import FreeEOSReference
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.eos import Composition
import time,traceback
out=io.create_directory('results','spherical_wind_probe');io.provenance(out,{'source':'Nakauchi2017 Figure2; prescribed sonic states; independent component'})
print(out,flush=True)
mesa=MesaMicrophysics(Composition(.7,.3,0),'backends/mesa/build_atmosphere_3/libqstars_mesa.dylib','/Users/abraiymbek/qsrd-runtime/mesa')
eos=FreeEOSReference('backends/freeeos/build_bridge_2/libqstars_freeeos.dylib',.7,mesa.radiation_constant)
rows=[]
for label,micro in [('ideal',IdealGasOpacity()),('freeeos',FreeEOSOpacity(eos,mesa))]:
 for rho in [9.65e-9,9.71e-9,9.85e-9,10.5e-9]:
  t=time.perf_counter();name=f'{label}_{rho}'
  try:
   w=DiffusionWind(micro,100*MSUN,50*RSUN,rho,2.5e5)
   f,s=w.integrate(inner_ratio=.5,outer_ratio=100,endpoint='paper_photosphere',samples=1201,rtol=1e-7)
   np.savez_compressed(out/(name+'.npz'),**f);s['mdot_msun_year']=w.mdot/MSUN*YEAR
  except Exception as e:s={'status':'numerical_or_domain_failure','reason':str(e),'traceback':traceback.format_exc()}
  s.update(name=name,seconds=time.perf_counter()-t);rows.append(s);io.json_write(out/'cases.json',rows)
  print(name, {k:s[k] for k in s if k in ['reason','terminations','endpoint_temperature','endpoint_radius','shell_mass_fraction','max_knudsen','endpoint_velocity_over_escape','mdot_msun_year','seconds','first_law_integral_relative_defect']},flush=True)
io.json_write(out/'provenance_check.json',{'source_snapshot_verified':io.finish_provenance(out)})
