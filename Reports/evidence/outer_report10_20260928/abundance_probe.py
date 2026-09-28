import os,json,numpy as np
from pathlib import Path
from qstars import io
from qstars.constants import MSUN,RSUN
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.freeeos_reference import FreeEOSReference
from qstars.spherical_wind import FreeEOSOpacity,DiffusionWind
class MassNumberAbundances:
    def __init__(self,eos):
        self.hydrogen=.7;self.radiation_constant=eos.radiation_constant
        self.factor=.7*1.007825031898+.3*4.002603254130/4
        self.eos=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7*1.007825031898/self.factor,self.radiation_constant)
    def raw(self,rho,T):
        x=self.eos.raw(rho*self.factor,T).copy()
        x[3:9]*=self.factor;x[9]-=np.log(self.factor)
        x[12]*=self.factor;x[15]*=self.factor;x[17]/=self.factor
        return x
out=io.create_directory('results','abundance_convention_probe');io.provenance(out,dict(experiment='mass-number abundance normalization; not claimed original paper EOS'));print(out,flush=True)
m=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
e=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,m.radiation_constant)
rows=[]
for tag,eos in [('isotope',e),('mass_number',MassNumberAbundances(e))]:
  try:
    w=DiffusionWind(FreeEOSOpacity(eos,m,'paper_gn93'),100*MSUN,50*RSUN,9.71e-9,2.5e5)
    f,s=w.integrate(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,samples=501,stalled_mach=1.5)
    np.savez_compressed(out/(tag+'.npz'),**f);rows.append(dict(tag=tag,summary=s))
    print(tag,s['endpoint_radius']/RSUN,s['endpoint_velocity_over_escape'],s['terminations'],flush=True)
  except Exception as exc:rows.append(dict(tag=tag,error=str(exc)));print(tag,str(exc),flush=True)
io.json_write(out/'cases.json',rows);io.finish_provenance(out)
