import os,json,time
from dataclasses import replace
import numpy as np
from qstars.constants import MSUN,RSUN
from qstars import io
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.freeeos_reference import FreeEOSReference
from qstars.spherical_wind import FreeEOSOpacity,DiffusionWind
out=io.create_directory('results','outer_opacity_probe');io.provenance(out,{'probe':'original_opacity_families','no_input_fitting':True});print(out,flush=True)
m=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
e=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,m.radiation_constant)
rows=[]
for family in ['default','paper_gs98','paper_gn93']:
    try:
        u=FreeEOSOpacity(e,m,family);w=DiffusionWind(u,100*MSUN,50*RSUN,9.71e-9,2.5e5)
        f,s=w.integrate(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,samples=501)
        s.update(family=family,photo_rsun=s['endpoint_radius']/RSUN)
        np.savez_compressed(out/(family+'.npz'),**f);rows.append(s)
        print(family,s['photo_rsun'],s['endpoint_velocity_over_escape'],flush=True)
    except Exception as exc:
        rows.append({'family':family,'error':str(exc)});print(family,repr(exc),flush=True)
io.json_write(out/'cases.json',rows);io.finish_provenance(out)
