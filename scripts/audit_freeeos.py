#!/usr/bin/env python3
"""Cross-check direct FreeEOS, MESA and gas-only thermodynamic identities."""
import argparse,ctypes as ct,json,time
from pathlib import Path
import numpy as np
from qstars import io
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.freeeos_reference import FreeEOSReference

p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--mesa-library',required=True);p.add_argument('--mesa-dir',required=True);a=p.parse_args()
out=io.create_directory('results','freeeos_audit');io.provenance(out,vars(a));start=time.perf_counter()
mesa=MesaMicrophysics(Composition(X=.7,Y=.3,Z=0),a.mesa_library,a.mesa_dir)
eos=FreeEOSReference(a.library,radiation_constant=mesa.radiation_constant)
rho,T=np.meshgrid(np.logspace(-22,-6,17),[4000.,5000.,6000.,8000.,1e4,2e4,4e4,1e5,1e6],indexing='ij');rho,T=rho.ravel(),T.ravel()
x=eos.raw(rho,T);q=eos.evaluate(rho,T);h=1e-4
defects=[];derivatives=[]
for axis in ('rho','T'):
    xp=eos.raw(rho*np.exp(h) if axis=='rho' else rho,T*np.exp(h) if axis=='T' else T)
    xm=eos.raw(rho*np.exp(-h) if axis=='rho' else rho,T*np.exp(-h) if axis=='T' else T)
    du=(xp[3]-xm[3])/(2*h);ds=(xp[6]-xm[6])/(2*h)
    # Screen the gas alone, so a dominant radiation term cannot hide a defect.
    rhs=du-q['pg']/rho if axis=='rho' else du
    defects.append(abs(T*ds-rhs)/(abs(du)+q['pg']/rho))
    col=4 if axis=='rho' else 5
    derivatives.append(abs(du-x[col])/(abs(du)+q['pg']/rho))
inverse=eos.density_from_gas_pressure(q['pg'],T)/rho-1
m=np.empty((21,rho.size),order='F');dr=np.empty_like(m);dt=np.empty_like(m);err=ct.c_int()
mesa.lib.qstars_eval(rho.size,np.ascontiguousarray(rho),np.ascontiguousarray(T),.7,m,dr,dt,ct.byref(err))
if err.value:raise ValueError(err.value)
k=mesa.opacity_external_eos(rho,T,q['free_e'],q['eta'])
names=('pg','u','s','free_e','mu');mv=(np.exp(m[0]),np.exp(m[1]),np.exp(m[2]),np.exp(m[9]),m[11])
overlap=(m[14]>.999999)&(m[18]<1e-6)
differences={name:float(np.max(abs(q[name][overlap]/v[overlap]-1))) for name,v in zip(names,mv)}
rows=[]
for r,t in ((1e-16,6000.),(1e-16,8000.),(1e-14,8000.),(1e-7,1e6)):
    i=int(np.argmin(abs(np.log(rho/r))+abs(np.log(T/t))))
    rows.append(dict(rho=float(rho[i]),T=float(T[i]),FreeEOS_free_e=float(q['free_e'][i]),MESA_free_e=float(np.exp(m[9,i])),
                     MESA_ideal_fraction=float(m[18,i]),native_opacity=float(m[19,i]),external_eos_opacity=float(k[i]),
                     direct_H_ionized=float(q['h_ionized'][i])))
report=dict(scope='Component audit, not a validated atmosphere or evolution',points=len(rho),
            gas_first_law_max=float(np.max(defects)),gas_derivative_max=float(np.max(derivatives)),
            pt_roundtrip_max=float(max(abs(inverse))),MESA_pure_FreeEOS_overlap_points=int(sum(overlap)),
            overlap_max_relative_discrepancy=differences,selected_points=rows,
            FreeEOS_build=eos.manifest,MESA_build=mesa.manifest,runtime_seconds=time.perf_counter()-start)
weak_rho=np.array([1e-9,1e-8,1e-7]);weak_T=np.full(3,1e6)
mq=mesa.evaluate(weak_rho,weak_T);fq=eos.evaluate(weak_rho,weak_T)
gas_screen=mesa.gas_pressure_limit(weak_rho,weak_T,mq)
report['radiation_dominated_overlap']=[dict(rho=float(weak_rho[i]),T=float(weak_T[i]),
    MESA_gas_pressure=float(mq['pg'][i]),direct_gas_pressure=float(fq['pg'][i]),
    MESA_classical_gas_defect=float(gas_screen['relative_error'][i]),coulomb_bound=float(gas_screen['coupling'][i]),
    pressure_difference_over_Pr=float((mq['pg'][i]-fq['pg'][i])/fq['pr'][i]),
    MESA_total_first_law_defect=float(mesa.consistency(weak_rho[i],weak_T[i]))) for i in range(3)]
report['interpretation']='Radiation-scaled offset is consistent with historical radiation subtraction; origin not proven from archived table-generation metadata. Do not fit it away. MESA loader explicitly floors table electron abundance at 1e-4.'
np.savez_compressed(out/'grid.npz',rho=rho,T=T,FreeEOS_raw=x,MESA_raw=m,gas_first_law_defect=defects,
                    gas_derivative_error=derivatives,pt_roundtrip=inverse,opacity_external_eos=k)
report['source_snapshot_verified']=io.finish_provenance(out);io.json_write(out/'audit.json',report)
print(out);print(json.dumps({k:v for k,v in report.items() if not k.endswith('_build')},indent=2))
