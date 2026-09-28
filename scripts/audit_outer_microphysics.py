#!/usr/bin/env python3
"""Audit additional low-density EOS queries and the native atmosphere boundary.

This bypasses the production Python query rectangle only for recorded library
audits. It does not silently authorize those states for envelope evolution.
"""
import argparse,ctypes as ct,json
from pathlib import Path
import numpy as np
from qstars import io
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics

p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--mesa-dir',required=True);a=p.parse_args()
out=io.create_directory('results','outer_microphysics_audit')
io.provenance(out,vars(a));eos=MesaMicrophysics(Composition(),a.library,a.mesa_dir)
def raw(rho,T):
    rho,T=np.broadcast_arrays(rho,T);r=np.ascontiguousarray(rho.ravel());t=np.ascontiguousarray(T.ravel())
    data=np.empty((21,r.size),order='F');dr=np.empty_like(data);dt=np.empty_like(data);err=ct.c_int()
    eos.lib.qstars_eval(r.size,r,t,eos.composition.X,data,dr,dt,ct.byref(err))
    if err.value:raise ValueError(f'MESA error {err.value}')
    return data
def first_law(rho,T):
    q=raw(rho,T);P=np.exp(q[0])+eos.radiation_constant*T**4/3;h=1e-4;errors=[]
    for density_axis in (True,False):
        a=raw(rho*np.exp(h) if density_axis else rho,T if density_axis else T*np.exp(h))
        b=raw(rho*np.exp(-h) if density_axis else rho,T if density_axis else T*np.exp(-h))
        du=(np.exp(a[1])-np.exp(b[1]))/(2*h);ds=(np.exp(a[2])-np.exp(b[2]))/(2*h)
        errors.append(abs(T*ds-du+(P/rho if density_axis else 0))/(abs(du)+P/rho))
    return q,np.maximum(*errors)
rho,T=np.meshgrid(np.logspace(-22,-12,11),[4000.,5000.,6000.,8000.,1e4,2e4,4e4],indexing='ij')
rho,T=rho.ravel(),T.ravel();q,defect=first_law(rho,T)
np.savez_compressed(out/'grid.npz',rho=rho,T=T,raw=q,first_law_defect=defect)
rows=[]
for teff,g in ((8000.,3.),(6000.,1.),(5777.,27400.)):
    row={'teff':teff,'gravity':g,'kappa':.34,'tau':2/3}
    for mode in ('fixed','varying_constant_opacity','varying'):
        try:row[mode]=eos.native_atmosphere(teff,g,mode=mode)
        except ValueError as e:row[mode]={'error':str(e)}
    row['constant_opacity_gas_pressure_ratio']=row['varying_constant_opacity']['gas_pressure']/row['fixed']['gas_pressure']
    rows.append(row)
report={'scope':'EOS query audit and MESA boundary-prescription diagnostic, not stellar evolution.',
        'query_points':len(rho),'max_first_law_defect':float(max(defect)),
        'failed_first_law_screen':int(sum(defect>.01)),
        'neutral_ideal_backstop_points':int(sum(q[18]>1e-6)),
        'native_atmospheres':rows,'build':eos.manifest,'v1_validated':False}
report['source_snapshot_verified']=io.finish_provenance(out)
io.json_write(out/'audit.json',report);print(out);print(json.dumps(report,indent=2))
