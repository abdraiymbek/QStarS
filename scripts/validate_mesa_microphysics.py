#!/usr/bin/env python3
"""Preserve direct API, derivative, blend, and first-law evidence."""
import argparse,os,subprocess,json
from pathlib import Path
import numpy as np
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.eos import Composition,IdealGasRadiation
from qstars import io
p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--mesa-dir',required=True)
p.add_argument('--parent',default='results');a=p.parse_args()
out=io.create_directory(a.parent,'mesa_microphysics_validation')
e=MesaMicrophysics(Composition(),a.library,a.mesa_dir)
rho,T=np.meshgrid(np.geomspace(1e-11,1e-3,9),np.geomspace(4000,2e7,17))
raw=e.raw(rho,T);defect=e.consistency(rho,T)
np.savez_compressed(out/'audit_grid.npz',rho=rho,T=T,raw=raw,first_law_defect=defect,
    fields=['lnPgas','lnE','lnS','Cv','Cp','gamma1','grad_ad','chiRho','chiT','lnfree_e','eta','mu',
            'HELM','OPAL_SCVH','FreeEOS','PC','Skye','CMS','ideal','opacity_radiative','opacity_effective'])
inputs=''.join(f'{r:.17e} {t:.17e} 0.75\n' for r,t in zip(rho.ravel(),T.ravel()))
env=os.environ.copy();env.update(MESA_DIR=str(Path(a.mesa_dir).resolve()),OMP_NUM_THREADS='1')
standalone=Path(a.library).resolve().parent/'check_microphysics'
r=subprocess.run([str(standalone)],input=inputs,text=True,capture_output=True,env=env)
(out/'standalone_output.txt').write_text(r.stdout);(out/'standalone_stderr.txt').write_text(r.stderr)
if r.returncode:raise RuntimeError('Standalone MESA check failed')
rows=np.array([[float(v) for v in line.split()] for line in r.stdout.splitlines() if len(line.split())==21]).T
relative=np.max(abs(rows-raw.reshape(21,-1))/np.maximum(abs(rows),1))
worst=np.unravel_index(np.argmax(defect),defect.shape)
summary={'scope':'H1=0.75, He4=0.25; direct MESA 24.08.1 EOS/opacity API',
 'standalone_relative_difference':float(relative),'interface_pass':bool(relative<1e-13),
 'grid_points':int(defect.size),'first_law_defect_max':float(defect.max()),'first_law_defect_p95':float(np.quantile(defect,.95)),
 'points_exceeding_one_percent':int(np.count_nonzero(defect>.01)),
 'worst_point':{'rho':float(rho[worst]),'T':float(T[worst]),'EOS_fractions':raw[12:19, *worst].tolist()},
 'eos_fraction_order':['HELM','OPAL_SCVH','FreeEOS','PC','Skye','CMS','ideal'],
 'production_gate':'CONDITIONAL for sampled supported states; full atmosphere/evolution gate FAIL',
 'build':e.manifest}
io.json_write(out/'validation.json',summary)
print(out);print(json.dumps(summary,indent=2))
