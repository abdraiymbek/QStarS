#!/usr/bin/env python3
"""Repair stale HMI naming and array stride in the native HI Rayleigh path."""
import argparse,difflib,json
from pathlib import Path
from qstars.cmfgen_backend import sha256
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
files={'new_main/mod_subs/comp_opac.f':None,'new_main/mod_subs/do_var_cont.f':None}
record=json.loads(Path('Reports/evidence/cool_20260929/rayleigh_source_defect.json').read_text())['evidence']
patch=[];hashes={}
for name in files:
    f=a.root/name
    if sha256(f)!=record[name]:raise ValueError('Unexpected native source '+name)
    old=f.read_text()
    new=old.replace("ION_ID(1) .EQ. 'HMI'","(ION_ID(1) .EQ. 'H0' .OR. ION_ID(1) .EQ. 'HMI')")
    if name.endswith('comp_opac.f'):
        new=new.replace('ATM(1)%NXzV_F,CONT_FREQ,ND)','ATM(ID)%NXzV_F,CONT_FREQ,ND)')
    else:
        new=new.replace('VJ(2,J,K)=VJ(2,J,K)', 'VJ(ATM(2)%EQXzV,J,K)=VJ(ATM(2)%EQXzV,J,K)')
        new=new.replace('VJ(1,J,K)=VJ(1,J,K)', 'VJ(ATM(1)%EQXzV,J,K)=VJ(ATM(1)%EQXzV,J,K)')
    if old==new:raise ValueError('Expected Rayleigh call site absent')
    f.write_text(new);hashes[name]=dict(original=record[name],patched=sha256(f))
    patch.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=name,tofile=name))
out=Path(__file__).parent
(out/'hminus_rayleigh.patch').write_text(''.join(patch))
(out/'hminus_rayleigh_patch_manifest.json').write_text(json.dumps(dict(source_hashes=hashes,
    physical_change='Use neutral HI, its actual array stride and equation index for Rayleigh opacity/linearization when H- is present',
    coefficients_changed=False),indent=2)+'\n')
