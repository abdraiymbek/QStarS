#!/usr/bin/env python3
"""Compact native evidence; do not copy linked atomic archives or large matrices."""
import argparse,gzip,json,re,shutil
from pathlib import Path
from qstars.cmfgen_backend import sha256

NAMES={'VADAT','IN_ITS','MODEL_SPEC','RVSIG_COL','RDINR','HYDRO_DEFAULTS','ADJUST_R_DEFAULTS',
       'IT_SPECIFIER','RVTJ','MOD_SUM','MODEL','OUTGEN','WARNINGS','GENCOOL','HYDRO','OBSFLUX',
       'CORRECTION_SUM','STEQ_VALS','POINT1','POINT2','SCRTEMP','QSTARS_EXPORT','EDDFACTOR_INFO',
       'native.log','run_manifest.json','run_result.json','atomic_manifest.json','reference_origin.json',
       'component_definition.json','continuation.json','changes.json','failure_attribution.json','termination_reason.txt'}

def archive(sources,destination):
    destination.mkdir(parents=True,exist_ok=False);files=[];external=[];omitted=[]
    for source in sources:
        for p in sorted(source.iterdir()):
            if p.is_symlink():
                external.append(dict(original=str(p),target=str(p.resolve()),sha256=sha256(p) if p.is_file() else None));continue
            if p.is_dir():
                if p.name=='initial_inputs':
                    for f in p.iterdir():
                        if f.is_file():
                            target=destination/source.name/'initial_inputs'/f.name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,target)
                            files.append(dict(original=str(f),archived=str(target.relative_to(destination)),sha256=sha256(f),bytes=f.stat().st_size,compressed=False))
                continue
            if not p.is_file():continue
            if p.name not in NAMES and not re.fullmatch(r'[A-Za-z0-9]+(?:OUT|_IN)',p.name):
                omitted.append(dict(original=str(p),bytes=p.stat().st_size,sha256=sha256(p),reason='large/intermediate native data retained locally; reproducible from pinned source and recorded inputs'));continue
            target=destination/source.name/(p.name+'.gz');target.parent.mkdir(parents=True,exist_ok=True)
            with p.open('rb') as src,gzip.open(target,'wb') as out:shutil.copyfileobj(src,out)
            files.append(dict(original=str(p),archived=str(target.relative_to(destination)),sha256=sha256(p),bytes=p.stat().st_size,compressed=True))
    manifest=dict(scope='Selected native input/output/restart evidence. Atomic links and large spectral/Jacobian files remain external, explicitly hashed. Not a standalone upstream distribution.',
                  files=files,external_atomic_links=external,omitted_local_files=omitted)
    (destination/'archive_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for row in files:
        p=destination/row['archived']
        import hashlib
        with (gzip.open(p,'rb') if row['compressed'] else p.open('rb')) as f:
            h=hashlib.sha256();count=0
            for block in iter(lambda:f.read(8*1024**2),b''):h.update(block);count+=len(block)
        if h.hexdigest()!=row['sha256'] or count!=row['bytes']:raise RuntimeError('Archive verification failed')
    (destination/'verification.json').write_text(json.dumps(dict(files_verified=len(files),failures=0),indent=2)+'\n')
    print(f'{len(files)} archived files verified; {len(omitted)} intermediate files explicitly omitted')
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--destination',type=Path,required=True);a.add_argument('sources',nargs='+',type=Path);archive(**vars(a.parse_args()))
