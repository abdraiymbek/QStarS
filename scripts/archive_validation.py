#!/usr/bin/env python3
"""Preserve selected evidence, keeping first/final profiles and all run budgets.

Original result directories are never modified. Upstream stellar source,
compiled libraries and initial stellar model binaries are not redistributed.
"""
import argparse,gzip,hashlib,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--destination',type=Path,required=True);p.add_argument('sources',nargs='+',type=Path);a=p.parse_args()
a.destination.mkdir(parents=True,exist_ok=False)
manifest=[];omitted=[]
for source in a.sources:
    if not source.is_dir():raise ValueError(f'Missing evidence directory: {source}')
    profiles={}
    for file in source.rglob('profile_*.npz'):profiles.setdefault(file.parent,[]).append(file)
    keep={file for files in profiles.values() for file in (min(files),max(files))}
    for file in sorted(source.rglob('*')):
        if not file.is_file():continue
        relative=file.relative_to(source)
        upstream='work' in relative.parts
        allowed_upstream=(file.name in ('history.data','inlist_project','inlist_project.original','history_columns.list','profile_columns.list'))
        if ((upstream and not allowed_upstream) or file.suffix in ('.dylib','.o','.mod','.smod','.a','.so')
            or file.name=='check_microphysics' or (file.name.startswith('profile_') and file.suffix=='.npz' and file not in keep)):
            omitted.append(str(file));continue
        target=a.destination/source.name/relative;target.parent.mkdir(parents=True,exist_ok=True)
        compress=file.name in ('events.jsonl','run.log','build.log','standalone.stdout','standalone.stderr')
        if compress:target=target.with_name(target.name+'.gz')
        data=file.read_bytes()
        if compress:
            with gzip.open(target,'wb') as stream:stream.write(data)
        else:target.write_bytes(data)
        manifest.append(dict(original=str(file),archived=str(target.relative_to(a.destination)),
                       original_sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),compressed=compress))
(a.destination/'archive_manifest.json').write_text(json.dumps(dict(files=manifest,omitted=omitted,
    scope='Selected complete histories, ledgers, checkpoints and source snapshots; first/final profiles. Original results remain intact.',
    omitted_upstream_assets='Source commit/input hashes identify separately installed MESA and cached QUEST; no upstream binaries or full stellar source redistributed.'),indent=2)+'\n')
print(a.destination,len(manifest),'files')
