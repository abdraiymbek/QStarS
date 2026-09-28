#!/usr/bin/env python3
"""Hash the installed microphysics datasets without modifying their contents."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
from qstars import io
p=argparse.ArgumentParser();p.add_argument('--mesa-dir',type=Path,required=True);p.add_argument('--parent',default='results');a=p.parse_args()
out=io.create_directory(a.parent,'mesa_data_inventory')
root=a.mesa_dir.resolve()/'data'
directories=['chem_data','eosDT_data','eosFreeEOS_data','eosCMS_data','ionization_data','kap_data']
start=time.perf_counter();files=[];total=0
for directory in directories:
    for file in sorted((root/directory).rglob('*')):
        if not file.is_file():continue
        digest=hashlib.sha256()
        with file.open('rb') as stream:
            for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
        stat=file.stat();total+=stat.st_size
        files.append(dict(path=str(file.relative_to(root)),bytes=stat.st_size,sha256=digest.hexdigest()))
    print(directory,len(files),total,flush=True)
payload=json.dumps(files,sort_keys=True,separators=(',',':')).encode()
with gzip.open(out/'data_files.json.gz','wb') as stream:stream.write(payload)
io.json_write(out/'inventory.json',dict(mesa_data_path=str(root),directories=directories,files=len(files),bytes=total,
    manifest_sha256=hashlib.sha256(payload).hexdigest(),runtime_seconds=time.perf_counter()-start,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope='Installed microphysics datasets including files present in these directories; not an assertion that every file is used by one model'))
print(out)
