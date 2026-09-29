#!/usr/bin/env python3
"""Acquire official archives with verified byte ranges and SHA256 provenance."""
import argparse, concurrent.futures, hashlib, json, os, time
from pathlib import Path
import requests

PINNED_SHA256 = {
 'source': 'fc81d2dcf2df6538058ee46eab27ab231e01152e02c12666906d583c4dae4a00',
 'atomic': '0dcb7d7284a1696493f9e1fcf7734eabdee3fa36ee29562d3763b9057658d16f',
 'models': '86a83075490924358b796d72a375a509226f2adeb42a146ab30c19cb648491fb',
}

FILES = {
 'source': ('cur_cmf_18jun25.tar.gz', 45060661),
 'atomic': ('atomic_data_21jun23.tar.gz', 1038599306),
 'models': ('three_recent_models.tar.gz', 748304324),
}

def acquire(label, directory, workers=8):
    name, size = FILES[label]
    url = 'https://sites.pitt.edu/~hillier/cmfgen_files/' + name
    path = directory / name
    part = path.with_suffix(path.suffix + '.partial')
    start = time.perf_counter()
    fd = os.open(part, os.O_CREAT | os.O_WRONLY | os.O_TRUNC, 0o644)
    os.ftruncate(fd, size)
    def fetch(offset):
        end = min(offset + 4*1024**2, size)-1
        for retry in range(3):
            try:
                r = requests.get(url, headers={'Range':f'bytes={offset}-{end}'}, timeout=60)
                r.raise_for_status()
                if r.status_code != 206 or r.headers.get('Content-Range') != f'bytes {offset}-{end}/{size}' or len(r.content) != end-offset+1:
                    raise RuntimeError('Unexpected or incomplete HTTP byte range')
                data = memoryview(r.content); position = offset
                while data:
                    n = os.pwrite(fd, data, position)
                    if n <= 0: raise OSError('Short archive write')
                    data = data[n:]; position += n
                return
            except Exception:
                if retry == 2: raise
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(fetch, range(0, size, 4*1024**2)))
        os.fsync(fd)
    finally: os.close(fd)
    digest = hashlib.sha256()
    with part.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2), b''): digest.update(b)
    if digest.hexdigest() != PINNED_SHA256[label]:
        raise RuntimeError('Archive SHA256 differs from inspected upstream release; retain partial for audit')
    part.replace(path)
    return dict(label=label, url=url, path=str(path), bytes=size,
                sha256=digest.hexdigest(), seconds=time.perf_counter()-start)

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('labels',nargs='+',choices=list(FILES))
    parser.add_argument('--directory',type=Path,default=Path('sources/cache'))
    parser.add_argument('--manifest',type=Path,required=True)
    a=parser.parse_args(); a.directory.mkdir(parents=True,exist_ok=True); rows=[]
    for label in a.labels:
        row=acquire(label,a.directory); rows.append(row)
        a.manifest.write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
