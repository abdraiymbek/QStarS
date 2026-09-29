#!/usr/bin/env python3
"""Record the capacity of a successfully built, specifically patched executable."""
import argparse,json,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
exe=a.root/'exe/cmfgen_dev.exe';module=a.root/'plane/mod_space_grid_v2.f';guard=a.root/'plane/define_grid_v2.f'
if 'RAY(2048)' not in module.read_text() or 'NP .GT. SIZE(RAY)' not in guard.read_text():raise ValueError('Capacity patch not present')
if exe.stat().st_mtime_ns<max(module.stat().st_mtime_ns,guard.stat().st_mtime_ns):raise ValueError('Executable predates source changes')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
record=dict(executable_sha256=sha(exe),max_rays=2048,explicit_native_capacity_guard=True,
 source_hashes={str(x.relative_to(a.root)):sha(x) for x in [module,guard,a.root/'new_main/cmfgen_sub.f']},
 scope='non-dust spherical transfer only; not production certification')
exe.with_suffix('.capabilities.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
