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
ray_manifest=Path(__file__).parent/'hminus_rayleigh_patch_manifest.json'
record['hminus_rayleigh_fixed']=False
if ray_manifest.exists():
    required=json.loads(ray_manifest.read_text())['source_hashes']
    record['hminus_rayleigh_fixed']=all(sha(a.root/name)==row['patched'] and
        exe.stat().st_mtime_ns>(a.root/name).stat().st_mtime_ns for name,row in required.items())
    record['source_hashes'].update({name:sha(a.root/name) for name in required})
ff_manifest=Path(__file__).parent/'hminus_ff_patch_manifest.json'
record['hminus_ff_fixed']=False
if ff_manifest.exists():
    required=json.loads(ff_manifest.read_text())['source_hashes']
    record['hminus_ff_fixed']=all(sha(a.root/name)==row['patched'] and
        exe.stat().st_mtime_ns>(a.root/name).stat().st_mtime_ns for name,row in required.items())
    record['source_hashes'].update({name:sha(a.root/name) for name in required})
for name in ['plane/mom_jrel_v9.f','plane/var_jrel_v5.f','new_main/mod_subs/comp_j_blank.f']:
    record['source_hashes'][name]=sha(a.root/name)
exe.with_suffix('.capabilities.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
