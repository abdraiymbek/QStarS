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
here=Path(__file__).parent
# Patch layers in application order. A later layer may overlay a file already
# patched by an earlier layer; a flag then stays true when the file matches the
# deepest applicable patched hash (each layer's base must equal the hash it
# replaces, so the chain is verified).
LAYER_FILES=['hminus_rayleigh_patch_manifest.json','hminus_ff_patch_manifest.json',
 'formal_ali_patch_manifest.json','ali_lin_patch_manifest.json',
 'jac_check_patch_manifest.json','krylov_patch_manifest.json']
layers=[]
for mf in LAYER_FILES:
    mp=here/mf
    if mp.exists():
        layers.append(json.loads(mp.read_text())['source_hashes'])
def baseof(row):
    return row.get('base',row.get('original'))
def chain_hashes(name,expected):
    # All hashes reachable from `expected` through later overlay layers; a
    # file at any link of a consistent chain counts as carrying the patch.
    out={expected}
    changed=True
    while changed:
        changed=False
        for lay in layers:
            row=lay.get(name)
            if row and baseof(row) in out and row['patched'] not in out:
                out.add(row['patched']);changed=True
    return out
def check(rows,flag):
    record[flag]=all(sha(a.root/name) in chain_hashes(name,row['patched']) and
        exe.stat().st_mtime_ns>(a.root/name).stat().st_mtime_ns for name,row in rows.items())
    record['source_hashes'].update({name:sha(a.root/name) for name in rows})
for mf,flag in [('hminus_rayleigh_patch_manifest.json','hminus_rayleigh_fixed'),
                ('hminus_ff_patch_manifest.json','hminus_ff_fixed'),
                ('formal_ali_patch_manifest.json','qstars_formal_ali'),
                ('ali_lin_patch_manifest.json','qstars_ali_lin'),
                ('jac_check_patch_manifest.json','qstars_jac_check'),
                ('krylov_patch_manifest.json','qstars_krylov')]:
    mp=here/mf
    record[flag]=False
    if mp.exists():
        check(json.loads(mp.read_text())['source_hashes'],flag)
for name in ['plane/mom_jrel_v9.f','plane/var_jrel_v5.f','new_main/mod_subs/comp_j_blank.f']:
    record['source_hashes'][name]=sha(a.root/name)
exe.with_suffix('.capabilities.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
