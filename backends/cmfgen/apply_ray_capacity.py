#!/usr/bin/env python3
"""Fix a verified upstream fixed-array overrun; no transfer equations changed."""
import argparse,hashlib,difflib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
patches=[];hashes={}
for name in ['plane/mod_space_grid_v2.f','plane/define_grid_v2.f']:
    f=a.root/name;old=f.read_text();new=old
    if name.endswith('mod_space_grid_v2.f'):
        assert old.count('TYPE (RAY_DATA) RAY(500)')==1
        new=old.replace('TYPE (RAY_DATA) RAY(500)','TYPE (RAY_DATA) RAY(2048)')
    else:
        marker='      LUER=ERROR_LU()'
        assert old.count(marker)==1
        new=old.replace(marker,marker+'''
! QStarS: reject unsupported grids before any characteristic-ray indexing.
      IF(NP .GT. SIZE(RAY))THEN
        WRITE(LUER,*)'QStarS ray capacity exceeded:',NP,SIZE(RAY)
        STOP 1
      END IF
''',1)
    f.write_text(new)
    hashes[name]=dict(original_sha256=hashlib.sha256(old.encode()).hexdigest(),patched_sha256=hashlib.sha256(new.encode()).hexdigest())
    patches.append(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=name,tofile=name)))
out=Path(__file__).parent
(out/'ray_capacity.patch').write_text(''.join(patches));(out/'ray_capacity_patch_manifest.json').write_text(json.dumps(dict(files=hashes,max_rays=2048,
 reason='665 radii + 15 core rays overran upstream RAY(500), corrupting the adjacent DTAU allocation descriptor',equations_changed=False,dust_capacity_changed=False),indent=2)+'\n')
