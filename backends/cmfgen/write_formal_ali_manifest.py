#!/usr/bin/env python3
"""Regenerate formal_ali.patch and its manifest from build_7 (base) and build_13 (patched)."""
import json,hashlib,subprocess
from pathlib import Path
here=Path(__file__).resolve().parent
FILES=['plane/mod_ray_mom_store.f','plane/solve_cmf_formal_v2.f','plane/cmf_formal_rel_v4.f','plane/var_jrel_v5.f',
       'plane/mom_jrel_v9.f','new_main/mod_subs/comp_j_blank.f','new_main/mod_subs/rd_control_variables.f',
       'newsubs/do_h0_ff.f','new_main/subs/solveba_v13.f']
base,new=here/'build_7/cur_cmf',here/'build_13/cur_cmf'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
patch=''.join(subprocess.run(['diff','-u','--label','a/'+f,'--label','b/'+f,str(base/f),str(new/f)],capture_output=True,text=True).stdout for f in FILES)
(here/'formal_ali.patch').write_text(patch)
m=dict(base_build='build_7 (Report 15 FF/stopping repairs, read-only traces)',
  source_hashes={f:dict(base=sha(base/f),patched=sha(new/f)) for f in FILES},patch_sha256=sha(here/'formal_ali.patch'),
  scope='QS_FRM_ALI transfer: J from the nonnegative linear short-characteristic relativistic CMF formal solution; '
        'coherent scattering by local ALI + Ng; moment operator retained only as the linearized response about the '
        'formal field (QS_FORMAL_MOM_STATE supplies its ratio state). H0 free-free: thread-local interpolation '
        'scalars and serialized table load. SOLVEBA: nonfinite Newton correction stops explicitly.')
(here/'formal_ali_patch_manifest.json').write_text(json.dumps(m,indent=2)+'\n');print(m['patch_sha256'])
