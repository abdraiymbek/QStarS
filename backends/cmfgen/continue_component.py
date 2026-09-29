#!/usr/bin/env python3
"""Continue in a new directory; retain all restart state and immutable prior trial."""
import argparse,json,shutil
from pathlib import Path
from prepare_reference import keyword
from qstars.cmfgen_backend import sha256

def prepare_continuation(source,destination,iterations=40):
    source=source.resolve()
    if not (source/'MOD_SUM').is_file():raise ValueError('Source has no final material summary')
    shutil.copytree(source,destination,symlinks=True,ignore=shutil.ignore_patterns(
        'run_manifest.json','run_result.json','native.log','QSTARS_SPECTRAL.bin','initial_inputs'))
    # These are numerical controls; composition, atomic data and flow are unchanged.
    keyword(destination/'IN_ITS','NUM_ITS',str(iterations))
    keyword(destination/'IN_ITS','DO_LAM_IT','F')
    keyword(destination/'IN_ITS','DO_LAM_AUTO','F')
    keyword(destination/'IN_ITS','DO_T_AUTO','F')
    keyword(destination/'VADAT','FIX_T','F')
    # Request FULL cycles; native acceleration/final-step reuse still applies.
    (destination/'IT_SPECIFIER').write_text('FULL 1\n')
    text=(destination/'VADAT').read_text()
    (destination/'VADAT').write_text(''.join(s for s in text.splitlines(True) if '[STOP_ON_WARN]' not in s))
    (destination/'continuation.json').write_text(json.dumps(dict(source=str(source),iterations=iterations,
        reason='resolve nonlinear corrections without changing physical inputs',
        source_state={n:sha256(source/n) for n in ['RVTJ','SCRTEMP','POINT1','POINT2','MODEL'] if (source/n).exists()},
        removed_unknown_keyword='STOP_ON_WARN was a no-op; native STOP_IF_MW remains enabled'),indent=2)+'\n')
    return destination
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('destination',type=Path)
    p.add_argument('--iterations',type=int,default=40);print(prepare_continuation(**vars(p.parse_args())))
