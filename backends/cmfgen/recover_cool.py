#!/usr/bin/env python3
"""Reproducible FF-repaired numerical initialization/refinement; never production."""
import argparse,json
from pathlib import Path
from prepare_cool import prepare
from continue_component import prepare_continuation
from prepare_reference import keyword
from run_component import run
from qstars.cmfgen_backend import sha256


def execute(source,destination,atomic,executable,factor=2):
    def trial(d):
        r=run(d,executable,2400)
        print(json.dumps(dict(directory=str(d),status=r['status'],error=r['error'],seconds=r['seconds'])),flush=True)
        if not r['status'].startswith('completed'):raise RuntimeError('Failed native stage preserved: '+str(d))
    bootstrap=destination.with_name(destination.name+'_bootstrap')
    prepare(source,bootstrap,atomic,case='T9',iterations=1,factor=factor)
    trial(bootstrap)
    fixed=destination.with_name(destination.name+'_fixedJ')
    prepare_continuation(bootstrap,fixed,100)
    for k,v in [('DO_LAM_IT','T'),('DO_LAM_AUTO','F'),('DO_T_AUTO','F')]:keyword(fixed/'IN_ITS',k,v)
    for k,v in [('USE_FIXED_J','T'),('FIX_T','T'),('MAX_LAM','10.0'),('DO_NG','F')]:keyword(fixed/'VADAT',k,v)
    (fixed/'IT_SPECIFIER').write_text('LAMBDA 1\n');trial(fixed)
    prepare_continuation(fixed,destination,60)
    for k,v in [('USE_FIXED_J','F'),('MAX_LIN','2.0'),('MAX_dT','0.05'),('DO_NG','F'),('N_FIX_BA','0'),('FIX_BA','0'),('NUM_LAM','0'),('LAM_VAL','1.0E99')]:keyword(destination/'VADAT',k,v)
    (destination/'recovery_control.json').write_text(json.dumps(dict(purpose='fully released T/J; force full linearization by numerical switching threshold; fixed-J stage only initialization',factor=factor,strict_all_level_stopping=True,checkpoint_sha256=sha256(fixed/'SCRTEMP'),physical_inputs_changed=False),indent=2)+'\n')
    trial(destination)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('destination',type=Path);p.add_argument('atomic',type=Path);p.add_argument('executable',type=Path);p.add_argument('--factor',type=int,default=2);execute(**vars(p.parse_args()))
