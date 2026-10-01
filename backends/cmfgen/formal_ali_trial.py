#!/usr/bin/env python3
"""Released cool continuation with formal-solution ALI transfer (QS_FRM_ALI); never production."""
import argparse,json
from pathlib import Path
from continue_component import prepare_continuation
from prepare_reference import keyword
from run_component import run
from qstars.cmfgen_backend import sha256

# Numerical controls only; identical to the Report 16 released trial except the transfer mode.
RELEASE=[('USE_FIXED_J','F'),('FIX_T','F'),('MAX_LIN','2.0'),('N_FIX_BA','0'),('FIX_BA','0'),
         ('LAM_VAL','1.0E99'),('NUM_LAM','0'),('COMP_F','T'),('MAX_dT','0.05'),
         ('CHK_J','NONE'),('CHK_H','NONE'),('CHK_XM','NONE')]


def prepare(source,destination,iterations,tol,max_it,extra=(),lambda_stage=False):
    prepare_continuation(source,destination,iterations)
    for k,v in RELEASE+[('QS_FRM_ALI','T'),('QS_ALI_TOL',f'{tol:.3E}'),('QS_ALI_MAX',str(max_it))]+list(extra):
        keyword(destination/'VADAT',k,v)
    if lambda_stage:
        # Populations relaxed with T held and J recomputed from the formal transfer each iteration.
        keyword(destination/'IN_ITS','DO_LAM_IT','T')
        keyword(destination/'VADAT','FIX_T','T')
        (destination/'IT_SPECIFIER').write_text('LAMBDA 1\n')
    (destination/'QSTARS_EXPORT').touch()
    (destination/'closure_control.json').write_text(json.dumps(dict(
        scope='Released cool H/He continuation. J is the nonnegative relativistic comoving-frame characteristic '
              'solution; coherent scattering converged by local ALI; moment operator used only as the linearized '
              'response about the formal field. No radiation/flux clipping or physics removal.',
        transfer='QS_FRM_ALI',ali_tolerance=tol,ali_max_formal_solutions=max_it,extra_controls=dict(extra),
        stage='lambda_fixed_T' if lambda_stage else 'full_release',
        physics_removed=False,production=False,
        input_source=str(source.resolve()),checkpoint_sha256=sha256(source/'SCRTEMP')),indent=2)+'\n')
    return destination


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('destination',type=Path)
    p.add_argument('executable',type=Path);p.add_argument('--iterations',type=int,default=40)
    p.add_argument('--tol',type=float,default=1e-10);p.add_argument('--max-it',type=int,default=2000)
    p.add_argument('--timeout',type=int,default=7200);p.add_argument('--threads',type=int,default=1)
    p.add_argument('--set',nargs=2,action='append',default=[],metavar=('KEY','VALUE'))
    p.add_argument('--lambda-stage',action='store_true')
    a=p.parse_args()
    prepare(a.source,a.destination,a.iterations,a.tol,a.max_it,a.set,a.lambda_stage)
    r=run(a.destination,a.executable,a.timeout,a.threads)
    print(json.dumps(dict(directory=str(a.destination),status=r['status'],error=r['error'],seconds=r['seconds'])))
