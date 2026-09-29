#!/usr/bin/env python3
"""Matched restart comparison of patched/export-on and unmodified executables."""
import argparse,json
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import read_rvtj,read_cooling,read_hhe_populations

def compare(patched,unmodified,output):
    comparisons={}
    for group,fun,args in [('material',read_rvtj,('RVTJ',)),('populations',read_hhe_populations,('RVTJ',)),
                           ('cooling',read_cooling,('GENCOOL',len(read_rvtj(patched/'RVTJ')['radius'])))]:
        a=fun(patched/args[0],*args[1:]);b=fun(unmodified/args[0],*args[1:])
        for k in a:comparisons[group+'_'+k]=float(np.max(abs(a[k]-b[k])/np.maximum(abs(b[k]),1e-100)))
    a=json.loads((patched/'run_result.json').read_text());b=json.loads((unmodified/'run_result.json').read_text())
    exclude={'QSTARS_EXPORT','changes.json','continuation.json'}
    ax={k:v for k,v in a['inputs'].items() if k not in exclude};bx={k:v for k,v in b['inputs'].items() if k not in exclude}
    same=ax==bx;completed=all(x['status']=='completed_component_not_certified' for x in (a,b))
    result=dict(status='PASS' if max(comparisons.values())<=1e-10 and same and completed else 'FAIL',limit=1e-10,
      comparisons=comparisons,shared_native_inputs_match=same,both_runs_completed=completed,
      scope='one Lambda iteration, identical completed H/He restart, patched/export enabled vs unmodified executable',
      timings_seconds={'patched_export_on':a['seconds'],'unmodified_export_off':b['seconds']},
      executable_sha256={'patched':a['executable_sha256'],'unmodified':b['executable_sha256']},physical_convergence_claim=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(result['status'])
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('patched',type=Path);a.add_argument('unmodified',type=Path);a.add_argument('output',type=Path);compare(**vars(a.parse_args()))
