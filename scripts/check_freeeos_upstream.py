#!/usr/bin/env python3
"""Report upstream byte comparisons and separately compare physical fields."""
import argparse,difflib,json,hashlib,re
from pathlib import Path
import numpy as np
from qstars import io
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);a=p.parse_args()
out=io.create_directory('results','freeeos_upstream_checks');io.provenance(out,{'build':str(a.build.resolve())})
rows=[]
for name in ('free_eos','exchange','fermi_dirac','rosenbrock','solve_linear'):
    ref=a.build/'source/test_FreeEOS'/f'{name}_test.stdout_standard'
    actual=a.build/'cmake/test_FreeEOS'/f'{name}_test.stdout'
    stderr=a.build/'cmake/test_FreeEOS'/f'{name}_test.stderr'
    expected=ref.read_text();got=actual.read_text()
    (out/f'{name}.stdout').write_text(got);(out/f'{name}.stderr').write_bytes(stderr.read_bytes())
    (out/f'{name}.diff').write_text(''.join(difflib.unified_diff(expected.splitlines(True),got.splitlines(True),fromfile='upstream',tofile='apple_silicon')))
    rows.append(dict(name=name,byte_equal=expected==got,reference_sha256=hashlib.sha256(ref.read_bytes()).hexdigest(),
                     actual_sha256=hashlib.sha256(actual.read_bytes()).hexdigest(),stderr_bytes=stderr.stat().st_size))
    if name=='free_eos':
        fields={}
        for label in ('degeneracy','pressure','density','energy','enthalpy','entropy','gamma1,gamma2,gamma3','h2rat,h2plusrat','fh2,fhe2,fhe3'):
            def values(s):return np.array([float(x) for x in re.search(r'^'+re.escape(label)+r' =\s*\n([^\n]+)',s,re.MULTILINE)[1].split()])
            x,y=values(expected),values(got)
            fields[label]=dict(reference=x.tolist(),actual=y.tolist(),
                max_scaled_absolute_error=float(max(abs(x-y)/np.maximum(abs(x),1e-12))))
        rows[-1]['first_state_physical_fields']=fields
report=dict(upstream_byte_tests=rows,byte_test_status='FAIL' if any(not r['byte_equal'] for r in rows) else 'PASS',
            scope='No claim of complete upstream suite success. First-state physical fields are quantified independently of platform-sensitive finite-difference diagnostics.',
            first_state_physical_fields_max_scaled_error=max(v['max_scaled_absolute_error'] for v in rows[0]['first_state_physical_fields'].values()))
report['source_snapshot_verified']=io.finish_provenance(out);io.json_write(out/'comparison.json',report)
print(out);print(report['first_state_physical_fields_max_scaled_error'])
