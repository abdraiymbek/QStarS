#!/usr/bin/env python3
"""Fresh-process Mac benchmarks; run after concurrent validation jobs finish."""
import json,platform,subprocess,sys,time
from pathlib import Path
from qstars import io
out=io.create_directory('results','performance_validation')
cases=[Path('configs')/(name+'.yaml') for name in
       ('closed_relaxation','hot_accreting','two_stream_pulse','bra_luminosity','dotan_launch','pressure_confined','hot_resolved_validation')]
mesa=Path('backends/mesa/build_verified/ten_year.yaml')
if mesa.exists():cases.append(mesa)
rows=[]
for file in cases:
    start=time.perf_counter()
    proc=subprocess.run([sys.executable,'-m','qstars','run',str(file),'--output',str(out)],text=True,capture_output=True)
    if not proc.stdout.strip():raise RuntimeError(proc.stderr)
    path=Path(proc.stdout.splitlines()[0]);report=json.loads((path/'termination.json').read_text())
    report.update(configuration=str(file),path=str(path),subprocess_wall_seconds=time.perf_counter()-start,exit_code=proc.returncode)
    rows.append(report);io.json_write(out/'performance.json',{'complete':False,'rows':rows})
    print(file.name,report['status'],report['runtime_seconds'],flush=True)
hardware={key:subprocess.check_output(['sysctl','-n',key],text=True).strip() for key in ('machdep.cpu.brand_string','hw.memsize')}
io.json_write(out/'performance.json',dict(complete=True,rows=rows,hardware=hardware,platform=platform.platform(),
    scope='Fresh sequential CPU processes; initialization and I/O included, plotting excluded. No GPU or cloud.',v1_validated=False))
print(out)
