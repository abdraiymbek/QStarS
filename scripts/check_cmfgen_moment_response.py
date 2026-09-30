#!/usr/bin/env python3
"""Actual native nonlocal MOM/VAR finite-difference fixture, fixed formal factors."""
import argparse, json, os, platform, shutil, subprocess, time
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import sha256


def check(root,output):
    root=root.resolve();output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).resolve().parents[1]/'backends/cmfgen/moment_response_probe.f90'
    shutil.copyfile(source,output/'probe.f90')
    cmd=['gfortran','-fopenmp','-ffree-line-length-none','-fcheck=all','-I'+str(root/'mod'),str(output/'probe.f90'),
         '-L'+str(root/'lib'),'-lsubs','-lmod','-lsubs']
    cmd+=['-framework','Accelerate'] if platform.system()=='Darwin' else ['-llapack','-lblas']
    cmd+=['-o',str(output/'probe')]
    env={**os.environ,'MACOSX_DEPLOYMENT_TARGET':'14.0','OMP_NUM_THREADS':'1'}
    start=time.perf_counter();build=subprocess.run(cmd,capture_output=True,text=True,cwd=output,env=env)
    (output/'build.log').write_text(build.stdout+build.stderr)
    if build.returncode:raise RuntimeError('Native fixture build failed; evidence retained')
    cases=[]
    def trial(mode,amp,label):
        inp=f'{mode} {amp}\n';r=subprocess.run([str(output/'probe')],input=inp,text=True,capture_output=True,cwd=output,env=env)
        (output/(label+'.txt')).write_text('INPUT '+inp+r.stdout+r.stderr)
        data=np.array([[float(v) for v in s.split()[1:]] for s in r.stdout.splitlines() if s.startswith('DATA ')])
        cases.append(dict(label=label,input=inp.strip(),exit_code=r.returncode,shape=list(data.shape)))
        if r.returncode or data.shape!=(5,3) or not np.isfinite(data).all():raise RuntimeError('Native fixture failure preserved: '+label)
        return data
    base=trial(0,0,'base');errors={}
    for mode,tag,column in [(1,'eta',1),(2,'chi',2)]:
        plus=trial(mode,1e-5,tag+'_plus');minus=trial(mode,-1e-5,tag+'_minus')
        fd=(plus[:,0]-minus[:,0])/2e-5;analytic=base[:,column]
        errors[tag]=dict(fd=fd.tolist(),analytic=analytic.tolist(),max_relative=float(max(abs(fd-analytic)/np.maximum(abs(fd),1e-300))))
    result=dict(cases=cases,errors=errors,target_relative=1e-6,passed=all(v['max_relative']<1e-6 for v in errors.values()),
        elapsed_seconds=time.perf_counter()-start,compiler_command=cmd,source_sha256=sha256(source),driver_sha256=sha256(__file__),
        binary_sha256=sha256(output/'probe'),native_sources={n:sha256(root/n) for n in ['plane/mom_jrel_v9.f','plane/var_jrel_v5.f','plane/edd_jrel_var_v3.f','subs/up_tx_tvx.f']},
        scope='Two frequencies, five depths, nonnegative angular fixture, coherent scattering, relativistic and radiation-advection terms. Actual native nonlocal eta/chi response with formal ratios held for this linearization control only. No full relaxed-factor, atomic, thermal, energy or wind certificate.')
    (output/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path)
    r=check(**vars(p.parse_args()));raise SystemExit(0 if r['passed'] else 2)
