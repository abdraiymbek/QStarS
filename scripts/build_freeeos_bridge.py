#!/usr/bin/env python3
"""Link a small C ABI against the independently built GPL FreeEOS library."""
import argparse,hashlib,json,os,shlex,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--freeeos-build',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();build=a.freeeos_build.resolve();out=a.output.resolve()
out.mkdir(parents=True,exist_ok=False)
meta=json.loads((build/'build.json').read_text())
if meta['return_codes']!=[0,0,0]:raise ValueError('FreeEOS build was not successful')
sdk=Path(meta['sdk_dir']);lib=build/'install/lib'
source=Path(__file__).resolve().parents[1]/'backends/freeeos/interface.f90'
command=['gfortran','-O2','-fPIC','-ffree-line-length-none','-dynamiclib',f'-I{lib}/fortran/modules/free_eos',f'-J{out}',str(source),str(lib/'libfree_eos.a')]
command += [str(sdk/f'math-slots/crmath/lib/lib{name}.a') for name in ('lapack','blas')]
command += [str(sdk/f'lib/lib{name}.a') for name in ('crmath','crlibm','scs')]
binary=out/'libqstars_freeeos.dylib';command+=['-o',str(binary)]
env=os.environ.copy();env.update(MESASDK_ROOT=str(sdk),MACOSX_DEPLOYMENT_TARGET='14.0',OMP_NUM_THREADS='1')
with (out/'build.log').open('w') as log:
    result=subprocess.run(['bash','-c','source '+shlex.quote(str(sdk/'bin/mesasdk_init.sh'))+' && '+shlex.join(command)],env=env,stdout=log,stderr=subprocess.STDOUT)
record=dict(freeeos_build=meta,command=command,return_code=result.returncode,bridge_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
if result.returncode==0:record['library_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest()
(out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
print(out);raise SystemExit(result.returncode)
