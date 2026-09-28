#!/usr/bin/env python3
"""Build pinned upstream FreeEOS in an isolated directory; never overwrite."""
import argparse,difflib,hashlib,json,os,platform,shlex,shutil,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--sdk-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--arm-kind-fix',action='store_true');a=p.parse_args()
a.source=a.source.resolve();a.sdk_dir=a.sdk_dir.resolve();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[1];cmake=root/'.venv/bin/cmake'
env=os.environ.copy();env.update(MESASDK_ROOT=str(a.sdk_dir),MACOSX_DEPLOYMENT_TARGET='14.0',OMP_NUM_THREADS='1',FC='gfortran')
prefix='source '+shlex.quote(str(a.sdk_dir/'bin/mesasdk_init.sh'))+' && '
probe=a.output/'precision_probe.f90'
probe.write_text('program probe\n print *, selected_real_kind(15,300),selected_real_kind(18,4900),selected_real_kind(32,4900)\nend program\n')
probe_exe=a.output/'precision_probe'
probe_cmd=['gfortran',str(probe),'-o',str(probe_exe)]
subprocess.run(['bash','-c',prefix+shlex.join(probe_cmd)],env=env,check=True)
precision=subprocess.run([str(probe_exe)],env=env,check=True,capture_output=True,text=True).stdout
(a.output/'precision_probe.txt').write_text(precision)
kinds=[int(k) for k in precision.split()]
original_source=a.source
if a.arm_kind_fix:
 if platform.machine() not in ('arm64','aarch64'):raise ValueError('This explicit portability repair is only for ARM')
 if not (kinds[0]!=kinds[1] and kinds[1]==kinds[2]>0):raise ValueError('Portability repair requires distinct double but identical extended/quadruple kinds')
 a.source=a.output/'source';shutil.copytree(original_source,a.source,symlinks=True)
 file=a.source/'src/mod_free_eos_debug.f90';before=file.read_text();after=before
 for rank in (1,2,3):
  line=f'     module procedure set_snan_ep_fp_kind_{rank}d\n'
  if after.count(line)!=1:raise ValueError('Upstream generic differs; re-audit required')
  after=after.replace(line,'')
 file.write_text(after)
 (a.output/'portability.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='upstream/src/mod_free_eos_debug.f90',tofile='arm/src/mod_free_eos_debug.f90')))
commands=[
 [str(cmake),'-S',str(a.source),'-B',str(a.output/'cmake'),'-DBUILD_SHARED_LIBS=OFF','-DCMAKE_POSITION_INDEPENDENT_CODE=ON','-DBUILD_TEST=ON','-DCMAKE_BUILD_TYPE=Release',f'-DCMAKE_INSTALL_PREFIX={a.output}/install',f'-DBLAS_LIBRARIES={a.sdk_dir}/math-slots/crmath/lib/libblas.a',f'-DLAPACK_LIBRARIES={a.sdk_dir}/math-slots/crmath/lib/liblapack.a;{a.sdk_dir}/math-slots/crmath/lib/libblas.a;{a.sdk_dir}/lib/libcrmath.a;{a.sdk_dir}/lib/libcrlibm.a;{a.sdk_dir}/lib/libscs.a'],
 [str(cmake),'--build',str(a.output/'cmake'),'--parallel','2'],
 [str(cmake),'--install',str(a.output/'cmake')]]
meta={'source':str(a.source),'original_source':str(original_source),'arm_kind_fix':a.arm_kind_fix,'source_file_sha256':{str(f.relative_to(a.source)):hashlib.sha256(f.read_bytes()).hexdigest() for f in a.source.rglob('*') if f.is_file()},'commands':commands,'sdk_dir':str(a.sdk_dir),'version':'3.0.0','return_codes':[]}
meta['precision_kinds']=dict(zip(('double','extended','quadruple'),kinds))
start=time.perf_counter()
for i,command in enumerate(commands):
 with (a.output/f'stage_{i}.log').open('w') as log:
  r=subprocess.run(['bash','-c',prefix+shlex.join(command)],env=env,stdout=log,stderr=subprocess.STDOUT)
 meta['return_codes'].append(r.returncode)
 print('stage',i,'return',r.returncode,flush=True)
 if r.returncode:break
meta['runtime_seconds']=time.perf_counter()-start
(a.output/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
print(a.output);raise SystemExit(meta['return_codes'][-1])
