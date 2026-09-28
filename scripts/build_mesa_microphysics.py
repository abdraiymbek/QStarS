#!/usr/bin/env python3
"""Build the QStarS C ABI against a separately installed MESA 24.08.1."""
import argparse,hashlib,json,os,platform,shlex,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--mesa-dir',required=True,type=Path);p.add_argument('--sdk-dir',required=True,type=Path)
p.add_argument('--output',type=Path,default=Path('backends/mesa/build'))
a=p.parse_args();root=Path(__file__).resolve().parents[1]
a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
source=root/'backends/mesa/microphysics.f90'
cmd=['make','-f',str(root/'backends/mesa/Makefile'),f'MESA_DIR={a.mesa_dir.resolve()}',f'DEST={a.output}',f'SOURCE={source}']
env=os.environ.copy();env.update(MESA_DIR=str(a.mesa_dir.resolve()),MESASDK_ROOT=str(a.sdk_dir.resolve()),OMP_NUM_THREADS='1')
if platform.system()=='Darwin':env.setdefault('MACOSX_DEPLOYMENT_TARGET','14.0')
shell='source '+shlex.quote(str(a.sdk_dir.resolve()/'bin/mesasdk_init.sh'))+' && '+shlex.join(cmd)
start=time.perf_counter()
r=subprocess.run(['bash','-c',shell],env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(a.output/'build.log').write_text(r.stdout)
def git(*args):return subprocess.run(['git','-C',str(a.mesa_dir),*args],capture_output=True,text=True).stdout.strip()
meta={'mesa_dir':str(a.mesa_dir.resolve()),'sdk_dir':str(a.sdk_dir.resolve()),'mesa_tag':git('describe','--tags','--always'),
      'mesa_commit':git('rev-parse','HEAD'),'mesa_status':git('status','--short'), 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'platform':platform.platform(),'deployment_target':env.get('MACOSX_DEPLOYMENT_TARGET'),
      'runtime_seconds':time.perf_counter()-start,'returncode':r.returncode}
(a.output/'mesa_build_changes.diff').write_text(git('diff'))
if r.returncode==0:
    meta['library_sha256']=hashlib.sha256((a.output/'libqstars_mesa.dylib').read_bytes()).hexdigest()
    import yaml
    sample={'name':'mesa_hot_two_stream','eos':'mesa','mesa':{'library_path':str(a.output/'libqstars_mesa.dylib'),'mesa_dir':str(a.mesa_dir.resolve())},
      'opacity':{'model':'mesa'},'numerics':{'zones':128,'mesh':'log_radius_cosine','solver_tolerance':1e-10,'initial_dt_years':.001,'max_dt_years':.001},
      'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},
      'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4},
      'feeding':{'rate_msun_per_year':.01,'thermodynamics':'surface'},'wind':{'model':'prescribed','rate_msun_per_year':.008},
      'evolution':{'end_time_years':.1}}
    (a.output/'example.yaml').write_text(yaml.safe_dump(sample,sort_keys=False))
(a.output/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
print(r.stdout);print(a.output)
raise SystemExit(r.returncode)
