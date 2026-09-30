#!/usr/bin/env python3
"""Run one isolated native component, retaining failures and input hashes."""
import argparse,json,os,platform,resource,signal,subprocess,time,shutil,re
from pathlib import Path
from qstars.cmfgen_backend import sha256,read_rvtj,read_spectral,LOADED_SOURCE_SHA256

RUNNER_SOURCE_SHA256 = sha256(__file__)

def check_grid_capacity(directory,executable):
    # Upstream MOD_SPACE_GRID_V2 has RAY(500), without a bounds check.
    # Only a matching build record may advertise the guarded larger array.
    limit=500;record=executable.with_suffix('.capabilities.json');info={}
    if record.is_file():
        info=json.loads(record.read_text())
        if info.get('executable_sha256')!=sha256(executable):raise RuntimeError('Build capability record does not match executable')
        if info.get('explicit_native_capacity_guard') is True and info.get('max_rays')==2048:limit=2048
    model=directory/'MODEL_SPEC'
    if model.exists():
        values=re.findall(r'(?m)^\s*(\d+)\s*\[NP\]',model.read_text())
        if len(values)!=1:raise RuntimeError('Missing/duplicate native ray count')
        if int(values[0])>limit:raise RuntimeError(f'Unsupported native ray count {values[0]} > verified capacity {limit}')
        # The current upstream release uses a stale HMI name in the Rayleigh
        # call and Jacobian, although the active negative ion is named H0.
        hminus=re.search(r'(?m)^\s*([1-9]\d*)\s*,\s*[1-9]\d*\s*,\s*[1-9]\d*\s*\[H0_ISF\]',model.read_text())
        vadat=directory/'VADAT'
        ray=vadat.exists() and re.search(r'(?m)^\s*T\s+\[INC_RAY\]',vadat.read_text())
        if hminus and ray and info.get('hminus_rayleigh_fixed') is not True:
            raise RuntimeError('H0 with Rayleigh scattering requires the verified HI selection/stride/Jacobian repair')
    return limit


def current_native_log(directory):
    """Exclude only the manifest-attested, byte-identical inherited log prefix."""
    path=directory/'OUTGEN'
    if not path.is_file():return ''
    data=path.read_bytes();initial=directory/'initial_inputs'/'OUTGEN'
    manifest=directory/'run_manifest.json'
    if initial.is_file() and manifest.is_file():
        expected=json.loads(manifest.read_text()).get('inputs',{}).get('OUTGEN')
        prefix=initial.read_bytes()
        if expected and sha256(initial)==expected and data.startswith(prefix):
            data=data[len(prefix):]
    return data.decode()


def check_native_finiteness(directory):
    """Zero nonlinear corrections can coexist with NaN luminosities upstream."""
    if re.search(r'(?<![A-Za-z])(?:NaN|[+-]?Infinity)(?![A-Za-z])',current_native_log(directory),re.I):
        raise RuntimeError('Nonfinite native diagnostic; zero corrections are not convergence')


def check_completion(directory,started_ns,export_required):
    check_native_finiteness(directory)
    log=current_native_log(directory)
    if re.search(r'Major Error|Major error|cant read BA File',log):
        raise RuntimeError('Native diagnostic reports a fatal error despite process status')
    if re.search(r'excessive iteration count|(?:RJ|Mean intensity) blowing up|You have poor flux accuracy',log,re.I):
        raise RuntimeError('Native radiation iteration failed; finite/stale output is not convergence')
    required=['RVTJ','MOD_SUM','OUTGEN']+(['QSTARS_SPECTRAL.bin'] if export_required else [])
    for name in required:
        p=directory/name
        if not p.is_file() or p.stat().st_mtime_ns<started_ns:
            raise RuntimeError(f'Native run did not produce fresh {name}')
    read_rvtj(directory/'RVTJ')
    if export_required:read_spectral(directory/'QSTARS_SPECTRAL.bin')


def run(directory,executable,timeout=3600,threads=1):
    directory=directory.resolve();executable=executable.resolve()
    if (directory/'run_manifest.json').exists():raise FileExistsError('Refuse to overwrite a recorded run')
    inputs={p.name:sha256(p) for p in directory.iterdir() if p.is_file() and p.name not in {'QSTARS_SPECTRAL.bin'}}
    manifest=dict(executable=str(executable),executable_sha256=sha256(executable),inputs=inputs,
        platform=platform.platform(),threads=threads,purpose='component_validation',production=False)
    initial=directory/'initial_inputs';initial.mkdir(exist_ok=False)
    for name in ['VADAT','IN_ITS','MODEL_SPEC','IT_SPECIFIER','RVSIG_COL','HYDRO_DEFAULTS','ADJUST_R_DEFAULTS','OUTGEN']:
        if (directory/name).is_file():shutil.copy2(directory/name,initial/name)
    manifest['adapter_sha256']=LOADED_SOURCE_SHA256
    manifest['runner_sha256']=RUNNER_SOURCE_SHA256
    (directory/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_STACKSIZE':'256M','VECLIB_MAXIMUM_THREADS':'1'}
    def stack_limit():
        soft,hard=resource.getrlimit(resource.RLIMIT_STACK)
        resource.setrlimit(resource.RLIMIT_STACK,(hard,hard))
    timer=['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v']
    started_ns=time.time_ns();start=time.perf_counter();status='numerical_failure';error=None
    with (directory/'native.log').open('wb') as log:
        try:
            manifest['verified_ray_capacity']=check_grid_capacity(directory,executable)
            p=subprocess.Popen(timer+[str(executable)],cwd=directory,env=env,stdin=subprocess.DEVNULL,
                stdout=log,stderr=subprocess.STDOUT,preexec_fn=stack_limit,start_new_session=True)
            try:
                while p.poll() is None:
                    remaining=timeout-(time.perf_counter()-start)
                    if remaining<=0:raise subprocess.TimeoutExpired(p.args,timeout)
                    try:p.wait(timeout=min(1.,remaining))
                    except subprocess.TimeoutExpired:check_native_finiteness(directory)
            except (subprocess.TimeoutExpired,RuntimeError) as exc:
                try:os.killpg(p.pid,signal.SIGTERM)
                except ProcessLookupError:pass
                try:p.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
                if isinstance(exc,subprocess.TimeoutExpired):
                    raise RuntimeError('Component runtime limit exceeded; process group stopped') from exc
                raise
            manifest['exit_code']=p.returncode
            if p.returncode:raise RuntimeError(f'Native exit code {p.returncode}')
            # A normal Fortran STOP can return zero for both fatal errors and success.
            check_completion(directory,started_ns,(directory/'QSTARS_EXPORT').exists())
            status='completed_component_not_certified'
        except Exception as exc:error=str(exc)
    manifest.update(status=status,error=error,seconds=time.perf_counter()-start,
        outputs={n:sha256(directory/n) for n in ['RVTJ','GENCOOL','HYDRO','OUTGEN','QSTARS_SPECTRAL.bin'] if (directory/n).is_file()})
    (directory/'run_result.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('directory',type=Path);a.add_argument('--executable',type=Path,required=True)
    a.add_argument('--timeout',type=float,default=3600);a.add_argument('--threads',type=int,default=1)
    r=run(**vars(a.parse_args()));print(json.dumps({k:r[k] for k in ['status','seconds','error']},indent=2))
    raise SystemExit(0 if r['status'].startswith('completed') else 2)
