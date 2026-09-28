#!/usr/bin/env python3
"""Isolated upstream MESA/QUEST evidence runs, never a QStarS physics claim.

Copies source and inlists, records every override and optional unit repair.
MESA's own energy error is retained; it is not the QStarS whole-system ledger.
"""
import argparse,hashlib,json,os,shlex,shutil,subprocess,time,resource
from pathlib import Path
from qstars import io
p=argparse.ArgumentParser()
p.add_argument('--family',choices=['mesa_smoke','quest_saturated','quest_bondi'],required=True)
p.add_argument('--mesa-dir',type=Path,required=True);p.add_argument('--sdk-dir',type=Path,required=True)
p.add_argument('--quest-dir',type=Path,default=Path('sources/cache/MESA-QUEST'))
p.add_argument('--startup-cgs-fix',action='store_true');p.add_argument('--max-models',type=int,default=100)
p.add_argument('--saved-seed',action='store_true')
p.add_argument('--max-wall-seconds',type=float,default=600.);p.add_argument('--parent',default='results')
a=p.parse_args();out=io.create_directory(a.parent,a.family+('_startup_cgs' if a.startup_cgs_fix else '_upstream')).resolve()
source=a.mesa_dir/'star/work' if a.family=='mesa_smoke' else a.quest_dir/'qs_satcon_winds_ApJ2026'
work=out/'work';work.mkdir()
for item in source.iterdir():
    if item.name in ('src','make'):
        shutil.copytree(item,work/item.name,ignore=shutil.ignore_patterns('*.o','*.mod','*.smod','.depend'))
    elif item.is_file() and (item.name.startswith('inlist') or item.name in ('mk','rn','re','clean','history_columns.list','profile_columns.list','initial_model.mod','Ki_tab.dat')):
        shutil.copy2(item,work/item.name)
original=(work/'inlist_project').read_text()
project=original.replace('pgstar_flag = .true.','pgstar_flag = .false.')
project=project.replace('&controls',f'&controls\n max_model_number = {a.max_models}\n')
changes=['Disabled PGSTAR (build has no graphics)',f'max_model_number={a.max_models}']
if a.saved_seed:
    if a.family=='mesa_smoke':raise ValueError('Saved QUEST seed requires QUEST family')
    project=project.replace('create_pre_main_sequence_model = .true.',"create_pre_main_sequence_model = .false.\n load_saved_model = .true.\n load_model_filename = 'initial_model.mod'")
    changes.append('Load upstream supplied initial_model.mod instead of rebuilding pre-main-sequence seed')
if a.family=='quest_bondi':
    project=project.replace('x_logical_ctrl(2) = .true.','x_logical_ctrl(2) = .false.')
    changes.append('Select published Bondi boundary using x_logical_ctrl(2)=false')
(work/'inlist_project').write_text(project)
if a.startup_cgs_fix:
    if a.family=='mesa_smoke':raise ValueError('No QUEST startup in ordinary MESA model')
    f=work/'src/run_star_extras.f90';s=f.read_text()
    old='s% M_center = new_core_mass !M_BH'
    if s.count(old)!=1:raise ValueError('Upstream startup expression differs; re-audit required')
    s=s.replace(old,'s% M_center = new_core_mass * Msun ! explicit QStarS audit: M_center is CGS')
    f.write_text(s);changes.append('Startup M_center solar-mass value converted to grams; no rate/structure coefficients changed')
headers=work/'history_columns.list'
if headers.exists():
    with headers.open('a') as f:
        f.write('\n! QStarS audit diagnostics\n total_energy\n rel_E_err\n')
def git(root,*args):return subprocess.run(['git','-C',str(root),*args],text=True,capture_output=True).stdout.strip()
meta={'family':a.family,'upstream_mesa_commit':git(a.mesa_dir,'rev-parse','HEAD'),
      'upstream_quest_commit':git(a.quest_dir,'rev-parse','HEAD') if a.family!='mesa_smoke' else None,
      'changes':changes,'status':'preparing','v1_validated':False,'max_wall_seconds':a.max_wall_seconds}
meta['input_sha256']={str(f.relative_to(work)):hashlib.sha256(f.read_bytes()).hexdigest() for f in work.rglob('*') if f.is_file()}
io.json_write(out/'run.json',meta)
(work/'inlist_project.original').write_text(original)
env=os.environ.copy();env.update(MESA_DIR=str(a.mesa_dir.resolve()),MESASDK_ROOT=str(a.sdk_dir.resolve()),OMP_NUM_THREADS='2',MACOSX_DEPLOYMENT_TARGET='14.0')
prefix='source '+shlex.quote(str(a.sdk_dir.resolve()/'bin/mesasdk_init.sh'))+' && '
start=time.perf_counter()
with (out/'build.log').open('w') as log:
    build=subprocess.run(['bash','-c',prefix+'./mk'],cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT)
meta['build_returncode']=build.returncode
if build.returncode:meta['status']='build_failure'
else:
    with (out/'run.log').open('w') as log:
        try:
            result=subprocess.run(['bash','-c',prefix+'exec ./star'],cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=a.max_wall_seconds)
            meta['returncode']=result.returncode;meta['status']='process_exited'
        except subprocess.TimeoutExpired:
            meta['status']='wall_budget_stop'
    # Exit code alone is not proof that the requested physical track completed.
    meta['last_log_lines']=(out/'run.log').read_text(errors='replace').splitlines()[-30:]
meta['wall_seconds']=time.perf_counter()-start
meta['child_peak_rss_MiB']=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024**2
io.json_write(out/'run.json',meta)
print(out);print(json.dumps(meta,indent=2))
