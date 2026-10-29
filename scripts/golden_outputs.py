"""Byte-level golden check of short QStarS outputs before and after a source change.

The repository stores no golden files: numeric bytes depend on the platform and
library build, so the check compares two runs made on the same machine.

    python scripts/golden_outputs.py snapshot BASE_SRC      # before editing: copy src/qstars
    python scripts/golden_outputs.py run BASE_SRC OUT_BEFORE
    python scripts/golden_outputs.py run src OUT_AFTER
    python scripts/golden_outputs.py compare OUT_BEFORE OUT_AFTER

`run` executes `qstars validate --quick` and 25 accepted steps of each example
configuration (2 steps of the 512-zone hot_resolved_validation), single-threaded
and niced, in fresh processes importing the given source tree. Run directories
carry UTC stamps and are matched by name. NPZ containers embed zip timestamps,
so their arrays are compared by dtype, shape and bytes. Provenance files
(metadata.json, source_snapshot.zip, working_tree.diff) and the wall-clock
fields of termination.json describe the process, not the model, and are skipped.
Never point OUT at results/.
"""
import argparse,json,os,re,shutil,subprocess,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CONFIGS=['hot_accreting','bra_luminosity','closed_relaxation','dotan_launch','pressure_confined','two_stream_pulse']
SKIP={'metadata.json','source_snapshot.zip','working_tree.diff'}

def run(src,out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,PYTHONPATH=str(Path(src).resolve()),OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    jobs=[('validation',['validate','--quick','--output',str(out/'validation')])]
    jobs+=[(c,['run',f'configs/{c}.yaml','--output',str(out/c),'--stop-after-steps','25']) for c in CONFIGS]
    jobs+=[('hot_resolved_validation',['run','configs/hot_resolved_validation.yaml','--output',str(out/'hot_resolved_validation'),'--stop-after-steps','2'])]
    for name,args in jobs:
        code=f'from qstars.cli import main; raise SystemExit(main({args!r}))'
        with open(out/f'{name}.log','w') as log:
            result=subprocess.run(['nice','-n','15',sys.executable,'-c',code],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        print(name,'exit',result.returncode)

def tree(root):
    root=Path(root)
    def key(p):return '/'.join(re.sub(r'^\d{8}T\d{6}\.\d{6}Z_','',x) for x in p.relative_to(root).parts)
    return {key(p):p for p in root.rglob('*') if p.is_file() and p.name not in SKIP and p.suffix!='.log'}

def compare(a,b):
    a,b=tree(a),tree(b);bad=[]
    if set(a)!=set(b):bad.append(('file sets differ',sorted(set(a)^set(b))))
    for k in sorted(set(a)&set(b)):
        x,y=a[k],b[k]
        if x.suffix=='.npz':
            with np.load(x) as f,np.load(y) as g:
                if sorted(f.files)!=sorted(g.files):bad.append((k,'arrays'));continue
                bad+=[(k,n) for n in f.files if (f[n].dtype,f[n].shape,f[n].tobytes())!=(g[n].dtype,g[n].shape,g[n].tobytes())]
        elif x.name=='termination.json':
            f,g=json.loads(x.read_text()),json.loads(y.read_text())
            for d in (f,g):d.pop('runtime_seconds',None);d.pop('process_peak_rss_MiB',None)
            if json.dumps(f)!=json.dumps(g):bad.append((k,'termination'))
        elif x.read_bytes()!=y.read_bytes():bad.append((k,'bytes'))
    print(f'compared {len(set(a)&set(b))} files; differences: {len(bad)}')
    for x in bad:print('  ',x)
    return 1 if bad else 0

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('snapshot');s.add_argument('destination')
    s=sub.add_parser('run');s.add_argument('src');s.add_argument('output')
    s=sub.add_parser('compare');s.add_argument('before');s.add_argument('after')
    args=p.parse_args()
    if args.command=='snapshot':
        shutil.copytree(ROOT/'src/qstars',Path(args.destination)/'qstars',ignore=shutil.ignore_patterns('__pycache__'))
    elif args.command=='run':run(args.src,args.output)
    else:raise SystemExit(compare(args.before,args.after))
