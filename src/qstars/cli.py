import argparse,copy,itertools,json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import yaml
from .config import load,resolve
from .evolution import run
from .io import load_checkpoint

def _sweep_one(args):
    c,parent=args
    p,r=run(c,parent)
    return dict(path=str(p),**r)

def main(argv=None):
    p=argparse.ArgumentParser(prog='qstars',description='QStarS research development framework; see validation gates in README.')
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run');r.add_argument('config');r.add_argument('--output',default='results')
    r.add_argument('--stop-after-steps',type=int)
    r=sub.add_parser('restart');r.add_argument('checkpoint');r.add_argument('--until-years',type=float);r.add_argument('--output',default='results')
    r=sub.add_parser('validate');r.add_argument('--output',default='results');r.add_argument('--quick',action='store_true')
    r.add_argument('--require-v1',action='store_true',help='Fail unless every v1 scientific gate is passed')
    r=sub.add_parser('plot');r.add_argument('directory')
    r=sub.add_parser('sweep');r.add_argument('config');r.add_argument('grid');r.add_argument('--workers',type=int,default=1);r.add_argument('--output',default='results')
    args=p.parse_args(argv)
    try:
        if args.command=='run':
            path,report=run(load(args.config),args.output,stop_after_steps=args.stop_after_steps)
        elif args.command=='restart':
            state,c=load_checkpoint(args.checkpoint)
            if args.until_years is not None:c['evolution']['end_time_years']=args.until_years
            path,report=run(resolve(c),args.output,restart=state)
        elif args.command=='validate':
            from .validation import campaign
            path,report=campaign(args.output,not args.quick)
            print(path/'validation.md')
            failures=[x['name'] for x in report['checks'] if not x['passed']]
            print(json.dumps(dict(v1_validated=False,failed_checks=failures),indent=2))
            return 2 if args.require_v1 or failures else 0
        elif args.command=='plot':
            from .plotting import plot_run
            print('\n'.join(plot_run(args.directory)));return 0
        else:
            c=load(args.config);grid=yaml.safe_load(Path(args.grid).read_text());configs=[]
            if args.workers<1:raise ValueError('workers must be positive')
            for i,values in enumerate(itertools.product(*grid.values())):
                model=copy.deepcopy(c)
                for key,value in zip(grid,values):
                    section,field=key.split('.',1)
                    if section not in model or field not in model[section]:raise ValueError(f'Unknown sweep axis {key}')
                    model[section][field]=value
                model['name']=f"{c['name']}_sweep_{i:04d}"
                configs.append((resolve(model),args.output))
            with ProcessPoolExecutor(max_workers=args.workers) as pool:reports=list(pool.map(_sweep_one,configs))
            print(json.dumps(reports,indent=2));return int(any(x['status']!='completed_requested_interval' for x in reports))
        print(path);print(json.dumps(report,indent=2))
        return 0 if report['status'] in ('completed_requested_interval','control_stop') else 1
    except (ValueError,OSError,KeyError) as exc:
        p.exit(2,f'qstars: {exc}\n')

if __name__=='__main__':raise SystemExit(main())
