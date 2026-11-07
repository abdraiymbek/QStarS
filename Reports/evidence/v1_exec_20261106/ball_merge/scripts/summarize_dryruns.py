"""Summaries of the goal-3 DRYRUNs: control check against the registered smoke run (history.csv), and per run the stop,
the rejection causes, the zones and the outer-front resolution. Writes data/dryrun_summary.json."""
from bm_common import *
import csv
def rows(path):
    with open(path) as f:return list(csv.DictReader(f))
reg=rows(RUN/'history.csv')
out={}
for d in sorted(BIG.glob('E*')):
    if not (d/'history.csv').exists():continue
    h=rows(d/'history.csv')
    steps=[r for r in h if r['event']=='step']
    ev=[json.loads(l) for l in open(d/'events.jsonl')]
    rej=[e for e in ev if e['event']=='rejected']
    acc=[e for e in ev if e['event']=='accepted']
    rep=json.loads((d/'DRYRUN_report.json').read_text()) if (d/'DRYRUN_report.json').exists() else None
    term=json.loads((d/'termination.json').read_text()) if (d/'termination.json').exists() else None
    kinds={}
    for e in rej:
        m=e['message'];key='controller' if 'delta_actual' in m else 'newton' if 'Quasi-static solve failed' in m else m[:60]
        kinds[key]=kinds.get(key,0)+1
    s=dict(running=rep is None,steps=len(steps),last_step=steps[-1]['step'] if steps else None,
           bh_mass=float(steps[-1]['bh_mass_msun']) if steps else None,zones_last=int(steps[-1]['zones']) if steps else None,
           zones_max=max(int(r['zones']) for r in steps) if steps else None,
           merges=sum(int(r['step_merges'] or 0) for r in steps),splits=sum(int(r['step_splits'] or 0) for r in steps),
           median_dt=float(np.median([float(r['step_dt']) for r in steps if r["step_dt"]])) if len(steps)>1 else None,
           rejections=len(rej),rejection_kinds=kinds,last_rejections=[(e['step'],round(e['dt_years'],4),e['message'][:140]) for e in rej[-10:]],
           report={k:(rep or {}).get(k) for k in ('status','reason','final_bh_mass_msun','final_step','wall_seconds')},
           termination=term and {k:term.get(k) for k in ('status','reason','final_bh_mass_msun','final_step')})
    import re
    evs=[r for r in steps[1:] if int(r['step_merges'] or 0) or int(r['step_splits'] or 0)]
    non=[r for r in steps[1:] if not (int(r['step_merges'] or 0) or int(r['step_splits'] or 0))]
    da=lambda rr:[float(r['metric_delta_actual']) for r in rr if r['metric_delta_actual']]
    s['event_steps']=dict(count=len(evs),max_delta_actual=max(da(evs),default=None),
                          over_delta_target=sum(x>.005 for x in da(evs)),max_delta_actual_nonevent=max(da(non),default=None))
    res=[float(m.group(1)) for e in rej for m in [re.search(r'residual ([0-9.eE+-]+)\)',e['message'])] if m]
    s['newton_rejection_residuals']=dict(min=min(res,default=None),max=max(res,default=None))
    s['passed_46_msun']=bool(steps) and float(steps[-1]['bh_mass_msun'])>46.0039
    if d.name.startswith('E0'):
        n=min(len(h),len(reg));cols=[c for c in h[0] if c in reg[0]]
        first=next((i for i in range(n) if any(h[i][c]!=reg[i][c] for c in cols)),None)
        s['control_vs_registered']=dict(rows_compared=n,registered_rows=len(reg),first_differing_row=first,
                                        identical=first is None)
    out[d.name]=s
    print(d.name,json.dumps(s,indent=0,default=str)[:2500],'\n')
jdump(out,DATA/'dryrun_summary.json')
