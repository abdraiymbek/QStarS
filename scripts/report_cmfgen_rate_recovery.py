#!/usr/bin/env python3
"""Report all retained recovery trials, excluding verified inherited log prefixes."""
import argparse,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qstars.cmfgen_backend import sha256,nonlinear_corrections


def trial_text(p,result):
    text=(p/'OUTGEN').read_text() if (p/'OUTGEN').exists() else ''
    prefix=p/'initial_inputs/OUTGEN'
    if not prefix.exists() and (p/'continuation.json').exists():
        prefix=Path(json.loads((p/'continuation.json').read_text())['source'])/'OUTGEN'
    if prefix.exists() and sha256(prefix)==result.get('inputs',{}).get('OUTGEN'):
        old=prefix.read_text()
        if text.startswith(old):return text[len(old):],dict(bytes=len(old),sha256=sha256(prefix),source=str(prefix))
    return text,None


def report(output):
    output.mkdir(parents=True,exist_ok=True);rows=[]
    for p in sorted(Path('results').glob('cmfgen_recovery_*')):
        if (p/'debugger_manifest.json').exists():
            d=json.loads((p/'debugger_manifest.json').read_text())
            rows.append(dict(case=p.name,status='debugger_native_crash',seconds=None,executable_sha256=d['executable_sha256'],error=d['native_outcome'],scope=d['scope'],hashes={'debugger_manifest.json':sha256(p/'debugger_manifest.json')}));continue
        if not (p/'run_result.json').exists():raise RuntimeError('Trial still active or lacks a preserved result: '+str(p))
        r=json.loads((p/'run_result.json').read_text());text,prefix=trial_text(p,r)
        history=nonlinear_corrections(text);log=(p/'native.log').read_text();model=(p/'MODEL_SPEC').read_text()
        nd=re.search(r'(?m)^\s*(\d+)\s+\[ND\]',model)
        current=re.findall(r'Current great iteration count is\s+(\d+)',text)
        row=dict(case=p.name,status=r['status'],error=r['error'],native_exit=r.get('exit_code'),seconds=r['seconds'],ND=int(nd[1]) if nd else None,
          correction_history=history,last_counter=int(current[-1]) if current else None,
          temperature_held_iterations=text.count('Temperature held fixed at all depths'),verified_inherited_log_prefix=prefix,
          radiation_instability=bool(re.search(r'blowing up|excessive iteration count',text)),executable_sha256=r['executable_sha256'],
          hashes={n:sha256(p/n) for n in ['run_manifest.json','run_result.json','OUTGEN','recovery_control.json','QSTARS_MOM_FAILURE','QSTARS_VAR_FAILURE','QSTARS_MOM_GROWTH','QSTARS_J_FAILURE'] if (p/n).exists()})
        for key,pattern in [('peak_footprint_bytes',r'(\d+)\s+peak memory footprint'),('maximum_resident_bytes',r'(\d+)\s+maximum resident set size')]:
            m=re.search(pattern,log);row[key]=int(m[1]) if m else None
        rows.append(row)
    result=dict(cases=rows,measured_native_wall_sum_seconds=sum(r['seconds'] or 0 for r in rows),timing_scope='sum of measured processes; concurrent cases, not elapsed campaign; debugger unmeasured',production_wind=False,driver_sha256=sha256(__file__))
    (output/'campaign_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    for ax,names,title in zip(axs,[['T9_fixedJ_ff','T9_refined_full_ff_fixedJ'],['T9_full_ff','T9_refined_full_ff']],['Fixed J and T: initialization only','Released J and T: failed closure']):
        for short in names:
            row=next(r for r in rows if r['case']=='cmfgen_recovery_'+short)
            h=row['correction_history'];y=[v['all_level_percent'] for v in h];x=[v['iteration'] for v in h]
            ax.semilogy(x,y,'o-',markersize=3,label=f"{row['ND']} radii")
        ax.axhline(.1,color='black',ls='--',label='0.1% target');ax.set(xlabel='Native iteration counter',ylabel='All-level correction (%)',title=title);ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(output/'rate_recovery.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(trials=len(rows),wall_sum=result['measured_native_wall_sum_seconds']),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);report(**vars(p.parse_args()))
