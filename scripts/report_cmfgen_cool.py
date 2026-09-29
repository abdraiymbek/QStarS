#!/usr/bin/env python3
"""Summarize retained cool trials without converting failed solves into data."""
import argparse,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qstars.cmfgen_backend import sha256

def report(output):
    rows=[]
    for p in sorted(Path('results').glob('cmfgen_cool_*')):
        if not (p/'run_manifest.json').exists():continue
        if not (p/'run_result.json').exists():raise RuntimeError('Trial is still active: '+str(p))
        result=json.loads((p/'run_result.json').read_text())
        text=(p/'OUTGEN').read_text() if (p/'OUTGEN').exists() else ''
        log=(p/'native.log').read_text()
        counters=[];corrections=[];current=None;lum=[]
        for line in text.splitlines():
            m=re.search(r'Current great iteration count is\s+(\d+)',line)
            if m:current=int(m[1])
            m=re.search(r'Maximm changes as returned by SOLVEBA_V13 is\s+(\S+)',line)
            if m and current is not None:
                value=float(m[1].replace('D','E'));counters.append(current);corrections.append(value if np.isfinite(value) else None)
            if 'Luminosity of star (d=1,ND)' in line:lum.append(line.strip())
        nonfinite=bool(re.search(r'(?<![A-Za-z])(?:NaN|Infinity)(?![A-Za-z])',text,re.I))
        termination=(p/'termination_reason.txt').read_text() if (p/'termination_reason.txt').exists() else None
        if nonfinite:classification='nonfinite_native_state; any zero correction is invalid'
        elif 'blowing up' in text or 'excessive iteration count' in text:classification='native_radiation_moment_instability'
        elif result['status']=='numerical_failure':classification='native_setup_or_iteration_failure'
        elif corrections and corrections[-1] is not None and corrections[-1]>.1:classification='completed_output_but_nonlinear_gate_failed'
        else:classification='component_snapshot_only; no production certificate'
        row=dict(case=p.name,directory=str(p),status=result['status'],classification=classification,error=result['error'],
            seconds=result['seconds'],native_exit=result.get('exit_code'),last_counter=current,
            iteration_counters=counters,correction_percent=corrections,nonfinite_native_diagnostic=nonfinite,
            luminosity_history_native=lum,termination=termination,executable_sha256=result['executable_sha256'],
            hashes={name:sha256(p/name) for name in ['run_manifest.json','run_result.json','OUTGEN','CORRECTION_LINK','CORRECTION_SUM'] if (p/name).exists()})
        for label,pattern in [('peak_footprint_bytes',r'(\d+)\s+peak memory footprint'),('maximum_resident_bytes',r'(\d+)\s+maximum resident set size')]:
            m=re.search(pattern,log);row[label]=int(m[1]) if m else None
        rows.append(row)
    result=dict(scope='cool prescribed-flow controls, including setup failures, instability and diagnostic-only cases',
        cases=rows,production_wind=False,published_evolution_runs=0,
        total_native_wall_seconds=sum(r['seconds'] for r in rows),
        total_wall_scope='sum of per-process wall times; cases overlap, not campaign elapsed time',
        driver_sha256=sha256(__file__))
    output.mkdir(parents=True,exist_ok=True)
    (output/'campaign_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    fig,ax=plt.subplots(1,2,figsize=(11.5,4.5),layout='constrained')
    chosen=['cmfgen_cool_T9_corrected','cmfgen_cool_T6_corrected','cmfgen_cool_T6_refined','cmfgen_cool_T9_newton']
    for row in rows:
        if row['case'] not in chosen:continue
        c=np.array([np.nan if x is None or x<=0 else x for x in row['correction_percent']])
        ax[0].semilogy(row['iteration_counters'],c,label=row['case'].replace('cmfgen_cool_',''))
    ax[0].axhline(.1,color='black',ls='--',lw=1,label='0.1% component target')
    ax[0].set(xlabel='Native iteration counter',ylabel='Maximum correction (%)',title='Cool nonlinear solves — failed trials retained')
    ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
    audit=output/'bootstrap_audit/cool_profiles.npz'
    if audit.exists():
        q=np.load(audit);r=q['radius']/q['radius'][0]
        ax[1].semilogx(r,q['material_temperature'],color='tab:red',label='Gas T (initialization only)')
        ax[1].axhspan(3000,12000,color='tab:red',alpha=.06)
        ax[1].set(xlabel='r / inner radius',ylabel='Gas temperature (K)',title='Recombination coverage ≠ converged closure')
        twin=ax[1].twinx();twin.semilogx(r,q['neutral_hydrogen_fraction'],color='tab:blue',label='Neutral H fraction')
        twin.set(ylabel='Neutral H fraction',ylim=(-.02,1.02))
        ax[1].legend(loc='upper right',fontsize=8);twin.legend(loc='center right',fontsize=8)
    fig.savefig(output/'cool_validation.png',dpi=180);fig.savefig(output/'cool_validation.pdf');plt.close(fig)
    print(f'{len(rows)} native trials recorded')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);report(**vars(p.parse_args()))
