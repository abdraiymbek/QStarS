#!/usr/bin/env python3
"""Summarize native component gates and performance without promoting production."""
import argparse,json,re,platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def analyze(evidence):
    audits={p.parent.name:json.loads(p.read_text()) for p in sorted(evidence.glob('*/validation.json'))}
    rows=[]
    for case in audits:
        d=Path('results')/('cmfgen_'+case if case.startswith(('hhe','zpup')) else 'cmfgen_'+case)
        if not (d/'native.log').exists():continue
        log=(d/'native.log').read_text()
        wall=re.search(r'([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys',log)
        rss=re.search(r'(\d+)\s+maximum resident set size',log)
        peak=re.search(r'(\d+)\s+peak memory footprint',log)
        rows.append(dict(case=case,wall_seconds=float(wall[1]) if wall else None,user_seconds=float(wall[2]) if wall else None,
             sys_seconds=float(wall[3]) if wall else None,max_resident_bytes=int(rss[1]) if rss else None,
             peak_footprint_bytes=int(peak[1]) if peak else None,threads=1,concurrent_work=True))
    summary=dict(outcome='BACKEND/PHYSICS BLOCKER',production_certified=False,full_evolution_replicas=0,qsrd_run=False,
        backend='CMFGEN 18-Jun-2025',case_checks={k:dict(passed=v['passed'],failed=v['failed']) for k,v in audits.items()},
        performance=rows,performance_scope='Mac M5 16 GiB, one OpenMP thread per process; concurrent work, no isolated scaling measurement',
        remaining=['independent complete atomic/thermal reference certification','material kinetic energy and derivatives consistent with populations and dissolution',
            'converged cold/recombining radius/frequency/atom-size campaign','global CMF force/source response and wind eigenvalue',
            'valid convective-envelope handoff and conservative escape','physical coupled evolution and complete published tracks'])
    from scipy.integrate import cumulative_trapezoid
    summary['interior_energy_diagnostics']={}
    for name in audits:
        q=np.load(evidence/name/'profiles.npz');r=q['radius'];L=q['luminosity_lab']
        exchange=4*np.pi*cumulative_trapezoid(r*r*q['energy_to_matter_lab'],r,initial=0)
        norm=max(abs(L[0]),abs(L[-1]));defect=abs(L-L[0]+exchange)/norm
        summary['interior_energy_diagnostics'][name]=dict(max_cumulative_relative=float(max(defect)),
            endpoint_relative=float(defect[-1]),scope='supplementary interior diagnostic; passing endpoint target does not bound every subvolume')
    (evidence/'campaign_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for name,color in [('hhe_2','#d27c2c'),('hhe_3','#286090'),('thermal_tight','#8b5da8'),('thermal_rel','#338865'),('thermal_equilibrium','#b04444')]:
        if name not in audits:continue
        d=Path('results')/('cmfgen_'+name)
        if not d.exists():continue
        text=(d/'OUTGEN').read_text();its=[];cor=[];last=None
        for line in text.splitlines():
            m=re.search(r'Current great iteration count is\s+(\d+)',line)
            if m:last=int(m[1])
            m=re.search(r'Maximm changes as returned by SOLVEBA_V13 is\s+([0-9.Ee+\-]+)',line)
            if m:its.append(last);cor.append(float(m[1])/100)
        # Continuation OUTGEN includes earlier iterations; plot only its new part.
        start=75 if name=='thermal_equilibrium' else (55 if name.startswith('thermal_') else (25 if name=='hhe_3' else 0))
        keep=np.array(its)>start
        axs[0,0].semilogy(np.array(its)[keep],np.array(cor)[keep],'.-',label=name,color=color)
        if not (evidence/name/'profiles.npz').exists():continue
        q=np.load(evidence/name/'profiles.npz');r=q['radius']/q['radius'][0]
        axs[0,1].semilogx(r,q['material_temperature'],label=name,color=color)
        neutral=q['population_HI'].sum(axis=0)/q['material_hydrogen_density']
        axs[1,0].loglog(r,neutral,label=name,color=color)
    axs[0,0].axhline(.001,color='black',ls='--',label='0.1% gate');axs[0,0].set(xlabel='Native iteration',ylabel='Maximum fractional correction',title='Nonlinear convergence')
    axs[0,1].set(xlabel='r / inner radius',ylabel='Gas temperature [K]',title='Internal material profile')
    axs[1,0].set(xlabel='r / inner radius',ylabel='Neutral hydrogen fraction',title='Hot H/He component — not a cool quasi-star')
    names=[];defects=[]
    for n,v in audits.items():
        entry=next((c for c in v['checks'] if c['name']=='global_radiation_energy'),None)
        if entry and entry['value'] is not None:names.append(n);defects.append(100*entry['value'])
    axs[1,1].barh(names,defects,color=['#b24b4b' if d>.1 else '#3a8470' for d in defects]);axs[1,1].set_xscale('log')
    axs[1,1].axvline(.1,color='black',ls='--');axs[1,1].set(xlabel='Independent radiation-energy defect [%]',title='0.1% gate; no luminosity renormalization')
    for ax in axs.flat:ax.grid(alpha=.2)
    for ax in axs.flat[:3]:ax.legend(fontsize=8)
    fig.suptitle('QStarS / CMFGEN component audit — production wind remains gated',fontsize=13)
    fig.savefig(evidence/'component_validation.png',dpi=170);fig.savefig(evidence/'component_validation.pdf');plt.close(fig)
    if (evidence/'mesh_665_fixed/profiles.npz').exists():
        reference=np.load(evidence/'mesh_665_fixed/profiles.npz')
        fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        for name in ['mesh_84','mesh_167','mesh_333','mesh_665_fixed']:
            q=np.load(evidence/name/'profiles.npz');r=q['radius'];x=r/r[0]
            force=np.interp(r,reference['radius'],reference['radiative_acceleration'])
            axs[0].semilogx(x,100*(q['radiative_acceleration']/force-1),label=name)
            L=q['luminosity_lab'];exchange=4*np.pi*cumulative_trapezoid(r*r*q['energy_to_matter_lab'],r,initial=0)
            defect=100*(L-L[0]+exchange)/max(abs(L[0]),abs(L[-1]))
            axs[1].semilogx(x,defect,label=name)
        axs[0].set(xlabel='r / inner radius',ylabel='Force difference from 665-point control [%]',title='Internal force refinement')
        axs[1].set(xlabel='r / inner radius',ylabel='Cumulative radiation-energy defect [%]',title='Independent luminosity/source balance')
        for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=8)
        fig.suptitle('Frozen incoming hot H/He state: spatial/angular transfer only')
        fig.savefig(evidence/'mesh_profiles.png',dpi=170);fig.savefig(evidence/'mesh_profiles.pdf');plt.close(fig)
    print(summary)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);analyze(**vars(p.parse_args()))
