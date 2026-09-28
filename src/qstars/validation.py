"""Reproducible evidence campaigns. Pending full-track benchmarks remain pending."""
import copy,csv,json,time
from pathlib import Path
import numpy as np
from . import io
from .config import resolve
from .evolution import run,physics
from .structure import initialize
from .reference import lane_emden,ball_limit,corrected_limit,saturated_inner,saturated_composite

def campaign(parent='results',convergence=True):
    out=io.create_directory(parent,'validation')
    checks=[]
    def record(name,actual,target,atol,source,level):
        checks.append(dict(name=name,actual=float(actual),reference=target,absolute_error=abs(float(actual)-target),
                      tolerance=atol,passed=bool(abs(actual-target)<=atol),source=source,level=level))
    l=lane_emden(1)
    record('Lane-Emden n=1 radius',l.t[-1],np.pi,1e-8,'analytic sin(x)/x','A')
    for n,ref in [(2,.127),(3,.102),(4,.089)]:
        d=ball_limit(n)
        record(f'Ball 2011 n={n} maximum BH fraction',d['bh_fraction'],ref,.0005,'ball2011 section 4.1','C-structure')
    for N,ref,tol in [(1,.017,.0001),(2,.0058,.0001),(1/3,.13,.005)]:
        d=corrected_limit(N)
        record(f'CB24 corrected Bondi N={N:g}',d['bh_fraction'],ref,tol,'coughlin2024 section 2, Figure 1 prose','C-structure')
    sat,_=saturated_inner()
    record('CB24 saturated entropy minimum radius',sat['x_entropy_min'],1.015,.001,'coughlin2024 section 3.1','C-structure')
    # Tight printed-value comparisons intentionally remain failing if inconsistent.
    record('CB24 printed Ki_max',sat['Ki_max'],.46554,.00005,'coughlin2024 section 3.1','C-open-discrepancy')
    record('CB24 printed BH/inner mass',sat['bh_inner_fraction'],.6321,.0001,'coughlin2024 section 3.1','C-open-discrepancy')
    _,inner=saturated_inner()
    m_at_one=float(inner.sol(0.)[0])
    record('CB24 printed cutoff Ki (x_tilde=1)',(1/m_at_one)**(5/3),.46554,.00001,'coughlin2024 printed cutoff interpretation','C-structure')
    record('CB24 printed cutoff BH/inner mass',1/m_at_one,.6321,.0001,'coughlin2024 printed cutoff interpretation','C-structure')
    printed=saturated_composite(1.)
    record('CB24 composite maximum at printed cutoff',printed['bh_fraction'],.618,.001,'coughlin2024 section 3.2','C-structure')
    curves=[saturated_composite(x) for x in np.linspace(.15,sat['x_entropy_min']*.999,30)]
    io.json_write(out/'saturated_sequence.json',curves)
    # Time integration/conservation with simultaneous mass streams.
    base=resolve({'name':'two_stream','engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4},
                  'feeding':{'rate_msun_per_year':.01,'thermodynamics':'surface'},
                  'wind':{'model':'prescribed','rate_msun_per_year':.008},
                  'evolution':{'end_time_years':.01},'numerics':{'zones':32,'max_dt_years':.0001,'initial_dt_years':.0001}})
    p,report=run(base,out)
    good=report['status']=='completed_requested_interval'
    checks.append(dict(name='Two-stream evolution completes',passed=good,actual=report['status'],level='A-integration'))
    if 'conservation' in report:
        for k in ('mass_residual_relative','energy_residual_relative'):
            record(k,report['conservation'][k],0.,1e-8,'independent ledger','A-integration')
    convergence_rows=[]
    if convergence:
        experiments=[('mesh',n,1e-4,1e-10) for n in (16,32,64,128)]
        experiments += [('timestep',32,dt,1e-10) for dt in (4e-4,2e-4,1e-4,5e-5)]
        experiments += [('solver',32,1e-4,tol) for tol in (1e-7,1e-9,1e-11)]
        experiments += [('confined_mesh',n,1e-4,1e-10) for n in (16,32,64,128)]
        for axis,n,dt,tol in experiments:
            c=copy.deepcopy(base);c['name']=f'{axis}_{n}_{dt}_{tol}'
            c['numerics'].update(zones=n,initial_dt_years=dt,max_dt_years=dt,solver_tolerance=tol)
            if axis=='confined_mesh':
                c['boundary']['pressure_cgs']=1e8
                c['boundary']['thermal']='luminosity'
                c['boundary']['luminosity_erg_s']=1e42
                c['numerics']['mesh']='core_refined'
            p,r=run(c,out)
            row=dict(axis=axis,zones=n,dt_years=dt,solver_tolerance=tol,path=str(p),**r)
            if (p/'history.csv').exists():
                h=np.atleast_1d(np.genfromtxt(p/'history.csv',delimiter=',',names=True))
                for key in ('radius_cm','luminosity_erg_s','bh_mass_msun','internal','gravity','kinetic','energy_residual_relative'):
                    row[key]=float(h[key][-1])
            convergence_rows.append(row)
            checks.append(dict(name=f"{axis} n={n}, dt={dt:g}, tol={tol:g} completes",
                passed=r['status']=='completed_requested_interval',actual=r['status'],level='B'))
        for axis in ('mesh','timestep','solver','confined_mesh'):
            group=[r for r in convergence_rows if r['axis']==axis]
            reference=group[-1]
            for row in group:
                comparable=(reference['status']==row['status']=='completed_requested_interval' and
                            abs(row['final_time_years']-reference['final_time_years'])<1e-12)
                row['comparable']=comparable
                row['relative_difference_to_finest']=({k:abs(row[k]-reference[k])/max(abs(reference[k]),1.)
                    for k in ('radius_cm','luminosity_erg_s','internal')} if comparable else {})
    summary=dict(checks=checks,convergence=convergence_rows,v1_validated=False,
        pending=['Full Ball 2011 Table 1 evolutionary replication','MESA-QUEST thermal tracks',
                 'Coupling stationary Fiacconi & Rossi winds to a conserved evolving envelope','Cantiello nonadiabatic stability',
                 'Production atmosphere and stellar evolution with audited MESA microphysics','Moving cavity/Bondi/saturated inner evolution',
                 'Long-duration quasi-static sequences and QSRD demonstration'])
    io.json_write(out/'validation.json',summary)
    lines=['# QStarS validation evidence','', '**QStarS v1 is not validated.** This campaign separates numerical identities from publication replication.','',
        '| Check | Result | Measured | Reference | Absolute error |','|---|---|---:|---:|---:|']
    for x in checks:
        lines.append(f"| {x['name']} | {'PASS' if x['passed'] else 'FAIL / investigate'} | {x.get('actual')} | {x.get('reference','—')} | {x.get('absolute_error','—')} |")
    lines+=['','## Convergence measurements','',
        '| Axis | Zones | dt [yr] | Solver tolerance | Status | Radius difference | Luminosity difference | Energy residual / U₀ |',
        '|---|---:|---:|---:|---|---:|---:|---:|']
    for x in convergence_rows:
        d=x.get('relative_difference_to_finest',{})
        lines.append(f"| {x['axis']} | {x['zones']} | {x['dt_years']:g} | {x['solver_tolerance']:g} | {x['status']} | {d.get('radius_cm',float('nan')):.3g} | {d.get('luminosity_erg_s',float('nan')):.3g} | {x.get('energy_residual_relative',float('nan')):.3g} |")
    lines+=['','The finest run is a comparison reference, not an exact solution. Small conservation residuals do not establish physical convergence.',
            '', '## Open gates','']+['- '+x for x in summary['pending']]
    (out/'validation.md').write_text('\n'.join(lines)+'\n')
    return out,summary
