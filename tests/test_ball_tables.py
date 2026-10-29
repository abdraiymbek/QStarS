"""Test U1: the printed Ball targets and the registered Gate 6 plan, as records.

Ball et al. (2011) Table 1 is checked against the boundary laws of Ball (2012) thesis
eqs 3.2 (r_B), 3.7 (M_cav), 3.15 (accretion, 16 pi) and 2.15 (Teff), using
ideal-gas-plus-radiation Gamma1 at the printed (rho0, T0). No QStarS model is built.
"""
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
import yaml
from qstars.eos import IdealGasRadiation,Composition
from qstars.constants import G,C,PI,MSUN,RSUN,LSUN,YEAR,SIGMA
from qstars.engine import rates

ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'benchmarks'/'ball2011_table1.csv'
THESIS=ROOT/'benchmarks'/'ball2012_thesis_table3_1.csv'
PLAN=ROOT/'benchmarks'/'ball2011_track_plan.json'
REGISTRATION=ROOT/'benchmarks'/'ball2011_track_plan.registration.json'
MANIFEST=ROOT/'benchmarks'/'manifest.yaml'
ETA=EPS=.1
EPS_PRIME=EPS/(1-EPS)
EOS=IdealGasRadiation(Composition(.7,.3,0.))
# Plan observable -> CSV column (both CSVs share the paper's column order).
COLUMNS={'t':'time_myr','Mdot_BH':'bh_growth_1e-4_msun_year','M_cav':'cavity_mass_msun',
         'L*':'luminosity_1e8_lsun','rho0':'rho_inner_cgs','T0':'T_inner_1e5_K',
         'Teff':'Teff_1e3_K','r_B':'r_inner_100_rsun','R*':'r_outer_1e4_rsun'}

def table(path):
    with open(path,newline='') as f: return list(csv.DictReader(f))

def implied(row):
    """Printed row -> boundary laws evaluated at the printed (M_BH, rho0, T0, r_B, L*, R*)."""
    mbh=float(row['bh_mass_msun'])*MSUN; rho=float(row['rho_inner_cgs'])
    q=EOS.evaluate(np.array([rho]),np.array([float(row['T_inner_1e5_K'])*1e5]))
    gamma1,cs=float(q['gamma1'][0]),float(q['cs'][0])
    r_inner=float(row['r_inner_100_rsun'])*100*RSUN
    radius=float(row['r_outer_1e4_rsun'])*1e4*RSUN
    luminosity=float(row['luminosity_1e8_lsun'])*1e8*LSUN
    eq3_15=16*PI*ETA/(EPS_PRIME*gamma1)*(G*mbh)**2*rho/(cs*C**2)
    return dict(mdot=eq3_15*YEAR/MSUN/1e-4,
                bondi=2*G*mbh/cs**2/(100*RSUN),
                cavity=8*PI/3*rho*r_inner**3/MSUN,
                teff=(luminosity/(4*PI*radius**2*SIGMA))**.25/1e3)

def log_ratio(x,row,column): return math.log(x/float(row[column]))

def check_row(row,label):
    q=implied(row)
    assert abs(log_ratio(q['bondi'],row,'r_inner_100_rsun'))<=.035,label
    if float(row['cavity_mass_msun'])==0:
        assert q['cavity']<.005,label  # printed 0.00
    else:
        assert abs(log_ratio(q['cavity'],row,'cavity_mass_msun'))<=.025,label
    assert abs(log_ratio(q['teff'],row,'Teff_1e3_K'))<=.005,label
    return q

def test_u1_table1_satisfies_thesis_boundary_laws_and_rejects_printed_4pi():
    rows=table(PAPER)
    assert [r['bh_mass_msun'] for r in rows]==['5','100','200','500','1000','1194']
    for row in rows:
        q=check_row(row,row['bh_mass_msun'])
        assert abs(log_ratio(q['mdot'],row,'bh_growth_1e-4_msun_year'))<=.02
        # Paper eq 15 prints 4 pi: a factor 4 below the 16 pi law, contradicted by the table.
        assert float(row['bh_growth_1e-4_msun_year'])/(q['mdot']/4)>=3.5

def test_u1_engine_16pi_mode_is_thesis_eq3_15_at_table_states():
    config=dict(model='quest2025_eq11',efficiency=EPS,coupling=ETA)
    for row in table(PAPER):
        state=SimpleNamespace(rho=np.array([float(row['rho_inner_cgs'])]),
                              T=np.array([float(row['T_inner_1e5_K'])*1e5]),
                              bh_mass=float(row['bh_mass_msun'])*MSUN,time=0.)
        r=rates(state,EOS,config)
        q=implied(row)
        assert r['bh_growth_rate']*YEAR/MSUN/1e-4==pytest.approx(q['mdot'],rel=1e-13)
        assert r['bondi_radius']/(100*RSUN)==pytest.approx(q['bondi'],rel=1e-13)
        printed=rates(state,EOS,dict(config,model='ball2011_printed'))
        assert printed['bh_growth_rate']==pytest.approx(r['bh_growth_rate']/4,rel=1e-13)

def test_thesis_table_same_checks_and_recorded_row0_rate_departure():
    rows=table(THESIS)
    assert [r['bh_mass_msun'] for r in rows]==['5','100','200','500','1000','1194']
    for row in rows:
        q=check_row(row,row['bh_mass_msun'])
        d=log_ratio(q['mdot'],row,'bh_growth_1e-4_msun_year')
        if row['bh_mass_msun']=='5':
            assert -.037<d<-.035  # recorded in plan sources.targets.consistency
        else:
            assert abs(d)<=.02

def plan():
    return json.loads(PLAN.read_text())

def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def test_plan_realization_spread_is_measured_from_both_tables():
    paper,thesis=table(PAPER),table(THESIS)
    assert [r['bh_mass_msun'] for r in paper]==[r['bh_mass_msun'] for r in thesis]
    columns=plan()['comparison']['columns']
    assert set(columns)==set(COLUMNS)
    for name,column in COLUMNS.items():
        spread=max(abs(math.log(float(a[column])/float(b[column])))
                   for a,b in zip(paper,thesis) if float(a[column])!=0)
        assert abs(spread-columns[name]['realization_spread'])<=5e-5,name

def test_plan_is_registered_and_consistent_with_table1():
    p=plan()
    registration=json.loads(REGISTRATION.read_text())
    assert registration['plan_sha256']==sha256(PLAN)
    # The tracked manifest anchors the untracked registration: a re-registration shows in git diff HEAD.
    entry={b['id']:b for b in yaml.safe_load(MANIFEST.read_text())['benchmarks']}['ball_stellar_tracks']
    assert entry['plan_sha256']==registration['plan_sha256']
    for name,digest in registration['target_files_sha256'].items():
        assert sha256(ROOT/name)==digest,name
    assert set(registration['target_files_sha256'])=={'benchmarks/ball2011_table1.csv',
                                                      'benchmarks/ball2012_thesis_table3_1.csv'}
    assert 'commit_rule' not in p['registration']
    cells=p['comparison']['cells']
    assert len(cells['row_0'])+4*len(cells['rows_1_to_4'])==cells['count']==44
    last=table(PAPER)[-1]
    assert p['endpoint']['ball']['M_BH_msun']==float(last['bh_mass_msun'])
    assert p['endpoint']['ball']['t_myr']==float(last['time_myr'])
    assert p['model']['initial_bh_mass_msun']==float(table(PAPER)[0]['bh_mass_msun'])

def printed_rows(text,first,last):
    """Rows of numbers between the cited 1-based lines; 'a x 10^-n' -> 'ae-n'."""
    out=[]
    for line in text.splitlines()[first-1:last]:
        line=re.sub(r'(\d+\.\d+)\s*×\s*10−(\d+)',r'\1e-\2',line)
        if re.fullmatch(r'[0-9.e\- ]+',line): out.append(line.split())
    return out

def cited_lines(entry):
    a,b=re.search(r'lines (\d+)-(\d+)',entry).groups()
    return int(a),int(b)

def test_csv_transcriptions_match_cited_source_lines():
    s=plan()['sources']
    paper_text=ROOT/s['paper']['cache']; thesis_text=ROOT/s['thesis']['cache']
    if not (paper_text.exists() and thesis_text.exists()):
        pytest.skip('source text cache absent (sources/cache is not distributed)')
    order=list(table(PAPER)[0])
    rows=printed_rows(paper_text.read_text(encoding='utf-8'),*cited_lines(s['paper']['table']))
    assert rows==[[r[k] for k in order] for r in table(PAPER)]
    rows=printed_rows(thesis_text.read_text(encoding='utf-8'),*cited_lines(s['thesis']['table']))
    top,bottom=rows[:6],rows[6:]
    assert len(bottom)==6 and all(a[:2]==b[:2] for a,b in zip(top,bottom))
    # Thesis blocks: (t, M_BH, Mdot, M_cav, T0, r0) and (t, M_BH, L*, rho0, Teff, R*).
    assert [[a[0],a[1],a[2],a[3],b[2],b[3],a[4],b[4],a[5],b[5]] for a,b in zip(top,bottom)]==\
           [[r[k] for k in order] for r in table(THESIS)]
