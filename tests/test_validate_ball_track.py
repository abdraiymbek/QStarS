"""Gate 6 Stage 5: scripts/validate_ball_track.py without MESA (no Ball configuration is run).

The registration gate (plan and target sha256, amendments), the registered
run configurations, the 44 comparison cells and their effective tolerances
against the printed table, the qualitative checks, the campaign thresholds
(pairs, Richardson, conservation), the verdict rules, and the restart
comparison on an analytic-cavity toy run. Stage 6 review: a relaxation that
ends in an error fails its seed only; relax.json is always written.
"""
import copy,importlib.util,json,math,shutil,sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('validate_ball_track',ROOT/'scripts/validate_ball_track.py')
vbt=importlib.util.module_from_spec(spec);spec.loader.exec_module(vbt)
PLAN=json.loads(vbt.PLAN.read_text())

def registered_copy(tmp_path):
    for name in ('benchmarks/ball2011_track_plan.json','benchmarks/ball2011_track_plan.registration.json',
                 'benchmarks/ball2011_table1.csv','benchmarks/ball2012_thesis_table3_1.csv'):
        (tmp_path/name).parent.mkdir(parents=True,exist_ok=True);shutil.copy(ROOT/name,tmp_path/name)
    return tmp_path/'benchmarks/ball2011_track_plan.json',tmp_path/'benchmarks/ball2011_track_plan.registration.json'

# ---------------------------------------------------------------- registration and amendments

def test_driver_refuses_unless_plan_and_targets_match_the_registration(tmp_path):
    record=vbt.verify_registration()
    assert record['plan_sha256']==json.loads(vbt.REGISTRATION.read_text())['plan_sha256']
    plan,registration=registered_copy(tmp_path)
    assert vbt.verify_registration(plan,registration,tmp_path)['plan_sha256']==record['plan_sha256']
    plan.write_bytes(plan.read_bytes()+b' ')
    with pytest.raises(vbt.Refused,match='plan.json sha256'):vbt.verify_registration(plan,registration,tmp_path)
    # The command line refuses before it creates any run directory.
    with pytest.raises(SystemExit,match='Refusing to run'):
        vbt.main(['--plan',str(plan),'--registration',str(registration),'--output',str(tmp_path/'out'),'compare'])
    assert not (tmp_path/'out').exists()
    shutil.copy(vbt.PLAN,plan)
    table=tmp_path/'benchmarks/ball2011_table1.csv';table.write_text(table.read_text().replace('8.71e-5','8.72e-5'))
    with pytest.raises(vbt.Refused,match='ball2011_table1.csv'):vbt.verify_registration(plan,registration,tmp_path)

FAILED_SMOKE='benchmarks/results/2026-10-08-ball-track/20261008T000000Z_smoke_C1_n200'

def failing_run(root,relative=FAILED_SMOKE,record='smoke.json',data=None):
    d=root/relative;d.mkdir(parents=True,exist_ok=True)
    (d/record).write_text(json.dumps(dict(passed=False) if data is None else data))
    return relative

def amendment(tmp_path,settings,plan_sha256=None,status='registered',run=FAILED_SMOKE):
    path=tmp_path/'ball2011_track_plan_amendment_1.json'
    path.write_text(json.dumps(dict(amendment_id='test',status=status,preserved_failing_run=run,settings=settings)))
    record=dict(amendment_sha256=vbt.sha256(path),plan_sha256=plan_sha256 or vbt.verify_registration()['plan_sha256'],
                registered_utc='2026-10-08T00:00:00Z')
    path.with_name(path.stem+'.registration.json').write_text(json.dumps(record))
    return path

def test_amendments_are_verified_and_limited_to_numerical_method_settings(tmp_path):
    sha=vbt.verify_registration()['plan_sha256'];failing_run(tmp_path)
    verify=lambda path:vbt.verify_amendment(path,sha,root=tmp_path)
    path=amendment(tmp_path,{'stepcontrol':{'metric':'mapped'},'quasistatic':{'base':'fixed_mass','base_mass_fraction':.05}})
    record,settings=verify(path)
    c=vbt.build_config(PLAN,'C0','lib','mesa',[(record,settings)])
    assert c['stepcontrol']['metric']=='mapped' and c['quasistatic']['base']=='fixed_mass' and c['replication']['amendments']==[str(path)]
    assert c['stepcontrol']['delta_target']==.005 and c['relaxation']['mesh_cycles']==0  # untouched
    assert record['preserved_failing_run']==[FAILED_SMOKE]
    path.write_text(path.read_text().replace('mapped','lagrangian'))
    with pytest.raises(vbt.Refused,match='differs from its registered sha256'):verify(path)
    with pytest.raises(vbt.Refused,match='another plan'):verify(amendment(tmp_path,{'stepcontrol':{'metric':'mapped'}},'0'*64))
    for bad in ({'stepcontrol':{'delta_target':.01}},{'cavity':{'profile':'none'}},{'quasistatic':{'tolerance':1e-8}},
                {'relaxation':{'tolerance':1e-8}},{}):
        with pytest.raises(vbt.Refused):verify(amendment(tmp_path,bad))
    lone=tmp_path/'lone.json';lone.write_text('{}')
    with pytest.raises(vbt.Refused,match='no registration record'):verify(lone)
    # The relaxation mesh cycles are off unless an amendment sets them, and they enter the row-0 signature.
    relax_record,relax_settings=verify(amendment(tmp_path,{'relaxation':{'mesh_cycles':40,'mesh_tolerance':1e-5}}))
    amended=vbt.build_config(PLAN,'C0','lib','mesa',[(relax_record,relax_settings)])
    assert amended['relaxation']['mesh_cycles']==40 and vbt.signature(amended)!=vbt.signature(vbt.build_config(PLAN,'C0','lib','mesa'))
    # Amendment 3: the seed surface and the seed radii (plan initial_model.seed, not an immutable section) are
    # amendable and enter the row-0 signature; nothing else of the replication section is.
    radii=[2000,6000,1000,12000,3464,1414,8485]
    seed_record,seed_settings=verify(amendment(tmp_path,{'replication':{'seed_surface':'grey','seed_radii_rsun':radii}}))
    seeded=vbt.build_config(PLAN,'C1','lib','mesa',[(seed_record,seed_settings)])
    plain=vbt.build_config(PLAN,'C1','lib','mesa')
    assert (plain['replication']['seed_surface'],plain['replication']['seed_radii_rsun'])==('eos',[2000,6000,1000,12000])
    assert (seeded['replication']['seed_surface'],seeded['replication']['seed_radii_rsun'])==('grey',radii)
    assert vbt.signature(seeded)!=vbt.signature(plain)
    for bad in ({'replication':{'total_mass_msun':1e5}},{'replication':{'policy':'production'}},{'replication':{'landings_msun':[50]}},
                {'replication':{'seed_surface':'grey','initial_bh_mass_msun':10.}},{'quasistatic':{'ptc':True}}):
        with pytest.raises(vbt.Refused):verify(amendment(tmp_path,bad))
    with pytest.raises(ValueError,match='seed_surface'):  # an unknown seed surface does not resolve
        vbt.build_config(PLAN,'C1','lib','mesa',[(seed_record,{'replication':{'seed_surface':'hot'}})])
    # Amendment 4: the bracket rule of photosphere.solve (registered: 'doubling'; 'scan' only when amended). It
    # changes which seeds can be built and the first step from a seed, so it enters the row-0 signature.
    bracket_record,bracket_settings=verify(amendment(tmp_path,{'quasistatic':{'photosphere_bracket':'scan'}}))
    scanned=vbt.build_config(PLAN,'C1','lib','mesa',[(bracket_record,bracket_settings)])
    assert plain['quasistatic']['photosphere_bracket']=='doubling' and scanned['quasistatic']['photosphere_bracket']=='scan'
    assert vbt.signature(scanned)!=vbt.signature(plain)
    assert {k:v for k,v in vbt._flat(scanned).items() if v!=vbt._flat(plain)[k]}=={'quasistatic.photosphere_bracket':'scan',
                                                                                  'replication.amendments':[str(bracket_record['path'])]}
    with pytest.raises(ValueError,match='photosphere_bracket'):
        vbt.build_config(PLAN,'C1','lib','mesa',[(bracket_record,{'quasistatic':{'photosphere_bracket':'bisect'}})])

def test_amendments_must_be_registered_and_cite_a_preserved_failing_run(tmp_path):
    """Plan registration.amendment_rule: no draft, and the cited run directory exists and records a failure."""
    sha=vbt.verify_registration()['plan_sha256'];settings={'stepcontrol':{'metric':'mapped'}}
    verify=lambda path:vbt.verify_amendment(path,sha,root=tmp_path)
    failing_run(tmp_path)
    placeholder='TO BE FILLED BY STAGE 7: benchmarks/results/<date>-ball-track/<smoke run, registered numerics>'
    with pytest.raises(vbt.Refused,match='draft'):verify(amendment(tmp_path,settings,status='DRAFT, not registered'))
    with pytest.raises(vbt.Refused,match='draft or has no status'):verify(amendment(tmp_path,settings,status=None))
    with pytest.raises(vbt.Refused,match='not a run directory'):verify(amendment(tmp_path,settings,run=placeholder))
    with pytest.raises(vbt.Refused,match='cites no preserved failing run'):verify(amendment(tmp_path,settings,run=None))
    outside='benchmarks/results/elsewhere/run';(tmp_path/outside).mkdir(parents=True);(tmp_path/outside/'smoke.json').write_text('{}')
    with pytest.raises(vbt.Refused,match='not a run directory'):verify(amendment(tmp_path,settings,run=outside))
    empty='benchmarks/results/2026-10-08-ball-track/empty';(tmp_path/empty).mkdir(parents=True)
    with pytest.raises(vbt.Refused,match='no smoke.json'):verify(amendment(tmp_path,settings,run=empty))
    passed=failing_run(tmp_path,'benchmarks/results/2026-10-08-ball-track/passed_smoke',data=dict(passed=True))
    with pytest.raises(vbt.Refused,match='records no failure'):verify(amendment(tmp_path,settings,run=passed))
    unique=failing_run(tmp_path,'benchmarks/results/2026-10-08-ball-track/relax_ok','relax.json',dict(status='converged_unique'))
    with pytest.raises(vbt.Refused,match='records no failure'):verify(amendment(tmp_path,settings,run=unique))
    blocked=failing_run(tmp_path,'benchmarks/results/2026-10-08-ball-track/relax_bad','relax.json',dict(status='not_unique'))
    assert verify(amendment(tmp_path,settings,run=[FAILED_SMOKE,blocked]))[0]['preserved_failing_run']==[FAILED_SMOKE,blocked]
    # An interrupted or requested stop (control_stop) shows no defect: it cannot back an amendment.
    for name,record,data in (('relax_stop','relax.json',dict(status='control_stop')),
                             ('smoke_stop','smoke.json',dict(status='control_stop',passed=False)),
                             ('track_stop','termination.json',dict(status='control_stop'))):
        stopped=failing_run(tmp_path,f'benchmarks/results/2026-10-08-ball-track/{name}',record,data)
        with pytest.raises(vbt.Refused,match='records no failure'):verify(amendment(tmp_path,settings,run=stopped))
    for name,record,data in (('relax_blocked','relax.json',dict(status='relaxation_blocker')),
                             ('track_failed','termination.json',dict(status='termination_trigger'))):
        assert vbt.preserved_failure(failing_run(tmp_path,f'benchmarks/results/2026-10-08-ball-track/{name}',record,data),tmp_path) is True
    # The unregistered draft of Amendment 1 as it stands (placeholder run) is refused.
    draft=dict(amendment_id='ball2011_track_plan_amendment_1',status='registered',preserved_failing_run=placeholder,settings=settings)
    path=tmp_path/'a.json';path.write_text(json.dumps(draft))
    path.with_name('a.registration.json').write_text(json.dumps(dict(amendment_sha256=vbt.sha256(path),plan_sha256=sha)))
    with pytest.raises(vbt.Refused,match='preserved_failing_run'):verify(path)

def test_amendment_drafts_are_refused_until_registered_and_stay_inside_the_amendable_settings():
    """The four drafts in benchmarks/: each is a DRAFT (refused by the driver even with a registration
    record), names the registered plan, changes only AMENDABLE settings that resolve, and the preserved
    failing runs Amendment 3 cites are the registered relax runs that ended relaxation_blocker. Amendments
    4, 2 and 1 cite runs Stage 7 has yet to make (placeholders the driver refuses)."""
    sha=vbt.verify_registration()['plan_sha256']
    for n,sections in ((1,{'quasistatic','stepcontrol'}),(2,{'relaxation'}),(3,{'replication'}),(4,{'quasistatic'})):
        path=ROOT/f'benchmarks/ball2011_track_plan_amendment_{n}.draft.json'
        draft=json.loads(path.read_text())
        assert draft['amendment_id']==f'ball2011_track_plan_amendment_{n}' and draft['plan_sha256']==sha
        assert draft['status'].startswith('DRAFT') and set(draft['settings'])==sections
        assert all(set(values)<=vbt.AMENDABLE[section] for section,values in draft['settings'].items())
        assert not path.with_name(path.stem+'.registration.json').exists()  # registering is the main session's step
        c=vbt.build_config(PLAN,'C1','lib','mesa',[(dict(path=str(path)),draft['settings'])])
        assert c['replication']['amendments']==[str(path)]
    three=json.loads((ROOT/'benchmarks/ball2011_track_plan_amendment_3.draft.json').read_text())
    assert three['settings']=={'replication':{'seed_surface':'grey','seed_radii_rsun':[2000,6000,1000,12000,3464,1414,8485]}}
    registered=PLAN['initial_model']['seed_radii_rsun']+PLAN['initial_model']['fallback_seed_radii_rsun']
    radii=three['settings']['replication']['seed_radii_rsun']
    assert radii[:4]==registered  # the registered radii first, in the registered order
    # The stated rule: geometric means of neighbouring registered radii, primary pair first, then the fallbacks'.
    assert radii[4:]==[round(math.sqrt(a*b)) for a,b in ((2000,6000),(1000,2000),(6000,12000))]
    assert three['changes'][0]['registered']==PLAN['initial_model']['seed']
    assert len(three['preserved_failing_run'])==3
    for run in three['preserved_failing_run']:
        if (ROOT/run).is_dir():  # results are not part of a source checkout
            assert vbt.preserved_failure(run) is True
            assert json.loads((ROOT/run/'relax.json').read_text())['status']=='relaxation_blocker'
    four=json.loads((ROOT/'benchmarks/ball2011_track_plan_amendment_4.draft.json').read_text())
    assert four['settings']=={'quasistatic':{'photosphere_bracket':'scan'}}
    for draft in (four,json.loads((ROOT/'benchmarks/ball2011_track_plan_amendment_2.draft.json').read_text()),
                  json.loads((ROOT/'benchmarks/ball2011_track_plan_amendment_1.draft.json').read_text())):
        assert isinstance(draft['preserved_failing_run'],str) and draft['preserved_failing_run'].startswith('TO BE FILLED BY STAGE 7')
        assert vbt.preserved_failure(draft['preserved_failing_run'])=='not a run directory under benchmarks/results/<date>-ball-track/'
    for draft in (three,four):
        for name in draft['tests_failed_before_pass_after']+draft.get('tests_asserting_the_registered_defect',[]):
            file,test=name.split(' ')[0].split('::')
            assert f'def {test}(' in (ROOT/file).read_text(),name

def test_registered_run_configurations():
    c=vbt.build_config(PLAN,'C0','lib','mesa')
    assert c['numerics']['zones']==400 and c['numerics']['conservation_tolerance']==1e-8
    assert c['engine']['model']=='ball2012_eq3_15' and c['engine']['coupling']==.1 and c['engine']['efficiency']==.1
    assert c['transport']=={'convection':'mlt_ml1','mixing_length_alpha':2.}
    assert c['composition']=={'X':.7,'Y':.3,'Z':0.} and c['eos']=='mesa' and c['boundary']['thermal']=='stars_photosphere'
    assert c['cavity']=={'active':True,'profile':'r-3/2','b':1.,'boundary_state':'face'}
    q=c['quasistatic'];assert (q['tolerance'],q['split'],q['merge'],q['base'],q['max_passes'])==(1e-12,2.,.5,'lagrangian',3)
    assert q['energy_per_step']==pytest.approx(5e-13) and q['mass_per_step']==1e-14
    sc=c['stepcontrol'];assert (sc['delta_target'],sc['dt_initial_years'],sc['dt_max_years'],sc['metric'])==(.005,1.,2e4,'lagrangian')
    assert c['opacity']['tables']=='stars2004_switch' and c['replication']['landings_msun']==[100,200,500,1000]
    # Resolved and cross-validated by qstars.config (Stage 6): quasi-static mechanics with the Bondi cavity.
    assert (c['numerics']['mechanics'],c['boundary']['inner'],c['replication']['plan'])==('quasi_static','bondi_cavity','ball2011_track_plan')
    assert c['replication']['seed_radii_rsun']==[2000,6000,1000,12000] and c['evolution']['max_steps']==20000
    assert c['validity']['min_temperature_K']==3000. and c['validity']['max_mach']==.2
    # Registered relaxation (plan initial_model): no STARS-Q mesh cycles without an amendment.
    r=c['relaxation'];assert r==vbt.RELAXATION and r['mesh_cycles']==0
    assert (r['dt_years'],r['dt_max_years'],r['converged_dt_years'],r['tolerance'],r['luminosity_tolerance'])==(.01,1e8,1e6,1e-10,1e-6)
    other=copy.deepcopy(PLAN);other['initial_model']['converged_when']=other['initial_model']['converged_when'].replace('1e-10','1e-9')
    with pytest.raises(vbt.Refused,match='converged_when'):vbt.build_config(other,'C0','lib','mesa')
    expect=dict(C1=('numerics','zones',200),C2=('numerics','zones',800),C3=('stepcontrol','delta_target',.01),
                C4=('stepcontrol','delta_target',.0025),C5=('quasistatic','tolerance',1e-10),C6=('cavity','boundary_state','cell'),
                C7=('quasistatic','merge',.67),S1=('opacity','tables','default'),
                S2=('transport','convection','mlt_thesis_lossless'),S4=('opacity','tables','stars2004_switch_strict'),
                V1=('cavity','profile','r-1/2'),V2=('cavity','profile','none'),V3=('replication','total_mass_msun',1e5))
    for run,(section,key,value) in expect.items():
        assert vbt.build_config(PLAN,run,'lib','mesa')[section][key]==value,run
    # Row-0 models: C3 and C4 change only the step controller and share C0's; every other run needs its own.
    sig=lambda r:vbt.signature(vbt.build_config(PLAN,r,'lib','mesa'))
    assert sig('C3')==sig('C0')==sig('C4') and len({sig(r) for r in ('C0','C1','C2','C5','C6','C7','S1','S2','S4','V1','V2','V3')})==12
    with pytest.raises(ValueError,match='Unknown registered run'):vbt.build_config(PLAN,'C9','lib','mesa')
    other=copy.deepcopy(PLAN);other['convergence_campaign']['runs']['C1']='N = 300'
    with pytest.raises(vbt.Refused,match='plan text of C1'):vbt.build_config(other,'C1','lib','mesa')
    options=vbt.track_options(PLAN,c)
    assert options['policy']=='published_assumption' and options['max_accepted_steps']==20000
    assert options['max_wall_seconds']==24*3600 and options['dt_floor']['steps']==500

# ---------------------------------------------------------------- cells, checks, verdict

TARGETS=vbt.table_rows(ROOT/'benchmarks/ball2011_table1.csv')

def printed_rows(scale=None):
    """History-like rows carrying the printed Table 1 values (scale: {(row, observable): factor})."""
    rows={}
    for i,t in enumerate(TARGETS[:5]):
        row=dict(event='start' if i==0 else 'landing',screen_warnings='')
        for name,(column,factor,history) in vbt.OBSERVABLES.items():
            row[history]=repr(float(t[column])*factor*(scale or {}).get((i,name),1.))
        rows[0 if i==0 else float(t['bh_mass_msun'])]=row
    return rows

def test_rounding_of_printed_cells():
    assert vbt.rounding('8.71e-5')==pytest.approx(.005/8.71,rel=1e-12)
    assert vbt.rounding('0.0166')==pytest.approx(.00005/.0166,rel=1e-12)
    assert vbt.rounding('1359')==pytest.approx(.5/1359,rel=1e-12) and vbt.rounding('3.83')==pytest.approx(.005/3.83,rel=1e-12)
    assert vbt.rounding('0.00')==math.inf

def test_the_printed_table_passes_all_44_cells_and_a_deviation_fails_one():
    out=vbt.cells(PLAN,TARGETS,printed_rows())
    assert len(out)==44==PLAN['comparison']['cells']['count'] and all(c['passed'] for c in out)
    assert all(c['deviation']==0 for c in out if c['kind']=='ratio')
    mcav=[c for c in out if c['kind']=='upper_limit'];assert len(mcav)==1 and mcav[0]['qstars']==0.
    L=next(c for c in out if c['row']==2 and c['observable']=='L*')
    assert L['effective_tolerance']==pytest.approx(.05+.0034+.005/3.40,rel=1e-12)
    assert L['strict_tolerance']==pytest.approx(.01+.0034+.005/3.40,rel=1e-12)
    out=vbt.cells(PLAN,TARGETS,printed_rows({(2,'L*'):math.exp(.06),(3,'T0'):math.exp(.03)}))
    failed=[(c['row'],c['observable']) for c in out if not c['passed']]
    assert failed==[(2,'L*')]
    T0=next(c for c in out if c['row']==3 and c['observable']=='T0');assert T0['passed'] and not T0['strict_passed']
    # A row missing from the run leaves its cells undecided, never passed.
    rows=printed_rows();del rows[1000.]
    out=vbt.cells(PLAN,TARGETS,rows);assert sum(c.get('missing',False) for c in out)==9 and sum(c['passed'] is None for c in out)==9
    # Report-only screen flags make a cell not production-valid.
    rows=printed_rows();rows[100.]['screen_convective_mach_violated']='True'
    cell=next(c for c in vbt.cells(PLAN,TARGETS,rows) if c['row']==1)
    assert cell['qstars_production_valid'] is False and cell['report_only_screens_fired']==['convective_mach']
    assert cell['passed'] and cell['enforced_screens_fired']==[]  # a report-only screen never fails a cell

def test_cells_of_a_row_outside_the_enforced_domain_never_pass():
    """Final review item 4: cells() looked only at the report-only flags, so a row 0 that violates an
    always-enforced screen (plan replication_mode.always_enforced; quasistatic.run ends such a run at its
    start row) could pass all its cells. Its cells are not evaluable (passed None) and name the screen:
    the plan's FAIL rule needs a cell outside tolerance, which such a row does not show (first review of
    the unblock: passed False made the verdict FAIL where it is INCONCLUSIVE), and it can never be PASS."""
    assert set(vbt._enforceable(PLAN))&set(vbt._report_only(PLAN))==set()
    assert {'radiation_knudsen_deep','free_streaming_deep','photosphere_gas_pressure','eos_temperature_domain','flow_mach',
            'compactness','mesa_first_law','mesa_ideal_backstop','mesa_dilute_gas','omitted_conduction','nonfinite_state'}<=set(vbt._enforceable(PLAN))
    assert 'outer_face_knudsen' not in vbt._enforceable(PLAN)  # replaced by the Eddington-relation check
    rows=printed_rows();rows[0]['screen_radiation_knudsen_deep_violated']='True';rows[0]['screen_radiation_knudsen_outer_violated']='True'
    out=vbt.cells(PLAN,TARGETS,rows)
    row0=[c for c in out if c['row']==0];rest=[c for c in out if c['row']>0]
    assert len(row0)==8 and all(c['passed'] is None and c['enforced_screens_fired']==['radiation_knudsen_deep'] for c in row0)
    assert all('radiation_knudsen_deep' in c['not_evaluable'] for c in row0)
    assert all(c['qstars_production_valid'] is False and c['report_only_screens_fired']==['radiation_knudsen_outer'] for c in row0)
    assert all(c['strict_passed'] is None for c in row0 if c['kind']=='ratio') and all(c['deviation']==0 for c in row0 if c['kind']=='ratio')
    assert all(c['passed'] and c['enforced_screens_fired']==[] and 'not_evaluable' not in c for c in rest)
    assert vbt.verdict(out,dict(outcome='PASS'),{},True,True,True)=='INCONCLUSIVE'  # before: FAIL
    assert vbt.verdict(out,dict(outcome='PASS'),{},False,True,True)=='INCONCLUSIVE'
    # A row outside tolerance elsewhere still fails the verdict when convergence is met.
    rows[100.]['radius_rsun']=repr(float(rows[100.]['radius_rsun'])*1.5)
    assert vbt.verdict(vbt.cells(PLAN,TARGETS,rows),dict(outcome='PASS'),{},True,True,True)=='FAIL'
    # endpoint() of a run that raised (track.json without a final state) is undecided, not an error.
    end=vbt.endpoint(PLAN,dict(status='numerical_failure',label='numerical_failure',reason='driver error: boom',
                               final_bh_mass_msun=None,final_time_years=None))
    assert end['outcome']=='INCONCLUSIVE' and 'driver error' in end['reason']

def test_fallback_f1_is_a_recorded_switch_outside_the_row0_signature(tmp_path):
    """Plan stop_rules.fallbacks_in_order[0]: --fallback F1 sets quasistatic.ptc for smoke, track and campaign.
    The relaxation never uses it, so it does not enter the row-0 signature; relax and compare refuse it."""
    assert vbt.FALLBACK_TEXTS['F1']==PLAN['stop_rules']['fallbacks_in_order'][0]
    plain=vbt.build_config(PLAN,'C1','lib','mesa')
    assert plain['quasistatic']['ptc'] is False and vbt.fallbacks_of(plain)==[]
    f1=vbt.apply_fallbacks(vbt.build_config(PLAN,'C1','lib','mesa'),['F1'])
    assert f1['quasistatic']['ptc'] is True and vbt.fallbacks_of(f1)==['F1'] and vbt.signature(f1)==vbt.signature(plain)
    assert {k:v for k,v in vbt._flat(f1).items() if v!=vbt._flat(plain)[k]}=={'quasistatic.ptc':True}
    with pytest.raises(vbt.Refused,match='unknown fallback F2'):vbt.apply_fallbacks(plain,['F2'])
    for command in (['relax'],['compare']):
        with pytest.raises(SystemExit) as exit:vbt.main(['--fallback','F1','--output',str(tmp_path/'out')]+command)
        assert exit.value.code==2 and not (tmp_path/'out').exists()

def test_verdict_rules():
    ok=[dict(passed=True)]*44;bad=ok[:-1]+[dict(passed=False)]
    end=lambda o:dict(outcome=o)
    q=lambda p:dict(Q1=dict(required=True,passed=p),Q6=dict(required=False))
    assert vbt.verdict(ok,end('PASS'),q(True),True,True,True)=='PASS'
    assert vbt.verdict(bad,end('PASS'),q(True),True,True,True)=='FAIL'
    assert vbt.verdict(bad,end('PASS'),q(True),False,True,True)=='INCONCLUSIVE'  # a failed cell without convergence
    assert vbt.verdict(ok,end('FAIL'),q(True),False,True,True)=='FAIL'
    assert vbt.verdict(ok,end('PASS'),q(False),True,True,True)=='FAIL'
    assert vbt.verdict(ok,end('INCONCLUSIVE'),q(True),True,True,True)=='INCONCLUSIVE'
    for flags in ((True,False,True),(True,True,False)):assert vbt.verdict(ok,end('PASS'),q(True),*flags)=='INCONCLUSIVE'
    assert vbt.verdict([dict(passed=None)]+ok[1:],end('PASS'),q(True),True,True,True)=='INCONCLUSIVE'
    # A required check that could not be evaluated is undecided, not failed.
    assert vbt.verdict(ok,end('PASS'),q(None),True,True,True)=='INCONCLUSIVE'
    assert vbt.verdict(ok,end('INCONCLUSIVE'),q(None),True,True,True)=='INCONCLUSIVE'
    # Opacity clamp above 0.01 Msun at a compared row invalidates the comparison; unreported blocks PASS.
    assert vbt.verdict(bad,end('PASS'),q(True),True,True,True,clamp_ok=False)=='INCONCLUSIVE'
    assert vbt.verdict(ok,end('PASS'),q(True),True,True,True,clamp_ok=None)=='INCONCLUSIVE'
    assert vbt.verdict(ok,end('PASS'),q(True),True,True,True,clamp_ok=True)=='PASS'

def stopped_c0(stop_msun=700.,failing=None):
    """A C0 that ended by numerical_failure at stop_msun: rows 0-3 and their profiles pass every check."""
    rows=printed_rows();rows={k:v for k,v in rows.items() if k<stop_msun}
    for k,row in rows.items():
        row.update(max_density_inversion='0.0' if k==0 else '0.01',stability_ratio='0.05',luminosity_over_law='1.0',
                   opacity_clamped_mass_msun='0.0')
    if failing:rows[failing[0]][failing[1]]=failing[2]
    history=[dict(event='start' if i==0 else 'step',mdot_bh_msun_per_year=repr(2.5e-4),stability_ratio='0.05',time_years=repr(1e5*i),
                  radius_rsun='1e4',bh_mass_msun=repr(5+i),effective_temperature_K='4000',screen_warnings='') for i in range(20)]
    profiles={m:dict(mass_above_outer_convective_boundary_msun=10.,convective_gamma1=1.34) for m in (100.,200.,500.,1000.) if m<stop_msun}
    track=dict(final_bh_mass_msun=stop_msun,final_time_years=3e6,label='numerical_failure')
    return rows,history,profiles,track

def test_qualitative_checks_are_three_valued_and_a_failed_run_is_inconclusive():
    """A C0 that stops at 700 Msun with rows 0-3 passing: missing epochs leave the required checks undecided."""
    rows,history,profiles,track=stopped_c0()
    q=vbt.qualitative(PLAN,rows,history,profiles)
    assert q['Q1']['passed'] is True and all(q[k]['passed'] is None for k in ('Q2','Q3','Q4','Q5'))
    assert q['Q2']['missing']==[1000.] and q['Q4']['missing']==[1000.] and q['Q5']['missing']==[1000.]
    end=vbt.endpoint(PLAN,track);assert end['outcome']=='INCONCLUSIVE'
    out=vbt.cells(PLAN,TARGETS,rows);assert all(c['passed'] for c in out if c['row']<4)
    assert vbt.verdict(out,end,q,False,True,True)=='INCONCLUSIVE'
    # An evaluated epoch that fails is a definite failure even when others are missing.
    rows,history,profiles,track=stopped_c0(failing=(200.,'max_density_inversion','0.0'))
    q=vbt.qualitative(PLAN,rows,history,profiles)
    assert q['Q2']['passed'] is False and vbt.verdict(vbt.cells(PLAN,TARGETS,rows),vbt.endpoint(PLAN,track),q,False,True,True)=='FAIL'
    profiles[200.]['convective_gamma1']=None  # no convective cells at row 2
    q=vbt.qualitative(PLAN,rows,history,profiles);assert q['Q5']['passed'] is False and q['Q5']['no_convective_cells']==[200.]
    assert vbt.qualitative(PLAN,{},[],{})['Q1']['passed'] is None

def test_opacity_clamp_limit_at_compared_rows():
    rows,history,_,_=stopped_c0(1e4)
    out=vbt.clamp_check(TARGETS,rows,history)
    assert out['valid'] is True and not out['S4_required'] and set(out['values'])=={0,100.,200.,500.,1000.}
    rows[500.]['opacity_clamped_mass_msun']='0.0101';history[3]['opacity_clamped_mass_msun']='1e-5'
    out=vbt.clamp_check(TARGETS,rows,history)
    assert out['valid'] is False and out['offending_rows']==[500.] and out['S4_required'] and 'S4' in out['follow_up']
    cell=next(c for c in vbt.cells(PLAN,TARGETS,rows) if c['row']==3);assert cell['opacity_clamped_mass_msun']==.0101
    del rows[1000.]
    out=vbt.clamp_check(TARGETS,rows);assert out['valid'] is False and out['unreported_rows']==[1000.]
    rows[500.]['opacity_clamped_mass_msun']='0.0'
    assert vbt.clamp_check(TARGETS,rows)['valid'] is None

def conservation_row(**change):
    row=dict(event='step',step_mass_residual_change='1e-16',mass_residual_relative='1e-15',energy_residual_relative='1e-12',
             engine_conversion_residual_erg='0.0',ledger_central_generated='1e45',screen_volume_geometry='1e-14',
             step_max_event_defect='1e30',ledger_initial_energy_scale='1e50',step_inflow_heating_over_bh_energy='1e-5')
    row.update(change);return row

def test_conservation_includes_the_per_event_reconstruction_defect():
    """Plan conservation.reconstruction_defect_per_event_relative = 1e-12 of the initial internal energy, every accepted step."""
    ok=vbt.conservation(PLAN,[dict(event='start'),conservation_row(),conservation_row(event='landing')])
    assert ok['passed'] and ok['accepted_steps']==2 and ok['worst']['reconstruction_event']==pytest.approx(1e-20)
    bad=vbt.conservation(PLAN,[conservation_row(),conservation_row(step_max_event_defect='-2e38')])
    assert not bad['passed'] and bad['violations']['reconstruction_event']==1 and bad['worst']['reconstruction_event']==pytest.approx(2e-12)
    missing=vbt.conservation(PLAN,[conservation_row(step_max_event_defect='')])
    assert not missing['passed'] and missing['missing']['reconstruction_event']==1
    flagged=vbt.conservation(PLAN,[conservation_row(step_inflow_heating_over_bh_energy='0.01')])
    assert flagged['passed'] and flagged['violations']['inflow_flagged']==1  # reported only

def test_relaxation_blocker_is_numerical_unless_gamma_ph_is_reached_on_a_converging_path():
    gamma='Photospheric gas pressure is not positive (Gamma_ph >= 1); no Ball surface exists'
    steps=lambda f:[dict(photosphere_eddington_factor=x) for x in f]
    record=lambda f,messages:dict(steps=steps(f),rejections=[dict(domain_message=m,domain_trials=3) for m in messages],
                                 final_rejection_streak=len(messages))
    rising=[.95,.97,.98,.985,.99,.995]
    exit_=vbt.relaxation_domain_exit(record(rising,[gamma]*21))
    assert exit_['domain_exit'] and exit_['final_rejections']==21 and exit_['domain_trials']==[3]*21
    # Any other domain message in the final streak, a factor that falls or stays below 0.99, or no streak: numerical.
    assert not vbt.relaxation_domain_exit(record(rising,[gamma]*20+['Supplied material lies outside the EOS temperature domain']))['domain_exit']
    assert not vbt.relaxation_domain_exit(record([.95,.97,.999,.98,.99,.995],[gamma]*21))['domain_exit']
    assert not vbt.relaxation_domain_exit(record([.5,.6,.7,.8,.9,.95],[gamma]*21))['domain_exit']
    assert not vbt.relaxation_domain_exit(record(rising,[]))['domain_exit']
    assert not vbt.relaxation_domain_exit(record(rising[:3],[gamma]*21))['domain_exit']
    assert not vbt.relaxation_domain_exit(record(rising,[None]*21))['domain_exit']

def test_endpoint_outcomes():
    """Plan endpoint: the fold must be confirmed at N = 400 (C0) and N = 800 (C2) and agree to 2%."""
    G=1.98847e33
    track=lambda m,t,label='model_domain_exit',confirmed=True,fold=None:dict(final_bh_mass_msun=m,final_time_years=t*1e6,label=label,
        fold=dict(fold_confirmed=confirmed,bh_mass=(fold or m)*G) if confirmed is not None else {})
    c2=track(1196,4.24)
    assert vbt.endpoint(PLAN,track(1194,4.23),c2)['outcome']=='PASS'
    out=vbt.endpoint(PLAN,track(1194*math.exp(.11),4.23),track(1194*math.exp(.11),4.23))
    assert out['outcome']=='FAIL' and out['fold_confirmed']
    assert vbt.endpoint(PLAN,track(1100,4.23,'numerical_failure',False),c2)['outcome']=='INCONCLUSIVE'
    assert vbt.endpoint(PLAN,track(1400,5.,'numerical_blocker (stop rule)',None),c2)['outcome']=='FAIL'
    assert vbt.endpoint(PLAN,None)['outcome']=='INCONCLUSIVE'
    # A C0 fold without N = 800 confirmation is a numerical_failure: no C2, a C2 without a fold, or folds 3% apart.
    for other,why in ((None,'no N = 800'),(track(1150,4.1,'numerical_failure',None),'C2 confirms no fold'),
                      (track(1194*math.exp(.03),4.4),'differ by')):
        out=vbt.endpoint(PLAN,track(1194,4.23),other)
        assert out['outcome']=='INCONCLUSIVE' and out['label']=='numerical_failure' and out['fold_confirmed_n400'] and why in out['reason']
    # ... and past 1.1 x 1194 Msun such a track still fails (no confirmed fold).
    assert vbt.endpoint(PLAN,track(1400,5.),None)['outcome']=='FAIL'
    out=vbt.endpoint(PLAN,track(1194,4.23),track(1194*math.exp(.019),4.3))
    assert out['outcome']=='PASS' and out['dln_fold_n400_n800']==pytest.approx(.019)

def test_campaign_pairs_and_richardson():
    base=printed_rows()
    def scaled(c):return {k:{**r,**{h:repr(float(r[h])*(1+c)) for _,(_,_,h) in vbt.OBSERVABLES.items()}} for k,r in base.items()}
    near=vbt.pair(PLAN,base,scaled(.005),{},{},.01);far=vbt.pair(PLAN,base,scaled(.02),{},{},.01)
    assert near['passed'] and not far['passed'] and len(near['observables'])==7+4*9
    # Second-order errors 4e-3, 1e-3, 2.5e-4 (h = 2, 1, 0.5): order 2, C0 error 1e-3 (limit min(0.01, tau/3)).
    out=vbt.richardson(PLAN,scaled(4e-3),scaled(1e-3),scaled(2.5e-4))
    first=out['1000.0:T0']
    assert first['order']==pytest.approx(2.,abs=5e-3) and first['error']==pytest.approx(1e-3,rel=5e-3) and first['passed']
    assert first['estimator']=='richardson' and first['expected_order_met']
    # Non-monotone differences: the spread |e1| + |e2| bounds the error (order and monotonicity are reported only).
    for c1,c0,c2,passed in ((1e-3,4e-3,2.5e-4,True),(4e-3,1.6e-2,1e-3,False)):
        out=vbt.richardson(PLAN,scaled(c1),scaled(c0),scaled(c2))
        spread=abs(math.log((1+c1)/(1+c0)))+abs(math.log((1+c0)/(1+c2)))
        assert not any(c['monotone'] for c in out.values()) and all(c['estimator']=='spread' for c in out.values())
        assert all(c['error']==pytest.approx(spread,rel=1e-9) and c['limit']==.01 and c['passed'] is passed for c in out.values())
    # Noise far below the limit (+2e-9, then -1e-9) passes every cell.
    out=vbt.richardson(PLAN,scaled(1e-9),scaled(-1e-9),scaled(0.))
    assert len(out)==43 and all(c['passed'] and c['estimator']=='spread' for c in out.values())

# ---------------------------------------------------------------- restart and compare on a toy run

def test_restart_check_and_compare_on_an_analytic_toy_run(tmp_path):
    sys.path.insert(0,str(ROOT/'tests'))
    from test_quasistatic import toy,EOS,OPACITY
    from test_ball_restart import POLICY
    from qstars import quasistatic
    from qstars.constants import MSUN
    s,c=toy(128,'face',eta=4e-4)
    s,_=quasistatic.relax(s,EOS,OPACITY,c,dict(mesh_cycles=0));s=quasistatic.zero_time(s,EOS)
    c['stepcontrol']=dict(metric='mapped',delta_target=.03,dt_initial_years=100.)
    run=lambda state,path,**k:quasistatic.run(state,EOS,OPACITY,c,tmp_path/path,policy=POLICY,landings=[1.05*MSUN],**k)
    run(s,'c0',stop_after_steps=22,restart_masses=[1.03*MSUN],segment_steps=5)
    state,_=quasistatic.load_checkpoint(tmp_path/'c0'/'restart_bh1.03msun.npz')
    run(state,'r1',stop_after_steps=5)
    check=vbt.restart_check(tmp_path/'c0',tmp_path/'r1','R1')
    assert check['passed'] and check['rows']==5 and check['crossed']['base_slice'] and check['crossed']['landing']
    # rejected_step counts only the rejections the segment itself replayed, and they equal C0's.
    replayed=[json.loads(x)['step'] for x in open(tmp_path/'r1'/'events.jsonl') if json.loads(x)['event']=='rejected']
    assert check['crossed']['rejected_step']==bool(replayed) and check['rejections_replayed']==len(replayed)
    assert check['rejections_equal_c0']
    # must_cross over the segments together, with the landing inside R2's segment.
    crossing=vbt.must_cross({'R1':check})
    assert crossing['unmet'][-1]=='landing (restart point 2)' and not crossing['passed']
    assert vbt.must_cross({'R2':check})['unmet']==[k for k in vbt.MUST_CROSS if not check['crossed'][k]]
    rows=(tmp_path/'r1'/'history.csv').read_text().splitlines();rows[-1]=rows[-1].replace(',step,',',landing,',1)
    (tmp_path/'r1'/'history.csv').write_text('\n'.join(rows)+'\n')
    assert not vbt.restart_check(tmp_path/'c0',tmp_path/'r1','R1')['passed']
    # compare() on a partial campaign: every product is written and the verdict cannot be PASS.
    (tmp_path/'c0'/'track.json').write_text(json.dumps(dict(final_bh_mass_msun=1.05,final_time_years=1e5,label='control_stop')))
    comparison,campaign=vbt.compare(PLAN,{'C0':str(tmp_path/'c0')},{'R1':str(tmp_path/'r1')},tmp_path,
                                    dict(plan_sha256=vbt.verify_registration()['plan_sha256']))
    assert (tmp_path/'comparison.json').exists() and (tmp_path/'campaign.json').exists()
    assert comparison['cell_count']==44 and campaign['verdict'] in ('FAIL','INCONCLUSIVE') and not campaign['restarts_bitwise']
    assert campaign['conservation']['C0']['passed'] and campaign['conservation']['C0']['accepted_steps']==22
    assert campaign['conservation']['C0']['missing']['reconstruction_event']==0
    assert 0<campaign['conservation']['C0']['worst']['reconstruction_event']<=1e-12  # split/merge events occurred
    assert 'R1' in campaign['restarts'] and 'restart_must_cross' in campaign
    # These library-made directories carry no driver provenance, so no common code is established: the
    # comparisons between runs do not hold and the verdict is INCONCLUSIVE whatever C0 shows on its own.
    sc=campaign['source_consistency'];assert set(sc['runs'])=={'C0','restart R1'} and not sc['one_package'] and not sc['all_snapshots_verified']
    assert not sc['one_code'] and not campaign['one_code'] and not campaign['c0_one_code'] and campaign['verdict']=='INCONCLUSIVE'
    assert not campaign['convergence_met'] and set(campaign['rejection_causes'])=={'C0','restart R1'}
    assert comparison['opacity_clamp']['valid'] is None  # the analytic opacity has no table-edge clamp
    assert comparison['endpoint']['outcome']=='INCONCLUSIVE' and 'radiation_knudsen_outer' in comparison['report_only_screens_fired']
    assert comparison['qualitative_checks']['Q1']['passed'] is False  # the toy accretes far below Ball's rate

def test_relax_and_track_paths_on_an_analytic_toy(tmp_path,monkeypatch):
    """The MESA-only subcommand paths with the Ball physics replaced by the analytic cavity toy: seed
    relaxation and uniqueness, a track that ends on a stop rule, a restart segment, and a track whose
    termination trigger fires and is classified by the frozen-entropy continuation."""
    sys.path.insert(0,str(ROOT/'tests'))
    from test_quasistatic import toy,EOS,OPACITY
    from test_ball_restart import POLICY
    from qstars import quasistatic
    from qstars.constants import MSUN
    _,base=toy(48,'cell')
    cycles=dict(mesh_cycles=40)  # as an amendment would set them; registered: none
    def build(plan,run_id,library='',mesa_dir='',amendments=()):
        c=copy.deepcopy(base);c['name']=f'toy_{run_id}'
        c['relaxation']=dict(vbt.RELAXATION,**cycles)
        c['stepcontrol']=dict(dt_initial_years=1e-9 if run_id=='C3' else 100.)
        c['replication']=dict(run=run_id,kind='convergence_campaign',plan_text='toy',policy='toy',opacity='toy',total_mass_msun=100.,
                              initial_bh_mass_msun=1.,landings_msun=[1.01],amendments=[])
        return c
    monkeypatch.setattr(vbt,'build_config',build);monkeypatch.setattr(vbt,'physics',lambda config:(EOS,OPACITY))
    monkeypatch.setattr(vbt,'track_options',lambda plan,config:dict(landings=[1.01*MSUN],policy=POLICY,max_wall_seconds=600.,
                                                                    max_accepted_steps=12,dt_floor=None))
    registration=vbt.verify_registration()
    config=build(PLAN,'C0')
    path=vbt.run_directory(tmp_path,'relax','C0',config,registration,[])
    out=vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.,100.],48,log=lambda *a:None)
    assert out['status']=='converged_unique' and len(out['seeds'])==2 and max(out['uniqueness']['max_abs_dln'].values())<=1e-6
    assert json.loads((path/'relax.json').read_text())['row0_checkpoint']=='row0_150rsun.npz'
    assert out['relaxation']['mesh_cycles']==40 and all(len(e['mesh_cycles'])>=1 for e in out['seeds'])
    row0=path/'row0_150rsun.npz'
    # relax.json vouches for its row-0 model: sha256, cell count and the screens of row 0 under the plan policy.
    assert out['row0_sha256']==vbt.sha256(row0)==out['seeds'][0]['row0_sha256'] and out['seeds'][0]['row0_zones']==48
    assert out['seeds'][0]['row0_screens']['policy']=='published_assumption'
    assert out['row0_enforced_screen_errors']==out['seeds'][0]['row0_screens']['enforced_screen_errors']
    # Final review item 3: a row-0 checkpoint is used only when its relax run ended converged_unique, names it
    # and recorded its sha256 (before: any t = 0 checkpoint with a matching configuration was accepted).
    def refuse(checkpoint,match):
        with pytest.raises(vbt.Refused,match=match):vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=checkpoint)
    assert vbt.row0_record(row0)[1] is None and vbt.row0_for(config,[row0]) is not None
    refuse(path/'row0_300rsun.npz',"names 'row0_150rsun.npz' as the row-0 checkpoint")  # the other seed: unique, not designated
    stray=tmp_path/'stray';stray.mkdir();shutil.copy(row0,stray/row0.name)
    refuse(stray/row0.name,'no relax.json next to it')
    record=json.loads((path/'relax.json').read_text())
    for change,match in ((dict(status='not_unique'),"ended 'not_unique', not 'converged_unique'"),
                         (dict(status='relaxation_blocker'),"ended 'relaxation_blocker'"),
                         (dict(row0_sha256='0'*64),'sha256 differs'),(dict(row0_sha256=None),'sha256 differs')):
        (stray/'relax.json').write_text(json.dumps(dict(record,**change)))
        refuse(stray/row0.name,match)
        assert vbt.row0_for(config,[stray/row0.name]) is None and vbt.row0_for(config,[stray/row0.name,row0]) is not None
    (stray/'relax.json').write_text(json.dumps(record))
    # First review of the unblock: a relax.json alone (hand-written, or a dry run's) used to vouch for a row 0.
    # The directory must be a relax run directory made under the registered plan with the same registered
    # amendments, and, for registered output, a registered run directory under benchmarks/results.
    reasons=[];assert vbt.row0_for(config,[stray/row0.name],reasons) is None and 'no run.json of kind relax' in reasons[0]
    for name in ('run.json','plan_registration.json','amendments.json'):shutil.copy(path/name,stray/name)
    scratch=vbt.row0_provenance(tmp_path,registration,[])
    assert scratch['results_root'] is None and vbt.row0_for(config,[stray/row0.name],None,scratch) is not None
    registered=vbt.row0_provenance(vbt.ROOT/'benchmarks/results/2026-10-14-ball-track',registration,[])
    assert registered['results_root']==(vbt.ROOT/'benchmarks/results').resolve()
    reasons=[];assert vbt.row0_for(config,[row0],reasons,registered) is None and 'not a registered run directory' in reasons[0]
    other=dict(registration,registration_sha256='0'*64)
    reasons=[];assert vbt.row0_for(config,[stray/row0.name],reasons,vbt.row0_provenance(tmp_path,other,[])) is None
    assert 'not made under the registered plan' in reasons[0]
    applied=[(dict(path='a.json',sha256='1'*64),{})]
    reasons=[];assert vbt.row0_for(config,[stray/row0.name],reasons,vbt.row0_provenance(tmp_path,registration,applied)) is None
    assert 'other registered amendments' in reasons[0]
    (stray/'run.json').write_text(json.dumps(dict(json.loads((path/'run.json').read_text()),kind='smoke')))
    reasons=[];assert vbt.row0_for(config,[stray/row0.name],reasons,scratch) is None and 'no run.json of kind relax' in reasons[0]
    reasons=[];assert vbt.row0_for(config,[path/'relaxed_150rsun.npz'],reasons) is None and 'not a t = 0 model' in reasons[0]
    track_dir,summary=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=row0)
    assert summary['status']=='stop_rule' and summary['label']=='numerical_blocker (stop rule)' and summary['final_step']==12
    assert summary['landed_msun']==[1.01] and (track_dir/'track.json').exists() and (track_dir/'plan_registration.json').exists()
    assert json.loads((track_dir/'run.json').read_text())['row0']==str(row0)
    # The relax run that vouched for the row-0 model is recorded with the sha256 of its relax.json.
    assert json.loads((track_dir/'run.json').read_text())['row0_relax']==dict(directory=str(path),
        relax_json_sha256=vbt.sha256(path/'relax.json'),checkpoint_sha256=vbt.sha256(row0))
    assert summary['rejection_causes']['total']==summary['rejections']
    # compare: the code of the track and of the relax run of its row-0 model. This pytest process has loaded
    # qstars modules outside the Ball closure, and relax_seeds called directly finishes no provenance, so
    # no code is established here (a driver process records one: test_source_consistency_*).
    closure=json.loads((track_dir/'code_closure.json').read_text())
    assert closure['modules']==list(vbt.BALL_MODULES) and closure['loaded_outside_closure'] and closure['missing']==[]
    sc=vbt.source_consistency({'C0':str(track_dir)})
    assert set(sc['runs'])=={'C0','C0 row0 relax'} and sc['one_package'] and sc['one_driver'] and not sc['one_code'] and not sc['c0_one_code']
    # (this test called relax_seeds directly; the relax subcommand finishes the provenance of its directory)
    assert sc['runs']['C0']['source_snapshot_verified'] is True and sc['runs']['C0 row0 relax']['source_snapshot_verified'] is None
    assert not sc['all_snapshots_verified']
    assert sc['runs']['C0 row0 relax']['path']==str(path) and sc['runs']['C0']['driver_sha256']==vbt.sha256(ROOT/'scripts/validate_ball_track.py')
    # A run whose configuration differs from the row-0 model's is refused (C6-like change).
    monkeypatch.setattr(vbt,'build_config',lambda *a,**k:{**build(*a,**k),'cavity':dict(boundary_state='face')})
    with pytest.raises(vbt.Refused,match='relax --run C6'):vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C6',row0=row0)
    monkeypatch.setattr(vbt,'build_config',build)
    # Restart segment from the run's last checkpoint (here: stop at the step limit, segment of 3 steps).
    segment_dir,seg=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,restart=track_dir/'checkpoint.npz',segment=3)
    assert seg['kind']=='restart' and seg['status']=='stop_rule'  # the run-wide accepted-step limit is already reached
    # The termination trigger (dt below t_dyn at once) is classified: no fold within 5% of this toy's M_BH.
    trig_dir,trig=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C3',row0=row0)
    assert trig['status']=='termination_trigger' and trig['label']=='numerical_failure' and 'no fold' in trig['label_reason']
    fold_out=json.loads((trig_dir/'fold.json').read_text())
    assert fold_out['runs'][0]['stop']=='lambda limit' and not fold_out['fold_confirmed']
    # Retries exhausted by enforced screens (production enforces the Saha screens this cool toy fails) are a
    # termination trigger too: the fold protocol labels the end, the screen exit is recorded next to it.
    monkeypatch.setattr(vbt,'track_options',lambda plan,config:dict(landings=[1.01*MSUN],policy='production',max_wall_seconds=600.,
                                                                    max_accepted_steps=12,dt_floor=None))
    # Final review item 4: the start row itself fails those screens, so the run ends there as a
    # model_domain_exit of row 0. No step is attempted and the fold protocol does not apply (no trigger fired).
    start_dir,start=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=row0)
    assert start['status']=='model_domain_exit' and start['label']==vbt.ROW0_EXIT_LABEL and not (start_dir/'fold.json').exists()
    assert start['attempts']==0 and any('Saha' in e for e in start['start_row_screen_errors']) and 'row 0' in start['label_reason']
    assert json.loads((start_dir/'track.json').read_text())['start_row_screen_errors']==start['start_row_screen_errors']
    assert vbt.endpoint(PLAN,start)['outcome']=='INCONCLUSIVE'
    # A restart segment whose restart row violates them: the exit of the restart row, not of row 0.
    seg_dir,seg_out=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,restart=track_dir/'checkpoint.npz',segment=3)
    assert seg_out['status']=='model_domain_exit' and seg_out['label']==vbt.RESTART_ROW_EXIT_LABEL and 'row 0' not in seg_out['label_reason']
    # First review of the unblock: a nonfinite start state is a numerical defect, never a domain exit.
    from qstars import validity
    check=validity.check
    def nonfinite_start(state,*args,**kwargs):
        errors,diag,conservation=check(state,*args,**kwargs)
        if state.steps==0:
            diag['screens']['nonfinite_state']=dict(diag['screens']['nonfinite_state'],violated=True,value=1.)
            errors=['Nonfinite model state']+errors
        return errors,diag,conservation
    monkeypatch.setattr(validity,'check',nonfinite_start)
    nan_dir,nan=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=row0)
    assert nan['status']=='numerical_failure' and nan['label']=='numerical_failure' and nan['attempts']==0
    assert nan['start_row_screen_errors'][0]=='Nonfinite model state' and 'nonfinite model state at the start row' in nan['label_reason']
    assert not (nan_dir/'fold.json').exists() and vbt.endpoint(PLAN,nan)['outcome']=='INCONCLUSIVE'
    monkeypatch.setattr(validity,'check',check)
    # With a start row inside the domain (here: the screens are applied to candidate steps only) the retries
    # are exhausted by enforced screens: the fold protocol labels the end and the screen exit is recorded.
    def candidates_only(state,*args,**kwargs):
        errors,diag,conservation=check(state,*args,**kwargs)
        return ([] if state.steps==0 else errors),diag,conservation
    monkeypatch.setattr(validity,'check',candidates_only)
    screen_dir,screen=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=row0)
    assert screen['status']=='model_domain_exit' and screen['label']=='numerical_failure' and (screen_dir/'fold.json').exists()
    assert 'all 13 consecutive retries' in screen['enforced_screen_exit']['reason'] and 'start_row_screen_errors' not in screen
    monkeypatch.setattr(validity,'check',check)
    monkeypatch.setattr(vbt,'track_options',lambda plan,config:dict(landings=[1.01*MSUN],policy=POLICY,max_wall_seconds=600.,
                                                                    max_accepted_steps=12,dt_floor=None))
    # Final review item 5: a failing run keeps its records. An error in the fold classification is a
    # numerical_failure of the protocol (never a fold); track.json and the provenance are written.
    from qstars import fold
    classify=fold.classify
    def broken_fold(*args,**kwargs):raise ValueError('Jacobian is singular at the reference state')
    monkeypatch.setattr(fold,'classify',broken_fold)
    bad_dir,bad=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C3',row0=row0)
    assert bad['status']=='termination_trigger' and bad['label']=='numerical_failure' and bad['fold']=={'fold_confirmed':False}
    assert bad['label_reason']=='fold classification failed: ValueError: Jacobian is singular at the reference state'
    assert json.loads((bad_dir/'track.json').read_text())['fold_error']==bad['fold_error'] and not (bad_dir/'fold.json').exists()
    assert 'Traceback' in json.loads((bad_dir/'fold_error.json').read_text())['traceback'] and (bad_dir/'termination.json').exists()
    assert 'source_snapshot_verified' in json.loads((bad_dir/'metadata.json').read_text())
    def interrupted_fold(*args,**kwargs):raise KeyboardInterrupt
    monkeypatch.setattr(fold,'classify',interrupted_fold)
    before=set(tmp_path.iterdir())
    with pytest.raises(KeyboardInterrupt):vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C3',row0=row0)
    left,=set(tmp_path.iterdir())-before
    assert json.loads((left/'track.json').read_text())['label']=='control_stop (track summary incomplete)'
    monkeypatch.setattr(fold,'classify',classify)
    # So is an error that escapes quasistatic.run (it records its own failures; this one is the driver's).
    run=quasistatic.run
    def broken_run(*args,**kwargs):raise KeyError('step_dt')
    monkeypatch.setattr(quasistatic,'run',broken_run)
    lost_dir,lost=vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=row0)
    assert lost['status']=='numerical_failure' and lost['label']=='numerical_failure' and "driver error: KeyError: 'step_dt'" in lost['reason']
    assert (lost_dir/'track.json').exists() and 'Traceback' in json.loads((lost_dir/'driver_error.json').read_text())['traceback']
    # The smoke gate writes smoke.json even when the run died before its start row (no history.csv).
    smoke_dir,smoke=vbt.smoke(None,PLAN,registration,[],'lib','mesa',tmp_path,row0)
    assert smoke['passed'] is False and smoke['history_written'] is False and smoke['status']=='numerical_failure'
    assert not (smoke_dir/'history.csv').exists() and json.loads((smoke_dir/'smoke.json').read_text())==json.loads(json.dumps(smoke))
    assert 'source_snapshot_verified' in json.loads((smoke_dir/'metadata.json').read_text())
    assert vbt.preserved_failure(str(smoke_dir.relative_to(tmp_path)),tmp_path)=='not a run directory under benchmarks/results/<date>-ball-track/'
    monkeypatch.setattr(quasistatic,'run',run)
    # A smoke run that works: M_BH = 1.01 Msun (this toy's first landing) reached; with --fallback F1 the
    # switch is recorded and, since every Newton solve converges, the history is bitwise the same.
    ok_dir,ok=vbt.smoke(None,PLAN,registration,[],'lib','mesa',tmp_path,row0)
    assert ok['reached_bh_mass_100'] and ok['history_written'] and ok['fallbacks']==[] and ok['start_row_screen_errors']==[]
    f1_dir,f1=vbt.smoke(None,PLAN,registration,[],'lib','mesa',tmp_path,row0,fallbacks=['F1'])
    assert f1['fallbacks']==['F1'] and json.loads((f1_dir/'run.json').read_text())['fallbacks']==['F1']
    assert (f1_dir/'history.csv').read_bytes()==(ok_dir/'history.csv').read_bytes()
    assert {k:v for k,v in f1.items() if k!='fallbacks'}=={k:v for k,v in ok.items() if k!='fallbacks'}
    import yaml
    assert yaml.safe_load((f1_dir/'config.yaml').read_text())['quasistatic']['ptc'] is True
    with pytest.raises(vbt.Refused,match='relax --run C1'):vbt.smoke(None,PLAN,registration,[],'lib','mesa',tmp_path,path/'row0_300rsun.npz')
    # A restart segment must be run with the fallbacks of its checkpoint.
    with pytest.raises(vbt.Refused,match="made with fallbacks \\['F1'\\], not \\[\\]"):
        vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,restart=f1_dir/'checkpoint.npz',segment=2)
    # Registered relaxation (no mesh cycles): each seed keeps its own mesh and row 0 is not unique.
    cycles.clear()
    plain=vbt.relax_seeds(PLAN,build(PLAN,'C0'),EOS,OPACITY,vbt.run_directory(tmp_path,'relax','C0',build(PLAN,'C0'),registration,[]),
                          [150.,300.],48,log=lambda *a:None)
    assert plain['relaxation']['mesh_cycles']==0 and plain['status']=='not_unique'
    assert all(e['mesh_cycles']==[] for e in plain['seeds']) and max(plain['uniqueness']['max_abs_dln'].values())>1e-6
    # Such a run still writes its row-0 checkpoints, and no track may start from them (item 3).
    plain_dir=max((d for d in tmp_path.iterdir() if d.name.endswith('relax_C0_n48')),key=lambda d:d.name)
    assert json.loads((plain_dir/'relax.json').read_text())['status']=='not_unique' and (plain_dir/'row0_150rsun.npz').exists()
    with pytest.raises(vbt.Refused,match="ended 'not_unique', not 'converged_unique'"):
        vbt.track(None,PLAN,registration,[],'lib','mesa',tmp_path,run_id='C0',row0=plain_dir/'row0_150rsun.npz')
    with pytest.raises(vbt.Refused,match="ended 'not_unique'"):
        vbt.smoke(None,PLAN,registration,[],'lib','mesa',tmp_path,plain_dir/'row0_150rsun.npz')

def test_relaxation_errors_fail_the_seed_and_keep_the_record(tmp_path,monkeypatch,capsys):
    """An error inside a relaxation that is not a step rejection (the P0 predictor on an accepted state, the
    mesh cycles, an EOS inversion outside a step) is that seed's failure: the fallback seeds are still tried,
    relax.json records the blocker, and the preserved directory can back an amendment."""
    sys.path.insert(0,str(ROOT/'tests'))
    from test_quasistatic import toy,EOS,OPACITY
    from qstars import quasistatic
    _,base=toy(48,'cell')
    def build(plan,run_id,library='',mesa_dir='',amendments=()):
        c=copy.deepcopy(base);c['name']=f'toy_{run_id}';c['relaxation']=dict(vbt.RELAXATION)
        c['replication']=dict(run=run_id,kind='convergence_campaign',plan_text='toy',policy='toy',total_mass_msun=100.,
                              initial_bh_mass_msun=1.,seed_radii_rsun=[150.,300.],landings_msun=[1.01],amendments=[])
        return c
    registration=vbt.verify_registration();config=build(PLAN,'C0')
    parent=tmp_path/'benchmarks/results/2026-10-08-ball-track'
    relax=quasistatic.relax;calls=[]
    def first_fails(*args,**kwargs):
        calls.append(1)
        if len(calls)==1:raise ValueError('Equidistributed faces are not increasing')
        return relax(*args,**kwargs)
    monkeypatch.setattr(quasistatic,'relax',first_fails)
    path=vbt.run_directory(parent,'relax','C0',config,registration,[])
    out=vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.,200.],48,log=lambda *a:None)
    assert [e['status'] for e in out['seeds']]==['failed','converged','converged'] and len(calls)==3
    assert out['seeds'][0]['reason']=='ValueError: Equidistributed faces are not increasing'
    assert out['seeds'][0]['domain_exit_evidence']==dict(domain_exit=False,domain_trials=[])
    assert out['status'] in ('converged_unique','not_unique') and out['row0_checkpoint']=='row0_300rsun.npz'
    assert out['uniqueness']['seeds']==[300.,200.] and json.loads((path/'relax.json').read_text())['status']==out['status']
    assert not (path/'relaxed_150rsun.npz').exists() and (path/'row0_300rsun.npz').exists()
    # Every seed errors: a numerical relaxation blocker (never a domain exit), recorded and citable.
    def broken(*args,**kwargs):raise RuntimeError('MESA temperature inversion failed')
    monkeypatch.setattr(quasistatic,'relax',broken)
    path=vbt.run_directory(parent,'relax','C0',config,registration,[])
    out=vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.],48,log=lambda *a:None)
    assert out['status']=='relaxation_blocker' and out['classification']=='numerical relaxation blocker'
    assert [e['reason'] for e in out['seeds']]==['RuntimeError: MESA temperature inversion failed']*2
    assert out['domain_trials']=={'150.0':[],'300.0':[]} and json.loads((path/'relax.json').read_text())['status']=='relaxation_blocker'
    assert vbt.preserved_failure(str(path.relative_to(tmp_path)),tmp_path) is True
    # First review of the unblock: an exception outside those types escaped relax_seeds and relax.json was lost
    # (a KeyError inside quasistatic.relax; a ValueError while saving the profile of a converged seed). Each is
    # now that seed's failure with its traceback, never a converged seed, and the next seed is tried.
    def key_error(*args,**kwargs):raise KeyError('step_dt')
    monkeypatch.setattr(quasistatic,'relax',key_error)
    path=vbt.run_directory(parent,'relax','C0',config,registration,[])
    out=vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.],48,log=lambda *a:None)
    assert out['status']=='relaxation_blocker' and [e['status'] for e in out['seeds']]==['failed','failed']
    assert out['seeds'][0]['reason']=="KeyError: 'step_dt'" and out['seeds'][0]['error_stage']=='relaxation'
    assert 'Traceback' in out['seeds'][0]['traceback'] and out['seeds'][0]['domain_exit_evidence']['domain_exit'] is False
    assert json.loads((path/'relax.json').read_text())['classification']=='numerical relaxation blocker'
    monkeypatch.setattr(quasistatic,'relax',relax)
    save_profile=quasistatic.save_profile;saved=[]
    def first_profile_fails(path_,*args,**kwargs):
        saved.append(1)
        if len(saved)==1:raise ValueError('nonfinite profile')
        return save_profile(path_,*args,**kwargs)
    monkeypatch.setattr(quasistatic,'save_profile',first_profile_fails)
    path=vbt.run_directory(parent,'relax','C0',config,registration,[])
    out=vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.,200.],48,log=lambda *a:None)
    first=out['seeds'][0]
    assert first['status']=='failed' and first['relaxation_status']=='converged' and first['reason']=='ValueError: nonfinite profile'
    assert first['error_stage']=='post-processing of a converged relaxation' and 'row0_checkpoint' not in first and 'row0' not in first
    assert [e['status'] for e in out['seeds'][1:]]==['converged','converged'] and out['uniqueness']['seeds']==[300.,200.]
    assert out['row0_checkpoint']=='row0_300rsun.npz' and json.loads((path/'relax.json').read_text())['status']==out['status']
    monkeypatch.setattr(quasistatic,'save_profile',save_profile)
    # An interrupt inside relax_seeds still leaves relax.json (status control_stop: no result, and no failure
    # an amendment could cite), and is raised again.
    def interrupt(*args,**kwargs):raise KeyboardInterrupt
    monkeypatch.setattr(quasistatic,'relax',interrupt)
    path=vbt.run_directory(parent,'relax','C0',config,registration,[])
    with pytest.raises(KeyboardInterrupt):vbt.relax_seeds(PLAN,config,EOS,OPACITY,path,[150.,300.],48,log=lambda *a:None)
    stopped=json.loads((path/'relax.json').read_text())
    assert stopped['status']=='control_stop' and 'row0_checkpoint' not in stopped
    assert vbt.preserved_failure(str(path.relative_to(tmp_path)),tmp_path)=='relax.json records no failure'
    monkeypatch.setattr(quasistatic,'relax',broken)
    # The relax subcommand: exit status 1, the record written and the provenance finished, also when relax_seeds raises.
    monkeypatch.setattr(vbt,'common',lambda a:(registration,PLAN,[],Path('lib'),'mesa',tmp_path/'cli'))
    monkeypatch.setattr(vbt,'build_config',build);monkeypatch.setattr(vbt,'physics',lambda config:(EOS,OPACITY))
    assert vbt.main(['relax'])==1
    run,=(tmp_path/'cli').iterdir()
    assert json.loads((run/'relax.json').read_text())['status']=='relaxation_blocker'
    assert 'source_snapshot_verified' in json.loads((run/'metadata.json').read_text())
    def interrupted(*args,**kwargs):raise KeyboardInterrupt
    monkeypatch.setattr(vbt,'relax_seeds',interrupted)
    with pytest.raises(KeyboardInterrupt):vbt.main(['relax'])
    other,=[d for d in (tmp_path/'cli').iterdir() if d!=run]
    assert 'source_snapshot_verified' in json.loads((other/'metadata.json').read_text()) and not (other/'relax.json').exists()

def code_run(tmp_path,name,closure='c'*64,driver=b'driver',library='l'*64,verified=True,outside=(),row0=None):
    """A run directory with the records source_consistency reads (as run_directory and finish() write them)."""
    d=tmp_path/name;d.mkdir(parents=True)
    (d/'metadata.json').write_text(json.dumps(dict(process_package_sha256='p'*64,source_snapshot_verified=verified)))
    (d/'run.json').write_text(json.dumps(dict(kind='track',library_sha256=library,row0=row0)))
    (d/'driver.py').write_bytes(driver)
    (d/'code_closure.json').write_text(json.dumps(dict(modules=list(vbt.BALL_MODULES),sha256=closure,missing=[],
                                                       loaded_outside_closure=list(outside))))
    return d

def test_source_consistency_is_enforced_between_compared_runs(tmp_path):
    """First review of the unblock: the code of each run was reported, not enforced, so a campaign made across
    source edits could pass. Comparisons between runs now hold only between runs made with one established
    code: the Ball import closure digest, the driver and the MESA library, a verified snapshot and no qstars
    module loaded outside the closure. C0 must share it with the relax run of its row-0 model."""
    relax=code_run(tmp_path,'relax');row0=str(relax/'row0_6000rsun.npz')
    runs={k:str(code_run(tmp_path,k,row0=row0)) for k in ('C0','C2')}
    sc=vbt.source_consistency(runs)
    assert sc['one_code'] and sc['c0_one_code'] and vbt.same_code(sc,'C0','C2') and set(sc['runs'])=={'C0','C2','C0 row0 relax','C2 row0 relax'}
    for change in (dict(closure='d'*64),dict(driver=b'other driver'),dict(library='m'*64),dict(verified=False),dict(outside=['plotting'])):
        runs['C2']=str(code_run(tmp_path,f'C2_{len(list(tmp_path.iterdir()))}',row0=row0,**change))
        sc=vbt.source_consistency(runs)
        assert not sc['one_code'] and sc['c0_one_code'] and not vbt.same_code(sc,'C0','C2'),change
    runs['C0']=str(code_run(tmp_path,'C0_other',closure='e'*64,row0=row0))
    assert not vbt.source_consistency(runs)['c0_one_code']  # C0 and its row-0 model were made with different code
    assert not vbt.source_consistency({'C0':str(tmp_path/'missing')})['one_code']
    # The closure of this driver: the modules a driver process imports (checked in a fresh process).
    assert 'quasistatic' in vbt.BALL_MODULES and 'fold' in vbt.BALL_MODULES and 'cmfgen_wind' not in vbt.BALL_MODULES
    import subprocess
    script=("import sys,json;sys.path.insert(0,'scripts');sys.path.insert(0,'src');import validate_ball_track as v;"
            "v.import_ball_modules();print(json.dumps(v.code_closure()))")
    out=json.loads(subprocess.run([sys.executable,'-c',script],cwd=ROOT,capture_output=True,text=True,check=True).stdout)
    assert out['loaded_outside_closure']==[] and out['missing']==[] and out['modules']==list(vbt.BALL_MODULES)

def test_rejection_causes_are_counted_from_the_events(tmp_path):
    events=[dict(event='rejected',kind='numerical',message='delta_actual 0.0188 exceeds 2 x delta_target (profile_dlnr)'),
            dict(event='rejected',kind='numerical',message='delta_actual 0.012 exceeds 2 x delta_target (profile_dlnr)'),
            dict(event='rejected',kind='numerical',domain_exit=False,
                 message='Quasi-static solve failed (pass 1, residual 1.91e-12): Newton line search or iteration limit'),
            dict(event='rejected',kind='numerical',domain_exit=False,
                 message='Quasi-static solve failed (pass 2, residual 0.404835): Newton line search or iteration limit'),
            dict(event='rejected',kind='numerical',domain_exit=True,
                 message='Quasi-static solve failed (pass 1, residual inf): Starting point outside the model domain'),
            dict(event='rejected',kind='screen',message='Enforced screens: Newtonian compactness limit exceeded'),
            dict(event='rejected',kind='numerical',message='Per-step conservation budget exceeded: energy 1e-12 (limit 5e-13), mass 0 (limit 1e-14)'),
            dict(event='rejected',kind='numerical',domain_exit=True,message='P3 predictor left the model domain (pass 1): Gamma_ph >= 1'),
            dict(event='rejected',kind='numerical',message='P2 base remesh failed (pass 1): EOS inversion'),
            dict(event='accepted',step=1)]
    (tmp_path/'events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in events))
    out=vbt.rejection_causes(tmp_path)
    assert out['total']==9 and out['controller_metric']=={'profile_dlnr':2} and out['newton_failures_at_tolerance_floor']==1
    assert out['causes']==dict(controller_metric=2,newton_failure=2,newton_failure_domain=1,enforced_screen=1,conservation_budget=1,
                               domain_exit_before_solve=1,error_before_solve=1)
    assert out['newton_failure_residuals'][:2]==[1.91e-12,0.404835] and math.isinf(out['newton_failure_residuals'][2])
    assert vbt.rejection_causes(tmp_path/'none')['total']==0

def test_mapped_plan_texts():
    vbt.check_mapped_texts(PLAN)
    other=copy.deepcopy(PLAN);other['qualitative_checks'][0]['test']='every accepted step lies in [1e-4, 5e-4] Msun/yr'
    with pytest.raises(vbt.Refused,match='Q1'):vbt.check_mapped_texts(other)
