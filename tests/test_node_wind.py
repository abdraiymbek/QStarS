"""Conservative evolution coupling to a synthetic CMFGEN-format wind node.

No CMFGEN run: the node is built from the initial 16-zone ideal-gas state so
its identities hold exactly. It is not a physical wind and stays unvalidated.
"""
import copy,hashlib,json
import numpy as np
import pytest
from qstars.config import resolve,load
from qstars.constants import A,G,PI,SIGMA,MSUN,YEAR
from qstars.evolution import physics,attempt_step,run
from qstars.structure import initialize
from qstars import ledger
from qstars.engine import rates
from qstars.mass_transport import exchange,DomainError
from qstars.numerics import StepRejected
from qstars.io import load_checkpoint
from qstars.winds import evaluate,rate
from qstars.wind_node import (validate,nodes_of,json_schema,interface_check,EVOLUTION_GAUGES,
                              VALIDATION,PROVENANCE,BASE,STATE,OUTPUT)

L0,MDOT,VINF,DELTA=1e40,1e-4*MSUN/YEAR,1e8,.01
BASECONF={'name':'node_wind','numerics':{'zones':16},'engine':{'model':'off'},'transport':{'convection':'off'},
          'boundary':{'thermal':'luminosity','luminosity_erg_s':L0},'evolution':{'end_time_years':.003}}

def make_node(s,e,node_id='synthetic',heat=1.,mdot=MDOT,vinf=VINF,gauge=EVOLUTION_GAUGES['ideal_gas_radiation']):
    M,R,rho,T=float(s.enclosed[-1]),float(s.r[-1]),float(s.rho[-1]),float(s.T[-1])*heat
    q=e.evaluate(rho,T);E=A*T**4;area=4*PI*R**2
    emergent=L0-mdot*(.5*vinf**2+G*M/R-float(q['h']))
    face=dict(radius=R,enclosed_mass=M,density=rho,velocity=mdot/(area*rho),gas_pressure=float(q['pg']),
              gas_internal_energy=float(q['u'])-E/rho,radiation_energy=E,radiation_flux=L0/area,
              radiation_pressure_rr=E/3,potential=-G*M/R,energy_gauge=gauge,gas_temperature=T)
    return dict(schema='qstars_wind_node',schema_version=1,node_id=node_id,M_total_g=M,R_base_cm=R,L_base_erg_s=L0,
                rho_base=rho,T_base=T,P_gas_base=face['gas_pressure'],P_rad_base=E/3,mdot_g_s=mdot,v_inf_cm_s=vinf,
                L_emergent_erg_s=emergent,R_phot_cm=2*R,Teff_emergent_K=(emergent/(4*PI*SIGMA*(2*R)**2))**.25,
                energy_gauge=gauge,base_face=face,validation={k:False for k in VALIDATION},
                provenance=dict(source_run='synthetic test node; no CMFGEN run',**{k:'0'*64 for k in PROVENANCE}))

def write(path,data):
    path.write_text(json.dumps(data));return path

def node_config(path,**raw):
    c=copy.deepcopy(BASECONF)
    for k,v in raw.items():c.setdefault(k,{}).update(v)
    c['wind']=dict(model='selected_nodes',node_path=str(path),drift_tolerance=DELTA)
    return resolve(c)

@pytest.fixture
def model(tmp_path):
    c=resolve(BASECONF);e,op=physics(c);s,_=initialize(c,e);ledger.initialize(s,e)
    node=make_node(s,e);path=write(tmp_path/'wind_node.json',node)
    return s,e,op,node_config(path),node,path

def test_contract_identities_and_schema(model,tmp_path):
    s,e,op,c,node,_=model
    validate(node)
    schema=json_schema()['$defs']['node']
    assert set(schema['required'])=={'schema','schema_version','node_id','energy_gauge','base_face','validation',
                                     'provenance',*BASE,*STATE,*OUTPUT}
    def broken(**change):
        n=copy.deepcopy(node);n.update(change);return n
    face=dict(node['base_face'],density=node['rho_base']*(1+1e-6))
    for n in [broken(mdot_g_s=-1.),broken(extra=1),broken(base_face=face),broken(Teff_emergent_K=node['Teff_emergent_K']*1.001),
              broken(L_base_erg_s=2*L0),broken(energy_gauge='other'),broken(R_phot_cm=node['R_base_cm']/2),
              broken(validation=dict(node['validation'],escapes=1)),broken(provenance=dict(node['provenance'],vadat_sha256='abc'))]:
        with pytest.raises(ValueError):validate(n)
    n=copy.deepcopy(node);del n['v_inf_cm_s']
    with pytest.raises(ValueError):validate(n)
    with pytest.raises(ValueError,match='share base data'):
        nodes_of(dict(schema='qstars_wind_node_set',schema_version=1,nodes=[node,dict(node,node_id='twin')]))
    bad=write(tmp_path/'bad.json',dict(node,v_inf_cm_s=0.))
    with pytest.raises(ValueError):run(node_config(bad),tmp_path/'runs')
    assert not (tmp_path/'runs').exists()  # refused before any run directory exists

def test_config_rejects_unevaluated_controls_and_incomplete_node_wind(model,tmp_path):
    _,_,_,c,_,path=model
    pulse=dict(kind='pulse',base=1e6,amplitude=1e5,start_years=1,end_years=2)
    # Previously accepted at config time, then a TypeError mid-run outside the retry tuple.
    for raw in [{'wind':{'terminal_speed_cm_s':pulse}},{'boundary':{'pressure_cgs':pulse}},{'wind':{'W':pulse}},
                {'boundary':{'bath_temperature_K':pulse}},{'boundary':{'external_optical_depth':pulse}},
                {'wind':{'model':'selected_nodes','drift_tolerance':.01}},
                {'wind':{'model':'selected_nodes','node_path':str(path)}},
                {'wind':{'model':'selected_nodes','node_path':str(path),'drift_tolerance':1.}},
                {'wind':{'model':'selected_nodes','node_path':str(path),'drift_tolerance':.01,'terminal_speed_cm_s':1e8}},
                {'wind':{'model':'selected_nodes','node_path':str(path),'drift_tolerance':.01,'rate_msun_per_year':.1}}]:
        with pytest.raises(ValueError):resolve(raw)
    assert resolve({'boundary':{'pressure_cgs':1e8}})['boundary']['pressure_cgs']==1e8
    (tmp_path/'run.yaml').write_text('wind:\n  model: selected_nodes\n  node_path: wind_node.json\n  drift_tolerance: 0.01\n')
    assert load(tmp_path/'run.yaml')['wind']['node_path']==str(path.resolve())

def test_single_node_values_and_drift_guard(model):
    s,e,op,c,node,_=model
    w=evaluate(s,e,op,c['wind'],L0)
    assert (w['rate'],w['terminal_speed'],w['emergent_luminosity'])==(node['mdot_g_s'],node['v_inf_cm_s'],node['L_emergent_erg_s'])
    assert w['max_abs_dln_drift']<1e-15 and rate(s,e,op,c['wind'],L0)==node['mdot_g_s']
    assert evaluate(s,e,op,c['wind'],L0*np.exp(.99*DELTA))['rate']==node['mdot_g_s']
    hot=s.copy();hot.T[-1]*=np.exp(1.1*DELTA)
    for state,L in [(hot,L0),(s,L0*np.exp(-1.1*DELTA)),(s,0.)]:
        with pytest.raises(DomainError):evaluate(state,e,op,c['wind'],L)

def test_node_sequence_loglinear_inside_hull_only(model,tmp_path):
    s,e,op,c,node,_=model
    far=make_node(s,e,'hot',heat=np.exp(.1),mdot=4*MDOT,vinf=2*VINF)
    c['wind']['node_path']=str(write(tmp_path/'set.json',dict(schema='qstars_wind_node_set',schema_version=1,nodes=[node,far])))
    def at(T=1.,rho=1.):
        x=s.copy();x.T[-1]*=T;x.dm[-1]*=rho  # rho=dm/volume
        return evaluate(x,e,op,c['wind'],L0)
    w=at(np.exp(.05))
    assert w['rate']==pytest.approx(2*MDOT,rel=1e-12) and w['terminal_speed']==pytest.approx(2**.5*VINF,rel=1e-12)
    assert w['emergent_effective_temperature']==pytest.approx((w['emergent_luminosity']/(4*PI*SIGMA*w['photospheric_radius']**2))**.25,rel=1e-12)
    assert at(np.exp(.1+.5*DELTA))['rate']==4*MDOT  # vertex value held within delta, never extrapolated
    for T,rho in [(np.exp(.1+2*DELTA),1.),(np.exp(-2*DELTA),1.),(np.exp(.05),np.exp(2*DELTA))]:
        with pytest.raises(DomainError):at(T,rho)

def test_step_ledger_node_terminal_speed_and_diagnostics(model):
    s,e,op,c,node,_=model;dt=.001*YEAR;x=s
    for _ in range(3):
        old=x;x,d=attempt_step(old,dt,e,op,c);r=ledger.residuals(x,e)
        assert abs(r['mass_residual_relative'])<1e-12 and abs(r['energy_residual_relative'])<1e-12
        assert x.ledger['wind_mass']-old.ledger['wind_mass']==pytest.approx(dt*MDOT,rel=1e-12)
        _,stream=exchange(old,dt,e,c,rates(old,e,c['engine']),MDOT,terminal_speed=VINF)
        assert d['wind_cost']==stream['wind_cost']>0
        assert x.ledger['wind_mechanical_energy']-old.ledger['wind_mechanical_energy']==pytest.approx(d['wind_cost'],rel=1e-12)
        assert d['wind_node']=='synthetic' and d['wind_terminal_speed_cm_s']==VINF
        # Envelope phi/h are discrete outer-cell values; the node uses -GM/R at the face.
        assert abs(d['interface_energy_defect'])<1e-3 and abs(d['wind_power_defect'])<1e-3

def test_retry_is_transactional(model):
    s,e,op,c,node,_=model;before=s.copy()
    with pytest.raises(StepRejected):attempt_step(s,1e6*YEAR,e,op,c)  # >10% of the boundary cell
    ramp=copy.deepcopy(c);ramp['boundary']['luminosity_erg_s']=dict(kind='table',time_years=[0,.001],values=[L0,1.02*L0])
    with pytest.raises(DomainError,match='neighbourhood'):attempt_step(s,.001*YEAR,e,op,ramp)
    np.testing.assert_array_equal(s.dm,before.dm);np.testing.assert_array_equal(s.T,before.T);assert s.ledger==before.ledger
    a,_=attempt_step(s,1e-4*YEAR,e,op,ramp);b,_=attempt_step(before,1e-4*YEAR,e,op,ramp)
    np.testing.assert_array_equal(a.T,b.T);np.testing.assert_array_equal(a.r,b.r);assert a.ledger==b.ledger

def test_restart_bitwise_with_hashed_node_copy(model,tmp_path):
    s,e,op,c,node,path=model
    p,report=run(c,tmp_path,stop_after_steps=2);restarted,cc=load_checkpoint(p/'checkpoint.npz')
    assert (p/'wind_node.json').read_bytes()==path.read_bytes()
    assert json.loads((p/'external_inputs.json').read_text())['wind_node_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    assert cc['wind']['node_path']==str((p/'wind_node.json').resolve())
    p2,r2=run(cc,tmp_path,restart=restarted);p3,r3=run(c,tmp_path)
    a,_=load_checkpoint(p2/'checkpoint.npz');b,_=load_checkpoint(p3/'checkpoint.npz')
    for k in ('T','r','dm','v'):np.testing.assert_array_equal(getattr(a,k),getattr(b,k))
    assert a.ledger==b.ledger and a.time==b.time
    assert r2['status']==r3['status']=='completed_requested_interval'
    assert r3['wind_node']==dict(validated=False,node_flags_all_true=False,interface_conservation_passed=True)
    assert abs(r3['conservation']['energy_residual_relative'])<1e-12
    h=np.genfromtxt(p3/'history.csv',delimiter=',',names=True)
    assert np.all(np.isnan(h['photospheric_radius_cm']))  # the node's held R_phot is not this state's photosphere
    np.testing.assert_array_equal(h['wind_node_photospheric_radius_cm'],node['R_phot_cm'])
    np.testing.assert_array_equal(h['wind_terminal_speed_cm_s'],VINF)
    events=[json.loads(x) for x in (p3/'events.jsonl').read_text().splitlines()]
    assert all('interface_energy_defect' in x for x in events if x['event']=='accepted')

def test_self_declared_flags_never_validate(model,tmp_path):
    s,e,op,c,node,_=model
    claimed=write(tmp_path/'claimed.json',dict(node,validation={k:True for k in VALIDATION}))
    p,r=run(node_config(claimed,evolution={'end_time_years':.001}),tmp_path)
    assert r['status']=='completed_requested_interval'
    assert r['wind_node']==dict(validated=False,node_flags_all_true=True,interface_conservation_passed=True)
    assert json.loads((p/'wind_node_interface.json').read_text())[0]['node_flags_all_true']

def test_remap_leaving_neighbourhood_is_rejected(model,tmp_path):
    s,e,op,c,node,path=model
    # Halving cells keeps the outer rho but lowers its T by ~7.9% (> DELTA); steps drift ~1e-4.
    p,r=run(node_config(path,numerics={'refine_every_steps':1,'max_zones':32}),tmp_path)
    events=[json.loads(x) for x in (p/'events.jsonl').read_text().splitlines()]
    rejected=[x for x in events if x['event']=='remap_rejected']
    assert rejected and all(len(x['errors'])==1 and 'neighbourhood' in x['errors'][0] for x in rejected)
    assert not any(x['event']=='remapped' for x in events)
    assert r['status']=='completed_requested_interval' and len(rejected)==r['accepted_steps_this_run']>=2
    last,_=load_checkpoint(p/'checkpoint.npz');assert len(last.dm)==16

def test_drift_and_outside_node_are_domain_exits(model,tmp_path):
    s,e,op,c,node,path=model
    ramp=node_config(path,boundary={'luminosity_erg_s':dict(kind='table',time_years=[0,.005],values=[L0,1.1*L0])},
                     evolution={'end_time_years':.005})
    p,r=run(ramp,tmp_path)
    assert r['status']=='model_domain_exit' and 'neighbourhood' in r['reason'] and r['accepted_steps_this_run']>0
    last,cc=load_checkpoint(p/'checkpoint.npz')
    assert evaluate(last,e,op,cc['wind'],float(np.interp(last.time/YEAR,[0,.005],[L0,1.1*L0])))['max_abs_dln_drift']<=DELTA
    away=write(tmp_path/'away.json',make_node(s,e,'away',heat=np.exp(5*DELTA)))
    p,r=run(node_config(away),tmp_path)
    assert r['status']=='model_domain_exit' and 'neighbourhood' in r['reason'] and r['accepted_steps_this_run']==0

def test_energy_gauge_mismatch_is_rejected(model,tmp_path):
    s,e,op,c,node,_=model
    cmfgen=write(tmp_path/'cmfgen_gauge.json',make_node(s,e,gauge='neutral_H_He_ground'))
    p,r=run(node_config(cmfgen),tmp_path)
    assert r['status']=='model_domain_exit' and 'gauges must agree' in r['reason']
    assert not (p/'history.csv').exists()
    with pytest.raises(DomainError,match='undeclared'):interface_check((node,),e,'mesa')
    assert interface_check((node,),e,'ideal_gas_radiation')[0]['conservation_passed']
