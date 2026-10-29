"""Gate 6 Stage 6: the quasi-static Ball path in the configuration schema, checkpoints and the CLI.

U13: configuration cross-validation (config.resolve). Every other
configuration resolves exactly as before; a quasi_static configuration
gets the cavity, quasistatic, stepcontrol, relaxation and replication
sections with every default, the Bondi cavity goes with quasi_static
mechanics, the flux-unknown solve, no wind and zero feeding, and a key the
quasi-static path does not read is an error, not a silent no-op.
U16: checkpoint schema 7 (io.save_checkpoint/load_checkpoint) round trip of
cavity_mass, photosphere_density, kinetic_tracked and aux; evolution-path
states are still written as schema 6 with the same arrays; schemas 1-6 load
with the defaults of the later fields. Also: the Ball history columns
(photospheric_radius_cm is r[-1] with the Ball photosphere only), `qstars
run` and `qstars restart` of a quasi_static configuration against an
uninterrupted library run (bitwise), and benchmarks/ball2011_track.yaml
against the registered plan. Stage 6 review: evolution.run and `qstars plot`
refuse quasi-static configurations and run directories, and a setting of the
wrong type is a configuration error (exit 2), not a traceback.

No Ball configuration is run, with one opt-in exception: the amendment 3
test builds the registered C1 seeds on the compiled MESA bridge and takes
single relaxation steps from them (nothing is written, nothing is compared
with the published tables). The CLI toy is the analytic cavity of
test_quasistatic.py made hot and compact (100 Msun, M_BH = 0.05 Msun,
15 Rsun, 32 zones, cell boundary state, eta = 5e-5), so that the fully
ionized EOS passes every enforced screen of the published_assumption policy.
"""
import copy,importlib.util,json,math,os,sys
from pathlib import Path
import numpy as np
import pytest
import yaml
sys.path.insert(0,str(Path(__file__).parent))
from test_ball_restart import same_state,rows
from qstars import config as configuration,io,ledger,quasistatic,stepcontrol
from qstars.cli import main
from qstars.config import load,resolve
from qstars.constants import MSUN,RSUN
from qstars.evolution import physics
from qstars.structure import initialize

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('validate_ball_track',ROOT/'scripts/validate_ball_track.py')
vbt=importlib.util.module_from_spec(spec);spec.loader.exec_module(vbt)
PLAN=json.loads(vbt.PLAN.read_text())

QUASI=dict(numerics=dict(mechanics='quasi_static'),boundary=dict(inner='bondi_cavity',thermal='stars_photosphere'),
           transport=dict(convection='mlt_ml1'),engine=dict(model='ball2012_eq3_15'))
MESA=dict(eos='mesa',mesa=dict(library_path='lib',mesa_dir='mesa'),opacity=dict(model='mesa'))
LANDING,END=.0500001,.05000015
TOY=dict(name='quasi_toy',composition=dict(X=.7,Y=.3,Z=0.),engine=dict(model='ball2012_eq3_15',coupling=5e-5),
         transport=dict(convection='mlt_ml1'),boundary=dict(inner='bondi_cavity',thermal='stars_photosphere'),
         numerics=dict(zones=32,mechanics='quasi_static'),cavity=dict(boundary_state='cell'),
         stepcontrol=dict(metric='mapped',delta_target=.03,dt_initial_years=1e-3),relaxation=dict(mesh_cycles=0),
         replication=dict(total_mass_msun=100.,initial_bh_mass_msun=.05,seed_radii_rsun=[15.],landings_msun=[LANDING],
                          end_bh_mass_msun=END))

def quasi(**changes):
    """QUASI with sections updated (a non-mapping value replaces the section)."""
    raw=copy.deepcopy(QUASI)
    for section,values in changes.items():
        if isinstance(values,dict) and isinstance(raw.get(section),dict):raw[section].update(values)
        else:raw[section]=values
    return raw

def refused(raw,match):
    with pytest.raises(ValueError,match=match):resolve(raw)

# ---------------------------------------------------------------- U13 configuration

def test_other_configurations_resolve_exactly_as_before():
    assert resolve({})==configuration.DEFAULT and 'tables' not in resolve({})['opacity']
    for path in sorted((ROOT/'configs').glob('*.yaml')):
        if path.name=='sweep.yaml':continue  # a sweep grid, not a configuration
        c=load(path)
        assert set(c)==set(configuration.DEFAULT) and all(set(c[k])==set(v) for k,v in configuration.DEFAULT.items() if isinstance(v,dict)),path.name
    # The thesis eq 3.15 law is the quest2025_eq11 law under its source name: any path may select it.
    assert resolve(dict(engine=dict(model='ball2012_eq3_15')))['engine']['model']=='ball2012_eq3_15'
    # The selected-node wind of the evolution path is untouched, and never a quasi-static wind.
    node=dict(model='selected_nodes',node_path='node.json',drift_tolerance=.1)
    assert resolve(dict(wind=node))['wind']['model']=='selected_nodes'
    refused(quasi(wind=node),'requires wind.model off')

def test_quasi_static_configuration_holds_every_default():
    c=resolve(QUASI)
    assert c['cavity']==quasistatic.DEFAULTS['cavity'] and c['quasistatic']==quasistatic.DEFAULTS['quasistatic']
    assert c['stepcontrol']==stepcontrol.DEFAULTS and c['relaxation']==quasistatic.RELAX_DEFAULTS
    assert c['replication']==quasistatic.REPLICATION_DEFAULTS and c['opacity']['tables']=='default'
    assert {k:c[k] for k in configuration.DEFAULT if k not in ('numerics','boundary','transport','engine','opacity')}== \
           {k:v for k,v in configuration.DEFAULT.items() if k not in ('numerics','boundary','transport','engine','opacity')}
    # The step and the controller read exactly the resolved sections; the module defaults are not shared.
    s=quasistatic.settings(c)
    assert s['cavity']==c['cavity'] and s['quasistatic']==c['quasistatic'] and stepcontrol.settings(c)==c['stepcontrol']
    c['replication']['landings_msun'].append(2000.);assert quasistatic.REPLICATION_DEFAULTS['landings_msun']==[100.,200.,500.,1000.]
    # A resolved configuration is a fixed point (checkpoints store it; `qstars restart` resolves it again).
    c=resolve(quasi(fidelity='research',cavity=dict(profile='r-1/2',boundary_state='cell'),stepcontrol=dict(metric='mapped'),
                    quasistatic=dict(base='fixed_mass',split=1.5,merge=.67),relaxation=dict(mesh_cycles=3),
                    replication=dict(total_mass_msun=100.,initial_bh_mass_msun=1.,landings_msun=[1.2,2],end_bh_mass_msun=3.)))
    assert c['numerics']['zones']==96 and c['cavity']['profile']=='r-1/2' and c['quasistatic']['merge']==.67
    assert resolve(c)==c and resolve(json.loads(json.dumps(c)))==c and resolve(yaml.safe_load(yaml.safe_dump(c)))==c

def test_quasi_static_choices_exist_only_on_the_quasi_static_path():
    refused(dict(boundary=dict(inner='bondi_cavity')),'go together')
    refused(quasi(boundary=dict(inner='fixed_excision')),'go together')
    refused(dict(boundary=dict(thermal='stars_photosphere')),'stars_photosphere requires numerics.mechanics quasi_static')
    for closure in ('mlt_ml1','mlt_thesis_lossless'):
        refused(dict(transport=dict(convection=closure)),f'{closure} requires numerics.mechanics quasi_static')
        assert resolve(quasi(transport=dict(convection=closure)))['transport']['convection']==closure
    for section in configuration.QUASI_STATIC_SECTIONS:
        refused({section:{}},f'{section} require numerics.mechanics quasi_static')
        refused({section:{},'numerics':dict(mechanics='backward_euler')},'require numerics.mechanics quasi_static')
        refused(quasi(**{section:5}),f'{section} must be a mapping')
        refused(quasi(**{section:dict(typo=1)}),f'Unknown configuration key {section}.typo')
    refused(dict(opacity=dict(tables='stars2004_switch')),'opacity.tables require numerics.mechanics quasi_static')

def test_bondi_cavity_cross_checks():
    """Spec Stage 6: the cavity requires quasi_static, newton_flux, wind off and zero feeding."""
    refused(quasi(numerics=dict(nonlinear_solver='newton')),'requires numerics.nonlinear_solver newton_flux')
    refused(quasi(numerics=dict(nonlinear_solver='hybr')),'requires numerics.nonlinear_solver newton_flux')
    refused(quasi(wind=dict(model='prescribed',rate_msun_per_year=1e-3)),'requires wind.model off')
    refused(quasi(wind=dict(model='dotan2011')),'requires wind.model off')
    refused(quasi(feeding=dict(rate_msun_per_year=.01)),'requires zero feeding')
    refused(quasi(feeding=dict(rate_msun_per_year=dict(kind='pulse',base=0.,amplitude=1.,start_years=1.,end_years=2.))),'requires zero feeding')
    refused(quasi(engine=dict(model='off')),'needs a central engine')
    refused(quasi(cavity=dict(active=False)),'cavity.active false')
    # Time controls have no breakpoints on this path; the Ball photosphere sets its own outer pressure.
    pulse=dict(kind='pulse',base=.1,amplitude=.1,start_years=1.,end_years=2.)
    refused(quasi(engine=dict(coupling=pulse)),'engine.coupling: expected a number')
    refused(quasi(boundary=dict(luminosity_erg_s=dict(kind='table',time_years=[0,1],values=[1e40,2e40]))),'luminosity_erg_s: expected a number')
    for key in configuration.PHOTOSPHERE_UNUSED:
        refused(quasi(boundary={key:1.}),f'boundary.{key} is not read')
    assert resolve(quasi(boundary=dict(thermal='grey',pressure_cgs=1e3,external_optical_depth=.5)))['boundary']['pressure_cgs']==1e3
    refused(quasi(opacity=dict(model='table',table_path='table.npz')),'electron_scattering or mesa')

def test_quasi_static_settings_are_validated(tmp_path,capsys):
    bad=dict(cavity=[dict(profile='r-2'),dict(profile=['r-3/2']),dict(profile={'r-3/2':1}),dict(boundary_state='edge'),dict(b=0.),dict(active=1)],
             quasistatic=[dict(inertia='midpoint'),dict(momentum='centred'),dict(base='ale'),dict(base_mass_fraction=.5),
                          dict(base_mass_fraction='a'),dict(base_mass_fraction=True),dict(tolerance=0.),
                          dict(max_iterations=0),dict(polish_iterations=-1),dict(max_passes=1.5),dict(pass_fraction=1.),dict(remesh=1),
                          dict(split=1.),dict(merge=1.),dict(split=1.5,merge=.8),dict(energy_per_step=0.),dict(mass_per_step=float('nan')),
                          dict(ptc=1),dict(ptc='F1'),dict(photosphere_bracket='bisect'),dict(photosphere_bracket=True)],
             stepcontrol=[dict(metric='eulerian'),dict(delta_target=0.),dict(safety=1.5),dict(max_factor=.9),dict(reject_ratio=1.),
                          dict(hold_steps=-1),dict(retry_factor=1.),dict(max_retries=0),dict(max_retries=True),dict(dt_initial_years=0.),
                          dict(dt_initial_years=10.,dt_max_years=1.)],
             relaxation=[dict(dt_years=0.),dict(dt_years=1.,dt_max_years=.5),dict(converged_dt_years=1e9),dict(tolerance=-1.),
                         dict(max_steps=0),dict(max_retries=-1),dict(mesh_cycles=-1),dict(mesh_cycles=1.5),dict(mesh_tolerance=0.)],
             replication=[dict(policy='lenient'),dict(plan=1),dict(total_mass_msun=0.),dict(initial_bh_mass_msun=1e4),dict(seed_radii_rsun=[]),
                          dict(seed_radii_rsun=[2000.,-1.]),dict(seed_radii_rsun=2000.),dict(landings_msun=[200.,100.]),
                          dict(landings_msun=[5.,100.]),dict(end_bh_mass_msun=5.),dict(end_bh_mass_msun=-1.),dict(amendments=[1]),
                          dict(seed_surface='hot'),dict(seed_surface=True)])
    for section,cases in bad.items():
        for values in cases:
            with pytest.raises(ValueError,match='|'.join(values)):resolve(quasi(**{section:values}))  # the error names the key
    assert resolve(quasi(replication=dict(landings_msun=[],end_bh_mass_msun=6.)))['replication']['end_bh_mass_msun']==6.
    # Options that are off unless a registered amendment or fallback switches them on: the grey seed surface
    # (plan amendment 3; it is the grey relation of the Ball photosphere) and fallback F1.
    c=resolve(QUASI);assert (c['replication']['seed_surface'],c['quasistatic']['ptc'],c['quasistatic']['photosphere_bracket'])==('eos',False,'doubling')
    c=resolve(quasi(replication=dict(seed_surface='grey'),quasistatic=dict(ptc=True,photosphere_bracket='scan')))
    assert (c['replication']['seed_surface'],c['quasistatic']['ptc'],c['quasistatic']['photosphere_bracket'])==('grey',True,'scan') and resolve(c)==c
    refused(quasi(boundary=dict(thermal='grey'),replication=dict(seed_surface='grey')),'seed_surface grey .* needs boundary.thermal stars_photosphere')
    # A value of the wrong type is a configuration error too: `qstars run` exits 2 with the key, without a traceback.
    for section,values in (('quasistatic',dict(base_mass_fraction='a')),('cavity',dict(profile=['r-3/2']))):
        (tmp_path/'bad.yaml').write_text(yaml.safe_dump({**TOY,section:{**TOY.get(section,{}),**values}}))
        with pytest.raises(SystemExit) as exit:main(['run',str(tmp_path/'bad.yaml'),'--output',str(tmp_path/'never')])
        err=capsys.readouterr().err
        assert exit.value.code==2 and err.startswith('qstars: ') and next(iter(values)) in err and not (tmp_path/'never').exists()

def test_opacity_tables_need_the_mesa_bridge_and_a_declared_policy():
    assert configuration.OPACITY_TABLES==('default','stars2004_switch','stars2004_switch_strict')
    refused(quasi(opacity=dict(tables='opal')),'Unsupported choice opal')
    for tables in ('stars2004_switch','stars2004_switch_strict'):
        refused(quasi(opacity=dict(tables=tables)),'needs opacity.model mesa')
    mesa=lambda tables,policy:resolve({**quasi(replication=dict(policy=policy)),**MESA,'opacity':dict(model='mesa',tables=tables)})
    # The STARS nearest-table-value rule is a published-assumption rule (plan model.opacity.off_table).
    assert mesa('stars2004_switch','published_assumption')['opacity']['tables']=='stars2004_switch'
    with pytest.raises(ValueError,match='published-assumption replication rule'):mesa('stars2004_switch','production')
    assert mesa('stars2004_switch_strict','production')['replication']['policy']=='production'
    assert mesa('default','production')['opacity']=={'model':'mesa','table_path':'','tables':'default'}
    # Without an ABI-7 bridge the switch is an error before any run directory exists.
    c=resolve(QUASI);eos,opacity=quasistatic.physics(c)
    assert quasistatic.opacity_tables(c,eos,opacity) is opacity
    c['opacity']['tables']='stars2004_switch'
    with pytest.raises(ValueError,match='ABI 7'):quasistatic.opacity_tables(c,eos,opacity)

def test_keys_the_quasi_static_path_does_not_read_are_errors():
    cases=[('initial',dict(radius_cm=2e14)),('initial',dict(bh_mass_msun=5.)),('evolution',dict(end_time_years=1e6)),
           ('output',dict(profile_every_steps=5)),('numerics',dict(mesh='log_radius')),('numerics',dict(max_dt_years=1.)),
           ('numerics',dict(solver_tolerance=1e-12)),('numerics',dict(max_retries=3)),('numerics',dict(refine_every_steps=2))]
    for section,values in cases:
        refused(quasi(**{section:values}),f'{section}.{next(iter(values))} is not read by quasi_static mechanics')
    # A fidelity preset is the default of its keys; the keys the path reads stay configurable.
    c=resolve(quasi(fidelity='validation',numerics=dict(zones=48,conservation_tolerance=1e-8),evolution=dict(max_steps=50,max_wall_seconds=60.),
                    output=dict(checkpoint_every_steps=3),validity=dict(min_temperature_K=2e5)))
    assert c['numerics']['solver_tolerance']==1e-11 and c['numerics']['zones']==48 and c['evolution']['max_steps']==50

# ---------------------------------------------------------------- U16 checkpoint schema 7

@pytest.fixture(scope='module')
def seeded():
    """A quasi-static state with every schema-7 field off its default (no step is taken)."""
    c=resolve(TOY);eos,opacity=quasistatic.physics(c);rep=c['replication']
    s,_=quasistatic.seed(c,eos,opacity,32,rep['initial_bh_mass_msun']*MSUN,rep['total_mass_msun']*MSUN,15*RSUN)
    ledger.initialize(s,eos,cavity=True)
    s.photosphere_density=1/3e9;s.face_luminosity=np.linspace(1e38,1.5e38,33);s.radius_compensation=np.linspace(0.,1e-3,33)
    s.aux={'controller':dict(dt_next=1/3,hold=2,retries=1),'predictor':dict(time=[0.,.1],rest_rate=[1/7,2/7],cavity_mass=[1.,2.]),
           'driver':dict(floor_run=0,landed=[LANDING*MSUN],screen_retries=0)}
    return s,c

def arrays(path):
    with np.load(path,allow_pickle=False) as f:return {k:f[k].copy() for k in f.files}

SCHEMA_6=['schema_version','time','bh_mass','dm','r','v','T','cell_volume','radius_compensation','face_luminosity',
          'reference_mass_fractions','ledger','steps','next_dt','config']
SCHEMA_7=SCHEMA_6+['cavity_mass','photosphere_density','kinetic_tracked','aux']

def test_schema_7_round_trip_is_exact(seeded,tmp_path):
    s,c=seeded
    assert s.cavity_mass>0 and not s.kinetic_tracked and s.next_dt==0
    io.save_checkpoint(tmp_path/'a.npz',s,c)
    data=arrays(tmp_path/'a.npz')
    assert int(data['schema_version'])==7 and sorted(data)==sorted(SCHEMA_7) and not (tmp_path/'a.tmp').exists()
    assert float(data['cavity_mass'])==s.cavity_mass and float(data['photosphere_density'])==s.photosphere_density
    assert not bool(data['kinetic_tracked']) and json.loads(str(data['aux']))==s.aux
    u,cc=io.load_checkpoint(tmp_path/'a.npz')
    same_state(s,u);assert u.next_dt==0 and cc==c
    # quasistatic.save_checkpoint/load_checkpoint are this format, with the configuration as stored.
    quasistatic.save_checkpoint(tmp_path/'b.npz',s,c)
    b=arrays(tmp_path/'b.npz')
    assert sorted(b)==sorted(data) and all(b[k].dtype==data[k].dtype and b[k].tobytes()==data[k].tobytes() for k in data)
    toy=copy.deepcopy(c);toy['numerics']['mechanics']='midpoint'  # like the hand-built toys of test_quasistatic.py: not a resolved configuration
    quasistatic.save_checkpoint(tmp_path/'c.npz',s,toy)
    v,stored=quasistatic.load_checkpoint(tmp_path/'c.npz')
    same_state(s,v);assert stored==json.loads(json.dumps(toy))
    with pytest.raises(ValueError):io.load_checkpoint(tmp_path/'c.npz')  # `qstars restart` resolves it
    # Fields at their defaults survive as such.
    t=s.copy();t.photosphere_density=None;t.face_luminosity=None;t.aux={};t.kinetic_tracked=True;t.cavity_mass=0.;t.next_dt=1e5
    quasistatic.save_checkpoint(tmp_path/'d.npz',t,c)
    assert int(arrays(tmp_path/'d.npz')['schema_version'])==7
    w,_=io.load_checkpoint(tmp_path/'d.npz')
    assert w.photosphere_density is None and w.face_luminosity is None and w.aux=={} and w.kinetic_tracked and w.cavity_mass==0

def test_evolution_checkpoints_stay_schema_6_and_older_schemas_get_the_defaults(tmp_path):
    c=resolve({'numerics':{'zones':16},'engine':{'model':'off'},'boundary':{'thermal':'insulated'},'transport':{'convection':'off'}})
    eos,_=physics(c);s,_=initialize(c,eos);ledger.initialize(s,eos)
    assert (s.cavity_mass,s.photosphere_density,s.kinetic_tracked,s.aux)==(0.,None,True,{})
    io.save_checkpoint(tmp_path/'six.npz',s,c)
    data=arrays(tmp_path/'six.npz')
    assert int(data['schema_version'])==6 and list(data)==SCHEMA_6
    with pytest.raises(ValueError,match='schema 7'):io.save_checkpoint(tmp_path/'bad.npz',s,c,schema=6)
    removed={1:['cell_volume','face_luminosity','reference_mass_fractions','radius_compensation'],
             2:['face_luminosity','reference_mass_fractions','radius_compensation'],3:['reference_mass_fractions','radius_compensation'],
             4:['radius_compensation'],5:[],6:[]}
    for schema,missing in removed.items():
        old={k:v for k,v in data.items() if k not in missing};old['schema_version']=np.array(schema)
        if schema<6:old['ledger']=np.array(json.dumps({k:v for k,v in s.ledger.items() if k!='numerical_kinetic_heating'}))
        np.savez_compressed(tmp_path/f'old{schema}.npz',**old)
        u,cc=io.load_checkpoint(tmp_path/f'old{schema}.npz')
        assert (u.cavity_mass,u.photosphere_density,u.kinetic_tracked,u.aux)==(0.,None,True,{}),schema
        assert np.array_equal(u.r,s.r) and np.array_equal(u.T,s.T) and np.array_equal(u.dm,s.dm) and u.ledger==s.ledger and cc==c,schema
        assert np.allclose(u.volume,s.volume,rtol=1e-14,atol=0)
    # A state that carries a quasi-static field is never written without it.
    t=s.copy();t.aux={'note':1};io.save_checkpoint(tmp_path/'seven.npz',t,c)
    assert int(arrays(tmp_path/'seven.npz')['schema_version'])==7 and io.load_checkpoint(tmp_path/'seven.npz')[0].aux=={'note':1}
    io.save_checkpoint(tmp_path/'forced.npz',s,c,schema=7)
    u,_=io.load_checkpoint(tmp_path/'forced.npz')
    assert list(arrays(tmp_path/'forced.npz'))==SCHEMA_7 and (u.cavity_mass,u.photosphere_density,u.kinetic_tracked,u.aux)==(0.,None,True,{})

def test_checkpoint_schema_and_quasi_static_fields_are_validated(seeded,tmp_path):
    s,c=seeded
    io.save_checkpoint(tmp_path/'a.npz',s,c);data=arrays(tmp_path/'a.npz')
    def changed(name,**values):
        np.savez_compressed(tmp_path/name,**{**data,**{k:np.array(v) for k,v in values.items()}});return tmp_path/name
    with pytest.raises(ValueError,match='Unsupported checkpoint schema'):io.load_checkpoint(changed('eight.npz',schema_version=8))
    for values in (dict(cavity_mass=-1.),dict(cavity_mass=np.inf),dict(photosphere_density=0.),dict(photosphere_density=-1e-9),dict(aux='[1]')):
        with pytest.raises(ValueError,match='invalid cavity, photosphere or driver state'):io.load_checkpoint(changed('bad.npz',**values))
    # A t = 0 quasi-static model has next_dt 0 (its steps are controlled from aux); an evolution checkpoint never does.
    assert io.load_checkpoint(changed('zero.npz',next_dt=0.))[0].next_dt==0
    six={k:v for k,v in data.items() if k in SCHEMA_6};six['schema_version']=np.array(6)
    np.savez_compressed(tmp_path/'six.npz',**six)
    with pytest.raises(ValueError,match='step controls'):io.load_checkpoint(tmp_path/'six.npz')
    # The Stage 5 quasi-static format ('qstars_quasistatic' 1) holds the schema-7 fields and still loads.
    stage5=changed('stage5.npz',schema_version=1,format='qstars_quasistatic')
    u,stored=quasistatic.load_checkpoint(stage5);same_state(s,u);assert stored==c
    with pytest.raises(ValueError,match='Unsupported checkpoint schema'):io.load_checkpoint(changed('other.npz',schema_version=1,format='other'))
    with pytest.raises(ValueError,match='Unsupported checkpoint schema'):io.load_checkpoint(changed('other.npz',format='qstars_quasistatic'))

# ---------------------------------------------------------------- CLI dispatch and history columns

def directory(parent):
    found=[p for p in Path(parent).iterdir() if p.is_dir()]
    assert len(found)==1;return found[0]

@pytest.fixture(scope='module')
def cli_run(tmp_path_factory):
    """`qstars run` of the toy configuration, stopped after 6 accepted steps."""
    root=tmp_path_factory.mktemp('quasi_cli')
    (root/'toy.yaml').write_text(yaml.safe_dump(TOY))
    assert main(['run',str(root/'toy.yaml'),'--output',str(root/'first'),'--stop-after-steps','6'])==0
    return root,directory(root/'first')

def test_cli_runs_a_quasi_static_configuration(cli_run):
    root,path=cli_run
    report=json.loads((path/'termination.json').read_text())
    assert report['status']=='control_stop' and report['final_step']==6 and report['rejections']==0
    assert abs(report['conservation']['energy_residual_relative'])<1e-13 and abs(report['conservation']['mass_residual_relative'])<1e-14
    stored=yaml.safe_load((path/'config.yaml').read_text())
    assert stored==load(root/'toy.yaml') and stored['numerics']['mechanics']=='quasi_static' and stored['relaxation']['mesh_cycles']==0
    relax=json.loads((path/'relax.json').read_text())
    assert [x['status'] for x in relax['seeds']]==['converged'] and relax['relaxation']==stored['relaxation']
    for name in ('row0.npz','profile_row0.npz','checkpoint.npz','profile_final.npz','validity_warnings.json','screen_violations.json',
                 'events.jsonl','metadata.json','source_snapshot.zip'):
        assert (path/name).exists(),name
    assert json.loads((path/'validity_warnings.json').read_text())['policy']=='published_assumption'
    # t = 0 is the relaxed model; both checkpoints are schema 7 with the resolved configuration.
    row0,c=io.load_checkpoint(path/'row0.npz');state,cc=io.load_checkpoint(path/'checkpoint.npz')
    assert int(arrays(path/'row0.npz')['schema_version'])==int(arrays(path/'checkpoint.npz')['schema_version'])==7
    assert c==cc==stored and row0.time==0 and row0.steps==0 and row0.aux=={} and row0.bh_mass==.05*MSUN
    assert state.steps==6 and set(state.aux)=={'controller','predictor','driver'} and not state.kinetic_tracked and state.cavity_mass>0

def test_ball_history_columns(cli_run):
    root,path=cli_run
    state,c=io.load_checkpoint(path/'checkpoint.npz');eos,opacity=quasistatic.physics(c)
    history=rows(path)
    assert [r['step'] for r in history]==[str(i) for i in range(7)] and [r['event'] for r in history]==['start']+['step']*6
    for column in ('time_myr','bh_mass_msun','mdot_bh_msun_per_year','cavity_mass_msun','luminosity_lsun','base_density_cgs',
                   'base_temperature_K','effective_temperature_K','inner_radius_rsun','radius_rsun','photospheric_radius_cm',
                   'photosphere_density_cgs','screen_policy','mass_residual_relative','energy_residual_relative','ledger_cavity_swept_mass'):
        assert all(r[column]!='' for r in history),column
    # The outer face is the tau = 2/3 photosphere of the Ball surface: photospheric_radius_cm = r[-1].
    assert float(history[-1]['photospheric_radius_cm'])==state.r[-1]
    assert all(math.isclose(float(r['photospheric_radius_cm']),float(r['radius_rsun'])*RSUN,rel_tol=1e-15) for r in history)
    # io.history_row is the Ball row for a quasi_static configuration.
    row=io.history_row(state,eos,opacity,c)
    assert row==quasistatic.history_row(state,eos,opacity,c) and row['photospheric_radius_cm']==state.r[-1]
    assert list(row)==list(history[0])
    # Any other surface has no photosphere at the outer face, as on the evolution path (NaN).
    grey=copy.deepcopy(c);grey['boundary']['thermal']='grey'
    assert math.isnan(io.history_row(state,eos,opacity,grey)['photospheric_radius_cm'])
    d=resolve({'numerics':{'zones':16}});e,op=physics(d);s,_=initialize(d,e);ledger.initialize(s,e)
    legacy=io.history_row(s,e,op,d)
    assert math.isnan(legacy['photospheric_radius_cm']) and 'cavity_mass_msun' not in legacy and 'radius_cm' in legacy

def test_cli_restart_continues_bitwise(cli_run,tmp_path,capsys):
    root,path=cli_run
    assert main(['restart',str(path/'checkpoint.npz'),'--output',str(tmp_path/'resumed')])==0
    resumed=directory(tmp_path/'resumed')
    report=json.loads((resumed/'termination.json').read_text())
    assert report['status']=='completed' and report['restarted'] and report['final_bh_mass_msun']>=END
    assert report['landed_msun']==pytest.approx([LANDING],rel=1e-15)
    assert not (resumed/'row0.npz').exists() and not (resumed/'relax.json').exists()
    assert (resumed/f'landing_{LANDING:g}msun.npz').exists() and (resumed/f'profile_landing_{LANDING:g}msun.npz').exists()
    # The uninterrupted track from the same t = 0 model, through the library.
    row0,c=quasistatic.load_checkpoint(path/'row0.npz');eos,opacity=quasistatic.physics(c)
    _,full=quasistatic.run(row0,eos,opacity,c,tmp_path/'full',landings=[LANDING*MSUN],end_bh_mass=END*MSUN,policy='published_assumption',
                           max_accepted_steps=c['evolution']['max_steps'],checkpoint_every=c['output']['checkpoint_every_steps'])
    assert full['status']=='completed' and full['final_step']==report['final_step']>8
    a,_=quasistatic.load_checkpoint(tmp_path/'full'/'checkpoint.npz');b,cc=quasistatic.load_checkpoint(resumed/'checkpoint.npz')
    same_state(a,b);assert cc==c
    first,second,whole=rows(path),rows(resumed),rows(tmp_path/'full')
    assert first==whole[:7] and second[1:]==whole[7:] and second[0]['event']=='restart'
    assert [r['event'] for r in whole].count('landing')==1
    landed=next(r for r in whole if r['event']=='landing')
    assert abs(float(landed['bh_mass_msun'])/LANDING-1)<=1e-12 and float(landed['landed_target_msun'])==pytest.approx(LANDING,rel=1e-15)
    # A quasi-static track ends on black-hole mass, not on time.
    capsys.readouterr()
    with pytest.raises(SystemExit) as exit:main(['restart',str(path/'checkpoint.npz'),'--until-years','1','--output',str(tmp_path/'never')])
    assert exit.value.code==2 and 'ends on black-hole mass' in capsys.readouterr().err and not (tmp_path/'never').exists()

def test_evolution_run_and_plot_refuse_quasi_static_runs(cli_run,tmp_path,capsys):
    # evolution.run has no moving cavity: it would run the fixed-excision dynamic step on the `initial` model.
    from qstars.evolution import run
    with pytest.raises(ValueError,match='quasistatic.run_track'):run(resolve(TOY),tmp_path/'never',stop_after_steps=2)
    assert not (tmp_path/'never').exists()
    # The Ball history and profiles have other columns and keys than the evolution plots read.
    root,path=cli_run;capsys.readouterr()
    with pytest.raises(SystemExit) as exit:main(['plot',str(path)])
    assert exit.value.code==2 and 'does not support quasi-static runs' in capsys.readouterr().err and not list(path.glob('*.png'))

def test_run_track_refusals_and_a_failed_relaxation(tmp_path,monkeypatch):
    # The identity-test inertia option is not a track run; refused before any output.
    c=resolve({**TOY,'quasistatic':dict(inertia='backward_euler')})
    with pytest.raises(ValueError,match='inertia none'):quasistatic.run_track(c,tmp_path/'inertia')
    assert not (tmp_path/'inertia').exists()
    # An error inside a relaxation (not a step rejection) is that seed's failure; the run directory stays complete.
    def broken(*args,**kwargs):raise ValueError('Equidistributed faces are not increasing')
    monkeypatch.setattr(quasistatic,'relax',broken)
    path,report=quasistatic.run_track(resolve(TOY),tmp_path/'broken')
    monkeypatch.undo()
    assert report['status']=='numerical_failure' and json.loads((path/'termination.json').read_text())['status']=='numerical_failure'
    assert json.loads((path/'relax.json').read_text())['seeds']==[dict(seed_radius_rsun=15.,status='failed',
                                                                      reason='ValueError: Equidistributed faces are not increasing')]
    # No loaded polytrope exists above M_BH/M* = 0.1015: every seed fails and the run is preserved as a failure.
    raw=copy.deepcopy(TOY);raw['replication'].update(initial_bh_mass_msun=20.,seed_radii_rsun=[15.,30.],landings_msun=[],end_bh_mass_msun=0.)
    (tmp_path/'heavy.yaml').write_text(yaml.safe_dump(raw))
    assert main(['run',str(tmp_path/'heavy.yaml'),'--output',str(tmp_path/'heavy')])==1
    path=directory(tmp_path/'heavy')
    report=json.loads((path/'termination.json').read_text());relax=json.loads((path/'relax.json').read_text())
    assert report['status']=='numerical_failure' and 'no seed relaxed' in report['reason']
    assert [x['status'] for x in relax['seeds']]==['seed_failed','seed_failed'] and not (path/'row0.npz').exists()

# ---------------------------------------------------------------- benchmarks/ball2011_track.yaml

def test_track_yaml_is_the_registered_c0_configuration(tmp_path,capsys):
    record=vbt.check_track_config(PLAN)
    assert record['sha256']==vbt.sha256(vbt.TRACK_CONFIG) and Path(record['path'])==ROOT/'benchmarks/ball2011_track.yaml'
    c=load(vbt.TRACK_CONFIG)
    assert c==vbt.build_config(PLAN,'C0',c['mesa']['library_path'],c['mesa']['mesa_dir'])
    assert c['replication']['plan']==PLAN['plan_id']=='ball2011_track_plan' and c['name']=='ball2011_C0_n400'
    assert c['mesa']['library_path']==str(ROOT/'backends/mesa/build_ball_abi7/libqstars_mesa.dylib')
    # Every registered run is a resolved quasi_static configuration (a fixed point of resolve).
    for run in vbt.RUNS:
        r=vbt.build_config(PLAN,run,'lib','mesa')
        assert resolve(r)==r and r['numerics']['mechanics']=='quasi_static' and r['replication']['run']==run,run
    # The driver refuses a file that is not its C0 mapping, naming what differs.
    text=vbt.TRACK_CONFIG.read_text()
    for old,new,key in (('zones: 400','zones: 200','numerics.zones'),('delta_target: 0.005','delta_target: 0.01','stepcontrol.delta_target'),
                        ('mesh_cycles: 0','mesh_cycles: 40','relaxation.mesh_cycles'),('tables: stars2004_switch ','tables: default ','opacity.tables'),
                        ('run: C0','run: C1','replication.run')):
        assert text.count(old)==1,old
        (tmp_path/'track.yaml').write_text(text.replace(old,new))
        with pytest.raises(vbt.Refused,match=key.replace('.',r'\.')):vbt.check_track_config(PLAN,tmp_path/'track.yaml')
    for broken in (text.replace('mechanics: quasi_static','mechanics: midpoint'),text.replace('numerics:\n','numerics: 5\nextra:\n'),text+'\n  - ['):
        (tmp_path/'track.yaml').write_text(broken)
        with pytest.raises(vbt.Refused,match='not a valid qstars configuration'):vbt.check_track_config(PLAN,tmp_path/'track.yaml')
    with pytest.raises(vbt.Refused,match='not a valid qstars configuration'):vbt.check_track_config(PLAN,tmp_path/'missing.yaml')
    # `qstars run` never executes a registered run: only the driver verifies the plan registration.
    capsys.readouterr()
    with pytest.raises(SystemExit) as exit:main(['run',str(vbt.TRACK_CONFIG),'--output',str(tmp_path/'never')])
    assert exit.value.code==2 and 'scripts/validate_ball_track.py' in capsys.readouterr().err and not (tmp_path/'never').exists()
    with pytest.raises(ValueError,match='registered run of ball2011_track_plan'):quasistatic.run_track(c,tmp_path/'never')
    assert not (tmp_path/'never').exists()

def test_registered_opacity_tables_on_the_compiled_bridge():
    """Opt-in (QSTARS_MESA_LIBRARY of an ABI 7 build, MESA_DIR): opacity.tables of the registered runs. No model is built."""
    library,directory=os.environ.get('QSTARS_MESA_LIBRARY'),os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Set QSTARS_MESA_LIBRARY and MESA_DIR for compiled MESA tests')
    config=lambda run:vbt.build_config(PLAN,run,library,directory)
    eos,opacity=quasistatic.physics(config('C0'))
    if eos.abi_version<7:pytest.skip('ABI 7 STARS opacity required')
    rho,T=np.array([1e-9,1e-7]),np.array([5e3,2e5])
    for run,kind,clamp in (('C0',quasistatic.StarsOpacity,True),('S4',quasistatic.StarsOpacity,False),('S1',None,None)):
        for build in (quasistatic.physics,vbt.physics):  # `qstars run` mapping and the registered driver's
            eos,opacity=build(config(run))
            if kind is None:assert opacity==eos.opacity and np.array_equal(opacity(rho,T),eos.opacity(rho,T))
            else:
                assert isinstance(opacity,kind) and opacity.clamp is clamp
                assert np.array_equal(opacity(rho,T),eos.opacity_stars(rho,T,clamp=clamp))

def test_registered_seed_blocker_and_the_grey_seed_surface_on_the_compiled_bridge():
    """Opt-in (QSTARS_MESA_LIBRARY of an ABI 7 build, MESA_DIR). The fail-before/pass-after test of plan
    amendment 3 on the registered C1 configuration (N = 200); the preserved registered runs are
    benchmarks/results/2026-10-14-ball-track/*_relax_C{1,0,2}_n*. No relaxation is run, nothing is written
    and nothing is compared with the published tables.

    Registered seed (T from the EOS at the polytrope (rho, p)) at 6000 Rsun, the only registered radius
    that builds: its last cell implies a grey-relation L_N more than ten times the seed-law luminosity,
    no photosphere exists there (Gamma_ph >= 1 at every admissible density) and every relaxation step is
    rejected at the predictor, whatever dt is. Grey seed surface (amendment 3): the same seed with the last
    cell on the grey relation has L_N = the law luminosity and relaxation steps are accepted at 6000 Rsun.
    At the appended 3464 Rsun (the geometric mean of the registered primary radii) the registered
    photosphere solve (bracket 'doubling') reports a false Gamma_ph >= 1 at the seed-law luminosity, so
    that seed cannot be built with amendment 3 alone; with the bracket scan of amendment 4 it builds and
    its relaxation steps are accepted. The registered 2000 and 1000 Rsun seeds have Gamma_ph >= 1 at their
    own law luminosity at every admissible density (true; confirmed by the scan) and 12000 Rsun lies below
    the 3000 K MESA floor under either surface and either bracket.
    """
    from qstars import photosphere
    from qstars.constants import YEAR
    from qstars.mass_transport import DomainError
    from qstars.numerics import StepRejected
    from qstars.transport import luminosities
    library,directory=os.environ.get('QSTARS_MESA_LIBRARY'),os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Set QSTARS_MESA_LIBRARY and MESA_DIR for compiled MESA tests')
    config=vbt.build_config(PLAN,'C1',library,directory)
    eos,opacity=vbt.physics(config)
    assert config['replication']['seed_surface']=='eos' and config['replication']['seed_radii_rsun']==[2000,6000,1000,12000]
    assert config['quasistatic']['photosphere_bracket']=='doubling'
    scan=copy.deepcopy(config);scan['quasistatic']['photosphere_bracket']='scan'  # amendment 4
    build=lambda radius,surface,c=config:quasistatic.seed(c,eos,opacity,200,5*MSUN,1e4*MSUN,radius*RSUN,surface=surface)
    grey_bc=copy.deepcopy(config);grey_bc['transport']['convection']='off'
    outer=lambda state:float(luminosities(state,eos,opacity,grey_bc,0.)[0][-1])
    # Before: the registered seed.
    registered,info=build(6000.,'eos')
    assert outer(registered)>10*info['law_luminosity']
    with pytest.raises(DomainError,match='Gamma_ph >= 1.*at every admissible'):
        photosphere.solve(outer(registered),registered.r[-1],registered.enclosed[-1],eos,opacity,bracket='scan')
    ledger.initialize(registered,eos,cavity=True)
    for dt in (1e-2,1e-5,1e-8):
        with pytest.raises(StepRejected,match='P3 predictor left the model domain') as rejected:
            quasistatic.bondi_step(registered,dt*YEAR,eos,opacity,config,mode='relax')
        assert rejected.value.domain_exit and 'Gamma_ph >= 1' in rejected.value.domain_message
    # After: the grey seed surface. Only the last temperature differs.
    grey,ginfo=build(6000.,'grey')
    for field in ('r','dm','cell_volume'):assert np.array_equal(getattr(grey,field),getattr(registered,field))
    assert np.array_equal(grey.T[:-1],registered.T[:-1]) and ginfo['surface_effective_temperature']<grey.T[-1]<registered.T[-1]
    assert outer(grey)==pytest.approx(info['law_luminosity'],rel=1e-10)
    for radius,c in ((6000.,config),(6000.,scan),(3464.,scan)):
        state,_=build(radius,'grey',c);ledger.initialize(state,eos,cavity=True)
        state,d=quasistatic.bondi_step(state,1e-2*YEAR,eos,opacity,c,mode='relax')
        assert d['nonlinear_residual']<=1e-12 and 0<d['photosphere_eddington_factor']<1
        state,d=quasistatic.bondi_step(state,1e-2*YEAR,eos,opacity,c,mode='relax')
        assert d['nonlinear_residual']<=1e-12 and 0<d['photosphere_eddington_factor']<1
    # Amendment 4's defect: with the registered bracket the 3464 Rsun grey seed reports a false Gamma_ph >= 1
    # (the same seed builds with the scan, and its photosphere has Gamma_ph < 1).
    with pytest.raises(DomainError,match='\\(Gamma_ph >= 1\\); no Ball surface exists'):build(3464.,'grey')
    seeded,sinfo=build(3464.,'grey',scan)
    rho=photosphere.solve(sinfo['law_luminosity'],3464.*RSUN,1e4*MSUN,eos,opacity,bracket='scan')
    assert 0<photosphere.evaluate(rho,sinfo['law_luminosity'],3464.*RSUN,1e4*MSUN,eos,opacity)['eddington_factor'].real<1
    # Radii no seed surface can build: the photosphere of the seed-law luminosity does not exist.
    for surface in quasistatic.SEED_SURFACES:
        for radius in (2000.,1000.):
            with pytest.raises(DomainError,match='Gamma_ph >= 1.*at every admissible'):build(radius,surface,scan)
            with pytest.raises(DomainError,match='\\(Gamma_ph >= 1\\); no Ball surface exists'):build(radius,surface)
        for c in (config,scan):
            with pytest.raises(ValueError,match='audited query rectangle'):build(12000.,surface,c)
