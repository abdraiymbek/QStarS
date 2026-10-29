"""Gate 6 Stage 5: the transactional quasi-static driver (quasistatic.run, relax), spec test I6.

I6: an analytic-cavity run stopped at step 20 and resumed to step 40 from
its checkpoint reproduces the uninterrupted run bitwise (every history row
and the end checkpoint: r, T, dm, v, cell_volume, radius_compensation,
face_luminosity, bh_mass, cavity_mass, photosphere_density, aux, ledger)
across base slices, split/merge events, rejected attempts and a landing;
the plan's restart points (first step at or above a mass, the step before a
landing) give segments that end bitwise on the uninterrupted run's
segment_end checkpoints. Also: the checkpoint round trip, relaxation with
the STARS-Q mesh cycles, and every driver end state.

Analytic cavity toy (as in test_quasistatic.py): 100 Msun, M_BH = 1 Msun,
200 Rsun, 128 zones, face boundary state, eta = 4e-4, ideal gas plus
radiation, electron scattering, ML1, Ball photosphere; mapped metric with
delta_target = 0.03 so that the run crosses every event within 40 steps.
The fully ionized EOS fails the Saha screens of this cool toy, so the test
policy reports them instead of enforcing them.
"""
import copy,csv,json,math,sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).parent))
from test_quasistatic import toy,EOS,OPACITY
from qstars import quasistatic,stepcontrol,validity
from qstars.constants import MSUN,YEAR

POLICY=dict(name='toy',report_only=validity.POLICIES['published_assumption']['report_only']
            +('ideal_eos_hydrogen_saha','ideal_eos_helium_saha'),replaced=('outer_face_knudsen',))
LANDING=1.2*MSUN
# A flow-Mach limit (validity.max_mach) that the relaxed start row meets (it is at rest: |r - rbar|/dt of
# its last 1e6 yr relaxation step) and every evolution step violates: an enforced screen at every retry.
MACH_LIMIT=1e-12

@pytest.fixture(scope='module')
def relaxed():
    s,c=toy(128,'face',eta=4e-4)
    s,record=quasistatic.relax(s,EOS,OPACITY,c,dict(mesh_cycles=0))
    assert record['status']=='converged'
    c=copy.deepcopy(c);c['stepcontrol']=dict(metric='mapped',delta_target=.03,dt_initial_years=100.)
    return quasistatic.zero_time(s,EOS),c

def rows(path):
    with open(Path(path)/'history.csv',newline='') as f:return list(csv.DictReader(f))

def events(path,kind):
    return [json.loads(x) for x in open(Path(path)/'events.jsonl') if json.loads(x)['event']==kind]

def same_state(a,b):
    for field in ('r','T','dm','v','cell_volume','radius_compensation','face_luminosity','reference_mass_fractions'):
        x,y=getattr(a,field),getattr(b,field)
        assert x.dtype==y.dtype and np.array_equal(x,y),field
    for field in ('time','bh_mass','cavity_mass','photosphere_density','steps','kinetic_tracked'):
        assert getattr(a,field)==getattr(b,field),field
    assert json.dumps(a.aux,sort_keys=True)==json.dumps(b.aux,sort_keys=True)
    assert json.dumps(a.ledger,sort_keys=True)==json.dumps(b.ledger,sort_keys=True)

def test_checkpoint_round_trip_is_exact(relaxed,tmp_path):
    s,c=relaxed
    t=s.copy();t.aux={'controller':dict(dt_next=1/3,hold=2,retries=1),'predictor':dict(time=[0.,.1],rest_rate=[1/7,2/7],cavity_mass=[1.,2.])}
    quasistatic.save_checkpoint(tmp_path/'a.npz',t,c)
    u,cc=quasistatic.load_checkpoint(tmp_path/'a.npz')
    same_state(t,u);assert cc==json.loads(json.dumps(c))
    assert not (tmp_path/'a.tmp').exists()
    t.photosphere_density=None;t.face_luminosity=None;quasistatic.save_checkpoint(tmp_path/'b.npz',t,c)
    u,_=quasistatic.load_checkpoint(tmp_path/'b.npz')
    assert u.photosphere_density is None and u.face_luminosity is None
    bad=t.copy();bad.cell_volume=bad.cell_volume*(1+1e-6);quasistatic.save_checkpoint(tmp_path/'c.npz',bad,c)
    with pytest.raises(ValueError,match='geometry'):quasistatic.load_checkpoint(tmp_path/'c.npz')
    z=quasistatic.zero_time(s,EOS)
    assert z.time==0 and z.steps==0 and z.aux=={} and z.ledger['initial_mass']==s.bh_mass+s.cavity_mass+s.dm.sum()
    assert 'cavity_swept_mass' in z.ledger and np.array_equal(z.r,s.r)

def test_restart_is_bitwise_across_slice_split_merge_rejection_and_landing(relaxed,tmp_path):
    """Spec I6: stop at step 20, resume to 40."""
    s,c=relaxed
    run=lambda state,path,**k:quasistatic.run(state,EOS,OPACITY,c,tmp_path/path,landings=[LANDING],policy=POLICY,**k)
    full,report=run(s,'full',stop_after_steps=40,restart_masses=[1.03*MSUN],restart_landing=LANDING,segment_steps=4)
    assert report['status']=='control_stop' and report['final_step']==40 and report['landed_msun']==[1.2]
    first,_=run(s,'first',stop_after_steps=20)
    state,config=quasistatic.load_checkpoint(tmp_path/'first'/'checkpoint.npz')
    assert state.steps==20 and config==json.loads(json.dumps(c))
    resumed,report=run(state,'resumed',stop_after_steps=20)
    assert report['restarted'] and report['final_step']==40
    same_state(full,resumed)
    end,_=quasistatic.load_checkpoint(tmp_path/'full'/'checkpoint.npz');again,_=quasistatic.load_checkpoint(tmp_path/'resumed'/'checkpoint.npz')
    same_state(end,again);same_state(end,full)
    # Every history row of the resumed segment equals the uninterrupted run's, bitwise (CSV repr).
    a={r['step']:r for r in rows(tmp_path/'full')};b=rows(tmp_path/'resumed')
    assert [r['step'] for r in b]==[str(k) for k in range(20,41)]
    for r in b[1:]:assert r==a[r['step']]
    start=b[0];assert start['event']=='restart'
    assert {k:v for k,v in start.items() if not k.startswith(('step_','metric_','event'))}=={k:v for k,v in a['20'].items() if not k.startswith(('step_','metric_','event'))}
    # The resumed segment crossed every event the plan names.
    seg=b[1:]
    assert all(float(r['step_slice_mass'])>0 for r in seg)
    assert any(int(r['step_splits'])>0 for r in seg) and any(int(r['step_merges'])>0 for r in seg)
    assert len({r['zones'] for r in seg})>1
    assert [r['event'] for r in seg].count('landing')==1
    assert len(events(tmp_path/'resumed','rejected'))>=1
    assert [e['step'] for e in events(tmp_path/'resumed','rejected')]==[e['step'] for e in events(tmp_path/'full','rejected') if e['step']>=20]
    landed=[r for r in seg if r['event']=='landing'][0]
    assert abs(float(landed['bh_mass_msun'])/1.2-1)<=1e-12 and float(landed['landed_target_msun'])==1.2
    assert (tmp_path/'full'/'landing_1.2msun.npz').exists() and (tmp_path/'full'/'profile_landing_1.2msun.npz').exists()
    # Conservation at every accepted step (plan conservation: mass 1e-14 per step, energy 5e-13 enforced per step).
    for r in a.values():
        if r['event']!='start':
            assert abs(float(r['step_mass_residual_change']))<=1e-14 and abs(float(r['step_energy_residual_change']))<=5e-13
    # Plan restart points: each segment from restart_<label>.npz ends bitwise on the uninterrupted segment_end.
    marks={e['label']:e['step'] for e in events(tmp_path/'full','restart_point')}
    assert marks=={'bh1.03msun':14,'before_landing_1.2msun':36}
    for label,step in marks.items():
        state,_=quasistatic.load_checkpoint(tmp_path/'full'/f'restart_{label}.npz')
        assert state.steps==step
        seg_end,_=run(state,f'segment_{label}',stop_after_steps=4)
        reference,_=quasistatic.load_checkpoint(tmp_path/'full'/f'segment_end_{step+4:07d}.npz')
        same_state(seg_end,reference)

def test_driver_end_states(relaxed,tmp_path):
    s,c=relaxed
    # completed: the end mass is a landing mass, hit exactly.
    _,r=quasistatic.run(s,EOS,OPACITY,c,tmp_path/'a',landings=[1.01*MSUN],end_bh_mass=1.01*MSUN,policy=POLICY)
    assert r['status']=='completed' and abs(r['final_bh_mass_msun']/1.01-1)<=1e-12 and r['landed_msun']==[1.01]
    # stop_rule: the dt floor (here any dt below 1e9 yr counts) after three consecutive steps.
    _,r=quasistatic.run(s,EOS,OPACITY,c,tmp_path/'b',policy=POLICY,
                        dt_floor=dict(dt=1e9*YEAR,steps=3,before_bh_mass=1e9*MSUN))
    assert r['status']=='stop_rule' and r['final_step']==3 and 'numerical blocker' in r['reason']
    _,r=quasistatic.run(s,EOS,OPACITY,c,tmp_path/'c',policy=POLICY,max_accepted_steps=2)
    assert r['status']=='stop_rule' and r['final_step']==2
    # model_domain_exit by retries: the start row passes every enforced screen, and every candidate step
    # violates one (a flow-Mach limit far below any step's |r - rbar|/dt; the relaxed start row is at rest).
    cm=copy.deepcopy(c);cm['validity']['max_mach']=MACH_LIMIT
    _,r=quasistatic.run(s,EOS,OPACITY,cm,tmp_path/'d',policy=POLICY)
    assert r['status']=='model_domain_exit' and r['final_step']==0 and r['rejections']==13 and r['start_row_screen_errors']==[]
    assert r['last_rejection']['kind']=='screen' and 'radial Mach' in r['last_rejection']['message'] and 'all 13' in r['reason']
    held,_=quasistatic.load_checkpoint(tmp_path/'d'/'checkpoint.npz')
    assert held.aux['driver']['screen_retries']==13==held.aux['controller']['retries']
    assert json.loads((tmp_path/'d'/'validity_warnings.json').read_text())['policy']=='toy'
    assert events(tmp_path/'d','start')[0]['enforced_screen_errors']==[]
    # termination_trigger: dt below the dynamical time.
    cc=copy.deepcopy(c);cc['stepcontrol']['dt_initial_years']=stepcontrol.dynamical_time(s)/YEAR/2
    _,r=quasistatic.run(s,EOS,OPACITY,cc,tmp_path/'e',policy=POLICY)
    assert r['status']=='termination_trigger' and 'dynamical time' in r['reason'] and r['attempts']==0
    for d in 'abcde':
        assert json.loads((tmp_path/d/'termination.json').read_text())['start_row_screen_errors']==[]
        assert json.loads((tmp_path/d/'termination.json').read_text())['status']
        assert (tmp_path/d/'checkpoint.npz').exists() and (tmp_path/d/'screen_violations.json').exists()
    summary=json.loads((tmp_path/'a'/'screen_violations.json').read_text())['screens']
    assert summary['radiation_knudsen_outer']['rows']>=1 and summary['radiation_knudsen_outer']['first_bh_mass_msun']==1.

def test_start_row_outside_the_enforced_domain_is_a_model_domain_exit_of_row_0(relaxed,tmp_path):
    """Final review item 4 (plan replication_mode.always_enforced): a start row that violates an enforced
    screen ends the run there, reported as a model_domain_exit of row 0; no step is attempted from it.
    Before the fix the violation was only logged in the start event and the run went on (here: to 13
    screen-rejected retries, reported as retries, never as the start row). production enforces the
    ideal-gas Saha screens, which this cool toy fails at t = 0."""
    s,c=relaxed
    _,r=quasistatic.run(s,EOS,OPACITY,c,tmp_path/'p',policy='production')
    assert r['status']=='model_domain_exit' and r['attempts']==0 and r['rejections']==0 and r['final_step']==0
    assert any('Saha' in e for e in r['start_row_screen_errors'])
    assert 'start row' in r['reason'] and 'model_domain_exit of row 0' in r['reason'] and r['last_rejection'] is None
    assert events(tmp_path/'p','start')[0]['enforced_screen_errors']==r['start_row_screen_errors']
    assert not events(tmp_path/'p','rejected') and not events(tmp_path/'p','accepted')
    h=rows(tmp_path/'p');assert [x['event'] for x in h]==['start'] and h[0]['screen_ideal_eos_helium_saha_violated']=='True'
    assert json.loads((tmp_path/'p'/'termination.json').read_text())['start_row_screen_errors']==r['start_row_screen_errors']
    assert json.loads((tmp_path/'p'/'validity_warnings.json').read_text())['policy']=='production'
    held,_=quasistatic.load_checkpoint(tmp_path/'p'/'checkpoint.npz')
    assert held.steps==0 and held.aux['controller']['retries']==0  # the controller never ran
    # A restart from a state outside the enforced domain ends the same way, without the row-0 label.
    again,_=quasistatic.load_checkpoint(tmp_path/'p'/'checkpoint.npz')
    _,r=quasistatic.run(again,EOS,OPACITY,c,tmp_path/'q',policy='production')
    assert r['status']=='model_domain_exit' and r['restarted'] and 'row 0' not in r['reason'] and r['start_row_screen_errors']

def test_retry_exhaustion_is_a_screen_exit_only_when_every_retry_was_screened(relaxed,tmp_path,monkeypatch):
    """One numerical rejection among the 13 consecutive retries makes the end the termination trigger,
    even though the last 12 (including the final one) were enforced-screen rejections."""
    s,c=relaxed
    c=copy.deepcopy(c);c['validity']['max_mach']=MACH_LIMIT  # every candidate step fails the flow-Mach screen
    step=quasistatic.bondi_step;calls=[]
    def first_fails(*a,**k):
        calls.append(1)
        if len(calls)==1:raise quasistatic.rejection('Quasi-static solve failed (test)')
        return step(*a,**k)
    monkeypatch.setattr(quasistatic,'bondi_step',first_fails)
    _,r=quasistatic.run(s,EOS,OPACITY,c,tmp_path/'m',policy=POLICY)
    assert r['status']=='termination_trigger' and r['rejections']==13 and r['last_rejection']['kind']=='screen'
    assert '12 of 13 by enforced screens' in r['reason']
    kinds=[e['kind'] for e in events(tmp_path/'m','rejected')];assert kinds==['numerical']+['screen']*12

def test_landing_just_above_the_black_hole_mass_lands_once_and_keeps_the_controlled_dt(relaxed,tmp_path,monkeypatch):
    """A landing mass a few 1e-9 above M_BH lands once (no second landing row, one restart point) and its
    forced dt does not seed the controller: no termination trigger (t_dyn of this toy is about 0.01 yr)."""
    s,c=relaxed
    target=s.bh_mass*(1+3e-9);tdyn=stepcontrol.dynamical_time(s)/YEAR
    def check(path):
        h=rows(tmp_path/path);landing=[r for r in h if r['event']=='landing']
        assert [r['step'] for r in landing]==['1'] and float(landing[0]['landed_target_msun'])==target/MSUN
        assert float(landing[0]['step_dt'])<1. and all(float(r['step_dt'])>=1. for r in h[2:])
        state,_=quasistatic.load_checkpoint(tmp_path/path/'checkpoint.npz')
        assert state.aux['driver']['landed']==[target] and state.aux['controller']['dt_next']/YEAR>tdyn
        marks=events(tmp_path/path,'restart_point');assert [e['step'] for e in marks]==[0]
        mark,_=quasistatic.load_checkpoint(tmp_path/path/f'restart_before_landing_{target/MSUN:g}msun.npz')
        assert mark.steps==0 and mark.aux['driver']['landed']==[] and mark.aux['controller']['retries']==0
    run=lambda path:quasistatic.run(s,EOS,OPACITY,c,tmp_path/path,landings=[target],policy=POLICY,stop_after_steps=3,
                                    restart_landing=target,segment_steps=2)
    _,r=run('direct')
    assert r['status']=='control_stop' and r['final_step']==3 and r['landed_msun']==[target/MSUN];check('direct')
    # The same target reached by an evolve step that passes it (the overshoot path) is landed once too.
    monkeypatch.setattr(stepcontrol,'landing_target',lambda *a:None)
    _,r=run('overshoot')
    assert r['status']=='control_stop' and r['final_step']==3;check('overshoot')
    monkeypatch.undo()
    # A landing mass at or below the starting M_BH that was never landed is an argument error.
    with pytest.raises(ValueError,match='never landed'):
        quasistatic.run(s,EOS,OPACITY,c,tmp_path/'below',landings=[s.bh_mass],policy=POLICY,stop_after_steps=1)

def test_registered_numerics_end_in_the_retry_trigger(relaxed,tmp_path):
    """The Stage 4 blocker through the driver: Lagrangian metric and base (the registered numerics);
    a base re-split jump is rejected at every dt until the 12-retry termination trigger fires."""
    s,c=relaxed
    cc=copy.deepcopy(c);cc['stepcontrol']=dict(dt_initial_years=100.)
    _,r=quasistatic.run(s,EOS,OPACITY,cc,tmp_path/'r',policy=POLICY,max_accepted_steps=400)
    assert r['status']=='termination_trigger' and '12 consecutive retries' in r['reason']
    rejected=events(tmp_path/'r','rejected')
    assert len(rejected)>=13 and len({e['step'] for e in rejected[-13:]})==1

def test_relaxation_mesh_cycles_reach_the_stars_q_fixed_point():
    """Two seeds of different radius relax to the same discrete model once the mesh is a fixed point
    of its own STARS-Q equidistribution (plan initial_model.uniqueness, here on a 48-zone cell toy)."""
    out=[]
    for radius in (150.,300.):
        s,c=toy(48,'cell',radius=radius)
        s,record=quasistatic.relax(s,EOS,OPACITY,c,dict(mesh_tolerance=1e-6))
        last=record['cycles'][-1]
        assert record['status']=='converged' and min(last['max_face_shift'],last['target_change'])<=1e-6 and len(record['cycles'])>=5
        assert abs(record['final']['luminosity_balance'])<=1e-6
        out.append(s)
    a,b=out
    for key in ('inner_radius','base_density','base_temperature','luminosity','radius','effective_temperature'):
        x,y=stepcontrol.observables(a,EOS,'cell')[key],stepcontrol.observables(b,EOS,'cell')[key]
        assert abs(math.log(x/y))<=1e-8,key
    assert abs(a.cavity_mass/b.cavity_mass-1)<=1e-8 and np.max(abs(np.log(a.dm/b.dm)))<=1e-5
    # Without the cycles each seed keeps its own Lagrangian mesh: the same observables differ far more.
    plain=[]
    for radius in (150.,300.):
        s,c=toy(48,'cell',radius=radius);plain.append(quasistatic.relax(s,EOS,OPACITY,c,dict(mesh_cycles=0))[0])
    assert abs(math.log(plain[0].r[0]/plain[1].r[0]))>=1e-5
