import copy
import numpy as np
import pytest
from qstars.eos import IdealGasRadiation,Composition
from qstars.constants import A,G,PI,MSUN,YEAR,C
from qstars.config import resolve
from qstars.evolution import physics,attempt_step,run
from qstars.structure import initialize,hydrostatic_residual
from qstars import ledger
from qstars.mass_transport import exchange,DomainError
from qstars.mesh import remap
from qstars.io import save_checkpoint,load_checkpoint
from qstars.engine import rates
from qstars.reference import ball_limit,corrected_limit,lane_emden,hassan_radius,hassan_critical_mass

@pytest.fixture
def model():
    config=resolve({'numerics':{'zones':16},'engine':{'model':'off'},
                    'boundary':{'thermal':'insulated'},'transport':{'convection':'off'}})
    eos,op=physics(config); state,_=initialize(config,eos); ledger.initialize(state,eos)
    return state,eos,op,config

def test_eos_first_law_and_derivatives():
    e=IdealGasRadiation();rho=np.geomspace(1e-9,1e-3,15);T=np.geomspace(1e5,1e7,15)
    q=e.evaluate(rho,T);h=1e-5
    qp=e.evaluate(rho*np.exp(h),T);qm=e.evaluate(rho*np.exp(-h),T)
    # T ds = du + P d(1/rho), independent finite differences.
    np.testing.assert_allclose(T*(qp['s']-qm['s'])/(2*h),(qp['u']-qm['u'])/(2*h)-q['p']/rho,rtol=1e-8)
    qp=e.evaluate(rho,T*np.exp(h));qm=e.evaluate(rho,T*np.exp(-h))
    np.testing.assert_allclose((qp['u']-qm['u'])/(2*h*T),q['cv'],rtol=1e-8)
    np.testing.assert_allclose(e.temperature(rho,q['u']),T,rtol=1e-12)
    np.testing.assert_allclose(e.temperature(rho,q['p'],'p'),T,rtol=1e-12)
    assert np.all((q['grad_ad']>=.25)&(q['grad_ad']<=.4))

def test_hydrostatic_and_energy(model):
    s,e,op,c=model
    assert max(abs(hydrostatic_residual(s,e,0)))<1e-9
    out,d=attempt_step(s,.001*YEAR,e,op,c)
    assert abs(ledger.residuals(out,e)['energy_residual_relative'])<1e-11
    assert max(abs(out.r/s.r-1))<1e-6

def test_n3_polytrope_does_not_create_artificial_convection(model):
    from qstars.transport import luminosities
    s,e,op,c=model;c['transport']['convection']='mlt_lossless'
    _,d=luminosities(s,e,op,c,0.)
    np.testing.assert_array_equal(d['convective'],0.)

def test_radiative_flux_exact_linear_T_fourth(model):
    from qstars.transport import luminosities
    from qstars.constants import A,C
    s,e,op,c=model
    slope=-1e11
    s.T=(1e26+slope*(s.rc-s.rc[0]))**.25
    _,d=luminosities(s,e,op,c,0.)
    rho=.5*(s.rho[:-1]+s.rho[1:])
    expected=-4*PI*s.r[1:-1]**2*A*C/(3*op.kappa*rho)*slope
    np.testing.assert_allclose(d['radiative'][1:-1],expected,rtol=1e-11)

def test_opacity_table_no_extrapolation(tmp_path):
    from qstars.opacity import OpacityTable
    p=tmp_path/'opacity.npz'
    np.savez(p,logrho=[-10.,-5.],logT=[5.,7.],logkappa=np.log10([[.34,.34],[.34,.34]]),source='controlled constant table')
    op=OpacityTable(p);assert np.isclose(float(op(1e-7,1e6)),.34)
    with pytest.raises(ValueError):op(1e-12,1e6)

def test_wind_energy_cap_and_subeddington(model):
    from qstars.winds import rate
    from qstars.engine import eddington
    s,e,op,c=model;c['wind']['model']='dotan2011'
    le=eddington(s.enclosed[-1],op.kappa)
    assert rate(s,e,op,c['wind'],.9*le)==0.
    md=rate(s,e,op,c['wind'],10*le)
    assert 0<md<=.9*10*le*s.r[-1]/(G*s.enclosed[-1])

def test_two_stream_equal_rates_do_not_cancel(model):
    s,e,op,c=model
    c['feeding'].update(rate_msun_per_year=.1,temperature_K=2e6,density_cgs=1e-7)
    engine=rates(s,e,c['engine'])
    out,flux=exchange(s,.001*YEAR,e,c,engine,.1*MSUN/YEAR)
    np.testing.assert_array_equal(out.dm,s.dm)
    assert out.ledger['supplied_mass']>0 and out.ledger['wind_mass']>0
    assert out.energies(e)['internal']!=s.energies(e)['internal']
    assert abs(ledger.residuals(out,e)['energy_residual_relative'])<1e-12

def test_accretion_mass_defect_and_energy(model):
    s,e,op,c=model;c['engine'].update(model='prescribed',rest_mass_rate_msun_per_year=1e-4)
    out,d=attempt_step(s,1e-4*YEAR,e,op,c)
    assert out.bh_mass>s.bh_mass
    assert out.ledger['accreted_rest_mass']>out.bh_mass-s.bh_mass
    assert abs(ledger.residuals(out,e)['mass_residual_relative'])<1e-14
    assert abs(ledger.residuals(out,e)['energy_residual_relative'])<1e-10
    assert np.isclose(out.ledger['central_generated'],out.ledger['mass_defect']*C**2,rtol=1e-13)

def test_retry_is_transactional(model):
    s,e,op,c=model;before=s.copy()
    c['feeding']['rate_msun_per_year']=1e9
    with pytest.raises(DomainError): attempt_step(s,YEAR,e,op,c)
    np.testing.assert_array_equal(s.dm,before.dm);np.testing.assert_array_equal(s.T,before.T)
    assert s.ledger==before.ledger
    c['feeding']['rate_msun_per_year']=0.
    a,_=attempt_step(s,.001*YEAR,e,op,c);b,_=attempt_step(before,.001*YEAR,e,op,c)
    np.testing.assert_array_equal(a.T,b.T)

def test_remap_conserves_mass_volume_total_energy(model):
    s,e,op,c=model
    s.v[1:]=np.linspace(0,1000,len(s.dm))
    fine=remap(s,np.repeat(s.dm/2,2),e)
    coarse=remap(fine,s.dm,e)
    assert np.isclose(fine.dm.sum(),s.dm.sum(),rtol=1e-14)
    assert np.isclose(fine.volume.sum(),s.volume.sum(),rtol=1e-14)
    scale=s.energies(e)['internal']
    for x in (fine,coarse): assert abs(x.energies(e)['total']-s.energies(e)['total'])/scale<1e-13

def test_checkpoint_equivalence(model,tmp_path):
    s,e,op,c=model;s,_=attempt_step(s,.001*YEAR,e,op,c)
    p=tmp_path/'checkpoint.npz';save_checkpoint(p,s,c); restored,cc=load_checkpoint(p)
    a,_=attempt_step(s,.001*YEAR,e,op,c);b,_=attempt_step(restored,.001*YEAR,e,op,cc)
    np.testing.assert_array_equal(a.r,b.r);np.testing.assert_array_equal(a.T,b.T)
    assert a.ledger==b.ledger

def test_restart_full_driver(model,tmp_path):
    s,e,op,c=model;c['evolution']['end_time_years']=.005
    p,report=run(c,tmp_path,stop_after_steps=2);restarted,cc=load_checkpoint(p/'checkpoint.npz')
    p2,r2=run(cc,tmp_path,restart=restarted);p3,r3=run(c,tmp_path)
    a,_=load_checkpoint(p2/'checkpoint.npz');b,_=load_checkpoint(p3/'checkpoint.npz')
    np.testing.assert_array_equal(a.T,b.T); assert a.ledger==b.ledger
    assert r2['status']==r3['status']=='completed_requested_interval'

def test_lane_emden_and_published_limits():
    s=lane_emden(1);assert abs(s.t[-1]-np.pi)<1e-8;assert abs(s.y[1,-1]-np.pi)<1e-8
    assert abs(ball_limit()['bh_fraction']-.102)<.0005
    assert abs(corrected_limit()['bh_fraction']-.017)<.0001

def test_hassan_cubic_and_domain():
    rho,cs=1e-8,1e7;crit=hassan_critical_mass(rho,cs)
    for ratio in (.001,.1,.9,.9999):
        mass=crit*ratio;r=hassan_radius(mass,rho,cs)
        assert np.isclose(r,2*G*(mass+8*PI/3*rho*r**3)/cs**2,rtol=1e-10)
    with pytest.raises(ValueError): hassan_radius(crit*1.01,rho,cs)

def test_domain_exit_and_preserved_run(tmp_path):
    c=resolve({'name':'unsupported_cool_model','validity':{'min_temperature_K':1e7}})
    p,r=run(c,tmp_path)
    assert r['status']=='model_domain_exit'
    assert (p/'termination.json').exists() and (p/'checkpoint.npz').exists()

@pytest.mark.parametrize('raw',[{'typo':1},{'numerics':{'zones':0}},{'composition':{'X':2}},
    {'engine':{'efficiency':1}},{'feeding':{'rate_msun_per_year':-1}},
    {'feeding':{'rate_msun_per_year':{'kind':'table','time_years':[1,0],'values':[1,2]}}}])
def test_bad_config(raw):
    with pytest.raises(ValueError): physics(resolve(raw))
