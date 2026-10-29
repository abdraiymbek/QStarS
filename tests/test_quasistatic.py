"""Gate 6 Stage 4: the quasi-static Ball step (quasistatic.py), step control and screens.

I1 well-balanced relaxed equilibria; I2 quasi-static vs backward-Euler inertia
scales as (t_dyn/dt)^2; I3 with a frozen cavity and inertia the step equals
numerics.thermal_mechanical_step(backward_euler); I5 an analytic-cavity
evolution over 1e4 yr closes mass, energy, engine and geometry; U14 the
screen policies; plus the declared Jacobian pattern, the predictor-corrector
passes, landing, and the endpoint momentum form against the secant memory.
Stage 4 review: the fixed-mass base (fallback F2) and the mapped controller
metric that draft amendment 1 proposes, against the registered numerics;
domain-exit records and step rejections; the sink flow work booked as
inner-boundary work; the per-step conservation budget; screen gaps.
Stage 6 review: an EOS temperature inversion or photosphere solve that does
not converge (RuntimeError) is a step rejection, retried at dt/2.
'Analytic cavity': fully ionized ideal gas plus radiation, electron
scattering, ML1 convection, Ball cavity and photosphere boundaries, the
thesis eq 3.15 engine with a reduced coupling (100 Msun, M_BH = 1 Msun) so
that the envelope stays sub-Eddington, and the Ball 2011 row-0 geometry
(1e4 Msun, M_BH = 5 Msun, eta = 0.1) relaxed with the face boundary state.
"""
import copy,json,math
import numpy as np
import pytest
from qstars import cavity,ledger,numerics,photosphere,quasistatic,stepcontrol,validity
from qstars.config import resolve
from qstars.constants import A,C,G,H,KB,ME,MU,PI,MSUN,RSUN,YEAR
from qstars.eos import IdealGasRadiation,Composition
from qstars.engine import rates
from qstars.ledger import residuals
from qstars.mass_transport import DomainError
from qstars.nonlinear import pattern_jacobian
from qstars.numerics import StepRejected
from qstars.opacity import ElectronScattering

EOS=IdealGasRadiation(Composition(.7,.3,0.))
OPACITY=ElectronScattering(EOS.composition)

def toy(zones,mode,bh=1.,total=100.,radius=200.,eta=1.5e-4):
    c=resolve(dict(composition=dict(X=.7,Y=.3,Z=0.),numerics=dict(zones=zones)))
    c['engine'].update(model='ball2012_eq3_15',coupling=eta)
    c['transport']['convection']='mlt_ml1';c['boundary']['thermal']='stars_photosphere'
    c['validity']['min_temperature_K']=3000.
    c['cavity']=dict(boundary_state=mode)
    s,_=quasistatic.seed(c,EOS,OPACITY,zones,bh*MSUN,total*MSUN,radius*RSUN)
    ledger.initialize(s,EOS,cavity=True)
    return s,c

def relax(s,c,dt=.01*YEAR,tol=1e-10,max_steps=200,direct=False):
    """Spec 6 relaxation (M_BH frozen): dt doubles to 1e8 yr, halves on a rejection; converged
    when dt >= 1e6 yr and two consecutive steps change ln T and ln r by <= tol. direct: the
    change is also measured cell by cell between consecutive states (a fixed-mass base moves
    its faces between steps until its regrid reaches a fixed point)."""
    quiet=0
    for _ in range(max_steps):
        try:new,d=quasistatic.bondi_step(s,dt,EOS,OPACITY,c,mode='relax')
        except StepRejected:
            dt*=.5;continue
        change=max(d['lagrangian_max_dlnT'],d['lagrangian_max_dlnr'])
        if direct:
            change=max(change,np.max(abs(np.log(new.r/s.r))),np.max(abs(np.log(new.T/s.T)))) if len(new.dm)==len(s.dm) else np.inf
        s=new;quiet=quiet+1 if dt>=1e6*YEAR and change<=tol else 0
        if quiet>=2:return s,d
        dt=min(2*dt,1e8*YEAR)
    raise RuntimeError('relaxation did not converge')

@pytest.fixture(scope='module')
def stable_cell():
    s,c=toy(64,'cell');return relax(s,c)+(c,)

@pytest.fixture(scope='module')
def stable_face():
    s,c=toy(128,'face');return relax(s,c)+(c,)

@pytest.fixture(scope='module')
def ball_face():
    s,c=toy(128,'face',bh=5.,total=1e4,radius=3000.,eta=.1);return relax(s,c)+(c,)

@pytest.fixture(scope='module')
def ball_fixed(ball_face):
    """The Ball row-0 toy relaxed again on the fixed-mass base (fallback F2, 5% of the envelope)."""
    s,_,c=ball_face
    c=copy.deepcopy(c);c['quasistatic']=dict(base='fixed_mass')
    return relax(s,c,dt=1e6*YEAR,direct=True)+(c,)

def zero_time(s):
    """t = 0 at the relaxed model: ledger re-initialised, driver state cleared (spec 6.5)."""
    s=s.copy();s.time=0.;s.steps=0;s.aux={};ledger.initialize(s,EOS,cavity=True)
    return s

# ---------------------------------------------------------------- legacy screens (validity.check at HEAD, verbatim)

def legacy_check(state,eos,opacity,config):
    errors=[]
    q=eos.evaluate(state.rho,state.T)
    eos_defect=float(max(eos.consistency(state.rho,state.T))) if hasattr(eos,'consistency') else 0.
    if eos_defect>.01:errors.append('MESA first-law defect exceeds 1% audit limit')
    gas_defect=float(np.max(eos.gas_pressure_limit(state.rho,state.T,q)['screened_relative_error'])) if hasattr(eos,'gas_pressure_limit') else 0.
    if gas_defect>.01:errors.append('MESA gas pressure differs by more than 1% from the dilute classical limit')
    conduction=float(np.max(q.get('conductive_fraction',0.)))
    if conduction>.01:errors.append('Omitted electron conduction exceeds 1% of diffusive heat flux')
    v=config["validity"]
    if min(state.T)<v["min_temperature_K"]: errors.append("EOS temperature below configured domain")
    if max(state.T)>v["max_temperature_K"]: errors.append("EOS temperature above configured domain")
    mach=float(max(abs(.5*(state.v[:-1]+state.v[1:]))/q["cs"]))
    if mach>v["max_mach"]: errors.append("Quasi-hydrostatic domain exceeded: radial Mach number")
    compactness=float(max(2*G*state.enclosed/(state.r*C**2)))
    if compactness>v["max_compactness"]: errors.append("Newtonian compactness limit exceeded")
    ne=state.rho/MU*(eos.composition.X+eos.composition.Y/2+eos.composition.Z/2)
    phase=(2*PI*ME*KB*state.T/H**2)**1.5
    eV=1.602176634e-12
    hydrogen_ratio=phase/ne*np.exp(-13.5984*eV/(KB*state.T))
    helium_ratio=phase/ne*np.exp(-54.4178*eV/(KB*state.T))
    ideal=config['eos']=='ideal_gas_radiation'
    if ideal and eos.composition.X>0 and min(hydrogen_ratio/(1+hydrogen_ratio))<.99:
        errors.append('Fully ionized EOS fails hydrogen Saha screening')
    if ideal and eos.composition.Y>0 and min(helium_ratio/(1+helium_ratio))<.99:
        errors.append('Fully ionized EOS fails He II/III Saha screening')
    degeneracy=ne*(H/np.sqrt(2*PI*ME*KB*state.T))**3
    if ideal and max(degeneracy)>.1: errors.append("Classical EOS fails electron nondegeneracy check")
    geometric=4*PI/3*np.diff(state.r)*(state.r[1:]**2+state.r[1:]*state.r[:-1]+state.r[:-1]**2)
    volume_error=float(max(abs(geometric/state.volume-1)))
    if volume_error>1e-8:raise StepRejected('Stored cell volumes disagree with radial geometry')
    kap=opacity(state.rho,state.T)
    tau=kap*state.rho*np.diff(state.r)
    # Diffusion concerns the physical gradient scale, not the numerical cell
    # width. Requiring tau_cell>1 incorrectly rejects convergent refinement.
    rho_face=.5*(state.rho[:-1]+state.rho[1:]);kap_face=.5*(kap[:-1]+kap[1:])
    kn=4*abs(np.diff(np.log(state.T)))/(np.diff(state.rc)*rho_face*kap_face)
    if v["require_optically_thick"] and max(kn)>.1:
        errors.append('Radiation gradient scale is outside the diffusion limit (Kn>0.1)')
    from qstars.transport import luminosities
    total_luminosity,transport=luminosities(state,eos,opacity,config,0.)
    significant=abs(transport['convective'][1:-1])>.01*(abs(transport['radiative'][1:-1])+abs(transport['convective'][1:-1]))
    min_pe=float(min(transport['peclet'][significant])) if np.any(significant) else None
    if min_pe is not None and min_pe<100:
        errors.append('Lossless MLT outside efficient-convection domain (significant flux with Pe<100)')
    if max(transport['convective_mach'],default=0)>1:
        errors.append('Lossless MLT predicts supersonic convection')
    from qstars.constants import A
    outer_kn=3*abs(total_luminosity[-1])/(4*PI*state.r[-1]**2*C*A*state.T[-1]**4)
    if v['require_optically_thick'] and outer_kn>.1:
        errors.append('Outer radiation flux is outside the LTE diffusion limit')
    free_stream=4*PI*state.r[1:-1]**2*C*A*.5*(state.T[:-1]**4+state.T[1:]**4)
    if np.any(abs(transport['radiative'][1:-1])>free_stream):
        errors.append('Radiative diffusion exceeds free-streaming flux')
    conservation=residuals(state,eos)
    tol=config["numerics"]["conservation_tolerance"]
    if max(abs(conservation["mass_residual_relative"]),abs(conservation["energy_residual_relative"]))>tol:
        raise StepRejected("Conservation tolerance exceeded")
    return errors,dict(max_mach=mach,max_compactness=compactness,min_zone_optical_depth=float(min(tau)),
                       max_eos_first_law_defect=eos_defect,
                       max_eos_dilute_gas_pressure_defect=gas_defect,
                       max_omitted_conductive_fraction=conduction,
                       max_radiation_knudsen=float(max(kn)),
                       outer_radiation_knudsen=float(outer_kn),
                       minimum_significant_convective_peclet=min_pe,
                       max_convective_mach=float(max(transport['convective_mach'],default=0)),
                       ideal_eos_helium_screen=float(min(helium_ratio/(1+helium_ratio))),
                       max_volume_geometry_error=volume_error,
                       max_temperature_jump=float(max(abs(np.diff(np.log(state.T)))))),conservation


# ---------------------------------------------------------------- settings, pattern, relaxation

def test_settings_and_mode_validation(stable_cell):
    s,_,c=stable_cell
    for bad in (dict(quasistatic=dict(unknown=1)),dict(quasistatic=dict(momentum='midpoint')),
                dict(quasistatic=dict(inertia='midpoint')),dict(cavity=dict(boundary_state='edge')),
                dict(cavity=dict(profile='r-2')),dict(quasistatic=dict(base='eulerian')),
                dict(quasistatic=dict(base_mass_fraction=.5)),dict(quasistatic=dict(base_mass_fraction='a')),
                dict(quasistatic=dict(base_mass_fraction=True)),dict(cavity=dict(profile=['r-3/2'])),
                dict(cavity=dict(profile={'r-3/2':1}))):
        with pytest.raises(ValueError):quasistatic.settings({**c,**bad})
    assert quasistatic.settings(c)['quasistatic']['momentum']=='endpoint'
    with pytest.raises(ValueError,match='frozen cavity'):
        quasistatic.bondi_step(s,YEAR,EOS,OPACITY,{**c,'cavity':dict(active=False)},mode='evolve')
    with pytest.raises(ValueError,match='target'):quasistatic.bondi_step(s,YEAR,EOS,OPACITY,c,mode='land')
    with pytest.raises(ValueError):quasistatic.bondi_step(s,YEAR,EOS,OPACITY,c,mode='walk')
    t=s.copy();t.kinetic_tracked=True
    with pytest.raises(ValueError,match='kinetic'):quasistatic.bondi_step(t,YEAR,EOS,OPACITY,c,mode='relax')

@pytest.mark.parametrize('momentum',['secant','endpoint'])
@pytest.mark.parametrize('mode',['evolve','relax','land'])
def test_declared_pattern_matches_dense_complex_step_jacobian(momentum,mode):
    s,c=toy(10,'face');c['quasistatic']=dict(momentum=momentum)
    S=quasistatic.settings(c);dt=YEAR
    target=s.bh_mass*(1+1e-6) if mode=='land' else None
    pred=quasistatic.predict(s,dt,EOS,c,mode,target,S)
    P=quasistatic.Problem(s,s.copy(),dt,EOS,OPACITY,c,mode,target,pred,S)
    x=P.x0+np.random.default_rng(6).normal(0,1e-4,len(P.x0))
    J=pattern_jacobian(P.evaluate,x,P.pattern,strict=True).toarray()
    dense=np.empty_like(J)
    for k in range(len(x)):
        xx=x.astype(complex);xx[k]+=1e-25j;dense[:,k]=P.evaluate(xx).imag/1e-25
    assert np.max(abs(J-dense))<=1e-13*np.max(abs(dense))
    assert P.x0.shape==(3*len(s.dm)+3,)

def test_relaxed_models_are_thermal_equilibria(stable_cell,stable_face,ball_face):
    for s,d,c in (stable_cell,stable_face,ball_face):
        assert abs(d['luminosity']/d['bh_luminosity']-1)<=1e-6  # spec 6.3: L* = L_BH
        assert abs(d['inner_radius']/d['bondi_radius']-1)<=1e-11  # piston: r0 = r_B (b = 1)
        assert abs(s.cavity_mass/d['cavity_mass_law']-1)<=1e-10  # thesis 3.7
        assert 0<d['photosphere_eddington_factor']<1 and d['stability_ratio']>0
        assert np.all(np.diff(s.r)>0) and s.photosphere_density>0
    # The Ball 2011 row-0 geometry (paper Table 1: L 1.33e42, rho0 8.71e-5) with this analytic microphysics.
    s,d,_=ball_face
    assert d['luminosity']==pytest.approx(1.33e42,rel=.05) and d['base_density']==pytest.approx(8.71e-5,rel=.15)

# ---------------------------------------------------------------- I1

@pytest.mark.parametrize('name,momentum',[('stable_cell','endpoint'),('stable_cell','secant'),
                                          ('stable_face','endpoint'),('ball_face','endpoint')])
def test_relaxed_equilibrium_is_well_balanced(request,name,momentum):
    s,_,c=request.getfixturevalue(name)
    c=copy.deepcopy(c);c['quasistatic']=dict(momentum=momentum)
    new,d=quasistatic.bondi_step(s,1e6*YEAR,EOS,OPACITY,c,mode='relax')
    assert len(new.dm)==len(s.dm) and d['splits']==d['merges']==0
    assert np.max(abs(np.log(new.T/s.T)))<=1e-10 and np.max(abs(np.log(new.r/s.r)))<=1e-10
    assert abs(new.cavity_mass/s.cavity_mass-1)<=1e-10 and abs(new.photosphere_density/s.photosphere_density-1)<=1e-10
    assert np.max(abs(new.face_luminosity/s.face_luminosity-1))<=1e-10
    assert new.bh_mass==s.bh_mass and abs(d['delta_fraction'])<=1e-10

# ---------------------------------------------------------------- momentum forms

def test_endpoint_momentum_removes_the_secant_memory(stable_cell):
    """An accepted secant state is off static balance by ~2 dr/r; in this envelope
    (3 Gamma1 - 4 = 0.3) the next step moves it again by most of that, whatever dt.
    The endpoint form leaves no such memory."""
    s,_,c=stable_cell
    x=(np.cumsum(s.dm)-.5*s.dm)/s.dm.sum()
    hot=s.copy();hot.T=s.T*(1+.02*np.exp(-((x-.5)/.15)**2));hot.aux={}
    ratio={}
    for form in ('endpoint','secant'):
        cc=copy.deepcopy(c);cc['quasistatic']=dict(momentum=form)
        a,_=quasistatic.bondi_step(hot,1*YEAR,EOS,OPACITY,cc,mode='relax')
        b,_=quasistatic.bondi_step(a,1e-6*YEAR,EOS,OPACITY,cc,mode='relax')
        assert len(a.r)==len(hot.r)==len(b.r)
        ratio[form]=np.max(abs(np.log(b.r/a.r)))/np.max(abs(np.log(a.r/hot.r)))
    assert ratio['endpoint']<=1e-5 and ratio['secant']>=.5

# ---------------------------------------------------------------- I2, I3

def frozen(s,c,inertia):
    """Frozen cavity and mesh, a fixed photospheric pressure, the engine as a constant source."""
    ph=photosphere.evaluate(s.photosphere_density,s.face_luminosity[-1],s.r[-1],s.enclosed[-1],EOS,OPACITY)
    c=copy.deepcopy(c);c['engine'].update(model='prescribed',rest_mass_rate_msun_per_year=
        float(s.face_luminosity[0])/(c['engine']['efficiency']*C**2)*YEAR/MSUN)
    c['boundary'].update(thermal='grey',external_optical_depth=2/3,bath_temperature_K=0.,pressure_cgs=float(ph['pressure'].real))
    c['cavity']=dict(active=False);c['quasistatic']=dict(remesh=False,inertia='backward_euler' if inertia else 'none')
    return c

def perturbed(s,amplitude,inertia):
    x=np.linspace(0,1,len(s.dm));t=s.copy();t.T=s.T*(1+amplitude*np.sin(3*np.pi*x));t.aux={}
    t.kinetic_tracked=inertia;t.v=np.zeros_like(s.v);ledger.initialize(t,EOS,cavity=True)
    return t

def test_inertia_difference_scales_as_dynamical_time_squared(stable_cell):
    s,_,c=stable_cell
    tdyn=stepcontrol.dynamical_time(s);ratios=[]
    for k in (10,30,100):
        a,da=quasistatic.bondi_step(perturbed(s,.01,False),k*tdyn,EOS,OPACITY,frozen(s,c,False),mode='relax')
        b,db=quasistatic.bondi_step(perturbed(s,.01,True),k*tdyn,EOS,OPACITY,frozen(s,c,True),mode='relax')
        ya,yb=np.log(a.r/s.r),np.log(b.r/s.r)
        ratios.append(np.max(abs(yb-ya))/np.max(abs(ya)))
        assert abs(db['energy_residual_change'])<=1e-12 and db['numerical_kinetic_heating_erg']>0
        assert abs(da['energy_residual_change'])<=1e-12
    slopes=np.diff(np.log(ratios))/np.diff(np.log([10,30,100]))
    assert np.all(abs(slopes+2)<=.2),slopes

def test_frozen_cavity_inertia_step_equals_thermal_mechanical_step(stable_cell):
    s,_,c=stable_cell
    tdyn=stepcontrol.dynamical_time(s)
    qc=frozen(s,c,True);qc['quasistatic']['momentum']='secant'
    tc=copy.deepcopy(qc);tc['numerics'].update(mechanics='backward_euler',nonlinear_solver='newton_flux',
                                              solver_tolerance=1e-12,max_fractional_change=10.)
    old=perturbed(s,.01,True)
    reference,_=numerics.thermal_mechanical_step(old,3*tdyn,EOS,OPACITY,tc,rates(old,EOS,tc['engine'])['luminosity'])
    mine,_=quasistatic.bondi_step(old,3*tdyn,EOS,OPACITY,qc,mode='relax')
    for field in ('r','T','face_luminosity','cell_volume'):
        assert np.max(abs(getattr(mine,field)/getattr(reference,field)-1))<=1e-9,field
    assert np.max(abs(mine.v-reference.v))<=1e-9*np.max(abs(reference.v))

# ---------------------------------------------------------------- I5, landing, passes

@pytest.mark.parametrize('name',['stable_cell','stable_face'])
def test_analytic_cavity_run_conserves_mass_energy_and_geometry(request,name):
    s,_,c=request.getfixturevalue(name)
    s=zero_time(s);S=stepcontrol.settings(c);mode=c['cavity']['boundary_state']
    dt=stepcontrol.start(s.aux,S);before=stepcontrol.observables(s,EOS,mode);accepted=0
    while s.time<1e4*YEAR*(1-1e-12):
        dt=min(dt,1e4*YEAR-s.time)
        new,d=quasistatic.bondi_step(s,dt,EOS,OPACITY,c,mode='evolve')
        after=stepcontrol.observables(new,EOS,mode);m=stepcontrol.metrics(before,after,d,S)
        assert not stepcontrol.rejected(m['delta_actual'],S)
        geometric=4*PI/3*np.diff(new.r)*(new.r[1:]**2+new.r[1:]*new.r[:-1]+new.r[:-1]**2)
        assert np.max(abs(geometric/new.volume-1))<=1e-8
        assert abs(d['energy_residual_change'])<=1e-10 and abs(d['mass_residual_change'])<=1e-14
        l=new.ledger
        assert abs(l['mass_defect']*C**2-l['central_generated'])<=1e-12*l['central_generated']
        assert abs(new.cavity_mass/d['cavity_mass_law']-1)<=1e-10 and abs(d['delta_fraction'])<=.05
        assert new.bh_mass>s.bh_mass and d['passes']<=quasistatic.DEFAULTS['quasistatic']['max_passes']
        json.dumps(new.aux)
        s,before=new,after;accepted+=1
        dt=stepcontrol.accept(s.aux,d['dt'],m['delta_actual'],d['newton_iterations'],S)
    r=residuals(s,EOS)
    assert abs(r['mass_residual_relative'])<=1e-14*accepted and abs(r['energy_residual_relative'])<=1e-10
    assert accepted>=5 and s.ledger['accreted_rest_mass']>0 and s.ledger['cavity_swept_mass']>0
    # Landing: dt is solved so that M_BH hits the target exactly.
    target=s.bh_mass*(1+1e-4)
    landed,d=quasistatic.bondi_step(s,dt,EOS,OPACITY,c,mode='land',target_mass=target)
    assert abs(landed.bh_mass/target-1)<=1e-12 and d['dt']>0 and landed.time==s.time+d['dt']
    assert d['law_bh_growth_rate']*d['dt']==pytest.approx(target-s.bh_mass,rel=1e-9)
    assert abs(d['energy_residual_change'])<=1e-10 and abs(d['mass_residual_change'])<=1e-14

def test_predictor_corrector_passes_recover_a_bad_prediction(stable_cell):
    s,_,c=stable_cell
    s=zero_time(s)
    s,_=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')
    s,good=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')
    # Mispredict the swept mass by 20% of the base cell (beyond the 5% pass threshold).
    extra=.2*s.dm[0]/(10*YEAR*s.aux['predictor']['rest_rate'][-1])
    bad=s.copy();bad.aux['predictor']['rest_rate']=[(1+extra)*x for x in s.aux['predictor']['rest_rate']]
    new,d=quasistatic.bondi_step(bad,10*YEAR,EOS,OPACITY,c,mode='evolve')
    assert d['passes']>=2 and abs(d['delta_fraction'])<=.05
    assert abs(d['mass_residual_change'])<=1e-14 and abs(d['energy_residual_change'])<=1e-10
    ref,_=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')
    assert new.bh_mass==pytest.approx(ref.bh_mass,rel=1e-9)

# ---------------------------------------------------------------- step control

def test_step_controller_rules():
    S=stepcontrol.settings({});aux={}
    assert stepcontrol.start(aux,S)==YEAR and aux['controller']==dict(dt_next=YEAR,hold=0,retries=0)
    assert stepcontrol.accept(aux,100.,.0045,3,S)==pytest.approx(100.)
    assert stepcontrol.accept(aux,100.,1e-9,3,S)==pytest.approx(150.)
    assert stepcontrol.accept(aux,100.,1.,3,S)==pytest.approx(50.)
    assert stepcontrol.accept(aux,100.,1e-9,9,S)==pytest.approx(120.)
    assert stepcontrol.accept(aux,1e12*YEAR,1e-9,3,S)==S['dt_max_years']*YEAR
    assert not stepcontrol.rejected(.01,S) and stepcontrol.rejected(.0101,S)
    assert stepcontrol.reject(aux,100.,S)==50. and aux['controller']['hold']==3
    assert stepcontrol.accept(aux,50.,1e-9,3,S)==50.  # no growth for 3 steps after a rejection
    assert aux['controller']['hold']==2 and aux['controller']['retries']==0
    for _ in range(12):assert stepcontrol.reject(aux,1.,S) is not None
    assert stepcontrol.reject(aux,1.,S) is None
    json.dumps(aux)
    # A nonfinite metric (here a state without solved face luminosities) is an error, not a pass.
    obs=dict(inner_radius=1.,base_density=1.,base_temperature=1.,bh_mass=1.,cavity_mass=0.,central_mass=1.,
             luminosity=1.,radius=1.,effective_temperature=1.)
    diag=dict(lagrangian_max_dlnT=0.,lagrangian_max_dlnr=0.)
    with pytest.raises(ValueError,match='Nonfinite'):
        stepcontrol.metrics({**obs,'luminosity':np.nan,'effective_temperature':np.nan},obs,diag,S)
    assert stepcontrol.metrics(obs,obs,diag,S)['delta_actual']==0.
    with pytest.raises(ValueError,match='finite'):stepcontrol.accept(aux,100.,float('nan'),3,S)
    targets=[100*MSUN,200*MSUN]
    assert stepcontrol.landing_target(99*MSUN,MSUN/YEAR,2*YEAR,targets)==100*MSUN
    assert stepcontrol.landing_target(99*MSUN,MSUN/YEAR,.5*YEAR,targets) is None
    assert stepcontrol.landing_target(100*MSUN,MSUN/YEAR,200*YEAR,targets[1:])==200*MSUN  # 100 landed: not pending
    assert stepcontrol.landing_target(100*MSUN,MSUN/YEAR,200*YEAR,[]) is None
    with pytest.raises(ValueError):stepcontrol.settings({'stepcontrol':{'metric':'eulerian'}})
    with pytest.raises(ValueError):stepcontrol.settings({'stepcontrol':{'delta':1}})

def test_step_metrics_and_dynamical_time(stable_cell):
    s,_,c=stable_cell
    s=zero_time(s);new,d=quasistatic.bondi_step(s,100*YEAR,EOS,OPACITY,c,mode='evolve')
    before,after=stepcontrol.observables(s,EOS,'cell'),stepcontrol.observables(new,EOS,'cell')
    for metric in stepcontrol.METRICS:
        m=stepcontrol.metrics(before,after,d,stepcontrol.settings({'stepcontrol':{'metric':metric}}))
        assert m['profile_dlnr']==d[metric+'_max_dlnr'] and m['delta_actual']==max(v for k,v in m.items() if k!='delta_actual')
        assert m['dln_bh_mass']==pytest.approx(math.log(new.bh_mass/s.bh_mass),rel=1e-12)
    assert stepcontrol.dynamical_time(s)==pytest.approx(math.sqrt(s.r[-1]**3/(G*s.enclosed[-1])),rel=1e-15)

# ---------------------------------------------------------------- U14 screens

def test_policy_none_is_the_legacy_screen_set(stable_cell,ball_face):
    for s,_,c in (stable_cell,ball_face):
        for config in (c,{**c,'validity':{**c['validity'],'require_optically_thick':False}}):
            new=validity.check(s,EOS,OPACITY,config);old=legacy_check(s,EOS,OPACITY,config)
            assert new[0]==old[0] and new[2]==old[2]
            assert json.dumps(new[1],sort_keys=True)==json.dumps(old[1],sort_keys=True)
    with pytest.raises(ValueError,match='policy'):validity.check(s,EOS,OPACITY,c,face_luminosity=s.face_luminosity)

def test_screens_report_only_in_replication_and_enforced_in_production(ball_face,stable_cell):
    for s,_,c in (ball_face,stable_cell):
        production=validity.check(s,EOS,OPACITY,c,policy='production',face_luminosity=s.face_luminosity)
        replica=validity.check(s,EOS,OPACITY,c,policy='published_assumption',face_luminosity=s.face_luminosity)
        screens=replica[1]['screens']
        # The outer cells are optically thin enough to fire the tau < 10 Knudsen screen.
        assert screens['radiation_knudsen_outer']['violated'] and screens['radiation_knudsen_outer']['action']=='report_only'
        assert 'radiation_knudsen_outer' in replica[1]['screen_warnings']
        kn=validity.SCREENS['radiation_knudsen_outer'][1];outer=validity.SCREENS['outer_face_knudsen'][1]
        assert kn in production[0] and outer in production[0] and kn not in replica[0] and outer not in replica[0]
        assert screens['outer_face_knudsen']['action']=='not_applicable'
        assert screens['outer_face_eddington']['value']<=validity.EDDINGTON_RELATION_TOLERANCE
        # ML1 includes radiative losses: Pe < 100 does not apply in either policy.
        assert screens['convective_peclet']['action']=='not_applicable'
        assert production[1]['screens']['convective_peclet']['action']=='not_applicable'
        # Always enforced in both: the ideal-gas Saha screen fires on this cool fully ionized toy.
        saha=validity.SCREENS['ideal_eos_helium_saha'][1]
        assert saha in production[0] and saha in replica[0]
        # Deep screens hold and the conservation and geometry rejections are recorded.
        assert not screens['radiation_knudsen_deep']['violated'] and not screens['free_streaming_deep']['violated']
        assert screens['volume_geometry']['action']=='reject' and screens['conservation']['action']=='reject'
        assert replica[1]['photosphere_eddington_factor']==pytest.approx(
            photosphere.eddington_factor(s.face_luminosity[-1],s.enclosed[-1],OPACITY.kappa),rel=1e-12)
        assert screens['density_inversion']['action']=='diagnostic' and not screens['density_inversion']['violated']
    # A custom policy moves an enforced screen to report-only.
    s,_,c=stable_cell
    custom=validity.check(s,EOS,OPACITY,c,policy=dict(name='test',report_only=('ideal_eos_helium_saha',)),
                          face_luminosity=s.face_luminosity)
    assert validity.SCREENS['ideal_eos_helium_saha'][1] not in custom[0] and 'ideal_eos_helium_saha' in custom[1]['screen_warnings']
    # Always-enforced domain screens fire in both policies.
    hot={**c,'validity':{**c['validity'],'max_temperature_K':1e6}}
    for policy in ('production','published_assumption'):
        assert validity.SCREENS['eos_temperature_domain'][1] in validity.check(s,EOS,OPACITY,hot,policy=policy,face_luminosity=s.face_luminosity)[0]
    t=s.copy();t.photosphere_density=s.photosphere_density;t.face_luminosity=s.face_luminosity*3
    edd=validity.check(t,EOS,OPACITY,c,policy='published_assumption',face_luminosity=t.face_luminosity)
    assert edd[1]['screens']['photosphere_gas_pressure']['violated']==(edd[1]['photosphere_eddington_factor']>=1)
    with pytest.raises(ValueError):validity.check(s,EOS,OPACITY,c,policy='lenient')
    table=validity.describe('published_assumption',c,EOS,OPACITY)
    assert table['screens']['convective_mach']['action']=='report_only' and table['screens']['flow_mach']['action']=='enforced'
    assert table['screens']['opacity_clamp']['action']=='not_applicable'

def test_stars_opacity_wrapper_routes_the_clamp_flags():
    calls=[]
    class Bridge:
        def opacity_stars(self,rho,T,clamp=False,flags=False):
            calls.append((clamp,flags));k=np.full(np.broadcast(rho,T).shape,.5)
            return (k,np.ones(k.shape,dtype=np.int32)) if flags else k
    clamped=quasistatic.StarsOpacity(Bridge(),clamp=True);strict=quasistatic.StarsOpacity(Bridge())
    assert np.all(clamped(np.ones(3),np.ones(3))==.5) and calls[-1]==(True,False)
    assert np.all(clamped.clamp_flags(np.ones(3),np.ones(3))==1) and calls[-1]==(True,True)
    assert np.all(strict.clamp_flags(np.ones(2),np.ones(2))==0)

# ---------------------------------------------------------------- Stage 4 review: base mesh and controller

def test_base_levels_are_exact_and_invertible(ball_face):
    s,_,c=ball_face
    lv=quasistatic.BaseLevels(s,EOS,c['numerics']['zones'])
    mu=np.cumsum(s.dm);N=lv(mu)
    assert lv(0.)==0. and np.all(np.diff(N)>0)
    # Round trip through the inner half (near the surface mu carries ~1e-16 of the envelope mass).
    levels=np.linspace(.1,float(lv(.5*mu[-1])),57)
    assert np.max(abs(lv(lv.inverse(levels))/levels-1))<=1e-12
    # The count is the integral of 1/t with ln t linear between centres (fine quadrature per segment).
    quad=[0.]
    for k in range(40):
        x=np.linspace(lv.mu[k],lv.mu[k+1],4001);f=np.exp(-np.interp(x,lv.mu,lv.a))
        quad.append(quad[-1]+math.fsum(.5*(f[1:]+f[:-1])*np.diff(x)))
    assert np.max(abs(lv.N[1:41]/np.array(quad[1:])-1))<=1e-7

def test_fixed_mass_regrid_conserves_and_keeps_the_lagrangian_part(ball_face):
    s,_,c=ball_face;zones=c['numerics']['zones']
    s=zero_time(s)
    cut,_=cavity.slice_base(s,3*s.dm[0],EOS)  # a slice deeper than the base cell
    out,info=quasistatic.regrid_base(s,cut,EOS,zones,.05)
    k,j=info['region_cells'],info['region_replaced_cells']
    lv=quasistatic.BaseLevels(s,EOS,zones);mu=np.cumsum(out.dm)
    assert k>=10 and np.max(abs(lv(mu[:k-1])-np.arange(1,k)))<=1e-9 and .5<=lv(mu[k-1])-(k-1)<1.5
    assert mu[k-1]>=.05*mu[-1] and np.sum(np.cumsum(cut.dm)[:j-1]>=.05*mu[-1])==0
    # Above the region everything is the sliced state, bitwise.
    for field in ('dm','T','cell_volume'):assert np.array_equal(getattr(out,field)[k:],getattr(cut,field)[j:])
    for field in ('r','v','face_luminosity'):assert np.array_equal(getattr(out,field)[k:],getattr(cut,field)[j:])
    # Mass, volume and total energy of the region; the EOS-inversion defect is booked.
    assert abs(math.fsum(out.dm[:k])/math.fsum(cut.dm[:j])-1)<=1e-15 and out.dm.sum()==pytest.approx(cut.dm.sum(),rel=1e-15)
    assert abs(math.fsum(out.volume[:k])/math.fsum(cut.volume[:j])-1)<=1e-13
    E0,E1=cut.energies(EOS),out.energies(EOS);scale=s.energies(EOS)['internal']
    assert abs(E1['total']-E0['total']-info['energy_defect'])<=1e-13*scale and abs(info['energy_defect'])<=1e-13*scale
    assert out.ledger['remap_energy_defect']-cut.ledger['remap_energy_defect']==pytest.approx(info['energy_defect'],abs=1e-16*scale)
    assert info['potential_reallocation']==pytest.approx(E1['gravity']-E0['gravity'],rel=1e-12)
    # split_merge leaves the protected region alone (only its top cell may take a sliver from above).
    post,_=quasistatic.remesh_base(s,out,EOS,zones,protect=k)
    assert np.array_equal(post.dm[:k-1],out.dm[:k-1]) and np.array_equal(post.r[:k],out.r[:k])

def test_registered_base_resplit_jump_is_dt_independent_and_reported(ball_face,ball_fixed):
    """The review's second blocker, forced at step 1: lowering the split threshold to 1.5 re-splits a
    base cell (x = 1.58), and the face boundary state jumps by about -11% whatever dt is. The
    counterfactual solve on the pre-event mesh attributes the jump to the remesh."""
    s,_,c=ball_face;z=zero_time(s);before=stepcontrol.observables(z,EOS,'face')
    cc=copy.deepcopy(c);cc['quasistatic']=dict(split=1.5)
    jumps=[]
    for dt in (1e-2,1e-3):
        _,d=quasistatic.bondi_step(z,dt*YEAR,EOS,OPACITY,cc,mode='evolve')
        total=math.log(d['base_density']/before['base_density'])
        assert d['splits']>=1 and d['remesh_effect_solved'] and abs(total)>=.09
        assert abs(d['remesh_dln_base_density']-total)<=.02*abs(total)
        assert abs(d['remesh_dln_inner_radius']-math.log(d['inner_radius']/before['inner_radius']))<=.02*abs(total)
        jumps.append(total)
    assert abs(jumps[0]/jumps[1]-1)<=.05
    _,d=quasistatic.bondi_step(z,1e-3*YEAR,EOS,OPACITY,c,mode='evolve')
    assert d['splits']==d['merges']==0 and d['remesh_effect_solved'] is None and d['remesh_dln_base_density']==0.
    # Fixed-mass base with the same thresholds: the region is never re-split and the base follows dt.
    s,_,c=ball_fixed;z=zero_time(s);before=stepcontrol.observables(z,EOS,'face')
    cc=copy.deepcopy(c);cc['quasistatic']['split']=1.5
    _,d=quasistatic.bondi_step(z,1e-3*YEAR,EOS,OPACITY,cc,mode='evolve')
    assert abs(math.log(d['base_density']/before['base_density']))<=1e-5

def test_fixed_mass_base_is_well_balanced_and_continuous_in_dt(ball_fixed):
    s,d,c=ball_fixed
    assert d['base_region_cells']>=10 and d['base_region_mass']>=.05*s.dm.sum()
    new,dd=quasistatic.bondi_step(s,1e6*YEAR,EOS,OPACITY,c,mode='relax')
    assert len(new.dm)==len(s.dm) and dd['splits']==dd['merges']==0
    assert np.max(abs(np.log(new.T/s.T)))<=1e-10 and np.max(abs(np.log(new.r/s.r)))<=1e-10
    assert np.max(abs(new.face_luminosity/s.face_luminosity-1))<=1e-10 and abs(new.cavity_mass/s.cavity_mass-1)<=1e-10
    z=zero_time(s);before=stepcontrol.observables(z,EOS,'face')
    for dt in (1e-2,1e-4,1e-6):
        new,dd=quasistatic.bondi_step(z,dt*YEAR,EOS,OPACITY,c,mode='evolve')
        after=stepcontrol.observables(new,EOS,'face')
        for key in ('base_density','base_temperature','inner_radius'):
            assert abs(math.log(after[key]/before[key]))<=1e-5,(dt,key)
        assert dd['base_region_cells']==d['base_region_cells'] and abs(dd['regrid_energy_defect'])<=1e-15*abs(z.ledger['initial_energy_scale'])
        assert abs(dd['energy_residual_change'])<=5e-13 and abs(dd['mass_residual_change'])<=1e-14

def test_registered_controller_pins_dt_and_amendment_releases_it(ball_face,ball_fixed):
    """Evidence for draft amendment 1 (the plan's amendment rule asks for a test that fails with the
    registered numerics and passes with the amended ones). Registered: Lagrangian metric and base;
    the profile change at face 1 tracks swept/dm0, so dt stays near 0.03-0.06 yr. Amended: mapped
    metric and fixed-mass base; no rejection and dt above 20 yr within 14 steps."""
    def run(s,c,metric,steps):
        c=copy.deepcopy(c);c['stepcontrol']=dict(metric=metric)
        S=stepcontrol.settings(c);s=zero_time(s);dt=stepcontrol.start(s.aux,S)
        before=stepcontrol.observables(s,EOS,'face');dts=[];rejections=0
        while len(dts)<steps:
            try:
                new,d=quasistatic.bondi_step(s,dt,EOS,OPACITY,c,mode='evolve')
                after=stepcontrol.observables(new,EOS,'face');m=stepcontrol.metrics(before,after,d,S)
                if stepcontrol.rejected(m['delta_actual'],S):raise StepRejected('delta_actual')
            except StepRejected:
                rejections+=1;dt=stepcontrol.reject(s.aux,dt,S);assert dt is not None;continue
            assert abs(d['energy_residual_change'])<=5e-13 and abs(d['mass_residual_change'])<=1e-14
            dts.append(d['dt']/YEAR);s,before=new,after
            dt=stepcontrol.accept(s.aux,d['dt'],m['delta_actual'],d['newton_iterations'],S)
        return np.array(dts),rejections
    registered,_=run(ball_face[0],ball_face[2],'lagrangian',12)
    amended,rejections=run(ball_fixed[0],ball_fixed[2],'mapped',14)
    assert registered[-6:].max()<=.1
    assert rejections==0 and amended.max()>=20.

# ---------------------------------------------------------------- Stage 4 review: rejections, ledger labels, budget, screens

def test_newton_counts_trial_points_outside_the_model_domain():
    from scipy.sparse import csc_matrix
    def fun(x):
        if np.real(x[0])<-3:raise DomainError('Gamma_ph >= 1 at this trial point')
        return np.arctan(x)
    x,info=quasistatic.newton(fun,np.array([2.]),csc_matrix(np.ones((1,1))),1e-12)
    assert info['success'] and abs(x[0])<=1e-12 and info['domain_trials']==1 and 'Gamma_ph' in info['domain_message']
    _,info=quasistatic.newton(fun,np.array([-4.]),csc_matrix(np.ones((1,1))),1e-12)
    assert not info['success'] and info['domain_trials']==1 and 'outside the model domain' in info['message']

def test_domain_exits_and_pre_solve_failures_are_step_rejections(stable_cell,monkeypatch):
    s,_,c=stable_cell;s=zero_time(s)
    # P1: a predicted cavity change of -0.2 dm0 returns more gas than fits: retried, not a crash.
    t=s.copy();t.aux['predictor']=dict(time=[-10*YEAR,0.],rest_rate=[0.,0.],cavity_mass=[s.cavity_mass+.2*s.dm[0],s.cavity_mass])
    with pytest.raises(StepRejected,match='P1 base slice failed') as e:quasistatic.bondi_step(t,10*YEAR,EOS,OPACITY,c,mode='evolve')
    assert e.value.domain_exit is False and e.value.domain_trials==0
    # P2: a split/merge failure.
    def negative(*a,**k):raise ValueError('Conservative split/merge would require negative internal energy')
    monkeypatch.setattr(quasistatic,'split_merge',negative)
    with pytest.raises(StepRejected,match='P2 base remesh failed'):quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')
    monkeypatch.undo()
    # P3: the predictor's photosphere does not exist (seed without rho_ph): a marked domain exit.
    def no_surface(*a,**k):raise DomainError('Photospheric gas pressure is not positive (Gamma_ph >= 1); no Ball surface exists')
    t=s.copy();t.photosphere_density=None
    monkeypatch.setattr(photosphere,'solve',no_surface)
    with pytest.raises(StepRejected,match='P3 predictor left the model domain') as e:
        quasistatic.bondi_step(t,YEAR,EOS,OPACITY,c,mode='relax')
    assert e.value.domain_exit and 'Gamma_ph' in e.value.domain_message
    monkeypatch.undo()
    # A failed solve whose trial points left the domain carries the record.
    def failing(fun,x,*a,**k):
        return x,dict(success=False,residual=.3,message='Newton line search or iteration limit',domain_trials=4,
                      domain_message='Photospheric gas pressure is not positive (Gamma_ph >= 1)')
    monkeypatch.setattr(quasistatic,'newton',failing)
    with pytest.raises(StepRejected,match='4 trial points outside the model domain') as e:
        quasistatic.bondi_step(s,YEAR,EOS,OPACITY,c,mode='evolve')
    assert e.value.domain_exit and e.value.domain_trials==4 and 'Gamma_ph' in e.value.domain_message
    monkeypatch.undo()
    _,d=quasistatic.bondi_step(s,YEAR,EOS,OPACITY,c,mode='evolve')
    assert d['domain_trials']==0 and d['domain_message'] is None

class FlakyEOS:
    """EOS proxy whose temperature inversion fails `failures` times, as MesaMicrophysics.temperature
    and IdealGasRadiation.temperature do when their Newton iteration does not converge."""
    def __init__(self,eos,failures=1):self.eos,self.failures,self.failed=eos,failures,0
    def __getattr__(self,name):return getattr(self.eos,name)
    def temperature(self,*args,**kwargs):
        if self.failed<self.failures:
            self.failed+=1;raise RuntimeError('MESA temperature inversion failed')
        return self.eos.temperature(*args,**kwargs)

def test_inversion_failures_are_step_rejections_retried_at_half_the_step(stable_cell,tmp_path,monkeypatch):
    """A RuntimeError of P1-P3 or of the P5 sink is a rejection like any other pre-solve failure; before
    the Stage 6 review one failed inversion ended a track as numerical_failure and escaped relax()."""
    s,_,c=stable_cell;s=zero_time(s)
    # P5: on a plain slice the only inversion of the step is the cell-0 sink.
    flaky=FlakyEOS(EOS)
    with pytest.raises(StepRejected,match='P5 sink failed.*temperature inversion failed') as e:
        quasistatic.bondi_step(s,10*YEAR,flaky,OPACITY,c,mode='evolve')
    assert flaky.failed==1 and e.value.domain_exit is False and e.value.domain_trials==0
    new,_=quasistatic.bondi_step(s,10*YEAR,flaky,OPACITY,c,mode='evolve')  # the proxy is exact once it stops failing
    assert np.array_equal(new.T,quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')[0].T)
    # P2: the fixed-mass base regrid inverts the EOS in every remapped cell.
    fixed=copy.deepcopy(c);fixed['quasistatic']=dict(base='fixed_mass')
    with pytest.raises(StepRejected,match='P2 base remesh failed.*temperature inversion failed'):
        quasistatic.bondi_step(s,10*YEAR,FlakyEOS(EOS),OPACITY,fixed,mode='evolve')
    # P3: the predictor's photosphere solve (seed without rho_ph) does not converge.
    def stalled(*a,**k):raise RuntimeError('Photospheric density solve failed at Teff=5000 K')
    t=s.copy();t.photosphere_density=None
    monkeypatch.setattr(photosphere,'solve',stalled)
    with pytest.raises(StepRejected,match='P3 predictor failed.*solve failed') as e:quasistatic.bondi_step(t,YEAR,EOS,OPACITY,c,mode='relax')
    assert e.value.domain_exit is False
    monkeypatch.undo()
    # Track driver: the attempt is logged as a rejection and retried at dt/2; the track goes on.
    flaky=FlakyEOS(EOS)
    _,report=quasistatic.run(s,flaky,OPACITY,c,tmp_path/'track',policy=None,stop_after_steps=2)
    events=[json.loads(x) for x in open(tmp_path/'track'/'events.jsonl')]
    rejected=[x for x in events if x['event']=='rejected'];accepted=[x for x in events if x['event']=='accepted']
    assert report['status']=='control_stop' and report['rejections']==1 and report['accepted_steps_this_run']==2 and flaky.failed==1
    assert len(rejected)==1 and 'temperature inversion failed' in rejected[0]['message'] and rejected[0]['domain_exit'] is False
    assert accepted[0]['dt']==pytest.approx(.5*rejected[0]['dt_years']*YEAR,rel=1e-15)
    # Relaxation: the same failure is retried, not raised.
    seed,sc=toy(48,'cell');flaky=FlakyEOS(EOS)
    _,record=quasistatic.relax(seed,flaky,OPACITY,sc,dict(mesh_cycles=0))
    assert record['status']=='converged' and flaky.failed==1
    assert 'temperature inversion failed' in record['rejections'][0]['message'] and record['rejections'][0]['dt_years']==.01
    assert record['steps'][0]['dt_years']==.005

def test_sink_flow_work_is_booked_as_inner_boundary_work(stable_face):
    """The sink takes delta h_0 at fixed volume; a slice of delta would take it at fixed (rho, u) and
    the piston (p_0 dV_0 plus the face-0 potential term) the rest. With the flow work
    delta p_face/rho_0 booked to inner_boundary_work, a 3% overprediction (slice 3% more, the sink
    returns it) gives the same solution and nearly the same channel split: what remains is second
    order (the review's p_0/rho_0 alone left 58% of the shift on this coarse base, p_b/p_0 = 2.1)."""
    s,_,c=stable_face;s=zero_time(s)
    for _ in range(2):s,_=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')
    over=s.copy();h=over.aux['predictor'];h['cavity_mass']=[h['cavity_mass'][0]-.03*s.dm[0],h['cavity_mass'][1]]
    runs=[(t,)+quasistatic.bondi_step(t,10*YEAR,EOS,OPACITY,c,mode='evolve') for t in (s,over)]
    (_,well,dw),(_,mis,dm_)=runs
    assert dm_['delta_fraction']<-.02 and dm_['passes']==1 and dm_['sink_flow_work']<0<dw['sink_flow_work']+1e300
    for t,new,d in runs:
        l0,l1=t.ledger,new.ledger
        assert l1['inner_boundary_work']-l0['inner_boundary_work']==pytest.approx(d['piston_work']-d['sink_flow_work'],rel=1e-9)
        assert sum(l1[k]-l0[k] for k in ('advected_energy','cavity_swept_energy'))==pytest.approx(d['swept_energy'],rel=1e-9)
        assert abs(d['energy_residual_change'])<=5e-13
    assert mis.bh_mass==pytest.approx(well.bh_mass,rel=1e-12)
    assert abs(dm_['swept_minus_piston_over_bh_energy']-dw['swept_minus_piston_over_bh_energy'])<=1e-6
    after=abs(dm_['swept_energy']-dw['swept_energy'])
    before=abs((dm_['swept_energy']+dm_['sink_flow_work'])-(dw['swept_energy']+dw['sink_flow_work']))
    assert after<=.2*before

def test_per_step_conservation_budget_is_enforced_outside_relaxation(stable_cell):
    s,_,c=stable_cell;s=zero_time(s)
    _,d=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c,mode='evolve')
    assert abs(d['energy_residual_change'])<=quasistatic.DEFAULTS['quasistatic']['energy_per_step']
    impossible=copy.deepcopy(c);impossible['quasistatic']=dict(energy_per_step=-1.)
    with pytest.raises(StepRejected,match='conservation budget'):quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,impossible,mode='evolve')
    quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,impossible,mode='relax')  # relaxation ledger is reset at t = 0

def test_screen_gaps_photosphere_clamp_and_disabled_diffusion_screens(stable_cell):
    s,_,c=stable_cell
    class PhotosphereClamp:
        """Clamps only at the photospheric point (rho_ph, Teff)."""
        def __call__(self,rho,T):return OPACITY(rho,T)
        def clamp_flags(self,rho,T):return (np.asarray(rho)==s.photosphere_density).astype(np.int32)
    for policy,action in (('published_assumption','report_only'),('production','enforced')):
        errors,diag,_=validity.check(s,EOS,PhotosphereClamp(),c,policy=policy,face_luminosity=s.face_luminosity)
        screen=diag['screens']['opacity_clamp']
        assert diag['opacity_clamped_cells']==0 and diag['opacity_photosphere_clamped'] is True
        assert screen['violated'] and screen['action']==action
        assert (validity.SCREENS['opacity_clamp'][1] in errors)==(action=='enforced')
    assert validity.check(s,EOS,OPACITY,c,policy='production',face_luminosity=s.face_luminosity)[1]['opacity_photosphere_clamped'] is None
    off={**c,'validity':{**c['validity'],'require_optically_thick':False}}
    table=validity.describe('published_assumption',off,EOS,OPACITY)['screens']
    for key in ('radiation_knudsen_deep','radiation_knudsen_outer','outer_face_knudsen'):
        assert table[key]['action']=='not_applicable' and 'require_optically_thick' in table[key]['reason']
    assert validity.describe('published_assumption',c,EOS,OPACITY)['screens']['radiation_knudsen_deep']['action']=='enforced'
    errors,diag,_=validity.check(s,EOS,OPACITY,off,policy='production',face_luminosity=s.face_luminosity)
    assert diag['screens']['radiation_knudsen_outer']['action']=='not_applicable' and not diag['screens']['radiation_knudsen_outer']['violated']
    assert validity.SCREENS['radiation_knudsen_outer'][1] not in errors

# ---------------------------------------------------------------- Stage 7 unblock: seed surface (plan amendment 3)

def outer_luminosity(s,c):
    """L_N of the step's grey surface relation for the last cell of s."""
    from qstars.transport import luminosities
    bc=copy.deepcopy(c);bc['transport']['convection']='off'
    return float(luminosities(s,EOS,OPACITY,bc,0.)[0][-1])

@pytest.mark.parametrize('radius',[150.,300.])
def test_grey_seed_surface_changes_only_the_last_cell_temperature(radius):
    """Seed surface 'grey' (off by default; plan amendment 3). The registered seed gives the last cell the
    EOS temperature of the polytrope, whose grey-relation luminosity is unrelated to the one the seed's
    photospheric pressure was evaluated at (here 0.67 and 3.1 times it; 20 times on the MESA Ball seed,
    where no photosphere then exists). 'grey' puts that one cell on the grey relation at the seed
    luminosity; nothing else changes."""
    registered,c=toy(48,'cell',radius=radius)  # toy() builds the default seed
    args=(c,EOS,OPACITY,48,MSUN,100*MSUN,radius*RSUN)
    same=lambda a,b:all(np.array_equal(getattr(a,f),getattr(b,f)) for f in ('r','dm','T','v','cell_volume')) and a.cavity_mass==b.cavity_mass
    default,info=quasistatic.seed(*args)
    assert same(default,registered) and same(quasistatic.seed(*args,surface='eos')[0],registered) and 'surface' not in info
    grey,ginfo=quasistatic.seed(*args,surface='grey')
    for f in ('r','dm','cell_volume','v'):assert np.array_equal(getattr(grey,f),getattr(registered,f)),f
    assert np.array_equal(grey.T[:-1],registered.T[:-1]) and grey.T[-1]!=registered.T[-1] and grey.cavity_mass==registered.cavity_mass
    law=info['law_luminosity']
    assert abs(outer_luminosity(registered,c)/law-1)>.1 and outer_luminosity(grey,c)==pytest.approx(law,rel=1e-12)
    assert (ginfo['surface'],ginfo['surface_temperature'],ginfo['surface_polytrope_temperature'])==('grey',grey.T[-1],registered.T[-1])
    assert grey.T[-1]>=ginfo['surface_effective_temperature']==pytest.approx(photosphere.effective_temperature(law,radius*RSUN),rel=1e-15)
    assert {k:v for k,v in ginfo.items() if not k.startswith('surface')}==info
    # Its photosphere is the one whose pressure closed the polytrope.
    rho_ph=photosphere.solve(outer_luminosity(grey,c),grey.r[-1],grey.enclosed[-1],EOS,OPACITY)
    closed=photosphere.evaluate(rho_ph,outer_luminosity(grey,c),grey.r[-1],grey.enclosed[-1],EOS,OPACITY)['pressure']
    assert closed==pytest.approx(info['external_pressure'],rel=1e-10)
    # Through the configuration (replication.seed_surface), as the registered driver passes it.
    cc=copy.deepcopy(c);cc['replication']=dict(seed_surface='grey')
    assert same(quasistatic.seed(cc,EOS,OPACITY,48,MSUN,100*MSUN,radius*RSUN)[0],grey)
    assert same(quasistatic.seed(cc,EOS,OPACITY,48,MSUN,100*MSUN,radius*RSUN,surface='eos')[0],registered)
    with pytest.raises(ValueError,match='seed surface must be one of'):quasistatic.seed(*args,surface='hot')
    cg=copy.deepcopy(c);cg['boundary']['thermal']='grey'
    with pytest.raises(ValueError,match='needs boundary.thermal stars_photosphere'):
        quasistatic.seed(cg,EOS,OPACITY,48,MSUN,100*MSUN,radius*RSUN,surface='grey')

def test_grey_surface_temperature_is_the_lowest_root_of_the_grey_relation():
    L,R,rho,dm=1e39,2e13,1e-10,1e30
    T=quasistatic.surface_temperature(L,R,rho,dm,OPACITY,1e7)
    tau=.5*OPACITY.kappa*dm/(4*PI*R**2)
    assert T==pytest.approx((3*L*(tau+4/3)/(4*PI*R**2*A*C))**.25,rel=1e-13)  # constant kappa: closed form
    assert T>photosphere.effective_temperature(L,R)
    # Two roots (kappa steps up by 30 above 2e4 K, between them): the lowest one is returned.
    step=lambda rho,T:OPACITY.kappa*np.where(np.asarray(T)<2e4,1.,30.)
    high=(3*L*(30*tau+4/3)/(4*PI*R**2*A*C))**.25
    assert T<2e4<high and quasistatic.surface_temperature(L,R,rho,dm,step,1e7)==pytest.approx(T,rel=1e-13)
    with pytest.raises(ValueError,match='no grey-relation temperature'):
        quasistatic.surface_temperature(L,R,rho,dm,OPACITY,.5*(T+photosphere.effective_temperature(L,R)))

def test_photosphere_bracket_option_reaches_every_solve_and_is_off_by_default(monkeypatch):
    """quasistatic.photosphere_bracket (plan amendment 4; registered 'doubling') is the bracket rule of every
    photosphere.solve of the Ball path: the seed's outer pressure, the predictor of a step from a state
    without a photospheric density (the first step from a seed) and the frozen-entropy problem of fold.py
    for such a state. Where the registered rule brackets the root (here: electron scattering) both rules
    give the same seed and the same step, bitwise."""
    from qstars import fold
    solve=photosphere.solve;calls=[]
    def traced(*args,**kwargs):
        calls.append(kwargs.get('bracket'));return solve(*args,**kwargs)
    monkeypatch.setattr(photosphere,'solve',traced)
    out={}
    for bracket in ('doubling','scan'):
        _,c=toy(48,'cell',radius=300.)
        c['quasistatic']=dict(c.get('quasistatic') or {},photosphere_bracket=bracket)
        calls.clear()
        s,_=quasistatic.seed(c,EOS,OPACITY,48,MSUN,100*MSUN,300*RSUN)
        assert calls==[bracket] and s.photosphere_density is None
        ledger.initialize(s,EOS,cavity=True)
        dt=1e-2*YEAR
        for _ in range(30):  # the first relaxation steps from this seed, halving dt on a rejection as relax() does
            try:new,d=quasistatic.bondi_step(s,dt,EOS,OPACITY,c,mode='relax');break
            except StepRejected:dt*=.5
        else:raise AssertionError('no relaxation step from the toy seed')
        assert calls[1:] and set(calls[1:])=={bracket}
        calls.clear();fold.Frozen(s,EOS,OPACITY,c)
        assert calls==[bracket]
        out[bracket]=(s,new)
    (a,b),(x,y)=out['doubling'],out['scan']
    assert np.array_equal(a.T,x.T) and np.array_equal(a.r,x.r) and np.array_equal(b.T,y.T) and np.array_equal(b.r,y.r)
    assert b.photosphere_density==y.photosphere_density

def test_registered_and_grey_seeds_relax_to_the_same_row0():
    """A seed is only the initial guess of the relaxation: with the mesh at its STARS-Q fixed point the
    registered and the grey-surface seed of one radius give the same row 0 (48-zone cell toy)."""
    out=[]
    for surface in ('eos','grey'):
        _,c=toy(48,'cell',radius=300.)
        s,_=quasistatic.seed(c,EOS,OPACITY,48,MSUN,100*MSUN,300*RSUN,surface=surface)
        ledger.initialize(s,EOS,cavity=True)
        s,record=quasistatic.relax(s,EOS,OPACITY,c,dict(mesh_cycles=40))
        assert record['status']=='converged' and len(s.dm)==48
        out.append(s)
    a,b=out
    assert np.max(abs(np.log(a.r/b.r)))<=1e-6 and np.max(abs(np.log(a.T/b.T)))<=1e-6
    assert abs(math.log(a.face_luminosity[-1]/b.face_luminosity[-1]))<=1e-6 and abs(math.log(a.cavity_mass/b.cavity_mass))<=1e-6

# ---------------------------------------------------------------- Stage 7 unblock: fallback F1 (pseudo-transient continuation)

def test_pseudo_transient_continuation_reaches_a_root_newton_cannot_start_from():
    """Kelley-Keyes SER continuation on a two-row differential-algebraic toy: x0^2 = 2 carries the
    pseudo-time term, x1 = x0 is a constraint. The Jacobian is singular at the start, so the plain Newton
    solve fails there; the continuation hands over near the root and Newton then meets its tolerance."""
    from scipy.sparse import csc_matrix
    fun=lambda x:np.array([x[0]*x[0]-2.,x[1]-x[0]])
    pattern=csc_matrix(np.array([[1.,0.],[1.,1.]]));pairs=(np.array([0]),np.array([0]));start=np.zeros(2)
    _,plain=quasistatic.newton(fun,start,pattern,1e-12)
    assert not plain['success'] and 'Jacobian failure' in plain['message']
    x,info=quasistatic.pseudo_transient(fun,start,pattern,pairs,1e-12,**quasistatic.PTC)
    assert info['handover'] and info['shrinks']==0 and 0<info['iterations']<quasistatic.PTC['max_iterations']
    assert info['delta']>=quasistatic.PTC['handover'] or info['residual']<=1e-12
    x,final=quasistatic.newton(fun,x,pattern,1e-12)
    assert final['success'] and final['residual']<=1e-12 and np.allclose(x,math.sqrt(2.),rtol=1e-13,atol=0)
    # Failures are reported, never a handover: the iteration limit, and a start outside the model domain.
    _,short=quasistatic.pseudo_transient(fun,start,pattern,pairs,1e-12,max_iterations=2)
    assert not short['handover'] and 'iteration limit' in short['message'] and short['iterations']==2
    def outside(x):raise DomainError('Gamma_ph >= 1 at this trial point')
    _,out=quasistatic.pseudo_transient(outside,start,pattern,pairs,1e-12)
    assert not out['handover'] and out['domain_trials']==1 and 'outside the model domain' in out['message']
    # The pairing of the step: the first-law rows with ln T and the growth row with X, nothing else.
    rows,cols=quasistatic.pairing(5)
    assert list(rows)==[5,6,7,8,9,17] and list(cols)==[6,7,8,9,10,17]

def state_arrays_equal(a,b):
    return all(np.array_equal(getattr(a,f),getattr(b,f)) for f in ('r','T','dm','v','cell_volume','face_luminosity','radius_compensation')) \
        and (a.bh_mass,a.cavity_mass,a.photosphere_density,a.time)==(b.bh_mass,b.cavity_mass,b.photosphere_density,b.time) \
        and json.dumps(a.ledger,sort_keys=True)==json.dumps(b.ledger,sort_keys=True) and json.dumps(a.aux,sort_keys=True)==json.dumps(b.aux,sort_keys=True)

def test_fallback_f1_is_off_by_default_and_never_changes_a_converging_step(stable_cell):
    """quasistatic.ptc False (default) is the registered solve; with ptc True a step whose Newton solve
    converges is bitwise the same, and relaxation steps never use the fallback."""
    s,_,c=stable_cell;s=zero_time(s)
    assert quasistatic.settings(c)['quasistatic']['ptc'] is False and quasistatic.DEFAULTS['quasistatic']['ptc'] is False
    # The photosphere bracket is the registered one unless amended (plan amendment 4: 'scan').
    assert quasistatic.DEFAULTS['quasistatic']['photosphere_bracket']=='doubling'==quasistatic.settings(c)['quasistatic']['photosphere_bracket']
    with pytest.raises(ValueError,match='photosphere_bracket must be one of'):
        quasistatic.settings({**c,'quasistatic':dict(photosphere_bracket='bisect')})
    reference,d0=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,c)
    assert not any(k.startswith('ptc') for k in d0)
    for ptc in (False,True):
        cc=copy.deepcopy(c);cc['quasistatic']=dict(ptc=ptc)
        new,d=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,cc)
        assert state_arrays_equal(new,reference) and all(d[k]==d0[k] or (d[k]!=d[k] and d0[k]!=d0[k]) for k in d0)
        assert ({k:v for k,v in d.items() if k not in d0}=={'ptc_entered':False,'ptc_iterations':0}) if ptc else set(d)==set(d0)
        landed,dl=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,cc,mode='land',target_mass=reference.bh_mass)
        assert landed.bh_mass==pytest.approx(reference.bh_mass,rel=1e-14) and ('ptc_entered' in dl)==ptc
        relaxed,dr=quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,cc,mode='relax')
        assert not any(k.startswith('ptc') for k in dr)
    cc=copy.deepcopy(c);cc['quasistatic']=dict(ptc='yes')
    with pytest.raises(ValueError,match='ptc must be true or false'):quasistatic.settings(cc)

def test_fallback_f1_ends_on_the_registered_equations_and_tolerance(stable_cell,monkeypatch):
    """Fallback F1: when the damped Newton solve fails, the system is restarted from the predictor with
    pseudo-transient continuation and the plain Newton solve takes over. The accepted unknowns satisfy the
    registered residual (no pseudo-time term) to the registered 1e-12, and they are the solution the plain
    solve finds; without the switch the same failure is a step rejection, as registered. (A 100 yr step of
    the 64-zone toy. The continuation is not attracting at every step solution: at 10 yr on this toy it
    stalls where Newton converges in two iterations, which is then reported as a failed fallback.)"""
    s,_,c=stable_cell;s=zero_time(s)
    on=copy.deepcopy(c);on['quasistatic']=dict(ptc=True)
    reference,d0=quasistatic.bondi_step(s,100*YEAR,EOS,OPACITY,c)
    newton=quasistatic.newton;calls=[];failed=[]
    def first_fails(fun,x,*a,**k):
        calls.append(1)
        if len(calls)>1:return newton(fun,x,*a,**k)
        _,info=newton(fun,x,*a,**k);failed.append(info['iterations'])
        return x,dict(info,success=False,message='Newton line search or iteration limit (forced by the test)')
    monkeypatch.setattr(quasistatic,'newton',first_fails)
    with pytest.raises(StepRejected,match='Quasi-static solve failed') as rejected:quasistatic.bondi_step(s,100*YEAR,EOS,OPACITY,c)
    assert 'F1' not in str(rejected.value)
    # The solve itself: P0-P2 as in bondi_step, then Problem.solve with the fallback.
    S=quasistatic.settings(on);predicted=quasistatic.predict(s,100*YEAR,EOS,on,'evolve',None,S)
    cut,_=cavity.slice_base(s,predicted['swept'],EOS)
    post,_=quasistatic.remesh_base(s,cut,EOS,on['numerics']['zones'])
    problem=quasistatic.Problem(s,post,100*YEAR,EOS,OPACITY,on,'evolve',None,predicted,S)
    calls.clear();x,info=problem.solve()
    assert info['success'] and info['ptc_entered'] and info['ptc_handover'] and info['ptc_iterations']>=1 and len(calls)==2
    assert info['residual']<=1e-12 and np.max(abs(problem.evaluate(x)))<=1e-12  # the registered residual, pseudo-time removed
    assert 'forced by the test' in info['newton_message'] and info['nfev']>info['ptc_iterations']
    # Second review: a rescued step reports every nonlinear iteration it needed (the failed Newton solve, the
    # pseudo-time iterations, the final Newton solve), so the registered controller rule 'extra factor 0.8
    # if Newton needed more than 8 iterations' sees them; before, only the final Newton iterations (1 here).
    assert info['iterations']==failed[-1]+info['ptc_iterations']+info['final_newton_iterations']>info['final_newton_iterations']
    # Through the step: accepted, booked and conserved like any other step, and the same solution.
    calls.clear();new,d=quasistatic.bondi_step(s,100*YEAR,EOS,OPACITY,on)
    assert d['ptc_entered'] and d['ptc_handover'] and d['nonlinear_residual']<=1e-12 and d['ptc_message']=='Pseudo-time removed'
    assert d['newton_iterations']==failed[-1]+d['ptc_iterations']+d['final_newton_iterations']
    assert np.max(abs(np.log(new.r/reference.r)))<=1e-10 and np.max(abs(np.log(new.T/reference.T)))<=1e-10
    assert new.bh_mass==pytest.approx(reference.bh_mass,rel=1e-13) and new.cavity_mass==pytest.approx(reference.cavity_mass,rel=1e-9)
    assert abs(d['energy_residual_change'])<=5e-13 and abs(d['mass_residual_change'])<=1e-14
    # When the continuation does not rescue the solve either, the rejection says so.
    def always_fails(fun,x,*a,**k):
        _,info=newton(fun,x,*a,**k)
        return x,dict(info,success=False,message='Newton line search or iteration limit (forced by the test)')
    monkeypatch.setattr(quasistatic,'newton',always_fails)
    with pytest.raises(StepRejected,match='F1 pseudo-transient continuation handed over, Newton then failed'):
        quasistatic.bondi_step(s,100*YEAR,EOS,OPACITY,on)

def test_fallback_f1_that_does_not_converge_is_a_reported_rejection(stable_cell,ball_face,monkeypatch):
    """Preserved negative results. (1) The pseudo-time flow of the step equations is not attracting at
    every step solution: at 10 yr on the 64-zone toy the continuation stalls although Newton converges in
    two iterations. (2) On the Ball row-0 toy a 3 yr step of the registered base mesh has no Newton
    solution from the predictor and F1 does not find one either (the base slice takes most of the base
    cell; fallback F2 addresses that). Both end as ordinary step rejections that name the fallback."""
    s,_,c=stable_cell;s=zero_time(s)
    on=copy.deepcopy(c);on['quasistatic']=dict(ptc=True)
    newton=quasistatic.newton;calls=[]
    def first_fails(fun,x,*a,**k):
        calls.append(1)
        if len(calls)>1:return newton(fun,x,*a,**k)
        _,info=newton(fun,x,*a,**k)
        return x,dict(info,success=False,message='Newton line search or iteration limit (forced by the test)')
    monkeypatch.setattr(quasistatic,'newton',first_fails)
    with pytest.raises(StepRejected,match='F1 pseudo-transient continuation failed: Pseudo-transient iteration limit'):
        quasistatic.bondi_step(s,10*YEAR,EOS,OPACITY,on)
    assert len(calls)==1  # no handover: the plain Newton solve was not called again
    monkeypatch.undo()
    s,_,c=ball_face;s=zero_time(s)
    on=copy.deepcopy(c);on['quasistatic']=dict(ptc=True)
    with pytest.raises(StepRejected,match='Quasi-static solve failed') as plain:quasistatic.bondi_step(s,3*YEAR,EOS,OPACITY,c)
    with pytest.raises(StepRejected,match='F1 pseudo-transient continuation failed') as fallback:quasistatic.bondi_step(s,3*YEAR,EOS,OPACITY,on)
    assert str(fallback.value).startswith(str(plain.value))
