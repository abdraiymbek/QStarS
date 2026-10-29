"""Gate 6 Stage 3: core extensions with bitwise-neutral defaults (U2, U7, U9, U12, U15)."""
import copy,math
from fractions import Fraction
import numpy as np
import pytest
from qstars.config import resolve
from qstars.constants import G,PI,MSUN,YEAR,C
from qstars.eos import IdealGasRadiation,Composition
from qstars.evolution import physics
from qstars.state import State
from qstars import ledger,structure
from qstars.engine import rates,rates_from
from qstars.mesh import remap,split_merge,q_targets,stars_q,STARS_Q
from qstars.nonlinear import colouring,pattern_jacobian

RSUN=6.957e10

def ball_config(zones,mesh='core_refined',R=3030.):
    # Ball 2011 row-0 geometry (M_BH = 5, r0 = 1.66 Rsun) on an ideal-gas seed.
    return resolve(dict(composition=dict(X=.7,Y=.3,Z=0.),
        initial=dict(envelope_mass_msun=1e4-5,bh_mass_msun=5,radius_cm=R*RSUN,inner_radius_cm=1.66*RSUN,
                     temperature_perturbation=0.),numerics=dict(zones=zones,mesh=mesh)))

@pytest.fixture(scope='module')
def ball():
    c=ball_config(200);e=IdealGasRadiation(Composition(.7,.3,0.))
    s,_=structure.initialize(c,e);ledger.initialize(s,e,cavity=True)
    s.v=np.linspace(0.,1e4,len(s.v))
    s.face_luminosity=np.linspace(1e42,2e42,len(s.r))
    s.radius_compensation=np.linspace(-1e-6,1e-6,len(s.r))
    return s,e

@pytest.fixture
def small():
    c=resolve({'numerics':{'zones':16},'engine':{'model':'off'},'boundary':{'thermal':'insulated'},'transport':{'convection':'off'}})
    e,op=physics(c);s,_=structure.initialize(c,e);ledger.initialize(s,e)
    return s,e,c

# ---------------------------------------------------------------- U12 neutrality

def legacy_rates(state,eos,config):
    from qstars.controls import value
    q=eos.evaluate(state.rho,state.T)
    rho,p,cs,gamma=(float(q[k][0]) if k!='rho' else float(state.rho[0]) for k in ('rho','p','cs','gamma1'))
    mode=config['model'];eps=config['efficiency'];eta=value(config['coupling'],state.time/YEAR)
    if mode=='off':luminosity=0.
    elif mode=='bra2008':luminosity=4*PI*G**2*eta*state.bh_mass**2*rho**1.5/p**.5
    elif mode in ('ball2011_printed','quest2025_eq11'):
        luminosity=(4 if mode=='ball2011_printed' else 16)*PI*eta*(G*state.bh_mass)**2*rho/(gamma*cs)
    elif mode=='prescribed':
        luminosity=eps*value(config['rest_mass_rate_msun_per_year'],state.time/YEAR)*MSUN/YEAR*C**2
    rest=luminosity/(eps*C**2)
    return dict(luminosity=luminosity,rest_mass_rate=rest,bh_growth_rate=(1-eps)*rest,bondi_radius=2*G*state.bh_mass/cs**2)

def legacy_remap(state,new_dm,eos):
    from qstars.mesh import overlap_integrals
    new_dm=np.asarray(new_dm,dtype=float)
    a,b=np.r_[0.,np.cumsum(state.dm)],np.r_[0.,np.cumsum(new_dm)];b[-1]=a[-1]
    gravity=-G*.5*(state.enclosed[:-1]/state.r[:-1]+state.enclosed[1:]/state.r[1:])
    specific=eos.evaluate(state.rho,state.T)['u']+gravity+.25*(state.v[:-1]**2+state.v[1:]**2)
    E,V=overlap_integrals(state.dm,new_dm,np.column_stack([state.dm*specific,state.volume])).T
    out=state.copy();out.dm=new_dm.copy();out.face_luminosity=None;out.reference_mass_fractions=new_dm/new_dm.sum()
    out.cell_volume=V.copy();out.r=np.r_[state.r[0],np.cbrt(state.r[0]**3+3*np.cumsum(V)/(4*PI))];out.r[-1]=state.r[-1]
    out.radius_compensation=np.zeros_like(out.r);out.v=np.interp(b,a,state.v)
    internal=E/new_dm+G*.5*(out.enclosed[:-1]/out.r[:-1]+out.enclosed[1:]/out.r[1:])-.25*(out.v[:-1]**2+out.v[1:]**2)
    out.T=eos.temperature(out.rho,internal)
    ds=float(np.dot(out.dm,eos.evaluate(out.rho,out.T)['s'])-np.dot(state.dm,eos.evaluate(state.rho,state.T)['s']))
    if out.ledger:out.ledger['remap_entropy_change']+=ds
    return out

def test_state_and_ledger_defaults_are_bitwise_neutral(small):
    s,e,c=small
    assert s.cavity_mass==0. and s.photosphere_density is None and s.kinetic_tracked and s.aux=={}
    assert s.enclosed.tobytes()==(s.bh_mass+np.r_[0.,np.cumsum(s.dm)]).tobytes()
    s.v[1:]=np.linspace(0,1e3,len(s.dm))
    q=e.evaluate(s.rho,s.T)['u'];internal=float(np.dot(s.dm,q))
    gravity=float(-G*np.sum(s.weights*(s.bh_mass+np.r_[0.,np.cumsum(s.dm)])/s.r));kinetic=float(.5*np.sum(s.weights*s.v**2))
    assert s.energies(e)==dict(internal=internal,gravity=gravity,kinetic=kinetic,total=internal+gravity+kinetic)
    assert set(s.ledger)==set(ledger.FIELDS)|{'initial_mass','initial_energy','initial_energy_scale'}
    assert not set(ledger.CAVITY_FIELDS)&set(s.ledger)
    assert s.ledger['initial_mass']==s.bh_mass+s.dm.sum()
    # Copies are independent, including the driver aux dictionary.
    s.aux['controller']={'dt':1.};t=s.copy();t.aux['controller']['dt']=2.
    assert s.aux['controller']['dt']==1.

def test_ledger_residual_unchanged_without_cavity_and_closed_with_it(small):
    s,e,c=small
    s.ledger['radiated_energy']=3e40;s.ledger['advected_energy']=-1e39;s.bh_mass*=1+1e-9
    l=s.ledger;en=s.energies(e)
    reference=(en['total']-l['initial_energy']+l['radiated_energy']+l['wind_energy']+l['advected_energy']
               -l['supplied_energy']-l['central_generated']-l['boundary_work']-l['central_potential_work'])
    r=ledger.residuals(s,e)
    assert r['energy_residual_erg']==float(reference)
    assert r['mass_residual_g']==float(s.bh_mass+s.dm.sum()+l['wind_mass']+l['mass_defect']-l['supplied_mass']-l['initial_mass'])
    t=s.copy();ledger.initialize(t,e,cavity=True)
    assert set(ledger.CAVITY_FIELDS)<=set(t.ledger)
    base=ledger.residuals(t,e)['energy_residual_erg']
    t.ledger['cavity_swept_energy']=5e40;t.ledger['inner_boundary_work']=2e40
    t.ledger['relaxation_injected_energy']=1e40;t.ledger['remap_energy_defect']=4e39
    assert ledger.residuals(t,e)['energy_residual_erg']==pytest.approx(base+5e40-2e40-1e40-4e39,rel=1e-14,abs=1e26)
    # The central mass includes the cavity; swept mass is a flow diagnostic only.
    t.cavity_mass=1e30;t.dm[0]-=1e30;t.ledger['cavity_swept_mass']=1e30
    assert abs(ledger.residuals(t,e)['mass_residual_relative'])<1e-15

def test_quasi_static_energies_drop_kinetic(small):
    s,e,c=small;s.v[1:]=1e5
    tracked=s.energies(e);s.kinetic_tracked=False;quasi=s.energies(e)
    assert quasi['kinetic']==0. and quasi['total']==tracked['internal']+tracked['gravity']
    s.kinetic_tracked=True;s.cavity_mass=3e30
    np.testing.assert_array_equal(s.enclosed,s.bh_mass+3e30+np.r_[0.,np.cumsum(s.dm)])

def test_default_remap_path_matches_legacy_bitwise(small):
    s,e,c=small;s.v[1:]=np.linspace(0,1e3,len(s.dm))
    for dm in (np.repeat(s.dm/2,2),s.dm*(1+.01*np.sin(np.arange(len(s.dm))))*s.dm.sum()/np.sum(s.dm*(1+.01*np.sin(np.arange(len(s.dm)))))):
        a,b=remap(s,dm,e),legacy_remap(s,dm,e)
        for f in ('r','T','dm','v','cell_volume','radius_compensation','reference_mass_fractions'):
            assert getattr(a,f).tobytes()==getattr(b,f).tobytes(),f
        assert a.ledger==b.ledger and a.face_luminosity is None and a.cavity_mass==0.

# ---------------------------------------------------------------- U2 engine

def test_ball_eq3_15_alias_and_rates_bitwise(small):
    s,e,c=small
    for T0 in (1.,1.3,.8):
        s.T=s.T*T0;s.bh_mass*=1.7
        for mode in ('off','bra2008','ball2011_printed','quest2025_eq11','prescribed'):
            cfg=dict(c['engine'],model=mode,rest_mass_rate_msun_per_year=2e-4)
            assert rates(s,e,cfg)==legacy_rates(s,e,cfg)
        alias=rates(s,e,dict(c['engine'],model='ball2012_eq3_15'))
        quest=rates(s,e,dict(c['engine'],model='quest2025_eq11'))
        assert alias==quest and all(type(alias[k]) is type(quest[k]) for k in alias)
        printed=rates(s,e,dict(c['engine'],model='ball2011_printed'))
        assert alias['luminosity']==4*printed['luminosity']

def test_rates_from_complex_step_derivatives():
    rho,p,cs,gamma,M=2.3e-5,4.1e11,1.6e8,1.37,7*MSUN
    cfg=dict(model='ball2012_eq3_15',efficiency=.1,coupling=.1)
    base=rates_from(rho,p,cs,gamma,M,cfg,0.)
    L=base['luminosity']
    assert L==pytest.approx(16*PI*.1*(G*M)**2*rho/(gamma*cs),rel=1e-15)
    h=1e-30
    for i,expected in enumerate((L/rho,0.,-L/cs,-L/gamma,2*L/M)):
        x=[rho,p,cs,gamma,M];x[i]=x[i]+h*1j*abs(x[i])
        d=rates_from(*x,cfg,0.)
        assert d['luminosity'].imag/(h*abs([rho,p,cs,gamma,M][i]))==pytest.approx(expected,rel=1e-13,abs=0)
        assert d['rest_mass_rate'].imag==pytest.approx(d['luminosity'].imag/(.1*C**2),rel=1e-15)
        assert d['bh_growth_rate'].imag==pytest.approx(.9*d['rest_mass_rate'].imag,rel=1e-15)
    d=rates_from(rho,p,cs+1e-20j,gamma,M,cfg,0.)
    assert d['bondi_radius'].imag/1e-20==pytest.approx(-2*base['bondi_radius']/cs,rel=1e-13)
    bra=dict(cfg,model='bra2008');Lb=rates_from(rho,p,cs,gamma,M,bra,0.)['luminosity']
    assert rates_from(rho*(1+1e-20j),p,cs,gamma,M,bra,0.)['luminosity'].imag/(1e-20*rho)==pytest.approx(1.5*Lb/rho,rel=1e-13)
    assert rates_from(rho,p*(1+1e-20j),cs,gamma,M,bra,0.)['luminosity'].imag/(1e-20*p)==pytest.approx(-.5*Lb/p,rel=1e-13)
    with pytest.raises(ValueError):rates_from(rho,p,cs,gamma,M,dict(cfg,model='unknown'),0.)

# ---------------------------------------------------------------- U7 slice, remap, split, merge

def cell_energy(s,e):
    m=s.enclosed;q=e.evaluate(s.rho,s.T)
    return s.dm*(q['u']-G*.5*(m[:-1]/s.r[:-1]+m[1:]/s.r[1:])+.25*(s.v[:-1]**2+s.v[1:]**2))

def test_base_slice_is_exact_and_booked_by_functional_difference(ball):
    s,e=ball;E0=s.energies(e)['total'];scale=s.energies(e)['internal']
    edges=[Fraction(0)]
    for x in s.dm:edges.append(edges[-1]+Fraction(float(x)))
    cut=1.5e-2*MSUN;k=max(i for i in range(len(s.dm)) if edges[i]<=Fraction(cut))
    remaining=float(edges[k+1]-Fraction(cut))
    new=np.r_[remaining,s.dm[k+1:]]
    out,info=remap(s,new,e,offset=cut,copy_unchanged=True,return_slab=True)
    # Exact mass partition (one rounding of the slab mass) and volume to rounding.
    exact=sum(Fraction(float(x)) for x in s.dm)-sum(Fraction(float(x)) for x in new)
    assert info['mass']==float(exact) and abs(info['mass']-cut)<=1e-13*cut
    assert out.cavity_mass==s.cavity_mass+info['mass']
    assert math.fsum(out.volume)+info['volume']==pytest.approx(math.fsum(s.volume),rel=2e-16)
    # Untouched cells are bitwise copies; the cut cell keeps (rho, T).
    for f,g in (('r',slice(1,None)),('v',slice(1,None)),('face_luminosity',slice(1,None)),('radius_compensation',slice(1,None))):
        assert getattr(out,f)[g].tobytes()==getattr(s,f)[k+1:].tobytes(),f
    assert out.T[1:].tobytes()==s.T[k+1:].tobytes() and out.volume[1:].tobytes()==s.volume[k+1:].tobytes()
    assert out.T[0]==s.T[k] and out.rho[0]==pytest.approx(s.rho[k],rel=1e-15)
    assert out.radius_compensation[0]==0.
    np.testing.assert_allclose(out.enclosed[1:],s.enclosed[k+1:],rtol=1e-15)
    geometric=4*PI/3*np.diff(out.r)*(out.r[1:]**2+out.r[1:]*out.r[:-1]+out.r[:-1]**2)
    assert max(abs(geometric/out.volume-1))<1e-12
    assert info['copied_cells']==len(new)-1 and info['fixed_state_cells']==1 and info['reconstructed_cells']==0
    # Booking is the functional difference of State.energies.
    assert abs(out.energies(e)['total']+info['energy']-E0)<=1e-13*scale
    # Independent rational overlap quadrature of the slab.
    f=(Fraction(cut)-edges[k])/(edges[k+1]-edges[k]);Ec=cell_energy(s,e)
    overlap=sum(Fraction(float(x)) for x in Ec[:k])+f*Fraction(float(Ec[k]))
    assert info['overlap_energy']==pytest.approx(float(overlap),rel=1e-14)
    # Functional minus overlap = re-faced potential and kinetic of the fixed-(rho, u) cut cell.
    m=out.enclosed
    kept=out.dm[0]*(-G*.5*(s.enclosed[k]/s.r[k]+s.enclosed[k+1]/s.r[k+1])+.25*(s.v[k]**2+s.v[k+1]**2)
                    +G*.5*(m[0]/out.r[0]+m[1]/out.r[1])-.25*(out.v[0]**2+out.v[1]**2))
    assert abs(info['energy']-info['overlap_energy']-kept)<=1e-13*scale
    # Entropy: the slab carries exactly its overlap entropy; nothing is booked as remapping.
    assert abs(out.ledger['remap_entropy_change']-s.ledger['remap_entropy_change'])<=1e-13*abs(np.dot(s.dm,e.evaluate(s.rho,s.T)['s']))

def test_slice_conservative_reconstruction_and_whole_cell_slice(ball):
    s,e=ball;E0=s.energies(e)['total'];scale=s.energies(e)['internal']
    # Without copy_unchanged every cell is rebuilt; total energy is conserved up to EOS inversion.
    cut=2e-2*MSUN;dm=s.dm[8:].copy();dm[0]=float(sum(Fraction(float(x)) for x in s.dm[:9])-Fraction(cut))
    out,info=remap(s,dm,e,offset=cut,return_slab=True)
    assert info['reconstructed_cells']==len(dm) and info['copied_cells']==0
    assert abs(info['reconstruction_defect'])<=1e-13*scale
    assert abs(out.energies(e)['total']+info['energy']-E0)<=1e-13*scale
    # A slice along an old face copies every remaining cell.
    whole=float(sum(Fraction(float(x)) for x in s.dm[:5]))
    out,info=remap(s,s.dm[5:],e,offset=whole,copy_unchanged=True,return_slab=True)
    assert info['copied_cells']==len(s.dm)-5 and out.r.tobytes()==s.r[5:].tobytes() and out.T.tobytes()==s.T[5:].tobytes()
    Ec=cell_energy(s,e)
    # E(before) - E(after) differences two totals of order E_int: its precision is relative to E_int.
    assert abs(info['energy']-math.fsum(Ec[:5]))<=1e-13*scale
    with pytest.raises(ValueError):remap(s,s.dm,e,offset=1e30,copy_unchanged=True)

def test_zero_offset_is_bitwise_copy_and_returned_gas(ball):
    s,e=ball
    out,info=remap(s,s.dm,e,copy_unchanged=True,return_slab=True)
    for f in ('r','T','dm','v','cell_volume','face_luminosity','radius_compensation','reference_mass_fractions'):
        assert getattr(out,f).tobytes()==getattr(s,f).tobytes(),f
    assert out.ledger==s.ledger and out.cavity_mass==s.cavity_mass and info['energy']==0.
    # Returned gas (negative offset) carries the cell-0 state out of the cavity.
    s=s.copy();s.cavity_mass=1e-3*MSUN
    g=1e-6*MSUN;E0=s.energies(e)['total'];scale=s.energies(e)['internal']
    out,info=remap(s,np.r_[g,s.dm],e,offset=-g,copy_unchanged=True,return_slab=True)
    assert out.cavity_mass==s.cavity_mass-g>0 and info['mass']==-g and out.dm[0]==g
    assert out.T[0]==s.T[0] and out.rho[0]==pytest.approx(s.rho[0],rel=1e-15)
    for f in ('r','T','v','cell_volume','face_luminosity'):
        assert getattr(out,f)[1:].tobytes()==getattr(s,f).tobytes(),f
    np.testing.assert_allclose(out.enclosed[1:],s.enclosed,rtol=1e-15)
    assert abs(out.energies(e)['total']+info['energy']-E0)<=1e-13*scale
    assert abs(info['reconstruction_defect'])<=1e-13*scale and info['energy']<0
    # Returned gas merged into cell 0: conservative reconstruction of that one cell.
    out,info=remap(s,np.r_[g+s.dm[0],s.dm[1:]],e,offset=-g,copy_unchanged=True,return_slab=True)
    assert info['reconstructed_cells']==1 and out.r[1:].tobytes()==s.r[1:].tobytes() and out.T[1:].tobytes()==s.T[1:].tobytes()
    assert out.r[0]<s.r[0] and abs(out.energies(e)['total']+info['energy']-E0)<=1e-13*scale
    assert abs(info['reconstruction_defect'])<=1e-13*scale
    # The cavity cannot return more than it holds, and returned gas must fit inside r0.
    with pytest.raises(ValueError,match='exceeds the cavity'):
        remap(s,np.r_[2e-3*MSUN,s.dm],e,offset=-2e-3*MSUN,copy_unchanged=True)
    with pytest.raises(ValueError,match='exceeds the cavity'):
        remap(ball[0],np.r_[g,ball[0].dm],e,offset=-g,copy_unchanged=True)
    s.cavity_mass=2e3*MSUN
    with pytest.raises(ValueError,match='does not fit'):
        remap(s,np.r_[1e3*MSUN,s.dm],e,offset=-1e3*MSUN,copy_unchanged=True)

def test_copy_unchanged_rezone_conserves_every_changed_cell(ball):
    s,e=ball;E0=s.energies(e)['total'];scale=s.energies(e)['internal']
    before=ledger.residuals(s,e)['energy_residual_erg']
    # Splitting old cells without a slice: no cell keeps a fixed state, nothing is booked.
    dm=np.r_[np.repeat(s.dm[:20]/2,2),s.dm[20:]]
    out,info=remap(s,dm,e,copy_unchanged=True,return_slab=True)
    assert info['fixed_state_cells']==0 and info['reconstructed_cells']==40 and info['copied_cells']==180
    assert abs(out.energies(e)['total']-E0)<=1e-13*scale and abs(info['reconstruction_defect'])<=1e-13*scale
    assert abs(ledger.residuals(out,e)['energy_residual_erg']-before)<=1e-13*scale
    assert out.T[40:].tobytes()==s.T[20:].tobytes() and out.r[40:].tobytes()==s.r[20:].tobytes()
    # A slice that also re-splits the next cell: only the cut remainder keeps (rho, T),
    # so the swept booking is the overlap plus the cut cell's re-faced potential only.
    edges=[Fraction(0)]
    for x in s.dm:edges.append(edges[-1]+Fraction(float(x)))
    cut=1.5e-2*MSUN;k=max(i for i in range(len(s.dm)) if edges[i]<=Fraction(cut))
    new=np.r_[float(edges[k+1]-Fraction(cut)),s.dm[k+1]/2,s.dm[k+1]/2,s.dm[k+2:]]
    out,info=remap(s,new,e,offset=cut,copy_unchanged=True,return_slab=True)
    assert info['fixed_state_cells']==1 and info['reconstructed_cells']==2 and out.T[0]==s.T[k]
    m=out.enclosed
    kept=out.dm[0]*(-G*.5*(s.enclosed[k]/s.r[k]+s.enclosed[k+1]/s.r[k+1])+.25*(s.v[k]**2+s.v[k+1]**2)
                    +G*.5*(m[0]/out.r[0]+m[1]/out.r[1])-.25*(out.v[0]**2+out.v[1]**2))
    assert abs(info['energy']-info['overlap_energy']-kept)<=1e-13*scale
    assert abs(out.energies(e)['total']+info['energy']-E0)<=1e-13*scale

def test_remap_without_kinetic_energy(ball):
    s,e=ball;q=s.copy();q.kinetic_tracked=False
    dm=np.repeat(q.dm/2,2)
    out=remap(q,dm,e)
    assert not out.kinetic_tracked
    assert abs(out.energies(e)['total']-q.energies(e)['total'])<=1e-13*q.energies(e)['internal']

def face_potential(s):
    m=s.enclosed
    return -G*.5*(m[:-1]/s.r[:-1]+m[1:]/s.r[1:])+(.25*(s.v[:-1]**2+s.v[1:]**2) if s.kinetic_tracked else 0.)

def changed_groups(s,out):
    """(old cells, new cells) between consecutive faces that both meshes keep bitwise."""
    common=np.intersect1d(s.r,out.r);i=np.searchsorted(s.r,common);j=np.searchsorted(out.r,common)
    assert i[0]==j[0]==0 and i[-1]==len(s.dm) and j[-1]==len(out.dm)
    return [(slice(i[a],i[a+1]),slice(j[a],j[a+1])) for a in range(len(common)-1) if not i[a+1]-i[a]==j[a+1]-j[a]==1]

def check_events(s,out,info,e):
    """T of every split or merged cell and the entropy decomposition of the event set.

    Each group conserves u + potential; the potential change dPhi comes out of
    u per unit mass, so dS = mixing + sum M [s(u_mix - dPhi/M) - s(u_mix)]
    = mixing - sum dPhi/T - sum dPhi^2/(2 M T^2 cv) + third order.
    """
    q=e.evaluate(s.rho,s.T);old,new=face_potential(s),face_potential(out)
    first=second=mixing=reallocation=0.
    for a,b in changed_groups(s,out):
        assert a.stop-a.start==1 or b.stop-b.start==1  # pure split or pure merge
        M=math.fsum(s.dm[a]);U=math.fsum(s.dm[a]*q['u'][a]);E=math.fsum(s.dm[a]*(q['u'][a]+old[a]))
        Phi=math.fsum(out.dm[b]*new[b]);dPhi=Phi-math.fsum(s.dm[a]*old[a]);reallocation+=dPhi
        rho=out.rho[b]
        np.testing.assert_allclose(out.T[b],e.temperature(rho,np.full(len(rho),(E-Phi)/M)),rtol=1e-14,atol=0)
        T=e.temperature(rho[:1],U/M);f=e.evaluate(rho[:1],T)
        group_mixing=M*f['s'][0]-math.fsum(s.dm[a]*q['s'][a])
        assert group_mixing>=-1e-14*math.fsum(s.dm[a]*abs(q['s'][a]))
        mixing+=group_mixing*(a.stop-a.start>1)
        first-=dPhi/T[0];second-=.5*dPhi**2/(M*T[0]**2*f['cv'][0])
    assert info['merge_mixing_entropy']==pytest.approx(mixing,rel=1e-10,abs=1e-15*np.dot(s.dm,abs(q['s'])))
    assert info['potential_reallocation']==pytest.approx(reallocation,rel=1e-10)
    residual=info['entropy_change']-mixing-first
    assert abs(residual)<=2*abs(second) and abs(residual-second)<=.2*abs(second)
    return first

def test_split_merge_conserves_and_books(ball):
    s,e=ball;scale=s.energies(e)['internal']
    target=q_targets(s,e,400)
    before=ledger.residuals(s,e)['energy_residual_erg']
    out,info=split_merge(s,e,target)
    assert info['splits']>0 and info['merges']==0
    assert math.fsum(out.dm)==math.fsum(s.dm) and math.fsum(out.volume)==pytest.approx(math.fsum(s.volume),rel=2e-16)
    geometric=4*PI/3*np.diff(out.r)*(out.r[1:]**2+out.r[1:]*out.r[:-1]+out.r[:-1]**2)
    assert max(abs(geometric/out.volume-1))<1e-12
    assert abs(info['energy_defect'])<=1e-13*scale
    assert out.ledger['remap_energy_defect']==s.ledger['remap_energy_defect']+info['energy_defect']
    # Per event (one split cell): the EOS-inversion defect, within the plan's 1e-12 of E_int per event.
    assert 0<=info['max_event_defect']<=1e-13*scale
    assert abs(ledger.residuals(out,e)['energy_residual_erg']-before)<=1e-13*scale
    # Kept faces and cells are bitwise; split children share one (rho, T).
    assert out.r[0]==s.r[0] and out.r[-1]==s.r[-1] and set(s.r.tolist())<=set(out.r.tolist())
    kept=np.isin(out.r,s.r);old_index=np.searchsorted(s.r,out.r[kept])
    assert out.v[kept].tobytes()==s.v[old_index].tobytes() and out.face_luminosity[kept].tobytes()==s.face_luminosity[old_index].tobytes()
    assert out.radius_compensation[kept].tobytes()==s.radius_compensation[old_index].tobytes() and not np.any(out.radius_compensation[~kept])
    k=int(np.argmax(s.dm/target));lo=np.searchsorted(out.r,s.r[k]);hi=np.searchsorted(out.r,s.r[k+1])
    assert hi-lo>=2 and len(set(out.T[lo:hi]))==1 and len(set(out.rho[lo:hi]))==1
    assert np.all(out.dm[lo:hi]==s.dm[k]/(hi-lo))
    unchanged=[j for j in range(len(s.dm)) if s.dm[j]<=2*target[j]]
    rows=np.searchsorted(out.r,s.r[unchanged])
    assert out.T[rows].tobytes()==s.T[unchanged].tobytes()
    # Total-energy conservation moves the refined face potential into u: every
    # split cools here (m/r is convex in mass inside a uniform-density cell).
    first=check_events(s,out,info,e)
    assert info['potential_reallocation']>0 and info['entropy_change']<0 and first<0

def test_merge_is_conservative_and_entropy_is_mixing_minus_reallocation(ball):
    s,e=ball;scale=s.energies(e)['internal']
    out,info=split_merge(s,e,s.dm/.45)
    assert info['merges']==len(s.dm)//2 and info['splits']==0 and len(out.dm)==len(s.dm)//2
    assert abs(math.fsum(out.dm)-math.fsum(s.dm))<=info['merges']*np.spacing(2*s.dm.max())
    assert out.r.tobytes()==s.r[::2].tobytes() and out.volume==pytest.approx(s.volume[::2]+s.volume[1::2],rel=1e-15)
    assert abs(info['energy_defect'])<=1e-13*scale
    check_events(s,out,info,e)
    S=float(np.dot(s.dm,abs(e.evaluate(s.rho,s.T)['s'])))
    assert info['merge_mixing_entropy']>=0
    # Distinct neighbours: mixing at fixed summed internal energy and volume strictly produces entropy.
    t=s.copy();t.T=t.T*(1+.05*(-1)**np.arange(len(t.T)))
    merged,mixed=split_merge(t,e,t.dm/.45)
    assert mixed['merge_mixing_entropy']>1e-6*S
    check_events(t,merged,mixed,e)
    # Nothing to do returns an exact copy.
    same,none=split_merge(s,e,s.dm)
    assert none['splits']==none['merges']==0 and same.T.tobytes()==s.T.tobytes() and same.ledger==s.ledger
    assert none['potential_reallocation']==0. and none['max_event_defect']==0.
    assert 0<=info['max_event_defect']<=1e-13*scale

def test_slice_remainder_never_survives_split_merge(ball):
    s,e=ball;scale=s.energies(e)['internal']
    edges=[Fraction(0)]
    for x in s.dm:edges.append(edges[-1]+Fraction(float(x)))
    k=3;cut=float(edges[k+1]-Fraction(float(s.dm[k]))/10**6)
    sliced=remap(s,np.r_[float(edges[k+1]-Fraction(cut)),s.dm[k+1:]],e,offset=cut,copy_unchanged=True)
    mu=lambda st:np.cumsum(st.dm)-.5*st.dm  # mass above the new base
    target=q_targets(s,e,400,mass_above_base=mu(sliced))
    assert sliced.dm[0]<1e-3*target[0] and sliced.dm[1]>2*target[1]
    E0=sliced.energies(e)['total']
    out,info=split_merge(sliced,e,target)
    assert info['merges']==1 and out.dm[0]>=.5*q_targets(s,e,400,mass_above_base=mu(out))[0]
    groups=changed_groups(sliced,out)
    for a,b in groups:assert np.all(out.dm[b]>=.5*target[a].min())
    assert math.fsum(out.dm)==pytest.approx(math.fsum(sliced.dm),rel=1e-15) and out.cavity_mass==sliced.cavity_mass
    assert math.fsum(out.volume)==pytest.approx(math.fsum(sliced.volume),rel=2e-16)
    geometric=4*PI/3*np.diff(out.r)*(out.r[1:]**2+out.r[1:]*out.r[:-1]+out.r[:-1]**2)
    assert max(abs(geometric/out.volume-1))<1e-12
    assert abs(out.energies(e)['total']-E0)<=1e-13*scale and abs(info['energy_defect'])<=1e-13*scale
    # The remainder and cell 1 form one group cut into equal-mass parts with
    # overlap shares of volume and internal energy: parts inside cell 1 keep its
    # density, and dPhi comes out of u per unit mass of the group.
    a,b=groups[0]
    assert a==slice(0,2) and b.start==0 and b.stop>=2 and len(set(out.dm[b]))==1
    np.testing.assert_allclose(out.rho[1:b.stop],sliced.rho[1],rtol=1e-13)
    old,new=face_potential(sliced),face_potential(out);U=sliced.dm[:2]*e.evaluate(sliced.rho,sliced.T)['u'][:2]
    ends=np.cumsum(out.dm[b]);starts=ends-out.dm[b];top=np.cumsum(sliced.dm[:2]);bottom=top-sliced.dm[:2]
    w=np.clip(np.minimum(ends[:,None],top[None,:])-np.maximum(starts[:,None],bottom[None,:]),0,None)/sliced.dm[None,:2]
    share=w@U;M=math.fsum(sliced.dm[:2])
    dPhi=math.fsum(out.dm[b]*new[b])-math.fsum(sliced.dm[:2]*old[:2])
    u=share/out.dm[b]-(share.sum()-math.fsum(U)+dPhi)/M
    np.testing.assert_allclose(out.T[b],e.temperature(out.rho[b],u),rtol=1e-13)

def test_stars_q_matches_thesis_eq_2_19_and_targets(ball):
    s,e=ball;q=e.evaluate(s.rho,s.T);m=.5*(s.enclosed[:-1]+s.enclosed[1:]);M=s.enclosed[-1]
    ours=stars_q(q['p'],s.T,m,s.rc,M)
    c1=9.99;ct1=10**(10*c1);c=STARS_Q
    literal=(c['c4']*np.log(q['p'])+c['c5']*np.log((q['p']+c['c9'])/(q['p']+ct1))+c['c2']*np.log((q['p']+c['c10'])/(q['p']+ct1))
             +c['c7']*np.log(s.T/(s.T+c['ct10']))-np.log((m/M)**(2/3)/c['c6']+1)-c['c3']*np.log(s.rc**2/c['c8']+1))
    np.testing.assert_allclose(np.diff(ours),np.diff(literal),rtol=1e-12,atol=1e-12*np.ptp(ours))
    assert np.all(np.diff(ours)<0)
    assert (c['c2'],c['c3'],c['c4'],c['c5'],c['c6'],c['c7'],c['c8'],c['c9'],c['c10'],c['ct10'])==(0.,.05,.5,.15,.02,.45,1e18,1e15,3e19,2e4)
    for n in (200,400,800):
        t=q_targets(s,e,n)
        assert np.all(t>0) and np.sum(s.dm/t)==pytest.approx(n,rel=1e-12)
    mid=np.cumsum(s.dm)-.5*s.dm
    np.testing.assert_allclose(q_targets(s,e,400,mass_above_base=mid),q_targets(s,e,400),rtol=1e-13)

# ---------------------------------------------------------------- U9 coloured Jacobian

def test_pattern_jacobian_equals_dense_complex_step():
    N=16;nr,nT,nL=N+1,N,N-1;n=nr+nT+nL+1
    R=lambda j:j;T=lambda k:nr+k;L=lambda j:nr+nT+j-1;X=n-1
    def fun(x):
        r=np.exp(x[:nr]);t=np.exp(x[nr:nr+nT]);l=np.r_[np.exp(x[X])*2.,x[nr+nT:nr+nT+nL],r[-1]**2*t[-1]**4*1e-3];mass=1+np.exp(x[X])+np.cumsum(np.r_[0.,np.full(N,.1)])
        momentum=[(t[j-1]*r[j-1]**-3-(t[j] if j<N else .5)*r[j]**-3)*r[j]**2+mass[j]/r[j]**2*np.sin(r[min(j+1,N)]) for j in range(1,N+1)]
        energy=[t[k]**1.5*r[k+1]-np.log(t[k])*r[k]+l[k+1]-l[k] for k in range(N)]
        closure=[l[j]-(t[j-1]-t[j])*r[j]**2*np.exp(.3*x[X]) for j in range(1,N)]
        piston=[np.log(r[0])-np.log(mass[0]/(t[0]+.2*t[1]))]
        growth=[np.exp(x[X])-1.1-.01*t[0]*r[0]]
        return np.array(momentum+energy+closure+piston+growth)
    rows=[];cols=[]
    def add(i,js):rows.extend([i]*len(js));cols.extend(js)
    for j in range(1,N+1):add(j-1,[R(j-1),R(j),R(min(j+1,N)),T(j-1)]+([T(j)] if j<N else [])+[X])
    for k in range(N):add(N+k,[R(k),R(k+1),T(k)]+[L(i) for i in (k,k+1) if 1<=i<=N-1]+([X] if k==0 else [])+([R(N),T(N-1)] if k==N-1 else []))
    for j in range(1,N):add(2*N+j-1,[L(j),T(j-1),T(j),R(j),X])
    add(3*N-1,[R(0),T(0),T(1),X]);add(3*N,[X,T(0),R(0)])
    from scipy.sparse import csc_matrix
    pattern=csc_matrix((np.ones(len(rows)),(rows,cols)),shape=(n,n))
    x=np.random.default_rng(3).normal(0,.1,n)
    dense=np.zeros((n,n))
    for j in range(n):
        xx=x.astype(complex);xx[j]+=1e-25j;dense[:,j]=fun(xx).imag/1e-25
    J=pattern_jacobian(fun,x,pattern,strict=True).toarray()
    assert np.max(abs(J-dense))<=1e-13*np.max(abs(dense))
    assert not np.any(dense[pattern.toarray()==0])
    colours=colouring(pattern)
    assert colouring(pattern) is colours and colours.max()+1<=12 and not colours.flags.writeable
    for c in range(colours.max()+1):
        assert np.max((pattern[:,colours==c]!=0).sum(axis=1))<=1
    # An undeclared dependence is detected in strict mode and dropped otherwise.
    missing=pattern.tolil();missing[3*N,R(0)]=0;missing=missing.tocsc()
    with pytest.raises(ValueError,match='outside the declared pattern'):pattern_jacobian(fun,x,missing,strict=True)
    assert pattern_jacobian(fun,x,missing).nnz==pattern.nnz-1

# ---------------------------------------------------------------- U15 Ball seed geometry

@pytest.mark.parametrize('zones',[200,400])
def test_ball_geometry_initializes(zones):
    c=ball_config(zones);e=IdealGasRadiation(Composition(.7,.3,0.))
    s,err=structure.initialize(c,e)
    assert err<1e-10 and np.all(s.dm>0) and np.all(np.diff(s.r)>0)
    assert max(abs(structure.hydrostatic_residual(s,e,c['boundary']['pressure_cgs'])))<1e-9
    assert s.r[0]==c['initial']['inner_radius_cm'] and s.r[-1]==pytest.approx(c['initial']['radius_cm'],rel=1e-12)

def test_monotone_seed_table_untouched_and_small_meshes_bypass(monkeypatch):
    import scipy.interpolate
    seen=[];original=scipy.interpolate.PchipInterpolator
    def spy(x,y,*a,**k):
        seen.append(len(x));return original(x,y,*a,**k)
    monkeypatch.setattr(scipy.interpolate,'PchipInterpolator',spy)
    e=IdealGasRadiation(Composition(.7,.3,0.))
    c=resolve({'numerics':{'zones':256,'mesh':'log_radius_cosine'},'boundary':{'pressure_cgs':1e8}})
    structure.initialize(c,e)
    assert seen==[30001]  # strictly increasing table: every BVP node used
    seen.clear();structure.initialize(ball_config(200),e)
    assert seen and seen[0]<30001  # Ball geometry: overshooting surface nodes dropped
    def forbidden(*a,**k):raise AssertionError('96 zones must use the dense initializer')
    monkeypatch.setattr(structure,'sparse_polytrope',forbidden)
    s,err=structure.initialize(ball_config(96),e)
    assert err<1e-8 and len(s.dm)==96
