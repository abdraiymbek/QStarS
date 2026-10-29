"""Ball 2011 replication, Stage 2: ML1 and thesis-lossless mixing-length closures.

The MESA comparison needs QSTARS_MESA_LIBRARY (bridge ABI >= 5) and MESA_DIR and
skips otherwise. Every other test runs in the default suite with no compiled code;
the cubic is checked against an independent 60-digit decimal reference instead.
"""
import os
from decimal import Decimal,getcontext
import numpy as np
import pytest
from qstars import mlt,numerics,nonlinear,ledger,transport
from qstars.constants import A as ARAD,C,PI,SIGMA,YEAR
from qstars.config import resolve
from qstars.evolution import physics
from qstars.structure import initialize
from qstars.transport import luminosities,inverse_transport_residual

EPS=np.finfo(float).eps
NAMES=('grad_rad','grad_ad','cp','Q','p','rho','T','kappa','g')


def reference_efficiency(delta):
    """60-digit positive root of a0 G^3 + G^2 + G = delta; monotone Newton from an upper bound."""
    getcontext().prec=60
    d=delta if isinstance(delta,Decimal) else Decimal(float(delta));a0=Decimal(9)/4
    g=min(d,(d/a0)**(Decimal(1)/3))
    for _ in range(400):
        step=(((a0*g+1)*g+1)*g-d)/((3*a0*g+2)*g+1);g-=step
        if abs(step)<=Decimal(10)**-50*g:return g
    raise RuntimeError('reference cubic did not converge')


def face_states(rng,n,log_ratio=(-3.,3.),log_A=(-2.,4.)):
    """Synthetic convective faces with prescribed nabla_rad/nabla_ad - 1 and A."""
    s=dict(grad_ad=rng.uniform(.1,.4,n),cp=10**rng.uniform(7,9,n),Q=rng.uniform(1,8,n),
           p=10**rng.uniform(2,10,n),rho=10**rng.uniform(-12,-5,n),T=10**rng.uniform(3.7,6.5,n),
           g=10**rng.uniform(-1,5,n))
    s['grad_rad']=s['grad_ad']*(1+10**rng.uniform(*log_ratio,n))
    unit=mlt.ml1_coefficient(s['cp'],s['Q'],s['p'],s['rho'],s['T'],1.,s['g'])
    s['kappa']=10**rng.uniform(*log_A,n)/unit  # A is linear in kappa
    return s


def test_ml1_cubic_and_excess_against_60_digit_reference():
    rng=np.random.default_rng(2011)
    delta=np.r_[10**rng.uniform(-30,30,600),1e-300,1e-30,1e-16,1e-8,1.,1e8,1e100,1e200]
    reference=[reference_efficiency(d) for d in delta]
    gamma=mlt._ml1_efficiency(delta)
    np.testing.assert_allclose(gamma,[float(g) for g in reference],rtol=4*EPS,atol=0)
    # Y = G(G+1)/A^2 is the exact excess (1 - zeta)(nabla_rad - nabla_ad), zeta = G^3/B^3.
    s=face_states(rng,300,log_ratio=(-6.,4.),log_A=(-4.,6.))
    out=mlt.ml1(**s)
    A=out['A'];D=s['grad_rad']-s['grad_ad']
    getcontext().prec=60
    for k in range(len(A)):
        a,d=Decimal(float(A[k])),Decimal(float(D[k]))
        g=reference_efficiency(a*a*d)
        identity=g*(g+1)/(a*a);zeta=Decimal(9)/4*g**3/(a*a*d)
        assert abs(identity-(1-zeta)*d)<=Decimal(10)**-40*d
        assert abs(float(out['excess'][k])/float(identity)-1)<=16*EPS
    # Away from cancellation the MESA form (1-zeta) nabla_rad + zeta nabla_ad agrees.
    zeta=mlt.A0*out['efficiency']**3/(A**2*D)
    ok=(zeta>1e-3)&(zeta<.5)
    assert ok.sum()>30
    np.testing.assert_allclose(out['gradient'][ok],((1-zeta)*s['grad_rad']+zeta*s['grad_ad'])[ok],rtol=1e-13)


def test_ml1_limits_and_schwarzschild_boundary():
    base=dict(grad_ad=.25,cp=3e8,Q=4.,p=1e6,rho=1e-9,T=2e4,g=1e3)
    unit=float(mlt.ml1_coefficient(base['cp'],base['Q'],base['p'],base['rho'],base['T'],1.,base['g']))
    D=.1
    # Efficient limit: nabla -> nabla_ad, Y (a0 A)^(2/3) / D^(2/3) = 1 + O(1/B).
    previous=np.inf
    for A in (1e3,1e5,1e7):
        out=mlt.ml1(grad_rad=.25+D,kappa=A/unit,**base)
        B=(A*A*D/mlt.A0)**(1/3);error=abs(out['excess']*(mlt.A0*A)**(2/3)/D**(2/3)-1)
        assert error<1/B and error<previous;previous=error
    # Radiative (inefficient) limit: nabla -> nabla_rad with 1 - Y/D = a0 A^4 D^2 (1 + O(A^2 D)).
    out=mlt.ml1(grad_rad=.25+D,kappa=1e-2/unit,**base)
    assert abs((1-out['excess']/D)/(mlt.A0*1e-8*D**2)-1)<1e-3
    assert abs(out['velocity']/(2*np.sqrt(base['Q']*base['p']/(8*base['rho']))*1e-2*D)-1)<1e-3
    # Stable faces return nabla_rad exactly with zero velocity and efficiency, like MESA.
    grad_rad=np.array([-.3,0.,.1,.25])
    out=mlt.ml1(grad_rad=grad_rad,kappa=1./unit,**base)
    np.testing.assert_array_equal(out['gradient'],grad_rad)
    np.testing.assert_array_equal(out['excess'],grad_rad-.25)
    assert not np.any(out['velocity']) and not np.any(out['efficiency']) and not np.any(out['convective'])
    # C1 across the boundary: Y -> D and dY/dnabla_rad -> 1 from both sides.
    for A in (1e-2,1.,1e3):
        for side in (-1,1):
            x=.25+side*1e-13
            y=mlt.ml1(grad_rad=x+1e-30j,kappa=A/unit,**base)['excess']
            assert abs(y.real/(x-.25)-1)<1e-6 and abs(y.imag/1e-30-1)<1e-6


def test_thesis_lossless_limits_and_cubic():
    base=dict(grad_ad=.25,cp=3e8,p=1e6,rho=1e-9,T=2e4,g=1e3)
    unit=float(mlt.thesis_coefficient(base['cp'],base['p'],base['rho'],base['T'],1.,base['g']))
    grad_rad=.25+np.array([1e-12,1e-6,1e-2,1.,1e2]);D=grad_rad-.25
    for U in (1e-6,1.,1e6):
        out=mlt.thesis_lossless(grad_rad=grad_rad,kappa=U/unit,**base)
        Y=out['excess'];y=np.sqrt(Y)
        np.testing.assert_allclose(U*Y*y+Y,D,rtol=8*EPS)
        hp=base['p']/(base['rho']*base['g'])
        np.testing.assert_allclose(out['velocity'],2*np.sqrt(base['g']*hp*Y/2),rtol=4*EPS)
    efficient=mlt.thesis_lossless(grad_rad=.35,kappa=1e12/unit,**base)['excess']
    assert abs(efficient/(.1/1e12)**(2/3)-1)<1e-3
    lossless=mlt.thesis_lossless(grad_rad=.35,kappa=1e-12/unit,**base)['excess']
    assert abs(lossless/.1-1)<1e-5
    out=mlt.thesis_lossless(grad_rad=np.array([-.1,.25]),kappa=1./unit,**base)
    np.testing.assert_array_equal(out['gradient'],[-.1,.25]);assert not np.any(out['velocity'])


@pytest.mark.parametrize('model',mlt.MODELS)
def test_forward_luminosity_inverts_closure(model):
    rng=np.random.default_rng(7)
    s=face_states(rng,400,log_ratio=(-6.,4.),log_A=(-3.,5.))
    s['grad_rad'][:40]=s['grad_ad'][:40]*rng.uniform(-1,.9,40)  # radiative faces
    mass=10**rng.uniform(33,38,400)
    closure=mlt.closure(model,**s)
    forward=mlt.forward_luminosity(closure['excess'],**{k:s[k] for k in NAMES[1:]},mass=mass,model=model)
    L=mlt.conductance(s['kappa'],s['p'],s['T'],mass)*s['grad_rad']
    K=mlt.conductance(s['kappa'],s['p'],s['T'],mass);floor=8*EPS*s['grad_ad']
    assert np.all(abs(forward['grad_rad']-s['grad_rad'])<=1e-13*abs(s['grad_rad'])+floor)
    assert np.all(abs(forward['luminosity']-L)<=1e-13*abs(L)+K*floor)
    np.testing.assert_allclose(forward['radiative']+forward['convective'],forward['luminosity'],rtol=1e-15)
    np.testing.assert_allclose(forward['velocity'],closure['velocity'],rtol=1e-12)
    assert np.all(forward['convective'][:40]==0) and np.all(forward['convective'][40:]>0)
    back=mlt.closure(model,mlt.radiative_gradient(forward['luminosity'],s['kappa'],s['p'],s['T'],mass),
                     **{k:s[k] for k in NAMES[1:]})
    assert np.all(abs(back['excess']-closure['excess'])<=1e-12*abs(closure['excess'])+floor)
    with pytest.raises(ValueError):mlt.closure('mlt_unknown',**s)
    with pytest.raises(ValueError):mlt.forward_luminosity(closure['excess'],**{k:s[k] for k in NAMES[1:]},mass=mass,model='x')


@pytest.mark.parametrize('model',mlt.MODELS)
def test_complex_step_derivatives_against_finite_differences(model):
    rng=np.random.default_rng(11)
    s=face_states(rng,60,log_ratio=(-2.,3.),log_A=(-1.,3.))
    mass=10**rng.uniform(33,38,60)
    def closure(x):return mlt.closure(model,**x)
    def forward(x):
        m=x.pop('mass');Y=x.pop('grad_rad')
        return mlt.forward_luminosity(Y,**x,mass=m,model=model)
    Y=closure(s)['excess'].real
    cases=((closure,dict(s),('gradient','excess','velocity')),
           (forward,dict(s,grad_rad=Y,mass=mass),('luminosity','radiative','convective','velocity')))
    for fun,x,keys in cases:
        value=fun(dict(x))
        for name in x:
            z=dict(x);z[name]=x[name]*np.exp(1e-30j)
            plus=dict(x);plus[name]=x[name]*np.exp(1e-5)
            minus=dict(x);minus[name]=x[name]*np.exp(-1e-5)
            complex_step,up,down=fun(z),fun(plus),fun(minus)
            for key in keys:
                exact=complex_step[key].imag/1e-30;finite=(up[key]-down[key])/2e-5
                assert np.all(abs(exact-finite)<=1e-6*abs(exact)+1e-9*abs(value[key])),(name,key)


def confined(n=32):
    c=resolve({'numerics':{'zones':n,'mesh':'core_refined','solver_tolerance':1e-10},
               'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},
               'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4}})
    e,op=physics(c);s,_=initialize(c,e);ledger.initialize(s,e)
    return s,e,op,c


@pytest.mark.parametrize('model',mlt.MODELS)
def test_manufactured_inverse_on_a_mixed_profile(model):
    s,e,op,c=confined();c['transport']['convection']=model
    x=np.linspace(0,1,len(s.dm))
    for amplitude in (.002,.01,.05):
        st=s.copy();st.T=st.T*np.exp(amplitude*np.sin(9*x))
        L,d=luminosities(st,e,op,c,0.)
        convective=d['convective'][1:-1]>0
        assert 0<convective.sum()<len(convective)
        assert np.all(d['convective'][1:-1][~convective]==0) and np.all(d['convective_mach'][convective]>0)
        residual=inverse_transport_residual(st,e,op,c,L)
        q=e.evaluate(st.rho,st.T)
        np.testing.assert_allclose(residual,0.,atol=2e-14)
        assert np.all(abs(residual)<=1e-10*abs(e.thermal_contrast(st.rho,st.T,q)))
        np.testing.assert_array_equal(L[1:-1][~convective],d['radiative'][1:-1][~convective])
    st.T=st.T[::-1].copy()  # pressure no longer decreases outward
    with pytest.raises(ValueError):luminosities(st,e,op,c,0.)


@pytest.mark.parametrize('model',mlt.MODELS)
def test_colored_flux_jacobian_matches_full_complex_step(monkeypatch,model):
    s,e,op,c=confined(8);c['transport']['convection']=model
    c['numerics'].update(nonlinear_solver='newton_flux',mechanics='backward_euler')
    s.T=s.T*np.exp(.03*np.sin(9*np.linspace(0,1,8)))
    _,d=luminosities(s,e,op,c,0.);assert np.any(d['convective'][1:-1]>0)
    original=numerics.newton;checked=[]
    def audit(fun,x,tolerance,**kwargs):
        jac=kwargs['jacobian'](fun,x).toarray();dense=np.empty_like(jac)
        for k in range(len(x)):
            xx=x.astype(complex);xx[k]+=1e-25j;dense[:,k]=fun(xx).imag/1e-25
        np.testing.assert_allclose(jac,dense,rtol=2e-12,atol=1e-12)
        checked.append(True)
        return original(fun,x,tolerance,**kwargs)
    monkeypatch.setattr(numerics,'newton',audit)
    try:numerics.thermal_mechanical_step(s,1e-4*YEAR,e,op,c,1e42)
    except numerics.StepRejected:pass
    assert checked


def test_stars_photosphere_thermal_branch():
    s,e,op,c=confined(16)
    for complex_step in (False,True):
        st=s.copy()
        if complex_step:st.T=st.T.astype(complex);st.T[-1]*=np.exp(1e-25j)
        c['boundary'].update(thermal='grey',external_optical_depth=2/3,bath_temperature_K=0.)
        grey,_=luminosities(st,e,op,c,0.)
        c['boundary'].update(thermal='stars_photosphere',external_optical_depth=0.,bath_temperature_K=5e3)
        stars,_=luminosities(st,e,op,c,0.)
        assert stars[-1].tobytes()==grey[-1].tobytes()
    # Eddington grey: T_c^4 = 3/4 Teff^4 (tau_half + 4/3), L = 4 pi R^2 sigma Teff^4.
    teff4=stars[-1].real/(4*PI*s.r[-1]**2*SIGMA)
    tau_half=.5*op(s.rho,s.T)[-1]*s.dm[-1]/(4*PI*s.r[-1]**2)
    assert abs(.75*teff4*(tau_half+4/3)/s.T[-1]**4-1)<1e-14


def test_default_closures_never_reach_mlt(monkeypatch):
    s,e,op,c=confined(16)
    expected={}
    for convection in ('mlt_lossless','off'):
        c['transport']['convection']=convection
        L,_=luminosities(s,e,op,c,0.)
        expected[convection]=(L,inverse_transport_residual(s,e,op,c,L))
    def forbidden(*args,**kwargs):raise AssertionError('default transport called mlt.py')
    for name in ('closure','forward_luminosity','radiative_gradient','_mixing_length_faces'):
        monkeypatch.setattr(transport,name,forbidden)
    for convection in ('mlt_lossless','off'):
        c['transport']['convection']=convection
        L,_=luminosities(s,e,op,c,0.)
        assert L.tobytes()==expected[convection][0].tobytes()
        assert inverse_transport_residual(s,e,op,c,L).tobytes()==expected[convection][1].tobytes()


@pytest.fixture(scope='module')
def mesa():
    library=os.environ.get('QSTARS_MESA_LIBRARY');directory=os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Set QSTARS_MESA_LIBRARY and MESA_DIR for compiled MESA tests')
    from qstars.eos import Composition
    from qstars.mesa_microphysics import MesaMicrophysics
    m=MesaMicrophysics(Composition(.7,.3,0.),library,directory)
    if m.abi_version<5:pytest.skip('Native MESA MLT requires bridge ABI 5')
    return m


def test_ml1_matches_mesa_native(mesa):
    """1000 MESA EOS states, nabla_rad/nabla_ad in [1+1e-6, 1e4], omega in [1e-4, 1e4]."""
    rng=np.random.default_rng(1);n=1000;alpha=2.
    logT=rng.uniform(4.,7.,n);T=10**logT;rho=10**(rng.uniform(-6.,-1.,n)+3*logT-18)
    q=mesa.evaluate(rho,T)
    omega=10**rng.uniform(-4,4,n);omega[:2]=[1e-4,1e4]
    ratio=1+10**rng.uniform(-6,4,n);ratio[2:4]=[1+1e-6,1e4]
    g=10**rng.uniform(-1,5,n);kappa=omega*g/(alpha*q['p'])
    X=np.array([q['chiT'],q['chiRho'],q['cp'],g,rho,q['p'],T,kappa,ratio*q['grad_ad'],q['grad_ad']])
    def ours(x):
        chiT,chiRho,cp,gg,rr,p,TT,kap,gr,ga=x
        out=mlt.ml1(gr,ga,cp,chiT/chiRho,p,rr,TT,kap,gg,alpha)
        return np.array([out['gradient'],out['velocity'],out['efficiency']])
    value,native=ours(X),mesa.mlt(X,'ML1',alpha)
    # MESA's closed-form Gamma cancels at small Gamma (absolute error ~ eps*scale);
    # it even returns Gamma<=0, v=0 below ~1e-16. Ours is polished (reference test).
    delta=value[2]*(value[2]+1)+mlt.A0*value[2]**3
    f=-2+9*mlt.A0+27*mlt.A0**2*delta;f1=np.cbrt(f+np.sqrt(f*f+4*(3*mlt.A0-1)**3))
    scale=(np.cbrt(4.)*f1+abs(2*np.cbrt(2.)*(1-3*mlt.A0)/f1)+2)/(6*mlt.A0)
    well=8*EPS*scale<=1e-11*value[2]
    assert well.sum()>=300 and (~well).sum()>=300
    # So v (spec: 1e-10 on all 1e3 states) is closed against 60-digit references on
    # every state: Gamma from our A and D, and v from C&G 14.81/14.86b on the inputs.
    out=mlt.ml1(X[8],X[9],X[2],X[0]/X[1],X[5],X[4],X[6],X[7],X[3],alpha);D=X[8]-X[9]
    getcontext().prec=60
    for k in range(n):
        a,d=Decimal(float(out['A'][k])),Decimal(float(D[k]))
        assert abs(value[2][k]/float(reference_efficiency(a*a*d))-1)<=4*EPS
        chiT,chiRho,cp,gg,rr,p,TT,kap,gr,ga=(Decimal(float(x)) for x in X[:,k])
        Q,al=chiT/chiRho,Decimal(alpha)
        a=4*cp*(p*Q*rr/8).sqrt()*al*(al*p*kap/gg)/(24*Decimal(ARAD)*Decimal(C)*TT**3)
        v=al*(Q*p/(8*rr)).sqrt()*reference_efficiency(a*a*(gr-ga))/a
        assert abs(value[1][k]/float(v)-1)<=1e-10
    np.testing.assert_allclose(value[0],native[0],rtol=1e-10,atol=0)
    np.testing.assert_allclose(value[1:,well],native[1:3,well],rtol=1e-10,atol=0)
    # The residual difference is MESA's radiation constant (3e-11 relative): with A
    # matched through kappa, Gamma agrees to MESA's closed-form rounding everywhere.
    matched=X.copy();matched[7]=kappa*mesa.radiation_constant/ARAD
    native=mesa.mlt(matched,'ML1',alpha)
    assert np.all(abs(value[2]-native[2])<=16*EPS*scale+1e-12*value[2])
    np.testing.assert_allclose(value[0],native[0],rtol=1e-11,atol=0)
    # Complex-step derivatives against MESA's auto_diff on well-conditioned states;
    # the floor covers cancellation in derivatives that are O(Gamma) of the value.
    for k in range(10):
        z=matched.astype(complex);z[k]*=np.exp(1e-30j)
        mine=ours(np.where(np.arange(10)[:,None]==7,z*ARAD/mesa.radiation_constant,z)).imag[:2]/1e-30
        theirs=mesa.mlt(z,'ML1',alpha).imag[:2]/1e-30
        assert np.all((abs(mine-theirs)<=1e-9*abs(theirs)+1e-10*abs(value[:2]))[:,well]),k
    stable=X.copy();stable[8]=X[9]*np.linspace(-1,1,n)
    out,native=ours(stable),mesa.mlt(stable,'ML1',alpha)
    np.testing.assert_array_equal(out[0],native[0]);np.testing.assert_array_equal(out[1],native[1])
