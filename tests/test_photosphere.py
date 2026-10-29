"""Gate 6 Stage 4, U5: the Ball thesis 2.15-2.16 photosphere (photosphere.py).

Constant-kappa Eddington grey atmosphere: the target gas pressure is exact,
Gamma_ph >= 1 is a model domain exit, and the discrete outer cell (endpoint
momentum row against p_ext, grey relation for T_N-1) converges to the
spherical grey atmosphere with order >= 1. Analytic microphysics only.
"""
import numpy as np
import pytest
from scipy.integrate import quad,solve_ivp
from qstars import photosphere
from qstars.config import resolve
from qstars.constants import A,C,G,PI,MSUN,RSUN,SIGMA
from qstars.eos import IdealGasRadiation,Composition
from qstars.mass_transport import DomainError
from qstars.opacity import ElectronScattering
from qstars.state import State
from qstars.transport import luminosities

EOS=IdealGasRadiation(Composition(.7,.3,0.))
OPACITY=ElectronScattering(EOS.composition)
KAPPA=OPACITY.kappa

@pytest.mark.parametrize('gamma',[.05,.5,.9,.99])
def test_constant_kappa_grey_atmosphere_gas_pressure_is_exact(gamma):
    M,R=100*MSUN,150*RSUN
    L=gamma*4*PI*G*C*M/KAPPA
    g=G*M/R**2;teff=photosphere.effective_temperature(L,R)
    assert photosphere.eddington_factor(L,M,KAPPA)==pytest.approx(gamma,rel=1e-14)
    assert 4*PI*R**2*SIGMA*teff**4==pytest.approx(L,rel=1e-14)
    target=photosphere.gas_pressure(g,KAPPA,teff)
    # Thesis 2.16 in both printed forms.
    assert target==pytest.approx(2*g/(3*KAPPA)*(1-gamma),rel=1e-12)
    # dP_g/dtau = g/kappa - dP_rad/dtau with P_rad = a T^4/3 and T^4 = 3/4 Teff^4 (tau + 2/3).
    integral,_=quad(lambda tau:g/KAPPA-A*.75*teff**4/3,0.,2/3,epsabs=0,epsrel=1e-13)
    assert integral==pytest.approx(target,rel=1e-12)
    assert (.75*teff**4*(2/3+2/3))**.25==pytest.approx(teff,rel=1e-15)
    assert photosphere.external_pressure(target,teff)==target+A*teff**4/3

def test_photosphere_row_solution_and_complex_step():
    M,R,L=100*MSUN,150*RSUN,5e39
    rho=photosphere.solve(L,R,M,EOS,OPACITY)
    out=photosphere.evaluate(rho,L,R,M,EOS,OPACITY)
    assert abs(out['residual'])<=1e-13
    pg=EOS.evaluate(np.atleast_1d(rho),np.atleast_1d(out['teff']))['pg'][0]
    assert out['pressure']==pytest.approx(pg+A*out['teff']**4/3,rel=1e-15)
    assert out['gas_pressure']==pytest.approx(out['target'],rel=1e-12)
    assert photosphere.implied_column(R,out['pressure'],out['gravity'])==pytest.approx(4*PI*R**2*out['pressure']*R**2/(G*M),rel=1e-14)
    # Complex-step derivatives through every input agree with central differences.
    base=(rho,L,R,M)
    for k in range(4):
        h=1e-30*base[k];x=list(base);x[k]=x[k]+1j*h
        exact=photosphere.evaluate(*x,EOS,OPACITY)['residual'].imag/h
        step=1e-6*base[k];up=list(base);down=list(base);up[k]+=step;down[k]-=step
        fd=(photosphere.evaluate(*up,EOS,OPACITY)['residual']-photosphere.evaluate(*down,EOS,OPACITY)['residual'])/(2*step)
        assert exact==pytest.approx(fd,rel=1e-6,abs=1e-9*abs(exact)+1e-300)
    # The row increases with rho_ph and has a single root on a wide bracket.
    grid=rho*np.logspace(-3,3,61)
    values=[photosphere.evaluate(x,L,R,M,EOS,OPACITY)['residual'] for x in grid]
    assert np.all(np.diff(values)>0)

def test_eddington_limit_is_a_model_domain_exit():
    M,R=100*MSUN,150*RSUN
    edd=4*PI*G*C*M/KAPPA
    assert photosphere.solve(edd*(1-1e-6),R,M,EOS,OPACITY)>0
    for L in (edd*(1+1e-9),1.5*edd):
        with pytest.raises(DomainError,match='Gamma_ph >= 1'):photosphere.solve(L,R,M,EOS,OPACITY)
        with pytest.raises(DomainError,match='Gamma_ph >= 1'):photosphere.evaluate(1e-12,L,R,M,EOS,OPACITY)

def rising_opacity(rho1):
    """kappa_es (1 + (rho/rho1)^2): Gamma_ph depends on the trial photospheric density, as with real tables."""
    return lambda rho,T:KAPPA*(1+(rho/rho1)**2)+0*T

def test_a_trial_density_with_gamma_ph_above_one_is_not_a_missing_photosphere():
    """Final review item 2: kappa_ph depends on the trial rho_ph, so Gamma_ph >= 1 at a trial density says
    nothing about the root. Here the root has Gamma_ph = 0.92 while the third point of the bracket scan from
    the default guess (rho = 1e-12 e^3.5) has Gamma_ph > 1. The registered solve (bracket 'doubling', the
    default) lets that trial's DomainError out ("Gamma_ph >= 1; no Ball surface exists"): that is the
    registered behaviour and stays the default. The option bracket='scan' (plan amendment 4) uses that the
    row is continuous through Gamma_ph = 1, so such a trial is an upper bound of the bracket."""
    M,R=100*MSUN,150*RSUN;opacity=rising_opacity(1e-11)
    L=.7*4*PI*G*C*M/KAPPA
    scan=1e-12*np.exp([0.,.5,1.5,3.5])
    assert [photosphere.evaluate(x,L,R,M,EOS,opacity)['residual']<0 for x in scan[:3]]==[True]*3
    with pytest.raises(DomainError,match='Gamma_ph >= 1'):photosphere.evaluate(scan[3],L,R,M,EOS,opacity)
    # Registered default: the false 'Gamma_ph >= 1' (no surface reported although a root exists).
    assert photosphere.BRACKETS==('doubling','scan')
    for kwargs in (dict(),dict(bracket='doubling')):
        with pytest.raises(DomainError,match='\\(Gamma_ph >= 1\\); no Ball surface exists'):photosphere.solve(L,R,M,EOS,opacity,**kwargs)
    with pytest.raises(ValueError,match='bracket must be one of'):photosphere.solve(L,R,M,EOS,opacity,bracket='bisect')
    rho=photosphere.solve(L,R,M,EOS,opacity,bracket='scan')
    out=photosphere.evaluate(rho,L,R,M,EOS,opacity)
    assert abs(out['residual'])<=1e-13 and scan[2]<rho<scan[3] and out['eddington_factor']==pytest.approx(.9204,abs=1e-4)
    # The same root from guesses on either side, including one inside the Gamma_ph >= 1 range.
    for guess in (1e-16,1e-13,rho,1e-10,1e-6):
        assert photosphere.solve(L,R,M,EOS,opacity,guess=guess,bracket='scan')==pytest.approx(rho,rel=1e-12)
    # Where the doubling scan brackets the root, 'scan' gives exactly the registered root (same iterates).
    for L2 in (5e39,.3*4*PI*G*C*M/KAPPA):
        assert photosphere.solve(L2,R,M,EOS,OPACITY,bracket='scan')==photosphere.solve(L2,R,M,EOS,OPACITY)
    wide=rising_opacity(1e-8)
    for L2 in (.3*L,L):
        assert photosphere.solve(L2,R,M,EOS,wide,bracket='scan')==photosphere.solve(L2,R,M,EOS,wide)

def test_gamma_ph_above_one_at_every_admissible_density_is_the_domain_exit():
    """A true Gamma_ph >= 1 case with a density-dependent opacity: kappa_ph >= kappa_es at every density,
    so L = 1.2 L_Edd(kappa_es) has no Ball surface. With a bounded microphysics domain, a root that lies
    outside it is 'not bracketed' as long as some admissible density has Gamma_ph < 1."""
    M,R=100*MSUN,150*RSUN;opacity=rising_opacity(1e-11);edd=4*PI*G*C*M/KAPPA
    for guess in (None,1e-16,1e-8):
        with pytest.raises(DomainError,match='Gamma_ph >= 1.*at every admissible'):
            photosphere.solve(1.2*edd,R,M,EOS,opacity,guess=guess,bracket='scan')
        with pytest.raises(DomainError,match='Gamma_ph >= 1'):  # the registered solve reports it too
            photosphere.solve(1.2*edd,R,M,EOS,opacity,guess=guess)
    solve=lambda *a,**k:photosphere.solve(*a,**k,bracket='scan')
    class Bounded:
        """The analytic EOS with a density floor, like the audited MESA query rectangle."""
        def __init__(self,floor):self.floor=floor
        def evaluate(self,rho,T):
            if np.any(np.real(rho)<self.floor):raise ValueError('outside the query rectangle')
            return EOS.evaluate(rho,T)
    root=solve(.7*edd,R,M,EOS,opacity)
    assert solve(.7*edd,R,M,Bounded(root/30),opacity)==pytest.approx(root,rel=1e-12)
    # Gamma_ph = 0.7 (1 + (rho/1e-11)^2) reaches 1 at 1.17 x the root. A floor between the two leaves
    # admissible densities with Gamma_ph < 1 but no root (it lies below the domain): not bracketed.
    # A floor above both leaves Gamma_ph >= 1 at every admissible density.
    for guess in (root*10,1e-8):  # guesses inside the bounded domain (a guess outside it is the caller's error)
        with pytest.raises(DomainError,match='not bracketed inside the microphysics domain'):
            solve(.7*edd,R,M,Bounded(root*1.07),opacity,guess=guess)
        with pytest.raises(DomainError,match='Gamma_ph >= 1.*at every admissible'):
            solve(.7*edd,R,M,Bounded(root*3),opacity,guess=guess)
    # Bounded from above too: a root above the domain is not bracketed either.
    class Capped(Bounded):
        def evaluate(self,rho,T):
            if np.any(np.real(rho)>self.floor):raise ValueError('outside the query rectangle')
            return EOS.evaluate(rho,T)
    with pytest.raises(DomainError,match='not bracketed inside the microphysics domain'):
        solve(.7*edd,R,M,Capped(root/3),opacity)
    # A guess outside the microphysics domain is the caller's error under either rule (second review:
    # the label stays numerical, never a domain exit; no registered caller passes one).
    for bracket in photosphere.BRACKETS:
        with pytest.raises(ValueError,match='outside the query rectangle'):
            photosphere.solve(.7*edd,R,M,Bounded(root*3),opacity,guess=root,bracket=bracket)

def test_registered_photosphere_roots_on_the_compiled_bridge():
    """Opt-in (QSTARS_MESA_LIBRARY of an ABI 7 build, MESA_DIR). The two cases of the final review, 1e4 Msun
    with the registered STARS opacity switch: (2.195e8 Lsun, 3030 Rsun) and (3.48e8 Lsun, 3500 Rsun) have
    ordinary roots at Gamma_ph 0.74 and 0.92, where the registered solve (bracket 'doubling', the default)
    raises 'Gamma_ph >= 1' from the default guess and bracket='scan' (amendment 4) finds them;
    (6.2e8 Lsun, 2000 Rsun), the registered 2000 Rsun seed, has Gamma_ph >= 1 at every admissible
    density. No model is built."""
    import importlib.util,json,os
    from pathlib import Path
    from qstars import quasistatic
    from qstars.constants import LSUN
    library,directory=os.environ.get('QSTARS_MESA_LIBRARY'),os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Set QSTARS_MESA_LIBRARY and MESA_DIR for compiled MESA tests')
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('validate_ball_track',root/'scripts/validate_ball_track.py')
    vbt=importlib.util.module_from_spec(spec);spec.loader.exec_module(vbt)
    eos,opacity=quasistatic.physics(vbt.build_config(json.loads(vbt.PLAN.read_text()),'C1',library,directory))
    if eos.abi_version<7:pytest.skip('ABI 7 STARS opacity required')
    M=1e4*MSUN
    for L,R,gamma in ((2.195e8,3030.,.74),(3.48e8,3500.,.92)):
        with pytest.raises(DomainError,match='\\(Gamma_ph >= 1\\); no Ball surface exists'):
            photosphere.solve(L*LSUN,R*RSUN,M,eos,opacity)  # the registered solve: false Gamma_ph >= 1
        rho=photosphere.solve(L*LSUN,R*RSUN,M,eos,opacity,bracket='scan')
        out=photosphere.evaluate(rho,L*LSUN,R*RSUN,M,eos,opacity)
        assert abs(out['residual'].real)<=1e-13 and out['eddington_factor'].real==pytest.approx(gamma,abs=.005)
        assert photosphere.solve(L*LSUN,R*RSUN,M,eos,opacity,guess=1e-8,bracket='scan')==pytest.approx(rho,rel=1e-10)
    with pytest.raises(DomainError,match='Gamma_ph >= 1.*at every admissible'):
        photosphere.solve(6.2e8*LSUN,2000*RSUN,M,eos,opacity,bracket='scan')
    with pytest.raises(DomainError,match='\\(Gamma_ph >= 1\\); no Ball surface exists'):
        photosphere.solve(6.2e8*LSUN,2000*RSUN,M,eos,opacity)

def test_outer_cell_converges_to_spherical_grey_atmosphere():
    """The discrete outer boundary has order >= 1 (order 2 asymptotically).

    Face N at R with p_ext from the photosphere; the last cell of mass dm
    obeys the endpoint momentum row 4 pi R^2 (p_N-1 - p_ext)/(dm/2) = G M/R^2
    and the grey relation T^4 = 3/4 Teff^4 (tau_half + 4/3), which
    transport.luminosities (stars_photosphere) reproduces. The exact
    spherical constant-kappa grey atmosphere (dP/dm = -G m/(4 pi r^4),
    dtau/dm = -kappa/(4 pi r^2)) is integrated to the cell's mass midpoint.
    """
    M,R,L=100*MSUN,150*RSUN,5e39
    rho_ph=photosphere.solve(L,R,M,EOS,OPACITY)
    ph=photosphere.evaluate(rho_ph,L,R,M,EOS,OPACITY)
    teff,pph=float(ph['teff']),float(ph['pressure'])
    config=resolve({'transport':{'convection':'off'}});config['boundary']['thermal']='stars_photosphere'
    def exact(dm):
        def rhs(m,y):
            r,P,tau=y;T4=.75*teff**4*(tau+2/3);rho=(P-A*T4/3)/(EOS.R*T4**.25)
            return [1/(4*PI*r*r*rho),-G*m/(4*PI*r**4),-KAPPA/(4*PI*r*r)]
        sol=solve_ivp(rhs,[M,M-dm/2],[R,pph,2/3],rtol=1e-12,atol=[1e-6*R,1e-12*pph,1e-14],method='DOP853')
        r,P,tau=sol.y[:,-1]
        return P,(.75*teff**4*(tau+2/3))**.25
    errors=[]
    for k in range(5):
        dm=4e-6*MSUN/2**k
        p=pph+G*M*dm/(8*PI*R**4)
        T=(.75*teff**4*(KAPPA*dm/(8*PI*R**2)+4/3))**.25
        # The grey relation is the one the solver uses for L_N.
        r=np.array([R*(1-1e-3),R]);rho=(p-A*T**4/3)/(EOS.R*T)
        state=State(0.,M-dm,np.array([dm]),r,np.zeros(2),np.array([T]),
                    cell_volume=np.array([dm/rho]))
        assert luminosities(state,EOS,OPACITY,config,0.)[0][-1]==pytest.approx(L,rel=1e-13)
        P,Tx=exact(dm)
        errors.append((abs(p/P-1),abs(T/Tx-1)))
    errors=np.array(errors);orders=np.log2(errors[:-1]/errors[1:])
    assert np.all(orders>=1) and np.all(orders[-1]>=1.8)
    assert errors[-1].max()<1e-4
