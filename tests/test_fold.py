"""Gate 6 Stage 5: frozen-entropy fold continuation (fold.py), spec test I4.

I4: an isentropic loaded polytrope with the Ball cavity (r^-3/2, b = 1,
face boundary state) continued in M_BH at frozen entropy and fixed total
mass folds where reference.ball_limit(n) puts the maximum of M_BH/M. The
spec names Gamma = 4/3; for n = 3 exactly the frozen problem is singular
along the whole branch (the homology r -> a r, rho -> rho/a^3, T -> T/a is
an exact symmetry of every row), which test_gamma_four_thirds_is_singular_by_homology
shows, so the identity is tested at n = 2.5 and 2.9 (Gamma = 1.4, 1.345).
Also: the frozen residual vanishes on an accepted quasi-static state, its
declared pattern, the det-sign/pivot/singular-value monitors against dense
linear algebra, and the classification labels.
"""
import numpy as np
import pytest
from qstars import fold,ledger,quasistatic,reference
from qstars.config import resolve
from qstars.constants import G,PI,MSUN,RSUN,KB,MU,YEAR
from qstars.eos import Composition,IdealGasRadiation
from qstars.nonlinear import pattern_jacobian
from qstars.opacity import ElectronScattering
from qstars.state import State


class Polytropic:
    """Ideal gas with constant gamma: p = rho R T, u = R T/(gamma - 1); its isentropes are p = K rho^gamma."""
    def __init__(self,gamma,mu=.6):
        self.gamma=gamma;self.R=KB/(mu*MU);self.composition=Composition(.7,.3,0.)
    def evaluate(self,rho,T):
        rho,T=np.asarray(rho),np.asarray(T)
        if np.any(np.real(rho)<=0) or np.any(np.real(T)<=0):raise ValueError('EOS needs positive density and temperature')
        g,R=self.gamma,self.R
        p=R*rho*T;cv=R/(g-1)+0*p
        return dict(p=p,pg=p,pr=0*p,u=R*T/(g-1),s=cv*np.log(T)-R*np.log(rho),cv=cv,cp=g*cv,gamma1=g+0*p,
                    grad_ad=(g-1)/g+0*p,delta=1+0*p,cs=np.sqrt(g*p/rho),h=g*R*T/(g-1))
    def temperature(self,rho,value,quantity='u'):
        rho,value=np.broadcast_arrays(rho,value)
        return value*(self.gamma-1)/self.R if quantity=='u' else value/(rho*self.R)


def config():
    c=resolve(dict(numerics=dict(zones=8)))
    # A vanishing mass defect: the total mass is then fixed along the branch, as in reference.ball_limit.
    c['engine'].update(model='ball2012_eq3_15',efficiency=1e-9)
    c['boundary'].update(thermal='grey',pressure_cgs=0.)
    c['cavity']=dict(boundary_state='face')
    return c

def polytrope(n,fraction,zones,total=1e4*MSUN,radius=3000*RSUN):
    """Loaded polytrope (paper eqs 24-33, index n, Gamma = 1 + 1/n), xi - xi0 ~ 1 - cos(pi u/2): the continuum,
    not yet discretely hydrostatic; fold.solve makes it so at frozen entropy."""
    eos=Polytropic(1+1/n)
    sol=quasistatic.loaded_polytrope(fraction,n=n)
    xi0,xis,phi0,phis=sol['xi0'],sol['xi_s'],sol['phi0'],sol['phi_s']
    alpha=radius/xis;rho_c=total/(4*PI*alpha**3*phis);K=4*PI*G*alpha**2*rho_c**(1-1/n)/(n+1)
    u=np.linspace(0,1,zones+1);xi=xi0+(xis-xi0)*(1-np.cos(PI*u/2))
    phi=sol['solution'].sol(xi)[1];phi[0],phi[-1]=phi0,phis
    dm=total*np.diff(phi)/phis;r=alpha*xi
    V=4*PI/3*np.diff(r)*(r[1:]**2+r[1:]*r[:-1]+r[:-1]**2)
    rho=dm/V;T=eos.temperature(rho,K*rho**(1+1/n),'p')
    state=State(0.,fraction*total,dm,r,np.zeros(zones+1),T,cell_volume=V,cavity_mass=2/3*xi0**3/phis*total,
                kinetic_tracked=False)
    return state,eos

def discrete(n,fraction,zones):
    """The frozen-entropy (uniform K) discrete solution at the continuum's M_BH, as a State."""
    state,eos=polytrope(n,fraction,zones);c=config()
    P=fold.Frozen(state,eos,None,c)
    x,info=fold.solve(P,P.x0(),0.)
    assert info['success']
    return P.state(x,0.),eos,c

# ---------------------------------------------------------------- I4

@pytest.mark.parametrize('n',[2.5,2.9])
def test_frozen_entropy_continuation_folds_at_the_loaded_polytrope_limit(n):
    limit=reference.ball_limit(n)['bh_fraction']
    state,eos=polytrope(n,.9*limit,100);c=config()
    P=fold.Frozen(state,eos,None,c)
    x,info=fold.solve(P,P.x0(),0.)
    assert info['success']
    run=fold.Arc(P).run(x,0.,1,.5)
    f=run['fold'];pts=run['points'];b=f['bracket']
    assert run['stop']=='fold passed' and f['refined']
    # Spec I4: 0.5%; the 100-zone discrete fold sits within 1e-4 of the continuum maximum.
    assert abs(f['bh_fraction']/limit-1)<=5e-3 and abs(f['bh_fraction']/limit-1)<=1e-4
    # det J changes sign at the turning point of lambda, where t_lambda changes sign.
    assert pts[b[0]]['det_sign']!=pts[b[1]]['det_sign'] and pts[b[0]]['t_lambda']>0>pts[b[1]]['t_lambda']
    assert all(p['det_sign']==pts[0]['det_sign'] for p in pts[:b[0]+1])
    assert max(p['lam'] for p in pts)<=f['lam']+1e-9
    # The smallest singular value vanishes as the square root of the distance to the fold.
    exponent,used=fold.sqrt_law(f)
    assert used==4 and abs(exponent-.5)<=.03
    # Fixed total mass up to the (vanishing) defect; the swept gas left through the base.
    assert abs(f['total_mass']/pts[0]['total_mass']-1)<=1e-8 and f['swept']>0

def test_fold_location_converges_with_resolution():
    """Second order: the 40- and 100-zone folds differ from the continuum limit by a ratio near (100/40)^2."""
    limit=reference.ball_limit(2.5)['bh_fraction'];errors=[]
    for zones in (40,100):
        state,eos=polytrope(2.5,.9*limit,zones);P=fold.Frozen(state,eos,None,config())
        x,_=fold.solve(P,P.x0(),0.)
        errors.append(fold.Arc(P).run(x,0.,1,.5)['fold']['bh_fraction']/limit-1)
    assert errors[0]>errors[1]>0 and 4<=errors[0]/errors[1]<=9

def test_gamma_four_thirds_is_singular_by_homology():
    """n = 3: d ln r = 1, d ln T = -1 (rho ~ a^-3, p ~ a^-4) is an exact null vector of J at a solution,
    so the frozen problem has no isolated solutions and no fold to locate; at n = 2.5 it is not."""
    out={}
    for n in (3.,2.5):
        limit=reference.ball_limit(n)['bh_fraction']
        state,eos,c=discrete(n,.9*limit,40)
        P=fold.Frozen(state,eos,None,c)
        x,_=fold.solve(P,P.x0(),0.)
        J,_=P.jacobian(x,0.)
        h=np.zeros(2*P.n+3);h[:P.n+1]=1.;h[P.n+1:2*P.n+1]=-1.
        out[n]=(J@h,fold.monitors(J)['min_singular'],P.n)
    # n = 3: J h vanishes to roundoff and so does the smallest singular value.
    assert np.max(abs(out[3.][0]))<=1e-8 and out[3.][1]<=1e-8
    # n = 2.5: only the entropy rows break the symmetry, by (3 - n)/(n + 1) = 1/7 (s = c_v ln T - R ln rho).
    Jh,sigma,cells=out[2.5]
    assert np.allclose(Jh[cells:2*cells],1/7,rtol=1e-6) and np.max(abs(np.delete(Jh,range(cells,2*cells))))<=1e-8
    assert sigma>=1e-4

# ---------------------------------------------------------------- residual, pattern, monitors

def test_declared_pattern_matches_dense_jacobian():
    limit=reference.ball_limit(2.5)['bh_fraction']
    state,eos,c=discrete(2.5,.9*limit,12)
    P=fold.Frozen(state,eos,None,c)
    x=P.x0()+np.random.default_rng(4).normal(0,1e-4,2*P.n+3)
    J=pattern_jacobian(lambda z:P.evaluate(z,.01),x,P.pattern,strict=True).toarray()
    dense=np.empty_like(J)
    for k in range(len(x)):
        xx=x.astype(complex);xx[k]+=1e-25j;dense[:,k]=P.evaluate(xx,.01).imag/1e-25
    assert np.max(abs(J-dense))<=1e-13*np.max(abs(dense))

def test_monitors_match_dense_linear_algebra():
    rng=np.random.default_rng(11)
    for trial in range(6):
        A=rng.normal(size=(9,9))+np.diag(rng.normal(size=9))
        A[:,0]*=1e-3*(trial+1)  # an increasingly small singular value
        m=fold.monitors(A)
        sign,_=np.linalg.slogdet(A)
        assert m['det_sign']==int(sign)
        assert m['min_singular']==pytest.approx(np.linalg.svd(A,compute_uv=False)[-1],rel=1e-8)
        assert m['log_abs_det']==pytest.approx(np.linalg.slogdet(A)[1],rel=1e-10)
    assert fold.parity([1,0,2])==-1 and fold.parity([1,2,0])==1 and fold.parity(list(range(5)))==1

def test_frozen_residual_vanishes_on_an_accepted_quasistatic_state():
    """The dt -> 0 limit of the step: an accepted Ball-toy state (endpoint balance, piston, cavity mass,
    photosphere with the grey L_N) solves the frozen problem at lambda = 0 to its Newton tolerance."""
    eos=IdealGasRadiation(Composition(.7,.3,0.));opacity=ElectronScattering(eos.composition)
    c=resolve(dict(composition=dict(X=.7,Y=.3,Z=0.),numerics=dict(zones=48)))
    c['engine'].update(model='ball2012_eq3_15',coupling=1.5e-4)
    c['transport']['convection']='mlt_ml1';c['boundary']['thermal']='stars_photosphere'
    c['validity']['min_temperature_K']=3000.;c['cavity']=dict(boundary_state='cell')
    s,_=quasistatic.seed(c,eos,opacity,48,MSUN,100*MSUN,200*RSUN);ledger.initialize(s,eos,cavity=True)
    s,record=quasistatic.relax(s,eos,opacity,c,dict(mesh_cycles=0))
    assert record['status']=='converged'
    s=quasistatic.zero_time(s,eos)
    new,_=quasistatic.bondi_step(s,10*YEAR,eos,opacity,c,mode='evolve')
    P=fold.Frozen(new,eos,opacity,c)
    f=P.evaluate(P.x0(),0.)
    assert np.max(abs(f))<=1e-9
    # Frozen continuation is defined from it; with M_BH/M = 0.01 no fold lies within 5% (search range).
    out=fold.classify(new,eos,opacity,c,dict(search=.01,window=.005))
    assert out['label']=='numerical_failure' and out['fold'] is None and 'no fold' in out['reason']
    assert out['runs'][0]['stop']=='lambda limit' and out['runs'][0]['points'][-1]['lam']>=.01

# ---------------------------------------------------------------- classification

def test_classification_confirms_a_fold_near_the_terminal_mass():
    """A model 0.5% below its frozen-entropy fold is a model_domain_exit; 10% below, the fold lies
    outside the 2% window and the end stays a numerical_failure (never relabelled)."""
    limit=reference.ball_limit(2.5)['bh_fraction']
    near,eos,c=discrete(2.5,.995*limit,60)
    out=fold.classify(near,eos,None,c)
    f=out['fold']
    assert out['label']=='model_domain_exit' and out['fold_confirmed'] and f['within_window'] and f['sqrt_law']
    assert 0<f['relative_offset']<=.02 and f['det_sign_change'] and abs(f['bh_fraction']/limit-1)<=1e-3
    far,eos,c=discrete(2.5,.9*limit,60)
    out=fold.classify(far,eos,None,c,dict(search=.5))
    assert out['fold'] is not None and not out['fold']['within_window'] and out['fold']['sqrt_law']
    assert out['label']=='numerical_failure' and 'within_window' in out['reason'] and not out['fold_confirmed']
