"""Mixing-length closures for the Ball et al. 2011 (STARS) replication.

ml1: Boehm-Vitense 1958 with radiative losses of the convective element, in
the Cox & Giuli 14.79-14.86 form with the MESA 24.08.1 ML1 constants
(turb/private/mlt.f90: ff1=1/8, ff2=1/2, ff3=24, ff4=0, so a0=9/4). Mixing
length alpha*Hp with the hydrostatic Hp=p/(rho g), Schwarzschild criterion
(fixed composition), no velocity or flux cap (STARS has none). The cubic
a0 G^3 + G^2 + G = A^2 (grad_rad - grad_ad) (C&G 14.82) gives exactly
Y = nabla - nabla_ad = G(G+1)/A^2, so the excess carries no cancellation of
(1 - zeta). Radiative faces: Y = grad_rad - grad_ad, C1-continuous at the
Schwarzschild boundary.

thesis_lossless: Ball thesis eqs 2.6-2.9 with (l/Hp)^2, which printed
eq 2.10 omits, delta = 1, no radiative losses, no cap. Sensitivity run S2.

Face inputs are given, not derived: grad_rad = 3 kappa p L/(16 pi a c G m T^4).
Every function accepts complex-step perturbed inputs; the regime switch and
the closed-form seed use real parts only, and Newton polishing carries the
imaginary parts.
"""
import numpy as np
from .constants import A as ARAD,C,G,PI

MODELS=('mlt_ml1','mlt_thesis_lossless')
FF1,FF2,FF3,FF4=.125,.5,24.,0.
A0=3/16*FF2*FF3  # MESA a0 with ff4=0: exactly 9/4

def conductance(kappa,p,T,mass):
    """Luminosity per unit gradient: L = conductance*nabla_rad."""
    return 16*PI*ARAD*C*G*mass*T**4/(3*kappa*p)

def radiative_gradient(luminosity,kappa,p,T,mass):
    return luminosity/conductance(kappa,p,T,mass)

def ml1_coefficient(cp,Q,p,rho,T,kappa,g,alpha=2.):
    """C&G 14.81 efficiency scale A in MESA's operation order; Q = chiT/chiRho."""
    omega=alpha*p/(rho*g)*rho*kappa
    return 4*cp*np.sqrt(FF1*p*Q*rho)*(alpha*omega)/(FF3*ARAD*C*T**3)

def thesis_coefficient(cp,p,rho,T,kappa,g,alpha=2.):
    """U = L_con/(L/nabla_rad)/Y^(3/2) for thesis eqs 2.6-2.9 with delta=1."""
    hp=p/(rho*g)
    return 3*alpha**2*cp*kappa*rho**2*hp*np.sqrt(g*hp/2)/(4*ARAD*C*T**3)

def _ml1_seed(delta):
    """MESA's closed-form root of a0 G^3+G^2+G = delta (real delta only)."""
    f=-2+9*A0+27*A0*A0*delta
    big=f>1e100
    f0=np.where(big,f,np.sqrt(np.where(big,0.,f)**2+4*(3*A0-1)**3))
    f1=np.cbrt(f+f0)
    return (np.cbrt(4.)*f1+2*np.cbrt(2.)*(1-3*A0)/f1-2)/(6*A0)

def _ml1_efficiency(delta):
    # Two Newton polishes restore relative accuracy where the closed form
    # cancels (small delta) and carry complex steps exactly.
    gamma=_ml1_seed(np.real(delta)).astype(np.result_type(delta,float))
    for _ in range(2):
        gamma=gamma-(((A0*gamma+1)*gamma+1)*gamma-delta)/((3*A0*gamma+2)*gamma+1)
    return gamma

def _thesis_root(u,d):
    """Positive root of u y^3 + y^2 = d (d>0): monotone Newton from an upper bound."""
    ur,dr=np.real(u),np.real(d)
    y=np.minimum(np.sqrt(dr),np.cbrt(dr/ur))
    for _ in range(100):
        step=((ur*y+1)*y*y-dr)/((3*ur*y+2)*y)
        y=y-step
        if np.all(abs(step)<=1e-14*y):break
    else:raise ValueError('Thesis MLT cubic failed to converge')
    y=y.astype(np.result_type(u,d,float))
    for _ in range(2):y=y-((u*y+1)*y*y-d)/((3*u*y+2)*y)
    return y

def ml1(grad_rad,grad_ad,cp,Q,p,rho,T,kappa,g,alpha=2.):
    """ML1 gradient for given face nabla_rad; matches MESA mlt('ML1',alpha)."""
    grad_rad,grad_ad,cp,Q,p,rho,T,kappa,g=np.broadcast_arrays(grad_rad,grad_ad,cp,Q,p,rho,T,kappa,g)
    A=np.asarray(ml1_coefficient(cp,Q,p,rho,T,kappa,g,alpha))
    D=np.asarray(grad_rad-grad_ad)
    dtype=np.result_type(D,A,float);unstable=np.asarray(D.real>0)
    Y=np.array(D,dtype=dtype);gamma=np.zeros(D.shape,dtype);velocity=np.zeros(D.shape,dtype)
    if np.any(unstable):
        a=A[unstable];a2=a*a;G1=_ml1_efficiency(a2*D[unstable])
        gamma[unstable]=G1;Y[unstable]=G1*(G1+1)/a2
        # C&G 14.86b
        velocity[unstable]=alpha*np.sqrt(Q[unstable]*p[unstable]/(8*rho[unstable]))*G1/a
    # Radiative faces return nabla_rad itself, as MESA does.
    return dict(gradient=np.where(unstable,grad_ad+Y,grad_rad),excess=Y,efficiency=gamma,
                velocity=velocity,A=A,convective=unstable)

def thesis_lossless(grad_rad,grad_ad,cp,p,rho,T,kappa,g,alpha=2.):
    """Lossless thesis MLT: U y^3 + y^2 = grad_rad - grad_ad, Y = y^2, v = alpha y sqrt(g Hp/2)."""
    grad_rad,grad_ad,cp,p,rho,T,kappa,g=np.broadcast_arrays(grad_rad,grad_ad,cp,p,rho,T,kappa,g)
    U=np.asarray(thesis_coefficient(cp,p,rho,T,kappa,g,alpha))
    D=np.asarray(grad_rad-grad_ad)
    dtype=np.result_type(D,U,float);unstable=np.asarray(D.real>0)
    Y=np.array(D,dtype=dtype);velocity=np.zeros(D.shape,dtype)
    if np.any(unstable):
        y=_thesis_root(U[unstable],D[unstable]);Y[unstable]=y*y
        hp=p[unstable]/(rho[unstable]*g[unstable])
        velocity[unstable]=alpha*y*np.sqrt(g[unstable]*hp/2)
    return dict(gradient=np.where(unstable,grad_ad+Y,grad_rad),excess=Y,velocity=velocity,
                U=U,convective=unstable)

def closure(model,grad_rad,grad_ad,cp,Q,p,rho,T,kappa,g,alpha=2.):
    if model=='mlt_ml1':return ml1(grad_rad,grad_ad,cp,Q,p,rho,T,kappa,g,alpha)
    if model=='mlt_thesis_lossless':return thesis_lossless(grad_rad,grad_ad,cp,p,rho,T,kappa,g,alpha)
    raise ValueError(f'Unknown mixing-length closure {model}')

def forward_luminosity(excess,grad_ad,cp,Q,p,rho,T,kappa,g,mass,alpha=2.,model='mlt_ml1'):
    """Face luminosity carried by a given excess Y = nabla - nabla_ad (predictors).

    Explicit inverse of the closure, no root find. ML1: G = 2A^2Y/(1+sqrt(1+4A^2Y)),
    nabla_rad - nabla = a0 G^3/A^2. Thesis: nabla_rad - nabla = U Y^(3/2).
    Radiative part K nabla, convective part K(nabla_rad - nabla), K = conductance.
    """
    excess,grad_ad,cp,Q,p,rho,T,kappa,g,mass=np.broadcast_arrays(excess,grad_ad,cp,Q,p,rho,T,kappa,g,mass)
    K=np.asarray(conductance(kappa,p,T,mass))
    dtype=np.result_type(excess,grad_ad,cp,Q,p,rho,T,kappa,g,mass,float);unstable=np.asarray(excess.real>0)
    extra=np.zeros(excess.shape,dtype);velocity=np.zeros(excess.shape,dtype);Y=excess[unstable]
    if model=='mlt_ml1':
        a=np.asarray(ml1_coefficient(cp,Q,p,rho,T,kappa,g,alpha))[unstable];a2=a*a
        G1=2*a2*Y/(1+np.sqrt(1+4*a2*Y))
        extra[unstable]=A0*G1**3/a2
        velocity[unstable]=alpha*np.sqrt(Q[unstable]*p[unstable]/(8*rho[unstable]))*G1/a
    elif model=='mlt_thesis_lossless':
        u=np.asarray(thesis_coefficient(cp,p,rho,T,kappa,g,alpha))[unstable];y=np.sqrt(Y)
        extra[unstable]=u*Y*y
        hp=p[unstable]/(rho[unstable]*g[unstable])
        velocity[unstable]=alpha*y*np.sqrt(g[unstable]*hp/2)
    else:raise ValueError(f'Unknown mixing-length closure {model}')
    radiative=K*(grad_ad+excess)
    return dict(luminosity=radiative+K*extra,radiative=radiative,convective=K*extra,
                grad_rad=grad_ad+excess+extra,velocity=velocity,convective_face=unstable)
