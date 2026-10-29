"""Ball thesis 2.15-2.16 surface: outer face at tau = 2/3 with no mass above it.

L* = 4 pi R*^2 sigma Teff^4 (2.15). The photospheric gas pressure is
P_g = (2 g/3 kappa_ph)(1 - L*/L_Edd), L_Edd = 4 pi G c M*/kappa_ph (2.16),
which equals 2 g/(3 kappa_ph) - a Teff^4/6: the gas pressure at tau = 2/3 of
a constant-kappa Eddington grey atmosphere, T^4 = 3/4 Teff^4 (tau + 2/3),
dP/dtau = g/kappa. kappa_ph = kappa(rho_ph, Teff) with rho_ph an implicit
unknown. The outer face feels P_ph = P_g + a Teff^4/3. Gamma_ph >= 1
(P_g <= 0) is a model domain exit and is never clipped.
"""
import numpy as np
from .constants import A,C,G,PI,SIGMA
from .mass_transport import DomainError

def effective_temperature(luminosity,radius):
    return (luminosity/(4*PI*radius**2*SIGMA))**.25

def eddington_factor(luminosity,mass,kappa):
    """Gamma_ph = L*/L_Edd with L_Edd = 4 pi G c M*/kappa_ph."""
    return kappa*luminosity/(4*PI*G*C*mass)

def gas_pressure(gravity,kappa,teff):
    """Thesis 2.16: (2 g/3 kappa)(1 - Gamma) = 2 g/(3 kappa) - a Teff^4/6."""
    return 2*gravity/(3*kappa)-A*teff**4/6

def external_pressure(pg,teff):
    """Total pressure on the outer face: P_g + a Teff^4/3 (T = Teff at tau = 2/3)."""
    return pg+A*teff**4/3

def implied_column(radius,pressure,gravity):
    """Mass a hydrostatic atmosphere above the face would need, 4 pi R^2 P_ph/g (diagnostic N4)."""
    return 4*PI*radius**2*pressure/gravity

def _row(rho_ph,luminosity,radius,mass,eos,opacity):
    """Row (f) and the photospheric state without the domain check (see evaluate)."""
    teff=effective_temperature(luminosity,radius)
    rho=np.atleast_1d(rho_ph);T=np.atleast_1d(teff)
    q=eos.evaluate(rho,T);kappa=np.atleast_1d(opacity(rho,T))[0]
    gravity=G*mass/radius**2
    scale=2*gravity/(3*kappa)
    target=scale-A*teff**4/6
    pg=q['pg'][0]
    return dict(residual=(pg-target)/scale,teff=teff,kappa=kappa,gravity=gravity,gas_pressure=pg,
                target=target,pressure=external_pressure(pg,teff),
                eddington_factor=eddington_factor(luminosity,mass,kappa))

def evaluate(rho_ph,luminosity,radius,mass,eos,opacity):
    """Row (f) and the photospheric state at a trial (rho_ph, L_N, R, M_N); complex-step safe.

    residual = [P_g,EOS(rho_ph, Teff) - (2 g/3 kappa_ph - a Teff^4/6)]/(2 g/3 kappa_ph).
    Raises DomainError when the target gas pressure is not positive.
    """
    out=_row(rho_ph,luminosity,radius,mass,eos,opacity)
    if not np.real(out['target'])>0:
        raise DomainError('Photospheric gas pressure is not positive (Gamma_ph >= 1); no Ball surface exists')
    return out

SCAN_POINTS=400  # fallback scan of the admissible ln rho_ph interval (solve, bracket 'scan')
# Bracket rules of solve. 'doubling' is the registered solve and the default of every caller
# (quasistatic.DEFAULTS photosphere_bracket); 'scan' is an option (Ball plan amendment 4).
BRACKETS=('doubling','scan')

def solve(luminosity,radius,mass,eos,opacity,guess=None,tolerance=1e-14,max_iterations=200,bracket='doubling'):
    """rho_ph of row (f) for given L*, R*, M*: safeguarded Newton in ln rho on a bracket.

    The residual increases with rho wherever kappa does not fall faster than
    1/P_g; a sign change is bracketed first, so a non-monotone opacity can
    only select a root, never a non-root.

    bracket 'doubling' (default; the registered solve, unchanged): the
    bracket is sought by doubling steps in ln rho from the guess, and every
    trial density is evaluated with the domain check, so a trial with
    Gamma_ph >= 1 raises DomainError('... (Gamma_ph >= 1); no Ball surface
    exists') even when an ordinary root exists: kappa_ph, and with it
    Gamma_ph, depends on the trial density (a false Gamma_ph >= 1, e.g.
    (2.195e8 Lsun, 3030 Rsun) and (3.48e8 Lsun, 3500 Rsun) at 1e4 Msun on
    the STARS opacity switch).

    bracket 'scan' (option, off unless set): the row is continuous through
    Gamma_ph = 1 (residual = P_g/scale - (1 - Gamma_ph), positive wherever
    Gamma_ph >= 1), so a trial density with Gamma_ph >= 1 is an upper bound
    of the bracket, not a failure, and a root has P_g = target > 0, i.e.
    Gamma_ph < 1. When the doubling scan from the guess finds no sign
    change, the whole admissible interval (the microphysics domain reached
    from the guess) is scanned at SCAN_POINTS densities before a DomainError
    is raised: 'Gamma_ph >= 1' only when Gamma_ph >= 1 at every admissible
    density sampled (no Ball surface exists), otherwise 'not bracketed inside
    the microphysics domain'. Where the doubling scan finds the bracket the
    iterates and the root are those of 'doubling'. A guess outside the
    microphysics domain raises its ValueError under either rule.
    """
    if bracket not in BRACKETS:raise ValueError(f'photosphere bracket must be one of {BRACKETS}')
    if bracket=='doubling':return _solve_doubling(luminosity,radius,mass,eos,opacity,guess,tolerance,max_iterations)
    def f(lnrho):
        x=complex(lnrho,1e-30)
        out=_row(np.exp(x),luminosity,radius,mass,eos,opacity)
        return float(out['residual'].real),float(out['residual'].imag/1e-30)
    teff=float(effective_temperature(luminosity,radius))
    x=float(np.log(1e-12 if guess is None else guess))
    fx,_=f(x)
    lo=hi=x;flo=fhi=fx
    edges={}
    for sign in (-1,1):
        step=.5
        while (flo>0 if sign<0 else fhi<0) and step<64:
            trial=(lo if sign<0 else hi)+sign*step;step*=2
            try:value,_=f(trial)
            except ValueError:  # EOS or opacity domain edge: no root beyond it
                edges[sign]=trial;break
            if sign<0:lo,flo=trial,value
            else:hi,fhi=trial,value
    if not (flo<=0<=fhi):
        lo,hi,flo,fhi=_scan(f,lo,hi,edges,luminosity,radius,mass,eos,opacity)
    x=lo if abs(flo)<abs(fhi) else hi
    for _ in range(max_iterations):
        fx,dfx=f(x)
        if abs(fx)<=tolerance:break
        if fx<0:lo=x
        else:hi=x
        trial=x-fx/dfx if dfx>0 else .5*(lo+hi)
        x=trial if lo<trial<hi else .5*(lo+hi)
        if hi-lo<=1e-15*max(1.,abs(x)):break
    else:raise RuntimeError(f'Photospheric density solve failed at Teff={teff:g} K')
    rho=float(np.exp(x))
    evaluate(rho,luminosity,radius,mass,eos,opacity)  # a root within rounding of Gamma_ph = 1 is a domain exit
    return rho

def _solve_doubling(luminosity,radius,mass,eos,opacity,guess,tolerance,max_iterations):
    """The registered solve (bracket 'doubling'), byte for byte the pre-amendment routine."""
    def f(lnrho):
        x=complex(lnrho,1e-30)
        out=evaluate(np.exp(x),luminosity,radius,mass,eos,opacity)
        return float(out['residual'].real),float(out['residual'].imag/1e-30)
    teff=float(effective_temperature(luminosity,radius))
    x=float(np.log(1e-12 if guess is None else guess))
    fx,_=f(x)
    lo=hi=x;flo=fhi=fx
    for sign in (-1,1):
        step=.5
        while (flo>0 if sign<0 else fhi<0) and step<64:
            trial=(lo if sign<0 else hi)+sign*step;step*=2
            try:value,_=f(trial)
            except ValueError:break  # EOS or opacity domain edge: no root beyond it
            if sign<0:lo,flo=trial,value
            else:hi,fhi=trial,value
    if not (flo<=0<=fhi):raise DomainError('Photospheric density not bracketed inside the microphysics domain')
    x=lo if abs(flo)<abs(fhi) else hi
    for _ in range(max_iterations):
        fx,dfx=f(x)
        if abs(fx)<=tolerance:return float(np.exp(x))
        if fx<0:lo=x
        else:hi=x
        trial=x-fx/dfx if dfx>0 else .5*(lo+hi)
        x=trial if lo<trial<hi else .5*(lo+hi)
        if hi-lo<=1e-15*max(1.,abs(x)):return float(np.exp(x))
    raise RuntimeError(f'Photospheric density solve failed at Teff={teff:g} K')

def _scan(f,lo,hi,edges,luminosity,radius,mass,eos,opacity):
    """Bracket of the first sign change on the admissible interval, or the DomainError that no root exists.

    [lo, hi] are the lowest and highest densities (ln) the doubling scan
    evaluated; where it stopped at a microphysics domain edge (edges), the
    edge is located by bisection first.
    """
    def valid(x):
        try:f(x);return True
        except ValueError:return False
    for sign,bad in edges.items():
        good=lo if sign<0 else hi
        for _ in range(60):
            mid=.5*(good+bad)
            if valid(mid):good=mid
            else:bad=mid
        if sign<0:lo=good
        else:hi=good
    grid=np.linspace(lo,hi,SCAN_POINTS);values=[];gammas=[]
    for x in grid:
        try:
            out=_row(np.exp(x),luminosity,radius,mass,eos,opacity)
            values.append(float(np.real(out['residual'])));gammas.append(float(np.real(out['eddington_factor'])))
        except ValueError:values.append(np.nan);gammas.append(np.nan)
    for k in range(len(grid)-1):
        if values[k]<=0<=values[k+1]:return grid[k],grid[k+1],values[k],values[k+1]
    known=[g for g in gammas if np.isfinite(g)]
    if known and min(known)>=1:
        raise DomainError('Photospheric gas pressure is not positive (Gamma_ph >= 1) at every admissible photospheric density; '
                          'no Ball surface exists')
    raise DomainError('Photospheric density not bracketed inside the microphysics domain')
