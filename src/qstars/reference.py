"""Published reference equations, independent of the time-dependent solver.

These are explicitly structure/law benchmarks, never full-track replication.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar,brentq
from .constants import G,PI

def _surface(x,y): return y[0]
_surface.terminal=True
_surface.direction=-1

def lane_emden(n=3,rtol=1e-10):
    x0=1e-6
    def rhs(x,y): return [-y[1]/x**2,x**2*max(y[0],0.)**n]
    sol=solve_ivp(rhs,[x0,1e3],[1-x0*x0/6,x0**3/3],events=_surface,
                  rtol=rtol,atol=rtol*.01,max_step=.1,dense_output=True)
    if len(sol.t_events[0])!=1: raise RuntimeError("Lane-Emden surface not found")
    return sol

def ball_loaded(x0,n=3,rtol=1e-10):
    """Ball+2011 eqs 24-26,30,32,33; BH-only Bondi radius + cavity."""
    bh=x0/(2*n)
    def rhs(x,y): return [-y[1]/x**2,x**2*max(y[0],0)**n]
    sol=solve_ivp(rhs,[x0,1e4],[1.,bh+2*x0**3/3],events=_surface,
                  rtol=rtol,atol=rtol*.01,max_step=.2)
    if len(sol.t_events[0])!=1: raise RuntimeError("Loaded polytrope has no surface")
    return dict(bh_fraction=bh/sol.y[1,-1],inner_outer_ratio=x0/sol.t[-1],xi0=x0)

def ball_limit(n=3,rtol=1e-10):
    sol=minimize_scalar(lambda x:-ball_loaded(x,n,rtol)["bh_fraction"],bounds=(.1,2.),method="bounded")
    return ball_loaded(sol.x,n,rtol)

def corrected_bondi(K,N=1.,rtol=1e-10):
    """Coughlin & Begelman 2024 eqs 2-7; mass includes cavity."""
    if K <= 6*N*(2/3)**(1/3): raise ValueError("Negative inferred black-hole mass")
    def rhs(x,y): return [-y[1]/(K*x*x),x*x*max(y[0],0)**3]
    sol=solve_ivp(rhs,[1,1e5],[6*N/K,1.],events=_surface,rtol=rtol,atol=rtol*.01,max_step=5)
    if len(sol.t_events[0])!=1: raise RuntimeError("Corrected Bondi surface not found")
    return (1-2/3*(6*N/K)**3)/sol.y[1,-1]

def corrected_limit(N=1.,rtol=1e-10):
    sol=minimize_scalar(lambda K:-corrected_bondi(K,N,rtol),
                         bounds=(6*N*(2/3)**(1/3)*(1+1e-6),50*N),method="bounded")
    return dict(N=N,K=sol.x,bh_fraction=-sol.fun)

def saturated_inner(rtol=1e-10):
    """CB24 eqs 13-16, integrated in log radius to avoid singular powers.

    z=(rho*x^0.5)^(1/3); dm/dlnx=x^2.5*z^3,
    dz/dlnx=1.5z-mz^3, from the printed ODE, without fitted coefficients.
    """
    x0=1e-7
    z0=np.sqrt(1.5)
    m0=1+3*np.sqrt(3)/(5*np.sqrt(2))*x0**2.5
    def rhs(t,y):
        m,z=y
        return [np.exp(2.5*t)*z**3,1.5*z-m*z**3]
    def entropy_min(t,y): return 3*y[0]*y[1]**2-16/3
    entropy_min.terminal=True; entropy_min.direction=1
    sol=solve_ivp(rhs,[np.log(x0),np.log(3)],[m0,z0],events=entropy_min,
                  rtol=rtol,atol=rtol*.01,max_step=.025,dense_output=True)
    if not len(sol.t_events[0]): raise RuntimeError("Entropy minimum not located")
    x=np.exp(sol.t[-1]); m,z=sol.y[:,-1]
    return dict(x_entropy_min=x,bh_inner_fraction=1/m,Ki_max=(x/m)**(5/3)),sol

def saturated_composite(xmatch,rtol=1e-9):
    """CB24 saturated core matched to adiabatic envelope, ending at eq 21."""
    info,inner=saturated_inner(rtol)
    if not 1e-5<xmatch<info["x_entropy_min"]: raise ValueError("Invalid saturated matching radius")
    m,z=inner.sol(np.log(xmatch))
    Ki=(xmatch/m)**(5/3)
    dmi=xmatch**2.5*z**3/m
    Ko=4*Ki/dmi
    def rhs(x,y): return [-y[1]/(Ko*x*x),x*x*max(y[0],0)**3]
    def transition(x,y):
        if x <= 1.000001: return 1e-10
        return (x*x*max(y[0],0)**3/dmi)**1.5/x-1
    transition.terminal=True; transition.direction=-1
    sol=solve_ivp(rhs,[1,1e5],[dmi**(1/3),1.],events=transition,
                  rtol=rtol,atol=rtol*.01,max_step=.2)
    if not len(sol.t_events[0]): raise RuntimeError("Radiative transition not found")
    return dict(Ki=Ki,bh_fraction=1/(m*sol.y[1,-1]),outer_inner_radius=sol.t[-1],
                adiabatic_mass_fraction=1-1/sol.y[1,-1])

def hassan_critical_mass(rho,cs,N=1):
    return cs**3/(12*np.sqrt(N**3*G**3*PI*rho))

def hassan_radius(mbh,rho,cs,N=1):
    """Small positive root continuous with zero BH mass, using stable brentq.

    This is equivalent to the paper's trigonometric Eq 6. Solving the cubic
    directly avoids cancellation in that expression as MBH/Mcrit tends to zero.
    """
    crit=hassan_critical_mass(rho,cs,N)
    if not 0<mbh<=crit: raise ValueError("No small positive corrected Bondi root")
    turning=cs/np.sqrt(16*PI*N*G*rho)
    ratio=mbh/crit
    x=brentq(lambda x:(x**3-3*x+2*ratio),0.,1.,xtol=1e-14)
    return turning*x

def run_reference(config,output='results'):
    """Configurable published STATIC families; no thermal-evolution claim."""
    import csv
    from . import io
    allowed={'family','n','N','minimum','maximum','samples'}
    if set(config)-allowed:raise ValueError('Unknown reference configuration key')
    family=config['family'];samples=config.get('samples',80)
    if not isinstance(samples,int) or not 4<=samples<=10000:raise ValueError('Reference samples must be 4-10000')
    if family=='ball2011':
        n=config.get('n',3)
        if n not in (2,3,4):raise ValueError('Validated Ball reference indices are 2,3,4')
        x=np.linspace(config.get('minimum',.05),config.get('maximum',1.5),samples)
        if min(x)<=0:raise ValueError('xi0 must be positive')
        rows=[ball_loaded(v,n) for v in x]
    elif family=='corrected_bondi':
        N=config.get('N',1.)
        if N<=0:raise ValueError('Bondi factor must be positive')
        x=np.linspace(config.get('minimum',5.5*N),config.get('maximum',25*N),samples)
        rows=[dict(K=v,bh_fraction=corrected_bondi(v,N)) for v in x]
    elif family=='saturated_convection':
        x=np.linspace(config.get('minimum',.15),config.get('maximum',1.),samples)
        rows=[dict(xmatch=v,**saturated_composite(v)) for v in x]
    else:raise ValueError('Unknown static reference family')
    path=io.create_directory(output,'reference_'+family);io.provenance(path,config)
    with (path/'sequence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report=dict(status='completed_reference',family=family,samples=samples,v1_validated=False,
                scope='Dimensionless published structure equations only; no time evolution or photosphere.')
    io.json_write(path/'termination.json',report)
    return path,report
