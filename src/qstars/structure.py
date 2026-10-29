"""Hydrostatic initialization on the same staggered mesh as the evolution."""
import numpy as np
from scipy.optimize import least_squares
from .constants import PI,G,MSUN
from .state import State

def sparse_polytrope(dm,bh_mass,r0,R,pext,continuum):
    """Same discrete hydrostatic equations, with local face-radius unknowns.

    The previous cumulative-volume parameterization has a dense Jacobian.
    Radius variables instead give a tridiagonal block plus one K column.
    A continuum solution is only a starting guess; the discrete equations
    determine the returned volumes and pressure independently.
    """
    from scipy.interpolate import PchipInterpolator
    from scipy.sparse import csc_matrix
    from .nonlinear import newton
    n=len(dm);mass=np.r_[0,np.cumsum(dm)/dm.sum()];mass[-1]=1.
    z=np.linspace(np.log(r0/R),0,30001)
    qm=continuum.sol(z)[0];qm[0]=0.;qm[-1]=1.
    if np.any(np.diff(qm)<=0):
        # Within the BVP tolerance m/M can overshoot 1 (or stall) next to a
        # boundary. The seed only needs a monotone table: keep the interior
        # nodes strictly inside (0, 1) that rise above every earlier node.
        interior=np.flatnonzero((qm[1:-1]>0)&(qm[1:-1]<1))+1
        rising=qm[interior]>np.maximum.accumulate(np.r_[0.,qm[interior][:-1]])
        keep=np.r_[0,interior[rising],len(qm)-1]
        z,qm=z[keep],qm[keep]
    seed=PchipInterpolator(qm,np.exp(z)*R)(mass);seed[0]=r0;seed[-1]=R
    enclosed=bh_mass+np.cumsum(dm)
    weights=.5*np.r_[dm[:-1]+dm[1:],dm[-1]]
    def unpack(x):
        r=np.r_[r0,seed[1:]*np.exp(x[:n])]
        V=4*PI/3*np.diff(r)*(r[1:]**2+r[1:]*r[:-1]+r[:-1]**2)
        if np.any(V<=0):raise ValueError('Hydrostatic trial has crossing faces')
        rho=dm/V;p=continuum.K*np.exp(x[-1])*rho**(4/3)
        support=G*enclosed*weights/(4*PI*r[1:]**4)
        return r,V,rho,p,support
    def residual(x):
        r,V,rho,p,support=unpack(x)
        return np.r_[np.log(p/(np.r_[p[1:],pext]+support)),np.log(r[-1]/R)]
    def jacobian(fun,x):
        r,V,rho,p,support=unpack(x)
        den=np.r_[p[1:],pext]+support
        next_fraction=np.r_[p[1:]/den[:-1],0.]
        left=(4/3)*4*PI*r[:-1]**3/V
        right=-(4/3)*4*PI*r[1:]**3/V
        diag=right+4*support/den
        diag[:-1]-=next_fraction[:-1]*left[1:]
        rows=np.r_[np.arange(n),np.arange(1,n),np.arange(n-1),np.arange(n),n]
        cols=np.r_[np.arange(n),np.arange(n-1),np.arange(1,n),np.full(n,n),n-1]
        values=np.r_[diag,left[1:],-next_fraction[:-1]*right[1:],1-next_fraction,1.]
        return csc_matrix((values,(rows,cols)),shape=(n+1,n+1))
    sol=newton(residual,np.zeros(n+1),2e-10,jacobian=jacobian)
    err=float(max(abs(residual(sol.x))))
    if err>1e-8:raise RuntimeError(f'Sparse hydrostatic initialization failed: residual {err:g}')
    r,V,rho,p,_=unpack(sol.x)
    return r,V,rho,p,err

def initialize(config,eos):
    c,n = config["initial"],config["numerics"]["zones"]
    # Cosine mass mesh resolves both boundaries more finely than the bulk.
    q = .5*(1-np.cos(np.linspace(0,np.pi,n+1)))
    if config['numerics']['mesh']=='core_refined': q=q*q
    dm = c["envelope_mass_msun"]*MSUN*np.diff(q)
    continuum=None
    if config['numerics']['mesh'] in ('log_radius','log_radius_cosine','log_radius_surface'):
        continuum=continuum_polytrope(config,eos)
        z=np.linspace(0,1,n+1)
        if config['numerics']['mesh']=='log_radius_cosine':z=.5*(1-np.cos(np.pi*z))
        if config['numerics']['mesh']=='log_radius_surface':z=np.sin(.5*np.pi*z)
        q=continuum.sol(np.log(c['inner_radius_cm']/c['radius_cm'])*(1-z))[0]
        q[0],q[-1]=0.,1.
        dm=c['envelope_mass_msun']*MSUN*np.diff(q)
    r0,R = c["inner_radius_cm"],c["radius_cm"]
    r = (r0**3+(R**3-r0**3)*q)**(1/3)
    r[0] = r0
    state = State(0.,c["bh_mass_msun"]*MSUN,dm,r,np.zeros(n+1),np.ones(n))
    vscale = 4*PI/3*(R**3-r0**3)
    pext = config["boundary"]["pressure_cgs"]
    pscale = G*state.enclosed[-1]**2/(4*PI*R**4)
    rhoscale = dm.sum()/vscale
    Kscale = pscale/rhoscale**(4/3)
    def unpack(x):
        V = vscale*np.exp(x[:n])
        rr = np.r_[r0,(r0**3+3*np.cumsum(V)/(4*PI))**(1/3)]
        rho = dm/V
        p = Kscale*np.exp(x[-1])*rho**(4/3)
        return rr,rho,p
    def residual(x):
        rr,rho,p = unpack(x)
        support = G*state.enclosed[1:]*state.weights[1:]/(4*PI*rr[1:]**4)
        hydro = np.log(p/(np.r_[p[1:],pext]+support))
        return np.r_[hydro,np.log(rr[-1]/R)]
    x0 = np.r_[np.log(np.diff(q)),0.]
    if n>128:
        if continuum is None:continuum=continuum_polytrope(config,eos)
        state.r,state.cell_volume,rho,p,err=sparse_polytrope(dm,state.bh_mass,r0,R,pext,continuum)
    else:
        sol = least_squares(residual,x0,xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=1500)
        err = float(max(abs(residual(sol.x))))
        if not sol.success or err > 1e-8:
            raise RuntimeError(f"Hydrostatic initialization failed: residual {err:g}")
        state.r,rho,p = unpack(sol.x)
        state.cell_volume = dm/rho
    state.reference_mass_fractions=dm/dm.sum()
    state.T = eos.temperature(rho,p,"p")
    perturb = c["temperature_perturbation"]
    state.T *= 1+perturb*np.exp(-((q[:-1]+q[1:])/2-c["perturbation_mass_center"])**2/
                               (2*c["perturbation_mass_width"]**2))
    state.next_dt = config["numerics"]["initial_dt_years"]*365.25*86400
    return state,err


def continuum_polytrope(config,eos,tolerance=1e-9):
    """Independent continuum n=3 boundary-value solution in log radius.

    y=(m/Menv, h/(G Menv/R)), h=4 K rho^(1/3). Unknown ln K is
    determined by m(rin)=0, m(R)=Menv, P(R)=Pext. This is an independent
    spatial reference and generates a radially resolved Lagrangian mesh.
    """
    import copy
    from scipy.integrate import solve_bvp
    seed=copy.deepcopy(config)
    seed['numerics'].update(zones=64,mesh='core_refined')
    seed['initial']['temperature_perturbation']=0.
    s,_=initialize(seed,eos)
    R=config['initial']['radius_cm'];r0=config['initial']['inner_radius_cm']
    M=s.dm.sum(); b=s.bh_mass/M; vscale=G*M/R; rhoscale=M/(4*PI*R**3)
    K=np.mean(eos.evaluate(s.rho,s.T)['p']/s.rho**(4/3))
    k0=K*rhoscale**(1/3)/vscale
    z=np.linspace(np.log(r0/R),0.,150)
    m=np.interp(np.exp(z)*R,s.r,np.r_[0,np.cumsum(s.dm)])/M
    h=np.interp(np.exp(z)*R,s.rc,4*K*s.rho**(1/3))/vscale
    def fun(z,y,k):
        density=(y[1]/(4*np.exp(k[0])))**3
        return np.array([np.exp(3*z)*density,-(b+y[0])/np.exp(z)])
    def boundary(a,bb,k):
        hsurf=(256*np.exp(k[0])**3*config['boundary']['pressure_cgs']/(rhoscale*vscale))**.25
        return np.array([a[0],bb[0]-1.,bb[1]-hsurf])
    sol=solve_bvp(fun,boundary,z,np.array([m,h]),p=[np.log(k0)],tol=tolerance,max_nodes=10000)
    if not sol.success: raise RuntimeError('Continuum polytrope failed: '+sol.message)
    sol.K=float(np.exp(sol.p[0])*vscale/rhoscale**(1/3))
    sol.rho_scale=rhoscale;sol.velocity_squared_scale=vscale
    return sol

def hydrostatic_residual(state,eos,pext):
    p = eos.evaluate(state.rho,state.T)["p"]
    acceleration = 4*PI*state.r[1:]**2*(p-np.r_[p[1:],pext])/state.weights[1:]
    gravity = G*state.enclosed[1:]/state.r[1:]**2
    return (acceleration-gravity)/gravity
