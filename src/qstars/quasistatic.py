"""Quasi-static Ball et al. 2011 step: hydrostatic, backward-Euler thermal, moving cavity.

bondi_step (spec 3.2): P0 predicts the swept mass; P1 slices it off the base
(cavity.slice_base); P2 re-splits against the STARS Q targets of the accepted
state, mapped in mass above the base in the inner half of the envelope and
in mass below the surface in the outer half, repeating mesh.split_merge
until no cell changes (remesh_base); P3 solves one fully implicit Newton
system; P4 repeats P1-P3 from the accepted state while the implicit cell-0
sink delta exceeds 5% of the base cell (at most 3 passes); P5 books every
transfer. No radial inertia, as in the STARS hydrostatic equations (thesis
2.2); inertia='backward_euler' adds (v - vbar)/dt with v = (r - rbar)/dt for
identity tests and the probe. The predictor maps the accepted state the same
two-sided way: radii from its cumulative volume (each new cell keeps the old
mean density of its mass range), ln T log-linear, r0 scaled to the predicted
r_B, the last-cell T kept on the accepted L_N through the grey relation.

Unknowns (3N+3): y_j = ln(r_j/rbar_j) on faces 0..N; ln(T_k/Tbar_k) on
cells; L_j/L_scale on interior faces; ln rho_ph; delta/dmbar_0; and
X = ln(M_BH^{n+1}/M_BH^n), replaced by ln(dt/dt_guess) in landing steps,
which fix M_BH^{n+1} at the target. Rows: momentum on faces 1..N with
p_N = p_ext; first law per cell with secant volumes; transport closure
on interior faces; piston ln r0 = ln(r_B/b); cavity mass; photosphere;
black-hole growth M_BH^{n+1} - M_BH^n = dt Mdot_BH(face state, M_BH^{n+1}).
First-law rows are scaled by dm u + dt(|L_k| + |L_k+1|) of the predictor.

Momentum rows. 'secant' (spec 3.2a): A_j (p_j-1 - p_j)/w_j = G m_j/(r_j rbar_j),
the discrete gradient of the face potential, so the work identity is exact
without any correction. Its accepted states are, however, off static balance
by 2 dr/r, and in a near-neutral (Gamma1 -> 4/3) envelope that mismatch drives
the next step by about as much again whatever dt is: an equilibrium error
propagates as e(n+1) = e(n)/(1 + 2k/a) for net restoring stiffness k and
gravitational stiffness a. 'endpoint' (default): 4 pi r_j^2 (p_j-1 - p_j)/w_j
= G m_j/r_j^2 at the end-of-step radii, and the gravity quadrature defect
D_j = w_j dr_j G m_j [1/(r_j rbar_j) - A_j/(4 pi r_j^4)] (about
2 G w m (dr/r)^2/r >= 0, the backward-Euler p dV excess) is given back in
the first-law rows of the two adjacent cells (all of D_N in the last cell).
The work identity below then holds to the Newton residual, as for 'secant',
and an equilibrium is preserved bitwise. With inertia the endpoint row
weights (v - vbar)/dt by 4 pi r_j^2/A_j, so its work is w (v - vbar) v.

Masses are end-of-step values: m_0 = M_BH^{n+1} + M_cav^{n+1} and
m_j = m_0 + sum_{k<j} dm_k with dm_0 = dmbar_0 - delta, which equals the
sliced masses minus eps X_acc for j >= 1 (correction C3: one dense column).
M_cav^{n+1} = (transit cavity) + delta - X_acc makes mass exact at every
pass; the cavity row then states M_cav^{n+1} = (8 pi/3) rho_b r0^3.

Energy: E(n+1) - E(n) = -(swept energy) + G eps X_acc sum(w'/rbar)
+ dt L_0 - dt L_N - p_ext dV_N + p_0 dV_0 + G w_0 m_0 (1/rbar_0 - 1/r_0)
+ (split/merge and regrid EOS-inversion defects), every term booked (spec 4). The
identity holds to the Newton residual. The sink removes delta at fixed
volume where a slice would remove it at fixed (rho, u) and let the piston
fill its volume; to first order the difference is the flow work
delta p_face/rho_0, with p_face = p_0 + G m_0 w_0/(4 pi r_0^4) because the
secant piston works at p_0 plus the face-0 potential term. That flow work
is booked to inner_boundary_work rather than to the swept energy, so the
channel split depends on how well P0 predicted the swept mass only at
second order in delta (a relabelling; the sum and the ledger identity are
unchanged). swept_minus_piston_over_bh_energy does not depend on the
labelling at all and is the primary swept-gas diagnostic. Outside
relaxation a step whose ledger residual moves by more than energy_per_step
(of the initial internal energy; 1e-8 over the 2e4-step stop rule) or
mass_per_step (plan conservation.mass_per_step_relative) is rejected.

Base mesh. base='lagrangian' (the registered plan: slice, then hysteresis
re-splits) re-faces the base in discrete events whenever a draining base
cell crosses a split or merge threshold, and the discrete solution near
the coarse base, with the face-extrapolated boundary state, then jumps by
the change of its truncation error, whatever dt is.
base='fixed_mass' is fallback F2 (plan stop_rules; it needs a registered
amendment): after the slice, every cell below the first Lagrangian face at
or above base_mass_fraction of the envelope is remapped (E and V
conserved, as mesh.remap) onto faces at integer levels of the STARS-Q
target count of the accepted state in mass above the base (BaseLevels),
so the base cells change continuously with the state and gas drains
through them (donor-cell advection). Only the top cell of that region
absorbs the Lagrangian drift; split_merge leaves the region alone. A step
with split or merge events reports the boundary-state change those events
cause at that dt (remesh_dln_base_density, _base_temperature,
_inner_radius), from a counterfactual solve on the pre-event mesh, so a
driver can tell a mesh-event jump from evolution.

Failures. Every failure of P1-P3 and of the P5 sink that a shorter step
may cure raises StepRejected (RETRY_ERRORS; RuntimeError is an EOS
temperature inversion in a remap or the sink, or the predictor's
photosphere solve, that did not converge). Its attribute domain_exit is
True when the cause, or a Newton trial point, left the model domain
(DomainError, e.g. Gamma_ph >= 1), with domain_trials and domain_message,
so a driver can apply the plan's domain-exit rule. Argument and
configuration errors stay ValueError; an error at the accepted state
itself (P0) is no rejection and propagates.

Fallback F1 (plan stop_rules.fallbacks_in_order; quasistatic.ptc, off by
default). When the damped Newton solve of an evolve or landing step fails,
the same system is restarted from the predictor with pseudo-transient
continuation (pseudo_transient: Kelley & Keyes 1998, switched evolution
relaxation), and the plain Newton solve then takes over from its iterate.
The pseudo-time term is gone before anything is accepted: a step is
accepted only on the registered residual, tolerance and polishing of
newton(). With ptc False no code path changes; with ptc True a step whose
Newton solve converges is bitwise the same. Relaxation steps never use it.
A step it rescues reports as its Newton iterations the sum of the failed
Newton, pseudo-time and final Newton iterations (final_newton_iterations
separately), so the controller's slow-iteration factor applies to it. F1
acts only on a failed Newton solve, never on a controller rejection.

Photosphere bracket (quasistatic.photosphere_bracket, default 'doubling',
the registered photosphere.solve; 'scan' is an option): the bracket rule of
photosphere.solve where a state carries no photospheric density (the seed's
outer pressure and the predictor of its first step).

Driver (spec 5 and 6, after bondi_step): save_checkpoint/load_checkpoint
(io checkpoint schema 7: every field a bitwise restart needs), zero_time
(t = 0 at the relaxed model), relax (seed relaxation with STARS-Q mesh
cycles), equidistribute, land (dt-unknown landing with the secant
fallback), history_row, save_profile and run, the transactional track
driver; physics (opacity.tables) and run_track, the `qstars run`/`restart`
entry of a quasi_static configuration. The end of a track is classified by
fold.py.
"""
import copy,json,math
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from .constants import G,PI,C,MSUN,A,YEAR,RSUN,LSUN
from .state import State
from .engine import rates_from
from .mesh import q_targets,split_merge,remap
from .nonlinear import colouring,pattern_jacobian
from .numerics import StepRejected
from .mass_transport import DomainError
from .transport import luminosities,inverse_transport_residual
from . import cavity,photosphere

DEFAULTS=dict(
    cavity=dict(active=True,profile='r-3/2',b=1.,boundary_state='face'),
    quasistatic=dict(inertia='none',momentum='endpoint',tolerance=1e-12,max_iterations=40,polish_iterations=1,
                     roundoff_residual=1e-10,roundoff_step=1e-13,max_passes=3,pass_fraction=.05,
                     remesh=True,split=2.,merge=.5,base='lagrangian',base_mass_fraction=.05,
                     energy_per_step=5e-13,mass_per_step=1e-14,ptc=False,photosphere_bracket='doubling'))
MODES=('evolve','relax','land')
BASES=('lagrangian','fixed_mass')
TRIAL_ERRORS=(ValueError,FloatingPointError,ZeroDivisionError,DomainError)
# Failures before or inside the solve that a shorter step may cure (DomainError is handled apart).
# RuntimeError: an EOS temperature inversion (eos.temperature in mesh.remap, split_merge and
# cavity.sink) or photosphere.solve that did not converge.
RETRY_ERRORS=(ValueError,FloatingPointError,ZeroDivisionError,RuntimeError)

def rejection(message,domain_trials=0,domain_message=None,domain_exit=None):
    """StepRejected carrying the domain-exit record (driver classification, plan initial_model.failure)."""
    exc=StepRejected(message)
    exc.domain_trials=int(domain_trials);exc.domain_message=domain_message
    exc.domain_exit=bool(domain_trials) if domain_exit is None else bool(domain_exit)
    return exc

def settings(config):
    """Quasi-static and cavity settings: module defaults updated by config['cavity'] and config['quasistatic']."""
    out={}
    for section,defaults in DEFAULTS.items():
        given=config.get(section) or {}
        unknown=set(given)-set(defaults)
        if unknown:raise ValueError(f'Unknown {section} settings {sorted(unknown)}')
        out[section]={**defaults,**given}
    if out['quasistatic']['inertia'] not in ('none','backward_euler'):raise ValueError('inertia must be none or backward_euler')
    if out['quasistatic']['momentum'] not in ('secant','endpoint'):raise ValueError('momentum must be secant or endpoint')
    if out['cavity']['boundary_state'] not in ('face','cell'):raise ValueError('boundary_state must be face or cell')
    profile,fraction=out['cavity']['profile'],out['quasistatic']['base_mass_fraction']
    if not isinstance(profile,str) or profile not in cavity.PROFILES:raise ValueError('Unknown cavity profile')
    if out['quasistatic']['base'] not in BASES:raise ValueError(f'base must be one of {BASES}')
    if not isinstance(out['quasistatic']['ptc'],bool):raise ValueError('ptc must be true or false')
    if out['quasistatic']['photosphere_bracket'] not in photosphere.BRACKETS:
        raise ValueError(f'photosphere_bracket must be one of {photosphere.BRACKETS}')
    if isinstance(fraction,bool) or not isinstance(fraction,(int,float)) or not 0<fraction<.5:
        raise ValueError('base_mass_fraction must lie in (0, 0.5)')
    return out


class StarsOpacity:
    """Opacity callable for the Ball replication: MESA bridge ABI 7 STARS switch.

    clamp=True applies the STARS nearest-table-value rule, which must be
    reported (validity screen opacity_clamp reads clamp_flags).
    """
    source='stars2004_switch'
    def __init__(self,eos,clamp=False):
        self.eos,self.clamp=eos,bool(clamp)
    def __call__(self,rho,T):
        return self.eos.opacity_stars(rho,T,clamp=self.clamp)
    def clamp_flags(self,rho,T):
        if not self.clamp:return np.zeros(np.broadcast(rho,T).shape,dtype=np.int32)
        return self.eos.opacity_stars(rho,T,clamp=True,flags=True)[1]


def loaded_polytrope(bh_fraction,n=3,b=1.,cavity_coefficient=2/3,rtol=1e-11):
    """Ball 2011 loaded polytrope (paper eqs 24-33) on the small-xi0 branch.

    theta(xi0) = 1 at r0 = r_B/b with Gamma = 1 + 1/n, so the black hole has
    phi = b xi0/(2n) and the cavity cavity_coefficient*xi0^3 (2/3 for
    rho ~ r^-3/2), phi = -xi^2 theta'. xi0 solves M_BH/M* = bh_fraction below
    the maximum (0.1015 for n = 3, reference.ball_limit). Returns the dense
    solution and xi0, xi_s, phi(xi0), phi_s.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq,minimize_scalar
    def rhs(x,y):return [-y[1]/x**2,x**2*max(y[0],0.)**n]
    def surface(x,y):return y[0]
    surface.terminal=True;surface.direction=-1
    def solve(x0,dense=False):
        sol=solve_ivp(rhs,[x0,1e4],[1.,b*x0/(2*n)+cavity_coefficient*x0**3],events=surface,rtol=rtol,
                      atol=rtol*1e-2,max_step=.05,dense_output=dense)
        if len(sol.t_events[0])!=1:raise RuntimeError('Loaded polytrope has no surface')
        return sol
    fraction=lambda x0:b*x0/(2*n)/solve(x0).y[1,-1]
    peak=minimize_scalar(lambda x:-fraction(x),bounds=(1e-3,2.),method='bounded',options=dict(xatol=1e-10)).x
    if not 0<bh_fraction<fraction(peak):raise ValueError('No loaded polytrope with this black-hole mass fraction')
    x0=brentq(lambda x:fraction(x)-bh_fraction,1e-9,peak,xtol=1e-15,rtol=1e-15)
    sol=solve(x0,dense=True)
    return dict(solution=sol,xi0=x0,xi_s=float(sol.t_events[0][0]),phi0=float(sol.y[1,0]),phi_s=float(sol.y[1,-1]))

SEED_SURFACES=('eos','grey')

def seed(config,eos,opacity,zones,bh_mass,total_mass,radius,s=None,fine=40,surface=None):
    """Spec 6 step 1: the n = 3 loaded polytrope scaled to R_seed, discretely hydrostatic on a STARS-Q mesh.

    Continuum: rho = rho_c theta^3, p = K rho^(4/3), K = pi G alpha^2 rho_c^(2/3),
    alpha = R_seed/xi_s; M_cav = (8 pi/3) rho_c r0^3. The outer pressure is
    the photospheric P_ph at the engine-law luminosity of the continuum base
    state (stars_photosphere), else boundary.pressure_cgs. Faces are equally
    spaced in the STARS mesh function Q of the continuum with p floored at
    P_ph, so the tail beyond it forms the last cell. The discrete hydrostatic
    n = 3 equations on these faces (structure.sparse_polytrope, central mass
    M_BH + M_cav, r0 and R_seed fixed, K free) then give r, rho, p, and T
    comes from the EOS at (rho, p). No ledger; kinetic_tracked=False.

    surface (default: replication.seed_surface of the configuration, else
    'eos'). 'eos' is the registered seed: every cell takes the EOS
    temperature at the polytrope (rho, p). The polytrope knows nothing of the
    photosphere, so its last cell, read through the grey surface relation of
    the step (T^4 = 3/4 Teff^4 (tau_half + 4/3)), implies a surface luminosity
    unrelated to the one P_ph was evaluated at; on MESA microphysics that
    luminosity is super-Eddington and no first step exists (plan amendment 3).
    'grey' (stars_photosphere only) changes the temperature of the last cell
    alone: it takes the grey-relation temperature at the seed luminosity (the
    engine-law luminosity P_ph was evaluated at), the lowest root above Teff
    (surface_temperature), so the seed's L_N equals that luminosity and its
    photosphere is the one whose pressure closes the polytrope. Radii,
    densities, masses and every other temperature are those of 'eos'. The
    seed is then no longer in pressure balance in its last cell; like its
    thermal imbalance everywhere else, that is what the relaxation removes.
    """
    from .mesh import stars_q
    from .structure import sparse_polytrope
    s=s or settings(config)
    if surface is None:surface=(config.get('replication') or {}).get('seed_surface','eos')
    if surface not in SEED_SURFACES:raise ValueError(f'seed surface must be one of {SEED_SURFACES}')
    if surface=='grey' and config['boundary']['thermal']!='stars_photosphere':
        raise ValueError('seed surface grey needs boundary.thermal stars_photosphere')
    cav=s['cavity']
    sol=loaded_polytrope(bh_mass/total_mass,b=cav['b'],cavity_coefficient=cavity.PROFILES[cav['profile']]/(4*PI))
    xi0,xis,phi0,phis=sol['xi0'],sol['xi_s'],sol['phi0'],sol['phi_s']
    alpha=radius/xis;rho_c=total_mass/(4*PI*alpha**3*phis);K=PI*G*alpha**2*rho_c**(2/3)
    r0=alpha*xi0;cavity_gas=float(cavity.cavity_mass(rho_c,r0,cav['profile']))
    base=eos.evaluate(np.atleast_1d(rho_c),eos.temperature(np.atleast_1d(rho_c),np.atleast_1d(K*rho_c**(4/3)),'p'))
    law=rates_from(rho_c,base['p'][0],base['cs'][0],base['gamma1'][0],bh_mass,config['engine'],0.)['luminosity']
    if config['boundary']['thermal']=='stars_photosphere':
        rho_ph=photosphere.solve(float(law),radius,total_mass,eos,opacity,bracket=s['quasistatic']['photosphere_bracket'])
        pext=float(photosphere.evaluate(rho_ph,float(law),radius,total_mass,eos,opacity)['pressure'].real)
    else:pext=float(config['boundary']['pressure_cgs'])
    xi=xi0*np.exp(np.linspace(0.,np.log(xis/xi0),fine*zones+1));xi[-1]=xis
    theta,phi=sol['solution'].sol(xi);theta=np.clip(theta,0.,None);theta[-1]=0.
    p=K*(rho_c*theta**3)**(4/3)
    above=p>=max(pext,1e-300)
    if above.sum()<2:raise ValueError('Seed radius leaves no envelope above the photospheric pressure')
    T=np.empty(len(xi));T[above]=eos.temperature(rho_c*theta[above]**3,p[above],'p');T[~above]=T[above][-1]
    Q=stars_q(np.maximum(p,p[above][-1]),T,total_mass*phi/phis,alpha*xi,total_mass)
    if np.any(np.diff(Q)>=0):raise ValueError('STARS mesh function is not monotone on the seed')
    xi_faces=np.interp(np.linspace(Q[0],Q[-1],zones+1)[::-1],Q[::-1],xi[::-1])[::-1]
    xi_faces[0],xi_faces[-1]=xi0,xis
    phi_faces=sol['solution'].sol(xi_faces)[1];phi_faces[0]=phi0;phi_faces[-1]=phis
    dm=total_mass*np.diff(phi_faces)/phis
    if np.any(dm<=0):raise ValueError('Seed mesh has nonpositive cells')
    class Continuum:  # interface of structure.continuum_polytrope used by sparse_polytrope
        pass
    continuum=Continuum();continuum.K=K
    continuum.sol=lambda z:np.array([(sol['solution'].sol(np.clip(xis*np.exp(z),xi0,xis))[1]-phi0)/(phis-phi0)])
    r,V,rho,p,error=sparse_polytrope(dm,bh_mass+cavity_gas,r0,radius,pext,continuum)
    state=State(0.,bh_mass,dm,r,np.zeros(zones+1),eos.temperature(rho,p,'p'),cell_volume=V,
                cavity_mass=cavity_gas,kinetic_tracked=False)
    state.reference_mass_fractions=dm/dm.sum()
    info=dict(xi0=xi0,xi_s=xis,alpha=alpha,rho_c=rho_c,K=K,external_pressure=pext,law_luminosity=float(law),
              hydrostatic_residual=error,envelope_mass=float(dm.sum()),
              mass_defect=float(total_mass-bh_mass-cavity_gas-dm.sum()))
    if surface=='grey':
        polytrope=float(state.T[-1])
        state.T[-1]=surface_temperature(float(law),radius,float(state.rho[-1]),float(dm[-1]),opacity,float(np.max(state.T)))
        info.update(surface='grey',surface_temperature=float(state.T[-1]),surface_polytrope_temperature=polytrope,
                    surface_effective_temperature=float(photosphere.effective_temperature(float(law),radius)))
    return state,info

def surface_temperature(luminosity,radius,rho,dm,opacity,upper,step=.01):
    """Last-cell temperature of the grey surface relation at a given L_N (seed surface 'grey').

    T^4 = 3 L (tau_half + 4/3)/(4 pi R^2 a c) with tau_half = kappa(rho, T) dm/(8 pi R^2):
    the relation the step uses for L_N (stars_photosphere: tau = 2/3 at the
    face, no bath). T = Teff solves it for tau_half = 0, so the relation asks
    for T >= Teff: ln T is marched up from ln Teff in steps of `step` to the
    first sign change, which brentq refines, so the lowest root is returned.
    Raises ValueError when there is none up to `upper` (seed(): the hottest
    polytrope temperature) or the opacity leaves its domain first.
    """
    from scipy.optimize import brentq
    area=4*PI*radius**2
    def g(lnT):
        kappa=float(np.atleast_1d(opacity(np.atleast_1d(rho),np.atleast_1d(np.exp(lnT))))[0])
        return lnT-.25*np.log(3*luminosity*(.5*kappa*dm/area+4/3)/(area*A*C))
    a=float(np.log(photosphere.effective_temperature(luminosity,radius)));fa=g(a)
    while fa<0 and a<np.log(upper):
        b=a+step;fb=g(b)
        if fb>=0:return float(np.exp(brentq(g,a,b,xtol=1e-15,rtol=8.9e-16)))
        a,fa=b,fb
    if fa==0:return float(np.exp(a))
    raise ValueError('Seed surface: no grey-relation temperature between Teff and the hottest seed temperature')

def predict(old,dt,eos,config,mode='evolve',target_mass=None,s=None):
    """P0: swept mass of the step, X_acc + dM_cav, from State.aux['predictor'].

    With two accepted points the end-of-step rest-mass rate and the cavity
    mass are extrapolated linearly in time (second-order local error, as the
    backward-Euler growth uses the end-of-step rate). Without history the law
    at the accepted face state is used, and the cavity change is the law
    cavity mass at r_B of that state minus the stored one.
    """
    s=s or settings(config)
    eps=config['engine']['efficiency'];cav=s['cavity']
    B=cavity.state_boundary(old,eos,cav['boundary_state'])
    law=rates_from(B['rho'],B['p'],B['cs'],B['gamma1'],old.bh_mass,config['engine'],old.time)
    history=old.aux.get('predictor')
    if history and history['time'][-1]==old.time and len(history['time'])==2 and history['time'][1]>history['time'][0]:
        span=history['time'][1]-history['time'][0]
        rate=max(0.,history['rest_rate'][1]+(history['rest_rate'][1]-history['rest_rate'][0])*dt/span)
        change=(history['cavity_mass'][1]-history['cavity_mass'][0])*dt/span
    else:
        rate=float(law['rest_mass_rate'])
        rB=cavity.bondi_radius(old.bh_mass,B['cs'],cav['b'])
        change=float(cavity.cavity_mass(B['rho'],rB,cav['profile'])-old.cavity_mass)
    if mode=='relax':accreted=0.
    elif mode=='land':accreted=(target_mass-old.bh_mass)/(1-eps)
    else:accreted=dt*rate
    if not cav['active']:accreted=change=0.
    return dict(swept=accreted+change,accreted=accreted,rest_rate=rate,growth_rate=(1-eps)*rate,
                cavity_change=change,boundary=B,law=law)


def pattern(n,wide=False):
    """Declared sparsity of the 3N+3 system (a superset; X is one dense column).

    wide: first-law rows also see the neighbouring cells, the outer pressure
    and (cell 1, through face 1) the sink correction: momentum='endpoint'
    moves the face work defect into them.
    """
    Y=lambda j:j;T=lambda k:n+1+k;L=lambda j:2*n+j;PH,D,X=3*n,3*n+1,3*n+2
    base=[Y(0),Y(1),Y(2),T(0),T(1),D,X];outer=[Y(n-1),Y(n),T(n-1)]
    rows=[];cols=[]
    def add(i,js):
        js=sorted(set(js));rows.extend([i]*len(js));cols.extend(js)
    for j in range(1,n+1):
        add(j-1,[Y(j-1),Y(j),T(j-1),X]+([Y(j+1),T(j)] if j<n else outer+[PH])+([D] if j==1 else []))
    for k in range(n):
        near=[Y(j) for j in range(max(0,k-1),min(n,k+2)+1)]+[T(j) for j in range(max(0,k-1),min(n-1,k+1)+1)] if wide else []
        add(n+k,[Y(k),Y(k+1),T(k),X]+([L(k)] if k>=1 else [])+([L(k+1)] if k+1<=n-1 else [])
            +(base if k==0 else [])+(outer if k==n-1 else [])+near+([PH] if wide and k==n-1 else [])
            +([D] if wide and k==1 else []))
    for j in range(1,n):
        add(2*n+j-1,[Y(j-1),Y(j),Y(j+1),T(j-1),T(j),L(j),X]+([D] if j==1 else []))
    add(3*n-1,base);add(3*n,base);add(3*n+1,[PH,X]+outer);add(3*n+2,base)
    return csc_matrix((np.ones(len(rows)),(rows,cols)),shape=(3*n+3,3*n+3))


def newton(fun,x,pattern_matrix,tolerance,max_iterations=40,polish_iterations=1,
           roundoff_residual=1e-10,roundoff_step=1e-13):
    """Damped Newton on a coloured complex-step Jacobian, factored by SuperLU.

    Same backtracking rule as nonlinear.newton (2-norm decrease, 25 halvings).
    After max|f| <= tolerance, polish_iterations further full steps are kept
    when they do not raise max|f|. A stalled line search is a roundoff exit
    (accepted, and counted) when max|f| <= roundoff_residual and the Newton
    step has max|dx| <= roundoff_step. nfev counts every residual evaluation.
    Residual evaluations outside the model domain (DomainError, e.g.
    Gamma_ph >= 1 at a trial point) fail like any trial; domain_trials counts
    them and domain_message keeps the last message.
    """
    colours=int(colouring(pattern_matrix).max())+1
    count=dict(nfev=0,domain=0,message=None)
    def call(z):
        count['nfev']+=1
        try:return fun(z)
        except DomainError as exc:
            count['domain']+=1;count['message']=str(exc);raise
    x=np.array(x,dtype=float);dx_norm=np.inf;polished=0
    info=dict(colours=colours,roundoff_exit=False)
    def result(**kw):
        return dict(info,nfev=count['nfev'],domain_trials=count['domain'],domain_message=count['message'],**kw)
    try:f=call(x)
    except TRIAL_ERRORS as exc:
        return x,result(success=False,iterations=0,residual=np.inf,polished=0,last_step=dx_norm,
                        message=f'Starting point outside the model domain: {exc}')
    for iteration in range(max_iterations+1):
        residual=float(np.max(abs(f))) if np.all(np.isfinite(f)) else np.inf
        if residual<=tolerance:
            for _ in range(polish_iterations):
                try:
                    jac=pattern_jacobian(call,x,pattern_matrix)
                    dx=splu(jac).solve(-f);trial=call(x+dx)
                except (RuntimeError,)+TRIAL_ERRORS:break
                if np.all(np.isfinite(trial)) and np.max(abs(trial))<=residual:
                    x,f,residual=x+dx,trial,float(np.max(abs(trial)));polished+=1;dx_norm=float(np.max(abs(dx)))
            return x,result(success=True,iterations=iteration,residual=residual,
                            polished=polished,last_step=dx_norm,message='Newton residual converged')
        if iteration==max_iterations:break
        try:
            jac=pattern_jacobian(call,x,pattern_matrix)
            dx=splu(jac).solve(-f)
        except (RuntimeError,)+TRIAL_ERRORS as exc:
            return x,result(success=False,iterations=iteration,residual=residual,polished=0,
                            last_step=dx_norm,message=f'Jacobian failure: {exc}')
        if not np.all(np.isfinite(dx)):break
        dx_norm=float(np.max(abs(dx)));norm=np.linalg.norm(f)
        for k in range(25):
            candidate=x+dx*2.**(-k)
            try:trial=call(candidate)
            except TRIAL_ERRORS:continue
            if np.all(np.isfinite(trial)) and np.linalg.norm(trial)<norm*(1-1e-4*2.**(-k)):
                x,f=candidate,trial;break
        else:
            if residual<=roundoff_residual and dx_norm<=roundoff_step:
                return x,result(success=True,iterations=iteration,residual=residual,polished=0,
                                last_step=dx_norm,roundoff_exit=True,message='Roundoff exit')
            break
    return x,result(success=False,iterations=iteration,
                    residual=float(np.max(abs(f))) if np.all(np.isfinite(f)) else np.inf,
                    polished=0,last_step=dx_norm,message='Newton line search or iteration limit')


PTC=dict(delta0=1.,max_iterations=60,handover=1e6,max_shrinks=20)

def pairing(n):
    """Fallback F1: the rows that carry the pseudo-time term and the unknown each is paired with.

    The step is an index-1 differential-algebraic system in time: only the
    first law of cell k (unknown ln T_k) and the black-hole growth row
    (unknown X, or ln dt in a landing) are evolution equations; hydrostatic
    balance, the transport closure, the piston, the cavity mass and the
    photosphere are constraints. Pseudo-time acts on the evolution equations
    alone and the constraints are solved at every iteration (Coffey, Kelley &
    Keyes 2003, SIAM J. Sci. Comput. 25, 553). Returns (rows, cols).
    """
    return np.r_[n+np.arange(n),3*n+2],np.r_[n+1+np.arange(n),3*n+2]

def pseudo_transient(fun,x,pattern_matrix,pairs,tolerance,delta0=1.,max_iterations=60,handover=1e6,max_shrinks=20):
    """Fallback F1: pseudo-transient continuation with switched evolution relaxation.

    Kelley & Keyes (1998, SIAM J. Numer. Anal. 35, 508): x_{k+1} = x_k -
    (V/delta_k + J(x_k))^-1 F(x_k) with delta_{k+1} = delta_k |F(x_k)|/|F(x_{k+1})|
    (2-norms). V is the fixed scaling matrix of the continuation: the entries
    of J at the starting point on the declared (row, unknown) pairs (unit
    where such an entry vanishes) and zero on every other row, so delta is
    dimensionless and the pseudo-time system is V dx/dtau = -F(x), a
    differential-algebraic system when the pairs cover only some rows. Steps
    are taken without a line search, as in the method; a step that cannot be
    evaluated (singular matrix, trial point outside the model domain,
    nonfinite residual) halves delta, at most max_shrinks times in all. The
    continuation ends, with handover True, when delta >= handover (the
    pseudo-time term is then below 1/handover of the paired entries) or
    max|F| <= tolerance: the caller must finish with the plain Newton solve,
    which alone decides convergence. Returns (x, info).
    """
    rows,cols=pairs
    count=dict(nfev=0,domain=0,message=None)
    def call(z):
        count['nfev']+=1
        try:return fun(z)
        except DomainError as exc:
            count['domain']+=1;count['message']=str(exc);raise
    x=np.array(x,dtype=float);delta=float(delta0);shrinks=0;V=None
    def result(**kw):
        return dict(nfev=count['nfev'],domain_trials=count['domain'],domain_message=count['message'],delta=delta,shrinks=shrinks,**kw)
    try:f=call(x)
    except TRIAL_ERRORS as exc:
        return x,result(handover=False,iterations=0,residual=np.inf,message=f'Starting point outside the model domain: {exc}')
    for iteration in range(max_iterations+1):
        residual=float(np.max(abs(f))) if np.all(np.isfinite(f)) else np.inf
        if residual<=tolerance or delta>=handover:
            return x,result(handover=True,iterations=iteration,residual=residual,message='Pseudo-time removed')
        if iteration==max_iterations:break
        try:jac=pattern_jacobian(call,x,pattern_matrix)
        except TRIAL_ERRORS as exc:
            return x,result(handover=False,iterations=iteration,residual=residual,message=f'Jacobian failure: {exc}')
        if V is None:
            paired=np.asarray(jac.tocsr()[rows,cols]).ravel()
            V=csc_matrix((np.where(np.isfinite(paired)&(paired!=0),paired,1.),(rows,cols)),shape=jac.shape)
        norm=np.linalg.norm(f)
        while True:
            try:
                dx=splu((jac+V/delta).tocsc()).solve(-f)
                trial=call(x+dx) if np.all(np.isfinite(dx)) else None
            except (RuntimeError,)+TRIAL_ERRORS:trial=None
            if trial is not None and np.all(np.isfinite(trial)):break
            delta*=.5;shrinks+=1
            if shrinks>max_shrinks:
                return x,result(handover=False,iterations=iteration,residual=residual,message='Pseudo-time step collapsed')
        x,f=x+dx,trial
        new=np.linalg.norm(f)
        delta=delta*norm/new if new>0 else float(handover)
    return x,result(handover=False,iterations=max_iterations,residual=residual,message='Pseudo-transient iteration limit')


class Problem:
    """P3 for one pass: residual, declared pattern, predictor and assembly from a post-remesh state."""
    def __init__(self,old,post,dt,eos,opacity,config,mode,target_mass,predicted,s):
        self.old,self.post,self.dt,self.eos,self.opacity,self.config=old,post,dt,eos,opacity,config
        self.mode,self.s,self.predicted,self.target_mass=mode,s,predicted,target_mass
        self.n=n=len(post.dm)
        if n<3:raise ValueError('Quasi-static step needs at least three cells')
        self.inertia=s['quasistatic']['inertia']=='backward_euler'
        self.endpoint=s['quasistatic']['momentum']=='endpoint'
        if post.kinetic_tracked!=self.inertia:
            raise ValueError('Quasi-static states carry no kinetic energy (kinetic_tracked=False); inertia needs it tracked')
        self.active=s['cavity']['active']
        self.photosphere=config['boundary']['thermal']=='stars_photosphere'
        self.eps=config['engine']['efficiency']
        self.M=old.bh_mass
        q=eos.evaluate(post.rho,post.T)
        self.q=q
        self.mbar=post.enclosed
        self.gbar=G*self.mbar/post.r**2
        self.carry=post.radius_compensation if post.radius_compensation is not None else np.zeros(n+1)
        self.storage=post.dm*q['u']
        self.contrast_scale=np.maximum(abs(np.diff(np.log(q['p'])))*.25,1e-3)
        self.bc_config=copy.deepcopy(config);self.bc_config['transport']['convection']='off'
        self.pext=float(config['boundary']['pressure_cgs'])
        if mode=='land':
            if not target_mass>self.M:raise ValueError('Landing target must exceed the current black-hole mass')
            self.X_land=float(np.log(target_mass/self.M))
        self.pattern=pattern(n,wide=self.endpoint)
        self.x0=self.initial()

    # ------------------------------------------------------------------ residual
    def evaluate(self,x,parts=False):
        n,post,eos=self.n,self.post,self.eos
        y,th,ell=x[:n+1],x[n+1:2*n+1],x[2*n+1:3*n]
        lnrho_ph,dhat,z=x[3*n],x[3*n+1],x[3*n+2]
        if self.mode=='land':X,dt=self.X_land,self.dt*np.exp(z)
        else:X,dt=z,self.dt
        time=self.old.time+np.real(dt)
        growth=self.M*np.expm1(X)
        accreted=growth/(1-self.eps) if self.mode!='relax' else 0.*growth
        delta=dhat*post.dm[0]
        dm=post.dm.astype(np.result_type(x,float));dm[0]=post.dm[0]-delta
        m0=self.mbar[0]+delta-self.eps*accreted
        dr=post.r*np.expm1(y)
        r=post.r+(dr-self.carry)
        if np.any(np.diff(r.real)<=0) or not r.real[0]>0:raise ValueError('Trial geometry has crossing radial faces')
        area=4*PI/3*(r*r+r*post.r+post.r*post.r)
        work=area*dr
        V=post.volume+np.diff(work)
        if np.any(V.real<=0):raise ValueError('Trial cell volume is nonpositive')
        T=post.T*np.exp(th)
        v=dr/dt
        st=State(time,m0,dm,r,v,T,cell_volume=V)
        m=st.enclosed
        rho=dm/V
        q=eos.evaluate(rho,T);p=q['p']
        B=cavity.boundary_state(rho[:2],T[:2],p[:2],r[0],m0,dm[0],eos,self.s['cavity']['boundary_state'])
        mass_new=self.M+growth if self.mode!='relax' else self.M+0.*growth
        law=rates_from(B['rho'],B['p'],B['cs'],B['gamma1'],mass_new,self.config['engine'],time)
        L0=self.eps*C**2*accreted/dt if self.mode!='relax' else law['luminosity']
        LN=luminosities(st,eos,self.opacity,self.bc_config,0.)[0][-1]
        if self.photosphere:
            ph=photosphere.evaluate(np.exp(lnrho_ph),LN,r[-1],m[-1],eos,self.opacity)
            pext=ph['pressure'];row_ph=ph['residual']
        else:
            ph=None;pext=self.pext;row_ph=lnrho_ph-self.x0[3*n]
        L=np.r_[L0,ell*self.Lscale,LN]
        w=.5*np.r_[dm[0],dm[:-1]+dm[1:],dm[-1]]
        dp=p-np.r_[p[1:],pext]
        force=area[1:]*dp/w[1:]
        secant=G*m[1:]/(r[1:]*post.r[1:])
        acceleration=(v[1:]-post.v[1:])/dt if self.inertia else 0.
        if self.endpoint:
            # Static balance at the end-of-step radii. The secant work identity then
            # leaves D_j = w_j dr_j (G m_j/(r_j rbar_j) - A_j dp_j/w_j + accel_j), which
            # on the row equals the gravity quadrature defect w dr G m (1/(r rbar)
            # - A/(4 pi r^4)), the backward-Euler p dV excess; the first-law rows of
            # the two adjacent cells give it back. Written with dp, the energy
            # identity holds at every iterate. Inertia is weighted by 4 pi r^2/A_j
            # so that its work is exactly w (v - vbar) v.
            shell=4*PI*r[1:]**2
            momentum=G*m[1:]/r[1:]**2-shell*dp/w[1:]+shell/area[1:]*acceleration
            defect=(secant-force+acceleration)*w[1:]*dr[1:]
            settling=np.zeros(n,dtype=defect.dtype)
            settling[:-1]+=.5*defect[:-1];settling[1:]+=.5*defect[:-1];settling[-1]+=defect[-1]
        else:
            momentum=secant-force+acceleration;settling=0.
        momentum=momentum/self.gbar[1:]
        internal=dm*(q['u']-self.q['u'])
        internal[0]=internal[0]+delta*(self.q['h'][0]-self.q['u'][0])
        heating=0.
        if self.inertia:
            dv=v-post.v;dv[0]=0.  # face 0 has no momentum row; its kinetic energy is piston work
            heating=.25*dm*(dv[:-1]**2+dv[1:]**2)
        energy=(internal+p*np.diff(work)+dt*np.diff(L)-heating+settling)/self.scale_energy
        closure=inverse_transport_residual(st,eos,self.opacity,self.config,L)/self.contrast_scale
        cav=self.s['cavity']
        if self.active:
            piston=np.log(r[0])-np.log(cavity.bondi_radius(mass_new,B['cs'],cav['b']))
            cavity_new=post.cavity_mass+delta-accreted
            row_cav=(cavity_new-cavity.cavity_mass(B['rho'],r[0],cav['profile']))/post.dm[0]
        else:
            piston,row_cav=y[0],dhat
        if self.mode=='relax':row_X=z
        else:row_X=(np.expm1(X)-dt*law['bh_growth_rate']/self.M)/self.X_scale
        f=np.r_[momentum,energy,closure,piston,row_cav,row_ph,row_X]
        if not parts:return f
        return f,dict(y=y,theta=th,dr=dr,r=r,area=area,work=work,V=V,T=T,v=v,dm=dm,m=m,w=w,p=p,q=q,
                      L=L,L0=L0,LN=LN,pext=pext,photosphere=ph,boundary=B,law=law,heating=heating,settling=settling,
                      delta=delta,accreted=accreted,growth=growth,dt=dt,time=time,mass_new=mass_new,
                      lnrho_ph=lnrho_ph)

    __call__=evaluate

    # ------------------------------------------------------------------ predictor
    def initial(self):
        """State n mapped in mass above the base (r, T log-linear) below half the envelope and
        in mass below the surface above it (so split children follow the old profile); r0
        scaled to the predicted r_B; face luminosities carried."""
        n,old,post,eos=self.n,self.old,self.post,self.eos
        x=np.zeros(3*n+3)
        half=.5*post.dm.sum()
        mu_old=np.r_[0.,np.cumsum(old.dm)];mu=np.r_[0.,np.cumsum(post.dm)]
        centre_old=mu_old[:-1]+.5*old.dm;centre=mu[:-1]+.5*post.dm
        def mapped(points,grid,values,total_old,total):
            # mass above the base in the inner half, mass below the surface in the outer half
            inner=np.interp(points,grid,values)
            outer=np.interp((total-points)[::-1],(total_old-grid)[::-1],values[::-1])[::-1]
            return np.where(points<=half,inner,outer)
        B=self.predicted['boundary']
        if self.mode=='evolve':growth=self.dt*self.predicted['growth_rate']
        elif self.mode=='land':growth=float(np.expm1(self.X_land))*self.M
        else:growth=0.
        r0=cavity.bondi_radius(self.M+growth,B['cs'],self.s['cavity']['b']) if self.active else old.r[0]
        # Radii from the accepted cumulative volume (exactly piecewise linear in mass
        # for uniform cells), so every new cell gets the old mean density of its mass range.
        volume=np.r_[0.,np.cumsum(old.volume)]
        inner=np.cbrt(r0**3+3*np.interp(mu,mu_old,volume)/(4*PI))
        above=volume[-1]-volume
        outer=np.cbrt(old.r[-1]**3-3*np.interp((mu[-1]-mu)[::-1],(mu_old[-1]-mu_old)[::-1],above[::-1])[::-1]/(4*PI))
        r=np.where(mu<=half,inner,outer);r[0]=r0;r[-1]=old.r[-1]
        if np.all(np.diff(r)>0):x[:n+1]=np.log(r/post.r)
        else:x[:n+1]=mapped(mu,mu_old,np.log(old.r),mu_old[-1],mu[-1])-np.log(post.r)
        x[n+1:2*n+1]=mapped(centre,centre_old,np.log(old.T),mu_old[-1],mu[-1])-np.log(post.T)
        x[2*n]=self.outer_temperature(old,x)-np.log(post.T[-1])
        L,_=luminosities(post,eos,self.opacity,self.bc_config,0.)
        if post.face_luminosity is not None and len(post.face_luminosity)==n+1:
            interior=np.array(post.face_luminosity[1:-1],dtype=float)
        else:
            law=self.predicted['law']['luminosity']
            interior=np.interp(mu[1:-1],[0.,mu[-1]],[float(np.real(law)),float(L[-1])])
        self.Lscale=max(float(np.max(abs(np.r_[interior,L[-1]]))),1.)
        x[2*n+1:3*n]=interior/self.Lscale
        # First-law rows are scaled by storage plus the flux terms they difference:
        # where dt exceeds a cell's thermal time, dt L_{k+1} - dt L_k cancels to
        # roundoff of dt L, far above any tolerance relative to dm u alone.
        if self.mode=='relax':inner=float(np.real(self.predicted['law']['luminosity']))
        elif self.mode=='land':inner=self.eps*C**2*self.predicted['accreted']/self.dt
        else:inner=self.eps*C**2*self.predicted['rest_rate']
        flux=abs(np.r_[inner,interior,float(L[-1])])
        self.scale_energy=self.storage+self.dt*(flux[:-1]+flux[1:])
        if self.photosphere:
            rho_ph=old.photosphere_density
            if rho_ph is None:
                rho_ph=photosphere.solve(float(L[-1]),float(post.r[-1]),float(post.enclosed[-1]),eos,self.opacity,
                                         bracket=self.s['quasistatic']['photosphere_bracket'])
            x[3*n]=np.log(rho_ph)
        elif post.photosphere_density is not None:x[3*n]=np.log(post.photosphere_density)
        if self.mode=='evolve':
            x[3*n+2]=np.log1p(self.dt*self.predicted['growth_rate']/self.M)
            self.X_scale=max(abs(float(np.expm1(x[3*n+2]))),1e-16)
        elif self.mode=='land':self.X_scale=abs(float(np.expm1(self.X_land)))
        else:self.X_scale=1.
        return x

    def outer_temperature(self,old,x):
        """ln T of the last cell that keeps the accepted L_N through the grey surface relation.

        A remeshed surface cell changes tau_half; at the mapped temperature the
        predicted L_N would jump (and Gamma_ph could pass 1 at the predictor).
        """
        n,post=self.n,self.post
        lnT=np.log(post.T[-1])+x[2*n]
        bc=self.config['boundary']
        if bc['thermal'] not in ('stars_photosphere','grey'):return lnT
        if old.face_luminosity is not None:L=float(old.face_luminosity[-1])
        else:L=float(luminosities(old,self.eos,self.opacity,self.bc_config,0.)[0][-1])
        external=2/3 if bc['thermal']=='stars_photosphere' else bc['external_optical_depth']
        bath=0. if bc['thermal']=='stars_photosphere' else bc['bath_temperature_K']
        r=post.r*np.exp(x[:n+1]);R=r[-1]
        rho=post.dm[-1]/(4*PI/3*(R**3-r[-2]**3))
        for _ in range(4):
            tau=.5*float(np.atleast_1d(self.opacity(np.atleast_1d(rho),np.atleast_1d(np.exp(lnT))))[0])*post.dm[-1]/(4*PI*R**2)+external
            T4=3*L*(tau+2/3)/(4*PI*R**2*A*C)+bath**4
            if not T4>0:return np.log(post.T[-1])+x[2*n]
            lnT=.25*np.log(T4)
        return lnT

    # ------------------------------------------------------------------ solve
    def solve(self):
        qs=self.s['quasistatic']
        solve=lambda x:newton(self.evaluate,x,self.pattern,qs['tolerance'],qs['max_iterations'],
                              qs['polish_iterations'],qs['roundoff_residual'],qs['roundoff_step'])
        x,info=solve(self.x0)
        if not qs['ptc'] or self.mode=='relax':return x,info
        if info['success']:return x,dict(info,ptc_entered=False,ptc_iterations=0)
        # Fallback F1: pseudo-transient continuation from the predictor, then the plain Newton solve.
        xp,ptc=pseudo_transient(self.evaluate,self.x0,self.pattern,pairing(self.n),qs['tolerance'],**PTC)
        final=None
        if ptc['handover']:x2,final=solve(xp)
        parts=[info,ptc]+([final] if final else [])
        extra=dict(nfev=sum(a['nfev'] for a in parts),domain_trials=sum(a['domain_trials'] for a in parts),
                   domain_message=next((a['domain_message'] for a in parts[::-1] if a['domain_message']),None),
                   ptc_entered=True,ptc_iterations=ptc['iterations'],ptc_delta=ptc['delta'],ptc_shrinks=ptc['shrinks'],
                   ptc_handover=ptc['handover'],ptc_message=ptc['message'],newton_message=info['message'],
                   newton_residual=info['residual'])
        if final and final['success']:
            # The registered controller rule 'extra factor 0.8 if Newton needed more than 8 iterations'
            # (stepcontrol.accept) sees every nonlinear iteration the step needed: the failed Newton
            # solve, the pseudo-time iterations and the final Newton solve.
            return x2,dict(final,**extra,final_newton_iterations=final['iterations'],
                           iterations=info['iterations']+ptc['iterations']+final['iterations'])
        why=(f"handed over, Newton then failed: {final['message']} (residual {final['residual']:g})" if final
             else f"failed: {ptc['message']}")
        return x,dict(info,**extra,message=f"{info['message']}; F1 pseudo-transient continuation {why}")

    def delta(self,x):
        return float(x[3*self.n+1]*self.post.dm[0])


def targets(old,cut,eos,zones):
    """P2 targets: STARS-Q cell masses of the accepted state, mapped onto the sliced state.

    ln dm* is interpolated in mass above the base in the inner half of the
    envelope, where the base structure is quasi-steady in that coordinate,
    and in mass below the surface in the outer half: a slice or returned gas
    shifts the base coordinate of every surface cell by the swept mass, far
    more than the thin surface cells hold.
    """
    lnt=np.log(q_targets(old,eos,zones))
    centre_old=np.cumsum(old.dm)-.5*old.dm;centre=np.cumsum(cut.dm)-.5*cut.dm
    total_old,total=old.dm.sum(),cut.dm.sum()
    inner=np.interp(centre,centre_old,lnt)
    outer=np.interp((total-centre)[::-1],(total_old-centre_old)[::-1],lnt[::-1])[::-1]
    return np.exp(np.where(centre<=.5*total,inner,outer))

class BaseLevels:
    """Target count below a mass above the base: N(mu) = int_0^mu dmu'/t(mu').

    t is the STARS-Q target cell mass of the accepted state (q_targets),
    ln t linear in mu between its cell centres and constant below the first
    and beyond the last, so N is exact (exponential segments), continuous
    and strictly increasing, and inverse() is closed form. Faces at integer
    N give cells of one target each that follow the state continuously.
    """
    def __init__(self,old,eos,zones):
        lnt=np.log(q_targets(old,eos,zones))
        self.mu=np.r_[0.,np.cumsum(old.dm)-.5*old.dm]
        self.a=np.r_[lnt[0],lnt]
        self.b=np.diff(self.a)/np.diff(self.mu)
        self.N=np.r_[0.,np.cumsum(self._count(self.a[:-1],self.b,np.diff(self.mu)))]

    @staticmethod
    def _count(a,b,x):
        """int_0^x exp(-(a + b s)) ds."""
        safe=np.where(b!=0,b,1.)
        return np.exp(-a)*np.where(b!=0,-np.expm1(-b*x)/safe,x)

    def __call__(self,mu):
        mu=np.asarray(mu,dtype=float)
        k=np.clip(np.searchsorted(self.mu,mu,side='right')-1,0,len(self.mu)-1)
        last=k==len(self.mu)-1;kk=np.where(last,len(self.b)-1,k)
        return self.N[k]+self._count(self.a[k],np.where(last,0.,self.b[kk]),mu-self.mu[k])

    def inverse(self,level):
        level=np.asarray(level,dtype=float)
        k=np.clip(np.searchsorted(self.N,level,side='right')-1,0,len(self.N)-1)
        last=k==len(self.N)-1;kk=np.where(last,len(self.b)-1,k)
        a=self.a[k];b=np.where(last,0.,self.b[kk]);dN=level-self.N[k]
        z=b*np.exp(a)*dN
        return self.mu[k]+np.where(b!=0,-np.log1p(-z)/np.where(b!=0,b,1.),np.exp(a)*dN)

def regrid_base(old,cut,eos,zones,fraction):
    """F2: remap the base region of the sliced state onto faces at fixed target levels.

    The region is every cell below the first face of `cut` at or above
    `fraction` of its envelope mass; that face stays Lagrangian (bitwise,
    with everything above it). Interior faces sit at N = 1, 2, ..., m of
    BaseLevels(old) with N(top) - m in [0.5, 1.5), so only the top region
    cell absorbs the drift of the Lagrangian face. Cells conserve volume and
    total energy (mesh.remap, overlap shares; T by EOS inversion); the
    inversion defect is booked to remap_energy_defect and the re-faced
    potential is reported as potential_reallocation (Amendment A1).
    Returns (state, info) with info['region_cells'].
    """
    levels=BaseLevels(old,eos,zones)
    mu=np.cumsum(cut.dm)
    j=int(np.searchsorted(mu,fraction*mu[-1]))
    if j>=len(mu)-1:raise ValueError('Fixed-mass base region would span the whole envelope')
    top=float(levels(mu[j]))
    count=max(int(math.floor(top-.5)),0)
    faces=levels.inverse(np.arange(1,count+1)) if count else np.zeros(0)
    if np.any(np.diff(np.r_[0.,faces,mu[j]])<=0):raise ValueError('Fixed-mass base faces are not increasing')
    region=np.diff(np.r_[0.,faces,mu[j]])
    region[-1]=math.fsum(cut.dm[:j+1])-math.fsum(region[:-1])
    new_dm=np.r_[region,cut.dm[j+1:]]
    before=cut.energies(eos);entropy=float(np.dot(cut.dm,eos.evaluate(cut.rho,cut.T)['s']))
    out,info=remap(cut,new_dm,eos,copy_unchanged=True,return_slab=True)
    if out.ledger:out.ledger['remap_energy_defect']=out.ledger.get('remap_energy_defect',0.)+info['reconstruction_defect']
    after=out.energies(eos)
    return out,dict(region_cells=len(region),region_replaced_cells=j+1,region_mass=float(mu[j]),
                    energy_defect=info['reconstruction_defect'],max_event_defect=abs(float(info['reconstruction_defect'])),
                    entropy_change=float(np.dot(out.dm,eos.evaluate(out.rho,out.T)['s']))-entropy,
                    potential_reallocation=after['gravity']-before['gravity'])

def remesh_base(old,cut,eos,zones,split=2.,merge=.5,max_passes=40,protect=0):
    """P2: split_merge against targets() until no cell changes.

    One pass cuts a cell into equal-mass parts against its own average
    target; after a slice the new base cell spans the steep base targets,
    so repeated passes (targets re-mapped at the new cell centres) restore
    the graded base. The first `protect` cells (a fixed-mass base region)
    are their own targets, so they are never split or merged themselves.
    Returns the state and the summed pass info.
    """
    total=dict(splits=0,merges=0,energy_defect=0.,entropy_change=0.,merge_mixing_entropy=0.,potential_reallocation=0.,
               passes=0,max_event_defect=0.)
    state=cut
    for _ in range(max_passes):
        goal=targets(old,state,eos,zones)
        if protect:goal[:protect]=state.dm[:protect]
        state,info=split_merge(state,eos,goal,split,merge)
        if not info['splits'] and not info['merges']:return state,total
        total['passes']+=1
        for key in ('splits','merges','energy_defect','entropy_change','merge_mixing_entropy','potential_reallocation'):
            total[key]+=info[key]
        total['max_event_defect']=max(total['max_event_defect'],info['max_event_defect'])
    raise StepRejected('Base re-split did not settle')

def bondi_step(old,dt,eos,opacity,config,mode='evolve',target_mass=None):
    """One quasi-static Ball step (P0-P5). Returns (state, diagnostics).

    mode: 'evolve' (rest mass swallowed, L_0 = eps c^2 X_acc/dt), 'relax'
    (M_BH frozen, L_0 from the engine law, booked as relaxation_injected_energy),
    or 'land' (M_BH^{n+1} = target_mass exactly; dt is solved, dt is its guess).
    Raises StepRejected (see rejection()) for every failure of P1-P3 or of
    the P5 sink (RETRY_ERRORS), more than max_passes passes, or a per-step
    conservation budget exceeded outside relaxation; ValueError only for
    invalid arguments or settings, or from P0 at an invalid accepted state.
    Nothing in `old` is modified.
    """
    if mode not in MODES:raise ValueError(f'Unknown quasi-static mode {mode}')
    if mode=='land' and target_mass is None:raise ValueError('Landing needs a target black-hole mass')
    if mode=='land' and not target_mass>old.bh_mass:raise ValueError('Landing target must exceed the current black-hole mass')
    s=settings(config);cav,qs=s['cavity'],s['quasistatic']
    if not cav['active'] and mode!='relax':raise ValueError('A frozen cavity sweeps no gas; use relax mode')
    if not dt>0:raise ValueError('Timestep must be positive')
    if old.kinetic_tracked!=(qs['inertia']=='backward_euler'):
        raise ValueError('Quasi-static states carry no kinetic energy (kinetic_tracked=False); inertia needs it tracked')
    predicted=predict(old,dt,eos,config,mode,target_mass,s)
    swept=predicted['swept']
    zones=config['numerics']['zones']
    fixed=qs['remesh'] and qs['base']=='fixed_mass'
    for passes in range(1,qs['max_passes']+1):
        stage='P1 base slice'
        try:
            if cav['active']:cut,cut_info=cavity.slice_base(old,swept,eos)
            else:cut,cut_info=old.copy(),dict(mass=0.,volume=0.,entropy=0.,energy=0.,copied_cells=len(old.dm),
                                              fixed_state_cells=0,reconstructed_cells=0,reconstruction_defect=0.)
            stage='P2 base remesh'
            regrid=None;premesh=cut
            if fixed:
                premesh,regrid=regrid_base(old,cut,eos,zones,qs['base_mass_fraction'])
                post,remesh=remesh_base(old,premesh,eos,zones,qs['split'],qs['merge'],protect=regrid['region_cells'])
            elif qs['remesh']:post,remesh=remesh_base(old,cut,eos,zones,qs['split'],qs['merge'])
            else:
                post,remesh=cut,dict(splits=0,merges=0,energy_defect=0.,entropy_change=0.,merge_mixing_entropy=0.,
                                     potential_reallocation=0.,max_event_defect=0.)
            stage='P3 predictor'
            problem=Problem(old,post,dt,eos,opacity,config,mode,target_mass,predicted,s)
        except StepRejected:raise
        except DomainError as exc:
            raise rejection(f'{stage} left the model domain (pass {passes}): {exc}',domain_message=str(exc),
                            domain_exit=True) from exc
        except RETRY_ERRORS as exc:
            raise rejection(f'{stage} failed (pass {passes}): {exc}') from exc
        x,solver=problem.solve()
        if not solver['success']:
            raise rejection(f"Quasi-static solve failed (pass {passes}, residual {solver['residual']:g}): {solver['message']}"
                            +(f" [{solver['domain_trials']} trial points outside the model domain: {solver['domain_message']}]"
                              if solver['domain_trials'] else ''),solver['domain_trials'],solver['domain_message'])
        delta=problem.delta(x)
        if abs(delta)<=qs['pass_fraction']*post.dm[0]:break
        swept+=delta
    else:
        raise rejection(f'Sink correction is {delta/post.dm[0]:.3g} of the base cell after {qs["max_passes"]} predictor passes')
    effect=_remesh_effect(problem,x,premesh) if remesh['splits'] or remesh['merges'] else None
    remesh=dict(remesh,regrid=regrid,effect=effect)
    new,out=_assemble(problem,x,solver,cut_info,remesh,passes,swept)
    if mode!='relax':
        energy,mass=abs(out['energy_residual_change']),abs(out['mass_residual_change'])
        if not (energy<=qs['energy_per_step'] and mass<=qs['mass_per_step']):
            raise rejection(f'Per-step conservation budget exceeded: energy {energy:.3g} (limit {qs["energy_per_step"]:g}), '
                            f'mass {mass:.3g} (limit {qs["mass_per_step"]:g})')
    return new,out

def _remesh_effect(problem,x,premesh):
    """Boundary-state change caused by P2's split/merge events at this dt.

    The accepted pass is compared with the same step (same prediction, one
    pass) solved on the mesh before those events: the sliced state, or the
    regridded state of a fixed-mass base, whose continuous regrid is part
    of the scheme. Only called on steps with a split or merge. A failed
    counterfactual solve is reported, not raised.
    """
    p=problem
    try:
        alt=Problem(p.old,premesh,p.dt,p.eos,p.opacity,p.config,p.mode,p.target_mass,p.predicted,p.s)
        xa,info=alt.solve()
    except (StepRejected,DomainError)+RETRY_ERRORS as exc:
        return dict(solved=False,message=f'{type(exc).__name__}: {exc}')
    if not info['success']:return dict(solved=False,message=info['message'])
    d=p.evaluate(x,parts=True)[1];da=alt.evaluate(xa,parts=True)[1]
    ratio=lambda key:float(np.log(np.real(d['boundary'][key])/np.real(da['boundary'][key])))
    return dict(solved=True,message=None,dln_base_density=ratio('rho'),dln_base_temperature=ratio('T'),
                dln_inner_radius=float(np.log(np.real(d['r'][0])/np.real(da['r'][0]))))


def _assemble(problem,x,solver,cut_info,remesh,passes,predicted_swept):
    """P5: the accepted state, its ledger bookings and the step diagnostics."""
    from . import ledger as ledger_module
    old,post,eos,mode=problem.old,problem.post,problem.eos,problem.mode
    _,d=problem.evaluate(x,parts=True)
    real=lambda a:np.real(a)
    dt=float(real(d['dt']));delta=float(real(d['delta']));accreted=float(real(d['accreted']))
    eps=problem.eps
    new=post.copy()
    new.time=old.time+dt;new.steps=old.steps+1
    new.dm=real(d['dm']).copy();new.r=real(d['r']).copy();new.cell_volume=real(d['V']).copy()
    new.T=real(d['T']).copy();new.v=real(d['v']).copy()
    new.radius_compensation=(new.r-post.r)-(real(d['dr'])-problem.carry)
    new.face_luminosity=real(d['L']).copy()
    new.bh_mass=old.bh_mass+float(real(d['growth'])) if mode!='relax' else old.bh_mass
    new.cavity_mass=post.cavity_mass+delta-accreted
    if problem.photosphere:new.photosphere_density=float(np.exp(real(d['lnrho_ph'])))
    new.reference_mass_fractions=new.dm/new.dm.sum()
    # Bookings (spec 4). 1: slice and sink, by functional difference. The sink takes
    # delta h_0 at fixed volume; its flow work delta p_face/rho_0 is what the piston
    # (p_0 dV_0 plus the face-0 potential term) would have done on a slice of delta,
    # so it goes to inner_boundary_work.
    try:
        sinked=cavity.sink(post,delta,eos)
        sink_energy=post.energies(eos)['total']-sinked.energies(eos)['total']
    except RETRY_ERRORS as exc:
        raise rejection(f'P5 sink failed (pass {passes}): {exc}') from exc
    p_face=float(problem.q['p'][0])+G*post.enclosed[0]*post.weights[0]/(4*PI*post.r[0]**4)
    flow=delta*p_face/post.rho[0]
    swept=cut_info['mass']+delta
    swept_energy=cut_info['energy']+sink_energy-flow
    swept_entropy=cut_info['entropy']+delta*float(problem.q['s'][0])
    l=new.ledger
    booked=cavity.bookings(swept_energy,swept_entropy,swept,accreted,eps)
    for key,value in booked.items():
        l[key]=l.get(key,0.)+value
    # 2: the eps share leaves the central mass at fixed rbar.
    defect=cavity.defect_work(sinked,eps*accreted)
    l['central_potential_work']+=defect
    # 3: geometry and heat at the new masses (secant work identity).
    w,m,dr=real(d['w']),real(d['m']),real(d['dr'])
    piston=cavity.piston_work(float(real(d['p'][0])),new.r[0],post.r[0],w[0],m[0],dr=dr[0])
    if problem.inertia:piston+=.5*w[0]*(new.v[0]**2-post.v[0]**2)
    l['inner_boundary_work']=l.get('inner_boundary_work',0.)+piston-flow
    pext=float(real(d['pext']));LN=float(real(d['LN']));L0=float(real(d['L0']))
    l['boundary_work']-=pext*float(real(d['work'][-1]))
    l['radiated_energy']+=dt*LN
    if mode=='relax':l['relaxation_injected_energy']=l.get('relaxation_injected_energy',0.)+dt*L0
    else:l['central_generated']+=eps*C**2*accreted
    heating=float(np.sum(real(d['heating']))) if problem.inertia else 0.
    if problem.inertia:l['numerical_kinetic_heating']+=heating
    # Diagnostic only (inside the cells, like numerical_kinetic_heating): the
    # backward-Euler p dV excess removed by the endpoint momentum form.
    settling=float(math.fsum(real(d['settling']))) if problem.endpoint else 0.
    if problem.endpoint:l['settling_work']=l.get('settling_work',0.)+settling
    # Driver state for the predictor (JSON, checkpointed with the model).
    B,law=d['boundary'],d['law']
    history=old.aux.get('predictor') if old.aux.get('predictor',{}).get('time',[None])[-1]==old.time else None
    point=dict(time=new.time,rest_rate=float(real(law['rest_mass_rate'])),cavity_mass=new.cavity_mass)
    keep=[dict(time=history['time'][-1],rest_rate=history['rest_rate'][-1],cavity_mass=history['cavity_mass'][-1])] if history else \
         [dict(time=old.time,rest_rate=float(real(problem.predicted['law']['rest_mass_rate'])),cavity_mass=old.cavity_mass)]
    keep.append(point)
    new.aux=copy.deepcopy(old.aux)
    new.aux['predictor']={k:[float(p_[k]) for p_ in keep] for k in ('time','rest_rate','cavity_mass')}
    # Diagnostics.
    q=d['q'];p=real(d['p']);ph=d['photosphere']
    n=problem.n
    mu_faces=np.r_[0.,np.cumsum(post.dm)];mu_cells=mu_faces[:-1]+.5*post.dm
    exclude=2*abs(swept)
    y,th=real(d['y']),real(d['theta'])
    reference=eos.evaluate(sinked.rho,sinked.T)['p']
    inflow=float(math.fsum(-.5*(p-reference)*np.diff(real(d['work']))))
    LBH_dt=L0*dt
    gamma1=real(q['gamma1'])
    regrid=remesh.get('regrid') or {}
    effect=remesh.get('effect')
    jump=effect or dict(solved=None,message=None,dln_base_density=0.,dln_base_temperature=0.,dln_inner_radius=0.)
    per_bh=lambda value:value/LBH_dt if LBH_dt and value is not None else None
    out=dict(dt=dt,mode=mode,passes=passes,newton_iterations=solver['iterations'],nfev=solver['nfev'],
             colours=solver['colours'],nonlinear_residual=solver['residual'],roundoff_exit=solver['roundoff_exit'],
             polished=solver['polished'],domain_trials=solver['domain_trials'],domain_message=solver['domain_message'],
             base_region_cells=regrid.get('region_cells',0),base_region_mass=regrid.get('region_mass',0.),
             regrid_energy_defect=regrid.get('energy_defect',0.),regrid_entropy_change=regrid.get('entropy_change',0.),
             regrid_potential_reallocation=regrid.get('potential_reallocation',0.),
             regrid_potential_reallocation_over_bh_energy=per_bh(regrid.get('potential_reallocation',0.)),
             # Boundary-state change caused by split/merge events (counterfactual solve on the
             # pre-event mesh; zero when there was none, None when that solve failed).
             remesh_effect_solved=jump['solved'],remesh_effect_message=jump['message'],
             remesh_dln_base_density=jump.get('dln_base_density'),
             remesh_dln_base_temperature=jump.get('dln_base_temperature'),
             remesh_dln_inner_radius=jump.get('dln_inner_radius'),
             predicted_swept_mass=float(problem.predicted['swept']),final_swept_prediction=float(predicted_swept),
             delta=delta,delta_fraction=delta/post.dm[0],swept_mass=swept,accreted_rest_mass=accreted,
             rest_mass_rate=accreted/dt,bh_luminosity=L0,luminosity=LN,radius=float(new.r[-1]),
             effective_temperature=float(photosphere.effective_temperature(LN,new.r[-1])),
             inner_radius=float(new.r[0]),bondi_radius=float(real(cavity.bondi_radius(real(d['mass_new']),B['cs'],problem.s['cavity']['b']))),
             base_density=float(real(B['rho'])),base_temperature=float(real(B['T'])),
             base_sound_speed=float(real(B['cs'])),base_gamma1=float(real(B['gamma1'])),
             law_luminosity=float(real(law['luminosity'])),law_bh_growth_rate=float(real(law['bh_growth_rate'])),
             cavity_mass=new.cavity_mass,
             cavity_mass_law=float(real(cavity.cavity_mass(B['rho'],new.r[0],problem.s['cavity']['profile']))),
             slice_mass=cut_info['mass'],slice_cells_copied=cut_info['copied_cells'],
             splits=remesh['splits'],merges=remesh['merges'],remesh_energy_defect=remesh['energy_defect'],
             # Largest |EOS-inversion defect| of one reconstruction event booked to remap_energy_defect
             # (a split/merge group, or the fixed-mass base regrid); plan reconstruction_defect_per_event_relative.
             max_event_defect=max(remesh['max_event_defect'],regrid.get('max_event_defect',0.)),
             remesh_entropy_change=remesh['entropy_change'],merge_mixing_entropy=remesh['merge_mixing_entropy'],
             potential_reallocation=remesh['potential_reallocation'],
             potential_reallocation_over_bh_energy=remesh['potential_reallocation']/LBH_dt if LBH_dt else None,
             swept_energy=swept_energy,sink_energy=sink_energy,sink_flow_work=flow,defect_potential_work=defect,
             piston_work=piston,inner_boundary_work=piston-flow,
             swept_energy_over_bh_energy=swept_energy/LBH_dt if LBH_dt else None,
             cavity_swept_energy_over_bh_energy=booked['cavity_swept_energy']/LBH_dt if LBH_dt else None,
             # Primary swept-gas diagnostic: independent of the channel labelling.
             swept_minus_piston_over_bh_energy=(swept_energy-(piston-flow))/LBH_dt if LBH_dt else None,
             inflow_heating=inflow,inflow_heating_over_bh_energy=inflow/LBH_dt if LBH_dt else None,
             numerical_kinetic_heating_erg=heating,settling_work=settling,
             settling_work_over_bh_energy=settling/LBH_dt if LBH_dt else None,
             lagrangian_max_dlnr=float(np.max(abs(y[mu_faces>=exclude]),initial=0.)),
             lagrangian_max_dlnT=float(np.max(abs(th[mu_cells>=exclude]),initial=0.)),
             mapped_max_dlnr=float(np.max(abs(y[1:]-problem.x0[1:n+1]))),
             mapped_max_dlnT=float(np.max(abs(th-problem.x0[n+1:2*n+1]))),
             stability_integral=float(math.fsum((3*gamma1-4)*p*new.cell_volume)),
             stability_ratio=float(math.fsum((3*gamma1-4)*p*new.cell_volume)/math.fsum(p*new.cell_volume)),
             zones=n)
    if 'ptc_entered' in solver:  # fallback F1 switched on (quasistatic.ptc): whether this step needed it
        out.update({k:solver[k] for k in ('ptc_entered','ptc_iterations','ptc_delta','ptc_shrinks','ptc_handover','ptc_message',
                                         'newton_message','newton_residual','final_newton_iterations') if k in solver})
    if ph is not None:
        out.update(photosphere_density=new.photosphere_density,photosphere_gas_pressure=float(real(ph['gas_pressure'])),
                   photosphere_eddington_factor=float(real(ph['eddington_factor'])),
                   photosphere_opacity=float(real(ph['kappa'])),external_pressure=pext,
                   implied_atmosphere_mass_msun=float(real(photosphere.implied_column(new.r[-1],pext,real(ph['gravity']))))/MSUN)
    # Plain floats (not numpy scalars): a checkpoint stores them as JSON, and Python's float
    # sum() rounds differently when numpy scalars are mixed in, so restarts stay bitwise.
    new.ledger={k:float(v) for k,v in new.ledger.items()}
    before,after=ledger_module.residuals(old,eos),ledger_module.residuals(new,eos)
    out.update(mass_residual_change=after['mass_residual_relative']-before['mass_residual_relative'],
               energy_residual_change=after['energy_residual_relative']-before['energy_residual_relative'],
               engine_conversion_residual_erg=after['engine_conversion_residual_erg'])
    return new,out


# ================================================================== driver (spec 5 and 6)
# Checkpoint, relaxation and the transactional track driver. Nothing here is
# used by bondi_step; default behaviour of every other module is untouched.

RELAX_DEFAULTS=dict(dt_years=.01,dt_max_years=1e8,converged_dt_years=1e6,tolerance=1e-10,luminosity_tolerance=1e-6,
                    max_steps=400,max_retries=20,mesh_cycles=40,mesh_tolerance=1e-5,cycle_dt_years=1e6)
# config 'replication' section: the published model a quasi-static track reproduces (defaults: the
# Ball et al. 2011 fiducial track) and the run's registration labels. plan names a registered plan
# (plan_id); such runs execute only through scripts/validate_ball_track.py. end_bh_mass_msun 0: none.
REPLICATION_DEFAULTS=dict(plan='',run='',kind='',plan_text='',policy='published_assumption',total_mass_msun=1e4,
                          initial_bh_mass_msun=5.,seed_radii_rsun=[2000.,6000.,1000.,12000.],seed_surface='eos',
                          landings_msun=[100.,200.,500.,1000.],end_bh_mass_msun=0.,amendments=[])
LEDGER_COLUMNS=None  # filled lazily from ledger.FIELDS + ledger.CAVITY_FIELDS

def _json(data):
    """JSON-safe copy: nonfinite floats become None, numpy scalars plain numbers."""
    if isinstance(data,dict):return {str(k):_json(v) for k,v in data.items()}
    if isinstance(data,(list,tuple)):return [_json(v) for v in data]
    if isinstance(data,(bool,np.bool_)):return bool(data)
    if isinstance(data,(int,np.integer)):return int(data)
    if isinstance(data,(float,np.floating)):return float(data) if math.isfinite(data) else None
    if isinstance(data,np.ndarray):return _json(data.tolist())
    if isinstance(data,complex):return _json(data.real)
    return data

def save_checkpoint(path,state,config):
    """Atomic checkpoint of a quasi-static state: io.save_checkpoint, schema 7.

    Holds every field a bitwise restart needs (plan restart.requirement):
    r, T, dm, v, cell_volume, radius_compensation, face_luminosity, bh_mass,
    cavity_mass, photosphere_density, kinetic_tracked, aux (controller,
    predictor and driver state), ledger, steps and the configuration.
    """
    from . import io
    io.save_checkpoint(path,state,config,schema=7)

def load_checkpoint(path):
    """(state, config) of a checkpoint; validates shapes, finiteness and geometry (io.load_checkpoint).

    The configuration is returned as stored (test toys and driver configurations
    are compared as saved); `qstars restart` resolves it.
    """
    from . import io
    return io.load_checkpoint(path,resolve_config=False)

def zero_time(state,eos):
    """Spec 6.5: t = 0 at the relaxed model; ledger re-initialised (cavity channels), driver state cleared."""
    from . import ledger
    out=state.copy();out.time=0.;out.steps=0;out.next_dt=0.;out.aux={}
    ledger.initialize(out,eos,cavity=True)
    return out

def law_at(state,eos,config,s=None):
    """Engine law (thesis eq 3.15) and boundary state of an accepted state."""
    s=s or settings(config)
    B=cavity.state_boundary(state,eos,s['cavity']['boundary_state'])
    law=rates_from(B['rho'],B['p'],B['cs'],B['gamma1'],state.bh_mass,config['engine'],state.time)
    return {k:float(np.real(v)) for k,v in law.items()},{k:float(np.real(v)) for k,v in B.items()}

def face_transport(state,eos,opacity,config,L=None):
    """Radiative/convective split of the face luminosities (the closure that solved them).

    Interior faces only; for the mlt.py closures the gradient follows from
    nabla_rad of the solved L, as in the replication screens.
    """
    from .transport import _mixing_length_faces
    from .mlt import MODELS,closure,conductance,radiative_gradient
    q=eos.evaluate(state.rho,state.T);kap=opacity(state.rho,state.T)
    if L is None:L=state.face_luminosity if state.face_luminosity is not None else luminosities(state,eos,opacity,config,0.)[0]
    L=np.asarray(L,dtype=float);model=config['transport']['convection']
    faces,mass=_mixing_length_faces(state,q,kap)
    nabla_rad=radiative_gradient(L[1:-1],faces['kappa'],faces['p'],faces['T'],mass)
    if model in MODELS:
        out=closure(model,nabla_rad,**faces,alpha=config['transport']['mixing_length_alpha'])
        radiative=conductance(faces['kappa'],faces['p'],faces['T'],mass)*out['gradient'];velocity=out['velocity']
    else:
        _,t=luminosities(state,eos,opacity,config,0.)
        radiative=t['radiative'][1:-1];velocity=np.zeros(len(radiative))
    convective=L[1:-1]-radiative
    return dict(luminosity=L,radiative=radiative,convective=convective,velocity=velocity,nabla_rad=nabla_rad,
                nabla_ad=faces['grad_ad'],convective_face=nabla_rad>faces['grad_ad'],q=q,opacity=kap)

def save_profile(path,state,eos,opacity,config):
    """Full profile of an accepted state (landings and run ends; plan numerics.landing.outputs)."""
    t=face_transport(state,eos,opacity,config);q=t['q']
    np.savez_compressed(path,time_seconds=state.time,bh_mass_g=state.bh_mass,cavity_mass_g=state.cavity_mass,
        photosphere_density_cgs=np.nan if state.photosphere_density is None else state.photosphere_density,
        dm_g=state.dm,r_edges_cm=state.r,enclosed_mass_edges_g=state.enclosed,density_cgs=state.rho,
        pressure_cgs=np.real(q['p']),temperature_K=state.T,entropy_cgs=np.real(q['s']),internal_energy_erg_g=np.real(q['u']),
        gamma1=np.real(q['gamma1']),grad_ad=np.real(q['grad_ad']),cp_cgs=np.real(q['cp']),opacity_cgs=np.real(t['opacity']),
        luminosity_edges_erg_s=t['luminosity'],radiative_interior_edges_erg_s=t['radiative'],
        convective_interior_edges_erg_s=t['convective'],convective_velocity_interior_edges_cm_s=t['velocity'],
        nabla_rad_interior_edges=t['nabla_rad'],nabla_ad_interior_edges=t['nabla_ad'],
        convective_interior_edges=t['convective_face'])

def structure_summary(state,eos,opacity,config):
    """Plan qualitative-check quantities of one state (Q2-Q5; thesis 3.2.1 stability integral).

    A cell is convective when either of its faces is (nabla_rad > nabla_ad of the solved L).
    """
    t=face_transport(state,eos,opacity,config);q=t['q']
    conv=np.r_[False,t['convective_face'],False]
    cells=conv[:-1]|conv[1:]
    p=np.real(q['p']);gamma1=np.real(q['gamma1'])
    outer=np.flatnonzero(conv)
    above=float(math.fsum(state.dm[outer[-1]:])) if len(outer) else float(math.fsum(state.dm))
    return dict(max_density_inversion=float(max(np.diff(state.rho)/state.rho[:-1],default=0.)),
                stability_integral=float(math.fsum((3*gamma1-4)*p*state.volume)),
                stability_ratio=float(math.fsum((3*gamma1-4)*p*state.volume)/math.fsum(p*state.volume)),
                mass_above_outer_convective_boundary_msun=above/MSUN,
                convective_mass_fraction=float(math.fsum(state.dm[cells])/math.fsum(state.dm)),
                convective_gamma1=float(math.fsum(state.dm[cells]*gamma1[cells])/math.fsum(state.dm[cells])) if np.any(cells) else None)

def _ledger_columns():
    global LEDGER_COLUMNS
    if LEDGER_COLUMNS is None:
        from .ledger import FIELDS,CAVITY_FIELDS
        LEDGER_COLUMNS=list(FIELDS)+list(CAVITY_FIELDS)+['settling_work','initial_mass','initial_energy','initial_energy_scale']
    return LEDGER_COLUMNS

STEP_COLUMNS=('dt','mode','passes','newton_iterations','nfev','colours','nonlinear_residual','roundoff_exit','polished',
              'domain_trials','delta_fraction','swept_mass','accreted_rest_mass','slice_mass','splits','merges',
              'base_region_cells','remesh_dln_base_density','remesh_dln_base_temperature','remesh_dln_inner_radius',
              'remesh_energy_defect','max_event_defect','remesh_entropy_change','merge_mixing_entropy','potential_reallocation_over_bh_energy',
              'regrid_energy_defect','regrid_potential_reallocation_over_bh_energy','swept_minus_piston_over_bh_energy',
              'cavity_swept_energy_over_bh_energy','inflow_heating_over_bh_energy','settling_work_over_bh_energy',
              'mass_residual_change','energy_residual_change','engine_conversion_residual_erg','lagrangian_max_dlnr',
              'lagrangian_max_dlnT','mapped_max_dlnr','mapped_max_dlnT','implied_atmosphere_mass_msun','landing_method')
METRIC_COLUMNS=('delta_actual','profile_dlnT','profile_dlnr','dln_inner_radius','dln_base_density','dln_base_temperature',
                'dln_bh_mass','dln_luminosity','dln_radius','dln_effective_temperature','cavity_change')

def history_row(state,eos,opacity,config,event='step',diag=None,metrics=None,screens=None,target=None,s=None):
    """One history row (fixed columns): Ball observables, step diagnostics, screens and the ledger.

    io.history_row returns this row for quasi_static configurations.

    Rows are a pure function of the state and the step record, so a
    restarted run writes bitwise identical rows. `screens` is the diagnostics
    dict of validity.check with a policy (None: no screen columns are filled).
    """
    from . import validity
    from .ledger import residuals
    s=s or settings(config)
    law,B=law_at(state,eos,config,s)
    L=float(state.face_luminosity[-1]) if state.face_luminosity is not None else float(np.real(luminosities(state,eos,opacity,config,0.)[0][-1]))
    q=eos.evaluate(state.rho,state.T)
    p=np.real(q['p']);gamma1=np.real(q['gamma1'])
    R=float(state.r[-1]);total=float(state.enclosed[-1])
    row=dict(step=state.steps,event=event,time_years=state.time/YEAR,time_myr=state.time/YEAR/1e6,
             bh_mass_msun=state.bh_mass/MSUN,landed_target_msun=None if target is None else target/MSUN,
             cavity_mass_msun=state.cavity_mass/MSUN,envelope_mass_msun=float(math.fsum(state.dm))/MSUN,total_mass_msun=total/MSUN,
             zones=len(state.dm),mdot_bh_msun_per_year=law['bh_growth_rate']*YEAR/MSUN,
             rest_rate_msun_per_year=law['rest_mass_rate']*YEAR/MSUN,law_luminosity_lsun=law['luminosity']/LSUN,
             luminosity_lsun=L/LSUN,luminosity_over_law=L/law['luminosity'] if law['luminosity'] else None,
             base_density_cgs=B['rho'],base_temperature_K=B['T'],base_gamma1=B['gamma1'],base_sound_speed_cm_s=B['cs'],
             inner_radius_rsun=float(state.r[0])/RSUN,bondi_radius_rsun=law['bondi_radius']/RSUN,radius_rsun=R/RSUN,
             # The outer face is the tau = 2/3 photosphere only with the Ball surface (thesis 2.15).
             photospheric_radius_cm=R if config['boundary']['thermal']=='stars_photosphere' else float('nan'),
             effective_temperature_K=float(photosphere.effective_temperature(L,R)),
             photosphere_density_cgs=state.photosphere_density,
             dynamical_time_years=math.sqrt(R**3/(G*total))/YEAR,
             stability_ratio=float(math.fsum((3*gamma1-4)*p*state.volume)/math.fsum(p*state.volume)),
             max_density_inversion=float(max(np.diff(state.rho)/state.rho[:-1],default=0.)))
    diag=diag or {};metrics=metrics or {}
    for k in STEP_COLUMNS:
        v=diag.get(k)
        row['step_'+k]=v/YEAR if k=='dt' and v is not None else v
    for k in METRIC_COLUMNS:row['metric_'+k]=metrics.get(k)
    sc=(screens or {}).get('screens',{})
    row['screen_policy']=(screens or {}).get('policy')
    for k in validity.SCREENS:
        entry=sc.get(k)
        row['screen_'+k]=None if entry is None else entry['value']
        row['screen_'+k+'_violated']=None if entry is None else entry['violated']
    row['photosphere_eddington_factor']=(screens or {}).get('photosphere_eddington_factor')
    row['opacity_clamped_mass_msun']=None if (screens or {}).get('opacity_clamped_mass_g') is None else screens['opacity_clamped_mass_g']/MSUN
    row['screen_warnings']=';'.join((screens or {}).get('screen_warnings',[]))
    row.update(residuals(state,eos))
    for k in _ledger_columns():row['ledger_'+k]=float(state.ledger.get(k,0.)) if state.ledger else None
    return row

def _sanitized_row(row):
    return {k:(None if isinstance(v,float) and not math.isfinite(v) else v) for k,v in row.items()}

def relax(state,eos,opacity,config,s=None,log=None):
    """Spec 6 relaxation of a seed (M_BH frozen; mode 'relax' of bondi_step).

    dt starts at dt_years and doubles to dt_max_years after every accepted
    step, halving on a rejection (at most max_retries in a row). Converged
    when two consecutive accepted steps with dt >= converged_dt_years change
    ln T and ln r by at most `tolerance`, measured cell by cell between the
    consecutive states (a fixed-mass base moves its faces between steps until
    its regrid reaches a fixed point; a step that changes the cell count is
    never quiet) and the Lagrangian change, and |L*/L_BH - 1| <= luminosity_tolerance.

    Mesh cycles (mesh_cycles > 0): the converged model is remapped onto
    `zones` cells equally spaced in the STARS-Q target count of itself
    (equidistribute) and relaxed again (from cycle_dt_years: the remap is a
    small conservative rezone of a relaxed model), until that remap moves no
    face by more than mesh_tolerance of its cell mass (max_face_shift) or the
    equidistributed faces move by less than that between consecutive cycles
    (target_change: with a fixed-mass base the regrid keeps its region on
    integer target levels, so the relaxed mesh is a fixed point of the cycle,
    not of the remap alone). The cycles contract linearly (by about 0.5 per
    cycle on the toys, slowest at the base). The relaxed mesh is then a fixed
    point of the STARS mesh function of the relaxed model, so seeds that
    start from different meshes relax to the same discrete model (plan
    initial_model.uniqueness); without it each seed keeps its own Lagrangian
    mesh and they differ by the truncation error.

    Returns (state, record); record['status'] is 'converged' or 'failed', with
    every step, rejection (its domain record) and cycle and, for a failure,
    the number of consecutive rejections it ended with (final_rejection_streak)
    and whether they were all domain exits (DomainError at a trial point or
    before the solve; the driver checks which domain).
    """
    rs={**RELAX_DEFAULTS,**(s or {})}
    zones=config['numerics']['zones']
    record=dict(settings=rs,steps=[],rejections=[],cycles=[],status='failed',reason=None)
    def note(kind,**data):
        data=_json(dict(kind=kind,**data))
        (record['steps'] if kind=='accepted' else record['rejections'] if kind=='rejected' else record['cycles']).append(data)
        if log:log(kind,data)
    previous=None
    for cycle in range(rs['mesh_cycles']+1):
        state,ok,streak=_relax_once(state,eos,opacity,config,rs,note,cycle,rs['dt_years'] if cycle==0 else rs['cycle_dt_years'])
        if not ok:
            last=record['rejections'][-streak:] if streak else []
            record.update(reason=f'relaxation failed in mesh cycle {cycle}',final_rejection_streak=streak,
                          final_rejections_all_domain_exits=bool(last) and all(r['domain_exit'] for r in last))
            return state,record
        if not rs['mesh_cycles']:break
        new,info=equidistribute(state,eos,zones)
        faces=info.pop('faces')
        info['target_change']=float(np.max(abs(faces-previous)/info.pop('local_mass'))) if previous is not None and len(faces)==len(previous) else np.inf
        note('mesh_cycle',cycle=cycle,**info)
        if min(info['max_face_shift'],info['target_change'])<=rs['mesh_tolerance']:break
        state,previous=new,faces
    else:
        record.update(reason=f"mesh not at its STARS-Q fixed point after {rs['mesh_cycles']} cycles")
        return state,record
    record.update(status='converged',reason=None,final=record['steps'][-1] if record['steps'] else None)
    return state,record

def _relax_once(state,eos,opacity,config,rs,note,cycle,dt_years):
    dt=dt_years*YEAR;quiet=0;retries=0
    for _ in range(rs['max_steps']):
        try:new,d=bondi_step(state,dt,eos,opacity,config,mode='relax')
        except StepRejected as exc:
            retries+=1
            note('rejected',cycle=cycle,dt_years=dt/YEAR,message=str(exc),domain_exit=getattr(exc,'domain_exit',False),
                 domain_trials=getattr(exc,'domain_trials',0),domain_message=getattr(exc,'domain_message',None))
            if retries>rs['max_retries']:return state,False,retries
            dt*=.5;continue
        retries=0
        same=len(new.dm)==len(state.dm)
        change=max(np.max(abs(np.log(new.r/state.r))),np.max(abs(np.log(new.T/state.T))),
                   d['lagrangian_max_dlnr'],d['lagrangian_max_dlnT']) if same else np.inf
        balance=abs(d['luminosity']/d['bh_luminosity']-1)
        quiet=quiet+1 if dt>=rs['converged_dt_years']*YEAR and change<=rs['tolerance'] and balance<=rs['luminosity_tolerance'] else 0
        note('accepted',cycle=cycle,dt_years=dt/YEAR,change=change,luminosity_balance=balance,zones=len(new.dm),
             passes=d['passes'],newton_iterations=d['newton_iterations'],splits=d['splits'],merges=d['merges'],
             photosphere_eddington_factor=d.get('photosphere_eddington_factor'),inner_radius=d['inner_radius'],
             base_density=d['base_density'],luminosity=d['luminosity'],radius=d['radius'],roundoff_exit=d['roundoff_exit'])
        state=new
        if quiet>=2:return state,True,0
        dt=min(2*dt,rs['dt_max_years']*YEAR)
    return state,False,retries

def equidistribute(state,eos,zones):
    """Conservative remap onto `zones` cells equally spaced in the STARS-Q target count of this state.

    Faces at levels k N_total/zones (k = 1..zones-1) of BaseLevels(state),
    the exact integral of 1/target in mass above the base; cells conserve E
    and V (mesh.remap, overlap shares), the EOS-inversion defect is booked to
    remap_energy_defect. max_face_shift is the largest face move in units of
    the smaller adjacent new cell mass (inf when the cell count changes);
    info also holds the new interior faces (mass above the base) and those
    cell masses (faces, local_mass).
    """
    levels=BaseLevels(state,eos,zones)
    mu=np.cumsum(state.dm);total=float(levels(mu[-1]))
    faces=levels.inverse(total*np.arange(1,zones)/zones)
    if np.any(np.diff(np.r_[0.,faces,mu[-1]])<=0):raise ValueError('Equidistributed faces are not increasing')
    new_dm=np.diff(np.r_[0.,faces,mu[-1]]);new_dm[-1]=math.fsum(state.dm)-math.fsum(new_dm[:-1])
    local=np.minimum(new_dm[:-1],new_dm[1:])
    shift=float(np.max(abs(faces-mu[:-1])/local)) if len(state.dm)==zones else np.inf
    out,info=remap(state,new_dm,eos,copy_unchanged=True,return_slab=True)
    if out.ledger:
        out.ledger['remap_energy_defect']=out.ledger.get('remap_energy_defect',0.)+info['reconstruction_defect']
        out.ledger={k:float(v) for k,v in out.ledger.items()}
    out.photosphere_density=state.photosphere_density
    return out,dict(max_face_shift=shift,levels=total,zones_before=len(state.dm),zones_after=zones,
                    energy_defect=info['reconstruction_defect'],faces=faces,local_mass=local)

def growth_rate(state,eos,config,s=None):
    """Black-hole growth rate of the engine law at an accepted state (landing prediction)."""
    return law_at(state,eos,config,s)[0]['bh_growth_rate']

def land(state,target,eos,opacity,config,tolerance=1e-6,max_secant=8,s=None):
    """Landing step (plan numerics.landing): M_BH^{n+1} = target exactly, ln dt solved.

    If the dt-unknown landing fails, the registered fallback is a secant on dt
    over evolve steps from the same accepted state, to |M_BH/target - 1| <=
    tolerance. Raises StepRejected when both fail.
    """
    guess=(target-state.bh_mass)/growth_rate(state,eos,config,s)
    try:
        new,d=bondi_step(state,guess,eos,opacity,config,mode='land',target_mass=target)
        d=dict(d,landing_method='solve')
    except StepRejected as first:
        def miss(dt):
            new,d=bondi_step(state,dt,eos,opacity,config,mode='evolve')
            return new,d,new.bh_mass/target-1
        dt0=guess;new,d,f0=miss(dt0);iterations=1
        if abs(f0)>tolerance:
            dt1=dt0*(target-state.bh_mass)/(new.bh_mass-state.bh_mass)
            for _ in range(max_secant):
                new,d,f1=miss(dt1);iterations+=1
                if abs(f1)<=tolerance:break
                dt0,f0,dt1=dt1,f1,dt1-f1*(dt1-dt0)/(f1-f0)
            else:raise rejection(f'Landing on {target/MSUN:g} Msun failed: dt solve ({first}); secant missed by {f1:.3g}')
        d=dict(d,landing_method='secant',landing_secant_steps=iterations,landing_solve_failure=str(first))
    if not abs(new.bh_mass/target-1)<=tolerance:
        raise rejection(f'Landing missed {target/MSUN:g} Msun by {new.bh_mass/target-1:.3g}')
    return new,d

def run(state,eos,opacity,config,path,landings=(),end_bh_mass=None,policy='published_assumption',
        restart_masses=(),restart_landing=None,segment_steps=None,stop_after_steps=None,max_wall_seconds=None,
        max_accepted_steps=20000,dt_floor=None,checkpoint_every=10):
    """Transactional quasi-static Ball track (spec 5; modelled on evolution.run). Returns (state, report).

    `state` is an accepted model (the t = 0 relaxed model, or a checkpoint for
    a restart); `path` an existing run directory. Each attempt: the next
    dt is the controller's (stepcontrol, State.aux['controller']); a step that
    would reach the next pending landing mass (one not in
    State.aux['driver']['landed']), predicted with the growth rate of the
    accepted state or found to reach it, becomes a landing step (land()); a
    candidate is rejected when bondi_step or the landing fails, delta_actual
    > 2 delta_target, or an enforced screen of `policy` fires (validity.check
    with the solved face luminosities). A landing's dt is forced by its
    target, so the controller continues from max(landing dt, attempted dt)
    and the dt < t_dyn trigger always compares t_dyn with a controlled step.
    A rejected candidate never touches the accepted state; only the
    controller record (retries, hold) and the driver's count of consecutive
    enforced-screen rejections in its aux change, and both are part of the
    checkpointed state, so a restart continues bitwise.

    Ends (report['status']):
      completed            end_bh_mass reached (landed when it is a landing mass);
      termination_trigger  dt < t_dyn (Ball's criterion) or 12 consecutive
                           retries (plan numerics.controller): classify with fold.py;
      model_domain_exit    12 consecutive retries, every one of them rejected by
                           an enforced screen (a termination trigger too: the
                           driver still applies the fold protocol); or the start
                           row itself violates an enforced screen of `policy`
                           (report['start_row_screen_errors'], no step attempted:
                           for a t = 0 start a model_domain_exit of row 0) other
                           than nonfinite_state, which is a numerical_failure
                           (start_row_screen_errors set, no step attempted);
      stop_rule            wall time, accepted-step limit or dt floor (plan
                           stop_rules): a numerical blocker;
      control_stop         stop_after_steps or an interrupt;
      numerical_failure    any other error.
    Writes events.jsonl, history.csv (history_row), checkpoint.npz (every
    checkpoint_every accepted steps and at the end), landing_<M>msun.npz and
    profile_landing_<M>msun.npz, restart_<label>.npz at restart_masses (first
    accepted step at or above each mass) and, when the restart_landing
    landing is accepted, of the accepted state before it (aux as accepted,
    so a segment replays the rejections in between), segment_end_<step>.npz
    segment_steps after each restart point, profile_final.npz,
    screen_violations.json and termination.json. Landing masses at or below
    the starting M_BH that were never landed are an argument error.
    """
    import time
    from pathlib import Path
    from . import stepcontrol,validity
    from . import io
    path=Path(path);path.mkdir(parents=True,exist_ok=True)
    s=settings(config);S=stepcontrol.settings(config);mode=s['cavity']['boundary_state']
    state=state.copy()
    restarted='controller' in state.aux
    if not restarted:stepcontrol.start(state.aux,S)
    state.aux.setdefault('driver',dict(floor_run=0,landed=[],screen_retries=0))
    landings=tuple(sorted(float(x) for x in landings))
    unreachable=[t/MSUN for t in landings if t not in state.aux['driver']['landed'] and not t>state.bh_mass]
    if unreachable:raise ValueError(f'Landing masses {unreachable} Msun lie at or below the starting black-hole mass and were never landed')
    pending=lambda:[t for t in landings if t not in state.aux['driver']['landed']]
    # aux of the last accepted state as accepted: rejections update state.aux in place.
    accepted_aux=copy.deepcopy(state.aux)
    clock=time.perf_counter();accepted=attempts=rejections=0
    status,reason='numerical_failure','initialization did not complete'
    last_rejection=None;segment_ends={};saved=set();summary={};start_errors=[]
    def log(event,**data):
        with (path/'events.jsonl').open('a') as f:f.write(json.dumps(_json(dict(event=event,**data)),allow_nan=False)+'\n')
    def screened(candidate):
        if policy is None:return [],None
        errors,diag,_=validity.check(candidate,eos,opacity,config,policy=policy,face_luminosity=candidate.face_luminosity)
        return errors,diag
    def record(row):
        row=_sanitized_row(row);io.append_history(path/'history.csv',row)
        for k in validity.SCREENS:
            if row.get('screen_'+k+'_violated'):
                e=summary.setdefault(k,dict(first_time_years=row['time_years'],first_bh_mass_msun=row['bh_mass_msun'],rows=0,extreme=None))
                e['rows']+=1;v=row.get('screen_'+k)
                if v is not None and (e['extreme'] is None or abs(v)>abs(e['extreme'])):e['extreme']=v
    try:
        if policy is not None and not (path/'validity_warnings.json').exists():
            io.json_write(path/'validity_warnings.json',_json(validity.describe(policy,config,eos,opacity)))
        errors,diag=screened(state)
        record(history_row(state,eos,opacity,config,event='restart' if restarted else 'start',screens=diag,s=s))
        log('start',step=state.steps,restarted=restarted,time_years=state.time/YEAR,bh_mass_msun=state.bh_mass/MSUN,
            enforced_screen_errors=errors)
        before=stepcontrol.observables(state,eos,mode)
        if errors:
            # The start row is outside the always-enforced domain: no step is attempted from it. A
            # nonfinite state is a numerical defect of the model, never a domain exit.
            start_errors=list(errors)
            if ((diag or {}).get('screens',{}).get('nonfinite_state') or {}).get('violated'):
                status,reason='numerical_failure','nonfinite model state at the start row: '+'; '.join(errors)
            else:
                status,reason='model_domain_exit',('enforced screens violated at the start row'
                                                   +(' (restart row)' if restarted else ' (row 0: model_domain_exit of row 0)')+': '+'; '.join(errors))
        while not start_errors:
            if end_bh_mass is not None and (state.bh_mass>=end_bh_mass or end_bh_mass in state.aux['driver']['landed']):
                status,reason='completed',f'black-hole mass {end_bh_mass/MSUN:g} Msun reached';break
            if stop_after_steps is not None and accepted>=stop_after_steps:
                status,reason='control_stop','requested checkpoint stop';break
            if state.steps>=max_accepted_steps:
                status,reason='stop_rule',f'maximum accepted steps {max_accepted_steps} reached (numerical blocker)';break
            if max_wall_seconds and time.perf_counter()-clock>=max_wall_seconds:
                status,reason='stop_rule',f'per-track wall time {max_wall_seconds/3600:g} h reached (numerical blocker)';break
            if dt_floor and state.aux['driver']['floor_run']>=dt_floor['steps']:
                status,reason='stop_rule',(f"dt < {dt_floor['dt']/YEAR:g} yr for {dt_floor['steps']} consecutive accepted steps "
                                           f"before M_BH = {dt_floor['before_bh_mass']/MSUN:g} Msun (numerical blocker)");break
            dt=state.aux['controller']['dt_next'];tdyn=stepcontrol.dynamical_time(state)
            if dt<tdyn:
                status,reason='termination_trigger',f'timestep {dt/YEAR:.4g} yr below the dynamical time {tdyn/YEAR:.4g} yr (Ball\'s criterion)'
                break
            waiting=pending()
            target=stepcontrol.landing_target(state.bh_mass,growth_rate(state,eos,config,s),dt,waiting)
            attempts+=1
            try:
                if target is not None:
                    candidate,d=land(state,target,eos,opacity,config,s=s)
                else:
                    candidate,d=bondi_step(state,dt,eos,opacity,config,mode='evolve')
                    if waiting and candidate.bh_mass>=min(waiting):
                        # The end-of-step growth outran the accepted-state prediction: land instead.
                        target=min(waiting);candidate,d=land(state,target,eos,opacity,config,s=s)
                after=stepcontrol.observables(candidate,eos,mode)
                m=stepcontrol.metrics(before,after,d,S)
                if stepcontrol.rejected(m['delta_actual'],S):
                    raise rejection(f"delta_actual {m['delta_actual']:.4g} exceeds {S['reject_ratio']:g} x delta_target "
                                    f"({max((v,k) for k,v in m.items() if k!='delta_actual')[1]})")
                errors,diag=screened(candidate)
                if errors:
                    exc=rejection('Enforced screens: '+'; '.join(errors));exc.kind='screen';raise exc
            except StepRejected as exc:
                rejections+=1
                kind=getattr(exc,'kind','numerical');driver=state.aux['driver']
                if kind=='screen':driver['screen_retries']=driver.get('screen_retries',0)+1
                last_rejection=dict(step=state.steps,dt_years=dt/YEAR,landing_target_msun=None if target is None else target/MSUN,
                                    kind=kind,message=str(exc),domain_exit=getattr(exc,'domain_exit',False),
                                    domain_trials=getattr(exc,'domain_trials',0),domain_message=getattr(exc,'domain_message',None))
                log('rejected',**last_rejection)
                if stepcontrol.reject(state.aux,dt,S) is None:
                    retries,screens=state.aux['controller']['retries'],driver.get('screen_retries',0)
                    if screens==retries:
                        status,reason='model_domain_exit',(f"enforced screens violated at all {retries} consecutive retries "
                                                           f"(last: {last_rejection['message']})")
                    else:
                        status,reason='termination_trigger',(f"{S['max_retries']} consecutive retries, {screens} of {retries} "
                                                             f"by enforced screens (last: {last_rejection['message']})")
                    break
                continue
            # A landing's dt is forced by its target: the controller continues from the attempted dt.
            controlled=max(d['dt'],dt) if target is not None else d['dt']
            stepcontrol.accept(candidate.aux,controlled,m['delta_actual'],d['newton_iterations'],S)
            driver=candidate.aux['driver']
            slow=dt_floor and controlled<dt_floor['dt'] and candidate.bh_mass<dt_floor['before_bh_mass']
            driver['floor_run']=driver['floor_run']+1 if slow else 0
            driver['screen_retries']=0
            if target is not None:
                driver['landed']=driver['landed']+[target]
                if restart_landing is not None and target==restart_landing:
                    # Plan restart point 2: the last accepted step before this landing, as accepted.
                    label=f'before_landing_{target/MSUN:g}msun'
                    mark=state.copy();mark.aux=copy.deepcopy(accepted_aux)
                    save_checkpoint(path/f'restart_{label}.npz',mark,config)
                    if segment_steps:segment_ends[state.steps+segment_steps]=label
                    log('restart_point',label=label,step=state.steps)
            state,before=candidate,after;accepted+=1
            accepted_aux=copy.deepcopy(state.aux)
            record(history_row(state,eos,opacity,config,event='landing' if target is not None else 'step',diag=d,metrics=m,
                               screens=diag,target=target,s=s))
            log('accepted',step=state.steps,landing_target_msun=None if target is None else target/MSUN,
                metrics=m,**{k:v for k,v in d.items() if not isinstance(v,(dict,list))})
            if target is not None:
                label=f'{target/MSUN:g}msun'
                save_checkpoint(path/f'landing_{label}.npz',state,config)
                save_profile(path/f'profile_landing_{label}.npz',state,eos,opacity,config)
                log('landed',target_msun=target/MSUN,step=state.steps,time_years=state.time/YEAR,method=d.get('landing_method'),
                    structure=structure_summary(state,eos,opacity,config))
            for mass in restart_masses:
                label=f'bh{mass/MSUN:g}msun'
                if label not in saved and not (path/f'restart_{label}.npz').exists() and state.bh_mass>=mass:
                    save_checkpoint(path/f'restart_{label}.npz',state,config);saved.add(label)
                    if segment_steps:segment_ends[state.steps+segment_steps]=label
                    log('restart_point',label=label,step=state.steps)
            if state.steps in segment_ends:
                save_checkpoint(path/f'segment_end_{state.steps:07d}.npz',state,config)
                log('segment_end',label=segment_ends[state.steps],step=state.steps)
            if checkpoint_every and state.steps%checkpoint_every==0:save_checkpoint(path/'checkpoint.npz',state,config)
    except KeyboardInterrupt:
        status,reason='control_stop','execution interrupted; last accepted state preserved'
    except (RuntimeError,ValueError,ArithmeticError,DomainError) as exc:
        status,reason='numerical_failure',f'{type(exc).__name__}: {exc}'
    finally:
        save_checkpoint(path/'checkpoint.npz',state,config)
        try:save_profile(path/'profile_final.npz',state,eos,opacity,config)
        except (ValueError,RuntimeError,ArithmeticError) as exc:log('final_profile_unavailable',error=str(exc))
        from .ledger import residuals
        report=dict(status=status,reason=reason,accepted_steps_this_run=accepted,attempts=attempts,rejections=rejections,
                    final_step=state.steps,final_time_years=state.time/YEAR,final_bh_mass_msun=state.bh_mass/MSUN,
                    final_dynamical_time_years=stepcontrol.dynamical_time(state)/YEAR,
                    final_dt_years=state.aux['controller']['dt_next']/YEAR,landed_msun=[x/MSUN for x in state.aux['driver']['landed']],
                    last_rejection=last_rejection,runtime_seconds=time.perf_counter()-clock,restarted=restarted,
                    start_row_screen_errors=start_errors,
                    conservation=residuals(state,eos) if state.ledger else None,v1_validated=False)
        log('termination',status=status,reason=reason)
        io.json_write(path/'screen_violations.json',_json(dict(policy=policy if isinstance(policy,(str,type(None))) else 'custom',
                                                                screens=summary)))
        io.json_write(path/'termination.json',_json(report))
    return state,report

def opacity_tables(config,eos,opacity):
    """The opacity of opacity.tables (config.OPACITY_TABLES).

    default: `opacity` (the opacity of opacity.model; for mesa, MESA's
    default handle); stars2004_switch: the MESA ABI-7 STARS switch with the
    STARS nearest-table-value rule, reported by the opacity_clamp screen;
    stars2004_switch_strict: the same switch, off-table states are errors.
    """
    tables=config['opacity'].get('tables','default')
    if tables=='default':return opacity
    if getattr(eos,'abi_version',0)<7:raise ValueError(f'opacity.tables {tables} needs a MESA bridge ABI 7 build')
    return StarsOpacity(eos,clamp=tables=='stars2004_switch')

def physics(config):
    """EOS and opacity of a quasi-static configuration (evolution.physics, then opacity_tables)."""
    from .evolution import physics as evolution_physics
    eos,opacity=evolution_physics(config)
    return eos,opacity_tables(config,eos,opacity)

def run_track(config,output_parent='results',restart=None,stop_after_steps=None):
    """`qstars run` and `qstars restart` of a quasi_static configuration; returns (path, report).

    Without restart, the seeds of the replication section (seed_radii_rsun, in
    order) are relaxed (config['relaxation']) until one converges; t = 0 is
    that relaxed model (zero_time; row0.npz, profile_row0.npz, relax.json).
    The track (run) then uses the replication landings, end mass and screen
    policy, evolution.max_steps and max_wall_seconds and
    output.checkpoint_every_steps. Seed uniqueness and the fold classification
    of a termination are the registered driver's (scripts/validate_ball_track.py):
    a configuration that names a registered plan (replication.plan) is refused
    here before any output exists, because only that driver verifies the
    plan registration. So is inertia backward_euler, the identity-test and
    probe option of bondi_step: seeds carry no kinetic energy.
    """
    from . import io,ledger
    rep=config['replication']
    if rep['plan']:
        raise ValueError(f"{config['name']} is a registered run of {rep['plan']}; execute it with scripts/validate_ball_track.py, "
                         "which verifies the plan registration")
    if settings(config)['quasistatic']['inertia']!='none':
        raise ValueError('quasistatic.inertia backward_euler is the identity-test and probe option; a track run needs inertia none')
    config=copy.deepcopy(config)
    eos,opacity=physics(config)  # invalid microphysics fails before a run directory exists
    path=io.create_directory(output_parent,config['name'])
    io.provenance(path,config)
    if hasattr(eos,'manifest'):io.json_write(path/'microphysics_build.json',eos.manifest)
    if restart is not None:state=restart.copy()
    else:
        seeds=[];state=None
        for radius in rep['seed_radii_rsun']:
            try:
                candidate,_=seed(config,eos,opacity,config['numerics']['zones'],rep['initial_bh_mass_msun']*MSUN,
                                 rep['total_mass_msun']*MSUN,radius*RSUN)
            except (ValueError,RuntimeError,DomainError) as exc:
                seeds.append(dict(seed_radius_rsun=radius,status='seed_failed',reason=f'{type(exc).__name__}: {exc}'));continue
            ledger.initialize(candidate,eos,cavity=True)
            try:candidate,record=relax(candidate,eos,opacity,config,config['relaxation'])
            except (ValueError,RuntimeError,ArithmeticError,DomainError) as exc:  # not a step rejection: relax retries those
                seeds.append(dict(seed_radius_rsun=radius,status='failed',reason=f'{type(exc).__name__}: {exc}'));continue
            seeds.append(dict(seed_radius_rsun=radius,status=record['status'],reason=record['reason'],accepted_steps=len(record['steps']),
                              rejections=len(record['rejections']),mesh_cycles=len(record['cycles'])))
            if record['status']=='converged':
                state=zero_time(candidate,eos);break
        io.json_write(path/'relax.json',_json(dict(relaxation=config['relaxation'],seeds=seeds)))
        if state is None:
            report=dict(status='numerical_failure',reason='no seed relaxed to a row-0 model (relaxation blocker; see relax.json)',
                        accepted_steps_this_run=0,v1_validated=False,source_snapshot_verified=io.finish_provenance(path))
            io.json_write(path/'termination.json',report)
            return path,report
        save_checkpoint(path/'row0.npz',state,config);save_profile(path/'profile_row0.npz',state,eos,opacity,config)
    state,report=run(state,eos,opacity,config,path,landings=[m*MSUN for m in rep['landings_msun']],
                     end_bh_mass=rep['end_bh_mass_msun']*MSUN or None,policy=rep['policy'],stop_after_steps=stop_after_steps,
                     max_wall_seconds=config['evolution']['max_wall_seconds'] or None,
                     max_accepted_steps=config['evolution']['max_steps'],checkpoint_every=config['output']['checkpoint_every_steps'])
    report['source_snapshot_verified']=io.finish_provenance(path)
    io.json_write(path/'termination.json',_json(report))
    return path,report
