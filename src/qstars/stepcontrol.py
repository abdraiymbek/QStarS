"""Step control for quasi-static Ball tracks (spec 5; plan numerics.controller).

delta_actual is the largest of: the profile change (metric 'lagrangian', the
registered plan: Lagrangian |dln T| per cell and |dln r| per face, excluding
gas within 2 x (the step's swept mass) of the base); |dln r0|, |dln rho_b|,
|dln T_b|; |dln M_BH|; |dM_cav|/M0; |dln L*|, |dln R*|, |dln Teff|. The base
condensation is controlled through these boundary observables, not through
the gas swept across r0 (correction C1).

metric 'mapped' (a proposed amendment, not the registered default) measures
the profile change against the step predictor instead: the accepted state
mapped in mass above the base in the inner half and in mass below the
surface in the outer half. With the Lagrangian metric, gas just above the
excluded band moves by about swept/M_cond per step near a base
condensation of mass M_cond, so whenever the swept mass is below M_cond the
controller settles at dt ~ delta_target M_cond/Mdot and never reaches the
regime the 2 x swept exclusion assumes. On the analytic 1e4 Msun, M_BH = 5
toy (128 zones, face state) it held profile_dlnr at 0.00455 for 58 steps
with dt falling from 0.044 to 0.026 yr; with the mapped metric and the
fixed-mass base (quasistatic base='fixed_mass', fallback F2) the same toy
runs at 40-58 yr steps without a rejection. Both need a registered
amendment (plan registration.amendment_rule) before any Ball run uses them.

dt_next = dt clip(0.9 delta_target/delta_actual, 0.5, 1.5), no growth for 3
steps after a rejection, x0.8 after more than 8 Newton iterations, at most
dt_max. A step is rejected when delta_actual > 2 delta_target, the Newton
solve fails or more than 3 predictor passes are needed; dt is then halved,
at most 12 consecutive times. The controller state lives in
State.aux['controller'] (plain JSON), so a restart continues bitwise.
A nonfinite metric or timestep is an error, never a silent pass.
"""
import math
import numpy as np
from .constants import G,YEAR

DEFAULTS=dict(delta_target=.005,safety=.9,min_factor=.5,max_factor=1.5,reject_ratio=2.,hold_steps=3,
              slow_iterations=8,slow_factor=.8,retry_factor=.5,max_retries=12,
              dt_initial_years=1.,dt_max_years=2e4,metric='lagrangian')
METRICS=('lagrangian','mapped')
OBSERVABLES=('inner_radius','base_density','base_temperature','bh_mass','luminosity','radius','effective_temperature')

def settings(config):
    given=config.get('stepcontrol') or {}
    unknown=set(given)-set(DEFAULTS)
    if unknown:raise ValueError(f'Unknown stepcontrol settings {sorted(unknown)}')
    out={**DEFAULTS,**given}
    if out['metric'] not in METRICS:raise ValueError(f'stepcontrol.metric must be one of {METRICS}')
    return out

def observables(state,eos,mode='face'):
    """Boundary and surface observables of an accepted state (face state as in cavity.state_boundary)."""
    from .photosphere import effective_temperature
    from .cavity import state_boundary
    B=state_boundary(state,eos,mode)
    L=float(state.face_luminosity[-1]) if state.face_luminosity is not None else np.nan
    return dict(inner_radius=float(state.r[0]),base_density=float(B['rho']),base_temperature=float(B['T']),
                bh_mass=float(state.bh_mass),cavity_mass=float(state.cavity_mass),
                central_mass=float(state.enclosed[0]),luminosity=L,radius=float(state.r[-1]),
                effective_temperature=float(effective_temperature(L,state.r[-1])))

def metrics(before,after,diag,s=None):
    """Components of delta_actual for one step; `before`/`after` from observables()."""
    metric=(s or DEFAULTS)['metric']
    out=dict(profile_dlnT=diag[metric+'_max_dlnT'],profile_dlnr=diag[metric+'_max_dlnr'])
    for key in OBSERVABLES:
        if key in before and key in after:out['dln_'+key]=abs(math.log(after[key]/before[key]))
    out['cavity_change']=abs(after['cavity_mass']-before['cavity_mass'])/after['central_mass']
    bad=sorted(k for k,v in out.items() if not math.isfinite(v))
    if bad:raise ValueError(f'Nonfinite step-control metrics {bad}')
    out['delta_actual']=max(out.values())
    return out

def start(aux,s,dt=None):
    """Initialise the controller state; returns the first timestep (s)."""
    dt=s['dt_initial_years']*YEAR if dt is None else dt
    aux['controller']=dict(dt_next=float(dt),hold=0,retries=0)
    return dt

def accept(aux,dt,delta_actual,iterations,s):
    """Record an accepted step; returns and stores the next timestep."""
    c=aux['controller']
    if not (math.isfinite(dt) and math.isfinite(delta_actual)):raise ValueError('Accepted step needs a finite dt and delta_actual')
    factor=min(max(s['safety']*s['delta_target']/max(delta_actual,1e-300),s['min_factor']),s['max_factor'])
    if c['hold']>0:factor=min(factor,1.)
    if iterations>s['slow_iterations']:factor*=s['slow_factor']
    dt_next=float(min(dt*factor,s['dt_max_years']*YEAR))
    if not (math.isfinite(dt_next) and dt_next>0):raise ValueError(f'Controller produced a nonfinite or nonpositive timestep {dt_next}')
    c.update(dt_next=dt_next,hold=max(c['hold']-1,0),retries=0)
    return c['dt_next']

def rejected(delta_actual,s):
    return delta_actual>s['reject_ratio']*s['delta_target']

def reject(aux,dt,s):
    """Record a rejected attempt; returns the retry timestep or None when retries are exhausted."""
    c=aux['controller']
    c.update(retries=c['retries']+1,hold=s['hold_steps'])
    if c['retries']>s['max_retries']:return None
    c['dt_next']=float(dt*s['retry_factor'])
    return c['dt_next']

def dynamical_time(state):
    """Ball's operational termination scale t_dyn = sqrt(R*^3/(G M*))."""
    return float(np.sqrt(state.r[-1]**3/(G*state.enclosed[-1])))

def landing_target(bh_mass,growth_rate,dt,pending):
    """The lowest pending landing mass if a step of dt would reach or pass it, else None (g).

    pending: the landing masses not landed yet (the driver keeps the landed
    ones in State.aux['driver']['landed']); a target is pending exactly when
    it has not been landed, with no tolerance on M_BH.
    """
    if not pending:return None
    target=min(pending)
    return target if bh_mass+growth_rate*dt>=target else None
