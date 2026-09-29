"""Stationary spherical LTE diffusion winds with a regular gas sonic point.

Nakauchi et al. (2017), equations 1--9, generalized to the derivatives and
ionization energy of a caller-supplied thermodynamic EOS. This component is
NOT an escaping-radiation closure or an enabled envelope mass-loss law.
The optional paper-photosphere endpoint deliberately extends diffusion beyond
its strict validity screen for benchmark comparisons; it never passes that
screen by construction. No MLT flux or mass-loss coefficient is introduced.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid, simpson
from .constants import G, C, KB, MU, A


@dataclass
class GasState:
    pg: float
    pr: float
    gas_u: float
    gas_s: float
    ct2: float
    b: float  # (dPg/dlnT)_rho / rho
    h: float
    s: float
    kappa: float
    ionized_h: float


class FreeEOSOpacity:
    """Direct FreeEOS gas derivatives plus matched MESA radiative opacity.

    This adapter does not invoke the historical MESA EOS table and its dilute
    gas-pressure cancellation. Constants and composition must match exactly.
    Bounds are an evaluation box, not evidence for LTE in an exposed wind.
    """
    def __init__(self, eos, mesa, opacity_family='default'):
        if eos.hydrogen != mesa.composition.X:
            raise ValueError('Wind EOS and opacity composition differ')
        if abs(eos.radiation_constant / mesa.radiation_constant - 1) > 1e-12:
            raise ValueError('Wind radiation constants differ')
        if opacity_family not in ('default','paper_gs98','paper_gn93'):
            raise ValueError('Unknown reference opacity family')
        self.eos, self.mesa = eos, mesa
        self.opacity_family=opacity_family
        self.arad = eos.radiation_constant

    def __call__(self, rho, temperature):
        if not (1e-22 <= rho <= 1e-4 and 3000 <= temperature <= 1e7):
            raise ValueError(f'Wind microphysics evaluation box exceeded: rho={rho}, T={temperature}')
        x = self.eos.raw(rho, temperature)
        pg = float(np.exp(x[0])); pr = self.arad * temperature**4 / 3
        kappa = float(self.mesa.opacity_external_eos(rho, temperature, x[15], x[16])
                      if self.opacity_family=='default' else
                      self.mesa.opacity_paper_reference(rho, temperature, self.opacity_family[6:]))
        return GasState(pg, pr, float(x[3]), float(x[6]), pg*x[1]/rho,
                        pg*x[2]/rho, float(x[3])+(pg+4*pr)/rho,
                        float(x[6])+4*pr/(rho*temperature), kappa, float(x[18]))


class IdealGasOpacity:
    """Fully ionized constant-opacity validation limit, never cool physics."""
    def __init__(self, mu=.62, kappa=.34, arad=A):
        if not all(np.isfinite(v) and v > 0 for v in (mu, kappa, arad)):
            raise ValueError('Invalid ideal reference constants')
        self.rgas = KB/(MU*mu); self.kappa=kappa; self.arad=arad

    def __call__(self, rho, temperature):
        pg=rho*self.rgas*temperature; pr=self.arad*temperature**4/3
        u=1.5*self.rgas*temperature
        s=self.rgas*(1.5*np.log(temperature)-np.log(rho))
        return GasState(pg,pr,u,s,self.rgas*temperature,self.rgas*temperature,
                        u+(pg+4*pr)/rho,s+4*pr/(rho*temperature),self.kappa,1.)


class WindBranchError(ValueError):
    """A local critical point is inadmissible; no global no-root claim."""


class WindMatchError(RuntimeError):
    def __init__(self,message,trials):
        super().__init__(message)
        self.trials=trials


def regular_sonic_slope(coefficients, temperature_slope, step=1e-4):
    """l'Hopital root for dlnv/dlnr=N/D with dlnT/dlnr supplied.

    coefficients([lnr,lnv,lnT]) returns nondimensional N,D near the critical
    point. Shared with the independent isothermal Parker verification.
    """
    h=step
    jac=np.column_stack([(coefficients(np.eye(3)[i]*h)-coefficients(-np.eye(3)[i]*h))/(2*h) for i in range(3)])
    n,d=jac;g=temperature_slope
    roots=np.roots([d[1],d[0]+d[2]*g-n[1],-n[0]-n[2]*g])
    physical=[float(s.real) for s in roots if abs(s.imag)<1e-8 and s.real>0
              and d[0]+d[1]*s.real+d[2]*g>0]
    if len(physical)!=1:
        raise WindBranchError(f'No unique accelerating subsonic-to-supersonic local branch; slopes={roots}')
    return physical[0], [dict(real=float(z.real),imag=float(z.imag)) for z in roots]


class DiffusionWind:
    def __init__(self, microphysics, mass, sonic_radius, sonic_density, sonic_temperature,
                 derivative_step=1e-4):
        if not all(np.isfinite(x) and x > 0 for x in
                   (mass,sonic_radius,sonic_density,sonic_temperature,derivative_step)):
            raise ValueError('Wind inputs must be positive and finite')
        self.micro=microphysics; self.mass=mass; self.rs=sonic_radius
        self.rhos=sonic_density; self.ts=sonic_temperature
        self.derivative_step=derivative_step
        q=microphysics(sonic_density,sonic_temperature)
        if q.ct2 <= 0 or q.kappa <= 0 or q.b <= 0:
            raise WindBranchError('Nonpositive sonic thermodynamic derivative or opacity')
        self.vs=np.sqrt(q.ct2)
        self.mdot=4*np.pi*self.rs**2*self.rhos*self.vs
        self.ls=4*np.pi*C*G*mass/q.kappa * (1-2*q.ct2*self.rs/(G*mass))/(1+q.b*self.rhos/(4*q.pr))
        if self.ls <= 0:
            raise WindBranchError('Specified sonic state requires nonpositive diffusion luminosity')
        self.energy=self.ls+self.mdot*(self.vs**2/2+q.h-G*mass/self.rs)
        self.sonic_q=q
        self.slope, self.slopes, self.sonic_temperature_slope = self._critical_slope()

    def state(self, x, z):
        r=self.rs*np.exp(x); v=self.vs*np.exp(z[0]); T=self.ts*np.exp(z[1])
        rho=self.rhos*np.exp(-2*x-z[0]); q=self.micro(rho,T)
        L=self.energy-self.mdot*(v*v/2+q.h-G*self.mass/r)
        g=-3*q.kappa*rho*L/(16*np.pi*self.micro.arad*C*r*T**4)
        gamma=q.kappa*L/(4*np.pi*C*G*self.mass)
        n=2*q.ct2-q.b*g+G*self.mass/r*(gamma-1)
        d=v*v-q.ct2
        return r,v,T,rho,q,L,g,n,d,gamma

    def _critical_slope(self):
        def coefficients(v):
            a=self.state(v[0],v[1:]); return np.array([a[7],a[8]])/self.vs**2
        g=self.state(0.,np.zeros(2))[6]
        slope,roots=regular_sonic_slope(coefficients,g,self.derivative_step)
        return slope,roots,g

    def boundary_flux(self,x,z,*,face_velocity=0.):
        """Lab-frame conserved fluxes across a possibly moving sphere.

        Reynolds transport subtracts face speed times the conserved density.
        Relative enthalpy advection needs an additional P*A*face_velocity work
        term. Omitting it would incorrectly remove pressure work at a comoving
        face. This interface alone does not implement a moving stellar update.
        """
        if not np.isfinite(face_velocity):raise ValueError('Face velocity must be finite')
        r,v,T,rho,q,L,*_=self.state(x,z)
        area=4*np.pi*r*r;mdot=area*rho*(v-face_velocity)
        moving_work=area*(q.pg+q.pr)*face_velocity
        total=L+mdot*(q.h+v*v/2-G*self.mass/r)+moving_work
        return dict(radius=r,density=rho,temperature=T,velocity=v,
                    face_velocity=face_velocity,mass=mdot,
                    momentum=mdot*v+area*(q.pg+q.pr),
                    radiative_energy=L,gas_internal_energy=mdot*q.gas_u,
                    gas_pressure_work=mdot*q.pg/rho,
                    radiation_internal_energy=mdot*3*q.pr/rho,
                    radiation_pressure_work=mdot*q.pr/rho,
                    moving_pressure_work=moving_work,
                    kinetic_energy=mdot*v*v/2,potential_energy=-mdot*G*self.mass/r,
                    entropy=mdot*q.s,total_energy=total)

    def rhs(self,x,z):
        r,v,T,rho,q,L,g,n,d,gamma=self.state(x,z)
        return np.array([n/d,g])

    def integrate(self, *, inner_ratio=.5, outer_ratio=100., rtol=1e-8,
                  sonic_offset=1e-6, max_log_step=.02, samples=1001,
                  endpoint='diffusion_limit', method='DOP853', stalled_mach=np.sqrt(1.01)):
        if not 0 < inner_ratio < 1 < outer_ratio or not 0 < rtol < 1e-3:
            raise ValueError('Invalid wind integration bounds or tolerance')
        if not 0 < sonic_offset < min(.01,-np.log(inner_ratio),np.log(outer_ratio)):
            raise ValueError('Invalid sonic offset')
        if endpoint not in ('diffusion_limit','paper_photosphere') or samples < 11 or max_log_step<=0:
            raise ValueError('Invalid endpoint or radial sampling')
        if not np.isfinite(stalled_mach) or stalled_mach<=1:
            raise ValueError('Stalled-branch endpoint must exceed gas Mach one')
        limit=.1 if endpoint=='diffusion_limit' else .75
        def radiation_limit(x,z):
            r,v,T,rho,q,L,*_=self.state(x,z)
            return 3*L/(4*np.pi*r*r*C*self.micro.arad*T**4)-limit
        def luminosity_zero(x,z):return self.state(x,z)[5]/self.ls
        def recross(x,z):
            a=self.state(x,z);return a[1]**2/a[4].ct2-stalled_mach**2
        def critical_floor(x,z):
            # A requested high-Mach comparison cutoff may never be reached.
            # It must not disable the independent second-sonic safety exit.
            a=self.state(x,z);return a[1]**2/a[4].ct2-1.01
        def too_cold(x,z):return self.ts*np.exp(z[1])-3100
        radiation_limit.direction=1;luminosity_zero.direction=-1
        recross.direction=-1;critical_floor.direction=-1;too_cold.direction=-1
        for f in (radiation_limit,luminosity_zero,recross,too_cold,critical_floor):f.terminal=True
        if radiation_limit(0,np.zeros(2))>=0:
            raise WindBranchError('Sonic state already outside selected radiation domain')
        parts=[]; terminations={}; raw=[]; retries=[]
        for direction,bound in ((-1,np.log(inner_ratio)),(1,np.log(outer_ratio))):
            x0=direction*sonic_offset
            z0=x0*np.array([self.slope,self.sonic_temperature_slope])
            events=(radiation_limit,luminosity_zero,recross,too_cold) if direction>0 else (luminosity_zero,)
            if direction>0 and stalled_mach!=np.sqrt(1.01):events=(*events,critical_floor)
            # An explicit RK trial can leave the EOS box although the accepted
            # solution stays inside. Retry from the identical sonic start with
            # a smaller step; never clamp a trial thermodynamic state.
            for attempt in range(4):
                try:
                    sol=solve_ivp(self.rhs,[x0,bound],z0,rtol=rtol,atol=rtol*.01,
                                  method=method,max_step=max_log_step/2**attempt,
                                  dense_output=True,events=events)
                    break
                except ValueError as exc:
                    retries.append(dict(direction=direction,attempt=attempt,
                                        max_log_step=max_log_step/2**attempt,reason=str(exc)))
                    if attempt==3:
                        error=RuntimeError('Integration retries exhausted; no physical no-root inference: '+str(exc))
                        error.integration_retries=retries
                        raise error from exc
            raw.append(sol)
            labels=[endpoint,'nonpositive_luminosity','approaching_second_sonic_point','microphysics_temperature_limit'] if direction>0 else ['nonpositive_luminosity']
            if direction>0 and stalled_mach!=np.sqrt(1.01):labels.append('approaching_second_sonic_point')
            reached=[labels[i] for i,v in enumerate(sol.t_events) if len(v)]
            terminations['outward' if direction>0 else 'inward']=reached[0] if reached else ('radial_bound' if sol.success else 'numerical_failure')
            # Retain the integrator's knots and subdivide them for independent
            # quadrature. A uniform output grid misses the sharp recombination
            # layer even when the ODE solution itself is converged.
            knots=np.concatenate([np.linspace(a,b,9)[:-1] for a,b in zip(sol.t[:-1],sol.t[1:])]) if len(sol.t)>1 else sol.t
            x=np.unique(np.r_[np.linspace(x0,sol.t[-1],samples),knots,sol.t[-1]])
            if direction<0:x=x[::-1]
            parts.append((x,sol.sol(x)))
        x=np.r_[parts[0][0][::-1],0.,parts[1][0]]
        z=np.c_[parts[0][1][:,::-1],[0.,0.],parts[1][1]]
        fields=self.profile(x,z)
        summary=self.summary(fields,terminations)
        summary.update(rtol=rtol,sonic_offset=sonic_offset,derivative_step=self.derivative_step,
                       max_log_step=max_log_step,samples=samples,endpoint=endpoint,method=method,
                       stalled_mach=stalled_mach,
                       solver_success=bool(all(s.success for s in raw)),
                       solver_messages=[s.message for s in raw],nfev=sum(s.nfev for s in raw),
                       integration_retries=retries)
        return fields,summary

    def profile(self,x,z):
        rows=[]
        for xi,zi in zip(x,z.T):
            r,v,T,rho,q,L,g,n,d,gamma=self.state(xi,zi)
            slope=self.slope if xi==0 else n/d
            rows.append([r,v,T,rho,q.pg,q.pr,q.kappa,L,self.mdot*q.gas_u,
                         self.mdot*q.pg/rho,self.mdot*4*q.pr/rho,self.mdot*v*v/2,
                         -self.mdot*G*self.mass/r,q.s,q.ionized_h,np.sqrt(q.ct2),
                         gamma,g,slope,q.h])
        names=('radius','velocity','temperature','density','gas_pressure','radiation_pressure',
               'opacity','radiative_luminosity','gas_internal_luminosity','gas_work_luminosity',
               'radiation_enthalpy_luminosity','kinetic_luminosity','potential_luminosity',
               'entropy','hydrogen_ionized','gas_sound_speed','eddington_factor',
               'temperature_log_slope','velocity_log_slope','enthalpy')
        f=dict(zip(names,np.array(rows).T));r=f['radius']
        f['log_radius']=x
        f['mass_luminosity']=4*np.pi*r*r*f['density']*f['velocity']
        f['gradient_knudsen']=3*f['radiative_luminosity']/(4*np.pi*r*r*C*self.micro.arad*f['temperature']**4)
        # Optical depth to the finite endpoint. No unspecified exterior added.
        tau=cumulative_trapezoid(f['opacity']*f['density'],r,initial=0)
        f['optical_depth_to_endpoint']=tau[-1]-tau
        f['shell_mass_from_inner']=cumulative_trapezoid(4*np.pi*r*r*f['density'],r,initial=0)
        f['mach']=f['velocity']/f['gas_sound_speed']
        f['escape_speed']=np.sqrt(2*G*self.mass/r)
        return f

    def summary(self,f,terminations):
        x=f['log_radius'];r=f['radius'];v=f['velocity'];rho=f['density'];T=f['temperature']
        # Independent differential/integral checks of momentum, diffusion and
        # first law; energy flux algebra alone is NOT a conservation test.
        a=f['velocity_log_slope'];b=f['temperature_log_slope']
        drho=rho*(-2-a)
        dp=np.empty(len(x))
        for i in range(len(x)):
            q=self.micro(rho[i],T[i]);dp[i]=q.ct2*drho[i]+rho[i]*q.b*b[i]+4*q.pr*b[i]
        momentum=v[-1]**2/2-v[0]**2/2+simpson(dp/rho+G*self.mass/r,x=x)
        diff=f['radiation_pressure'][-1]-f['radiation_pressure'][0]+simpson(f['opacity']*rho*f['radiative_luminosity']/(4*np.pi*C*r),x=x)
        # Integration by parts avoids numerical entropy derivatives.
        entropy_integral=T[-1]*f['entropy'][-1]-T[0]*f['entropy'][0]-simpson(f['entropy']*T*b,x=x)
        energy_defect=f['radiative_luminosity'][-1]-f['radiative_luminosity'][0]+self.mdot*entropy_integral
        channels=('radiative_luminosity','gas_internal_luminosity','gas_work_luminosity',
                  'radiation_enthalpy_luminosity','kinetic_luminosity','potential_luminosity')
        total=sum(f[k] for k in channels);norm=max(abs(total).max(),self.ls)
        return dict(family='nakauchi2017_lte_diffusion',mass=self.mass,sonic_radius=self.rs,
                    sonic_density=self.rhos,sonic_temperature=self.ts,sonic_velocity=self.vs,
                    sonic_luminosity=self.ls,mass_flux=self.mdot,total_energy_flux=self.energy,
                    sonic_slope=self.slope,critical_slopes=self.slopes,
                    sonic_regularity_residual=abs(self.state(0,np.zeros(2))[7])/self.vs**2,
                    mass_flux_relative_defect=float(np.max(abs(f['mass_luminosity']/self.mdot-1))),
                    energy_algebra_relative_defect=float(np.max(abs(total-self.energy))/norm),
                    momentum_integral_relative_defect=float(abs(momentum)/(G*self.mass/r[0])),
                    diffusion_integral_relative_defect=float(abs(diff)/max(f['radiation_pressure'])),
                    first_law_integral_relative_defect=float(abs(energy_defect)/norm),
                    shell_mass_fraction=float(f['shell_mass_from_inner'][-1]/self.mass),
                    point_mass_screen_passed=bool(f['shell_mass_from_inner'][-1]/self.mass<.01),
                    diffusion_screen_passed=bool(max(f['gradient_knudsen'])<=.1*(1+1e-8)),
                    spherical_gravity_change=float(1-(r[0]/r[-1])**2),
                    max_knudsen=float(max(f['gradient_knudsen'])),
                    max_velocity_over_escape=float(max(v/f['escape_speed'])),
                    max_velocity_over_c=float(max(v)/C),
                    nonrelativistic_screen_passed=bool(max(v)/C<.01),
                    endpoint_velocity_over_escape=float(v[-1]/f['escape_speed'][-1]),
                    endpoint_temperature=float(T[-1]),endpoint_radius=float(r[-1]),
                    terminations=terminations,escaping_wind_validated=False,
                    convection_assumption='Zero convective flux in wind; requires separate matching validity',
                    physical_impossibility_claimed=False,v1_validated=False)


def match_base(microphysics, *, mass, radius, density, temperature, luminosity,
               sonic_guess, rtol=1e-7, match_tolerance=1e-5, max_nfev=200,
               initial_trust_radius=.2, derivative_step=1e-4,
               sonic_offset=1e-7, max_log_step=.01):
    """Shoot sonic r,rho,T to supplied base rho,T,total outward energy flux.

    No mdot input: it follows from sonic regularity and the three base data.
    This is a component matching operation, not authorization to remove mass
    from an evolving envelope. A successful local match cannot certify escape.
    The three luminosity channels on the right of Nakauchi eq9 include gas/
    radiation enthalpy and kinetic energy; gravitational flux is explicit.
    """
    if not all(np.isfinite(v) and v>0 for v in (mass,radius,density,temperature,luminosity)):
        raise ValueError('Base data must be finite and positive')
    if (not 0 < rtol < 1e-3 or not 0 < match_tolerance < .01 or max_nfev < 1
            or not 0 < initial_trust_radius <= 2 or not 0 < sonic_offset < .001
            or not 0 < derivative_step < .01 or not 0 < max_log_step < 1):
        raise ValueError('Invalid matching controls')
    trace=[]
    scale=np.array([radius,density,temperature])
    def decode(y):
        return scale*np.array([1+np.exp(y[0]),np.exp(y[1]),np.exp(y[2])])
    def evaluate(y):
        rs,ds,ts=decode(y)
        wind=DiffusionWind(microphysics,mass,rs,ds,ts,derivative_step)
        eps=min(sonic_offset,.01*np.log(rs/radius))
        sol=solve_ivp(wind.rhs,[-eps,np.log(radius/rs)],-eps*np.array([wind.slope,wind.sonic_temperature_slope]),
                      method='DOP853',rtol=rtol,atol=rtol*.01,max_step=max_log_step)
        if not sol.success:
            raise RuntimeError('Base matching integration failed; no physical no-root inference: '+sol.message)
        r,v,T,rho,q,L,*_=wind.state(sol.t[-1],sol.y[:,-1])
        total=L+wind.mdot*(v*v/2+q.h)
        if total<=0:raise WindBranchError('Negative outward total base energy flux')
        value=[np.log(rho/density),np.log(T/temperature),np.log(total/luminosity)]
        return np.asarray(value)
    lower=np.array([np.log(1e-8),-30.,-10.]);upper=np.array([np.log(99),10.,10.])
    def residual(y, purpose):
        if len(trace)>=max_nfev:
            raise WindMatchError('Matching evaluation budget exhausted; no physical no-root inference',trace)
        row=dict(log_coordinates=y.tolist(),purpose=purpose)
        try:
            if np.any(y<lower) or np.any(y>upper):
                raise WindBranchError('Trial outside declared sonic search box')
            value=evaluate(y)
            if not np.all(np.isfinite(value)):
                raise WindBranchError('Nonfinite matching residual')
            row.update(status='admissible',residual=value.tolist(),
                       sonic_radius=float(decode(y)[0]),sonic_density=float(decode(y)[1]),
                       sonic_temperature=float(decode(y)[2]))
        except (ValueError,RuntimeError,FloatingPointError) as exc:
            row.update(status='rejected_local_trial',reason=str(exc));value=None
        trace.append(row)
        return value
    guess=np.asarray(sonic_guess)/scale
    if guess.shape!=(3,) or np.any(guess<=0) or guess[0]<=1:
        raise ValueError('Sonic guess must lie above base with positive rho,T')
    guess[0]-=1
    y=np.log(guess);value=residual(y,'initial');iterations=[];trust=initial_trust_radius
    if value is None:
        raise WindMatchError('Initial sonic branch is inadmissible; continuation needs an admissible anchor',trace)
    while max(abs(value))>match_tolerance:
        columns=[]
        for k in range(3):
            h=1e-4
            for attempt in range(10):
                delta=np.eye(3)[k]*h
                plus=residual(y+delta,'jacobian_plus');minus=residual(y-delta,'jacobian_minus')
                if plus is not None and minus is not None:
                    columns.append((plus-minus)/(2*h));break
                h*=.5
            else:
                raise WindMatchError('No admissible centered Jacobian stencil; local numerical/domain failure',trace)
        jac=np.column_stack(columns)
        step,_,rank,singular=np.linalg.lstsq(jac,-value,rcond=1e-12)
        if rank<3:
            raise WindMatchError('Matching Jacobian lost rank; branch continuation required',trace)
        step*=min(1.,trust/max(abs(step)))
        old_norm=float(np.linalg.norm(value))
        for backtrack in range(16):
            fraction=2.**-backtrack
            trial=residual(y+fraction*step,'line_search')
            if trial is not None and np.linalg.norm(trial)<old_norm*(1-1e-4*fraction):
                y+=fraction*step;value=trial
                iterations.append(dict(residual=value.tolist(),fraction=fraction,
                                       jacobian_condition=float(singular[0]/singular[-1]),trust=trust))
                trust=min(initial_trust_radius*4,trust*(1.5 if backtrack==0 else .5))
                break
        else:
            raise WindMatchError('Admissible line search stalled; no physical nonexistence claim',trace)
    defect=float(max(abs(value)))
    wind=DiffusionWind(microphysics,mass,*decode(y),derivative_step)
    return wind,dict(max_log_residual=defect,nfev=len(trace),trials=trace,iterations=iterations,
                     rejected_trials=sum(t['status']=='rejected_local_trial' for t in trace),
                     method='admissible_damped_newton',mass_flux_is_eigenvalue=True,
                     escape_not_certified=True,all_branches_enumerated=False)


def spherical_hydrostatic(microphysics, *, mass, radius, density, temperature,
                         luminosity, outer_ratio=2., self_gravity=False, rtol=1e-9,
                         samples=501):
    """Spherical radiative hydrostatic component, without an MLT assumption.

    Convective stability is not implied. Self-gravity is included when chosen;
    the point-mass option is retained for the analytical radiative-zero test.
    No photospheric boundary or zero external-pressure condition is imposed.
    """
    if not all(np.isfinite(v) and v>0 for v in (mass,radius,density,temperature,luminosity,outer_ratio,rtol)) or outer_ratio<=1:
        raise ValueError('Invalid hydrostatic reference inputs')
    def rhs(x,z):
        r=radius*np.exp(x);rho=density*np.exp(z[0]);T=temperature*np.exp(z[1])
        m=mass*np.exp(z[2]);q=microphysics(rho,T)
        g=-3*q.kappa*rho*luminosity/(16*np.pi*microphysics.arad*C*r*T**4)
        density_slope=(-G*m/r-(q.b+4*q.pr/rho)*g)/q.ct2
        return [density_slope,g,4*np.pi*r**3*rho/m if self_gravity else 0.]
    grid=np.linspace(0,np.log(outer_ratio),samples)
    sol=solve_ivp(rhs,[0,grid[-1]],[0.,0.,0.],t_eval=grid,method='DOP853',rtol=rtol,atol=rtol*.01,max_step=.02)
    if not sol.success:raise RuntimeError('Spherical hydrostatic integration failed: '+sol.message)
    return dict(radius=radius*np.exp(sol.t),density=density*np.exp(sol.y[0]),
                temperature=temperature*np.exp(sol.y[1]),enclosed_mass=mass*np.exp(sol.y[2]),
                luminosity=np.full(len(sol.t),luminosity))


def run_spherical_wind_reference(config, output='results'):
    """Reproducible CLI component run; a completed integration is not escape."""
    import os
    from . import io
    from .constants import MSUN
    allowed={'family','microphysics','mass_msun','sonic_radius_cm','sonic_density',
             'sonic_temperature','hydrogen','inner_ratio','outer_ratio','rtol',
             'sonic_offset','derivative_step','max_log_step','samples','endpoint','method'}
    if set(config)-allowed:raise ValueError('Unknown spherical wind configuration key')
    path=io.create_directory(output,'spherical_diffusion_wind');io.provenance(path,config)
    try:
        choice=config.get('microphysics','freeeos_mesa')
        if choice=='freeeos_mesa':
            from .freeeos_reference import FreeEOSReference
            from .mesa_microphysics import MesaMicrophysics
            from .eos import Composition
            X=config.get('hydrogen',.7)
            mesa=MesaMicrophysics(Composition(X,1-X,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
            eos=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],X,mesa.radiation_constant)
            micro=FreeEOSOpacity(eos,mesa)
            io.json_write(path/'microphysics.json',dict(freeeos=eos.manifest,mesa=mesa.manifest,
                          hydrogen=X,radiation_constant=micro.arad))
        elif choice=='ideal_constant_opacity_reference':micro=IdealGasOpacity()
        else:raise ValueError('Unrecognized wind microphysics family')
        wind=DiffusionWind(micro,config['mass_msun']*MSUN,config['sonic_radius_cm'],
                            config['sonic_density'],config['sonic_temperature'],
                            config.get('derivative_step',1e-4))
        keys={'inner_ratio','outer_ratio','rtol','sonic_offset','max_log_step','samples','endpoint','method'}
        fields,report=wind.integrate(**{k:v for k,v in config.items() if k in keys})
        np.savez_compressed(path/'wind.npz',**fields)
        report['status']='completed_reference' if report['solver_success'] else 'numerical_failure'
    except (ValueError,RuntimeError,KeyError,OSError) as exc:
        report=dict(status='component_failure',reason=str(exc),physical_impossibility_claimed=False,v1_validated=False)
    report['source_snapshot_verified']=io.finish_provenance(path)
    io.json_write(path/'termination.json',report)
    return path,report
