"""Fiacconi & Rossi (2016), equations 12–22, constant-opacity wind.

This is the published gas-pressure-free alternative, not a gas sonic-point
eigenvalue solver. The base Mach number is an EXPLICIT boundary prescription
(the paper adopts 0.8 +/- 0.05). It does not yet attach an atmosphere or wind
to QStarS's time-dependent envelope. No evolutionary mass loss is inferred.
"""
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.optimize import brentq


def integrate_wind(alpha,energy,beta,w_infinity,*,delta=1e-6,rtol=1e-8,method='Radau'):
    """Integrate inward using t=ln(r/R), y=u exp(2t), and positive log states."""
    if not all(np.isfinite(x) and x>0 for x in (alpha,energy,beta,w_infinity,delta,rtol)):
        raise ValueError('Wind parameters must be positive and finite')
    if not 1e-9<=delta<=1e-3 or not 1e-12<=rtol<=1e-4:
        raise ValueError('Wind boundary delta or integration tolerance outside supported range')
    gamma_inf=(energy-alpha*w_infinity)/(1+2*np.sqrt(w_infinity)/beta)
    w_start=w_infinity-(gamma_inf-1)*delta
    if gamma_inf<=0 or w_start<=0:raise ValueError('Nonpositive asymptotic wind boundary')
    if alpha*beta*delta/np.sqrt(w_start)>.01:
        raise ValueError('Asymptotic boundary is not optically thin; decrease delta')
    def rhs(t,z):
        w,y=np.exp(z)
        tau=alpha*beta*np.exp(-t)/np.sqrt(w)
        g=(1+1.5*tau)/(1+4.5*tau)
        gamma=energy-alpha*(w-np.exp(-t))-(1+g)*y*np.sqrt(w)/beta
        return [np.exp(-t)*(gamma-1)/w,2-2*gamma*(1+1.5*tau)/y]
    def stalled(t,z):return z[0]-np.log(1e-10)
    stalled.terminal=True
    result=solve_ivp(rhs,[np.log(1/delta),0.],np.log([w_start,gamma_inf]),
                     method=method,events=stalled,rtol=rtol,atol=rtol*.01,dense_output=True)
    result.w_infinity=w_infinity
    return result


def wind_fields(t,z,alpha,energy,beta):
    w,y=np.exp(z);radius=np.exp(t)
    tau=alpha*beta/(radius*np.sqrt(w));g=(1+1.5*tau)/(1+4.5*tau)
    advected=(1+g)*y*np.sqrt(w)/beta
    gamma=energy-alpha*(w-1/radius)-advected
    sound_squared=(2/3)*g*y*np.sqrt(w)/(alpha*beta)
    lte_sound_squared=(2/9)*y*np.sqrt(w)/(alpha*beta)
    return dict(radius_ratio=radius,w=w,u=y/radius**2,optical_depth_proxy=tau,
                luminosity=gamma,advected_luminosity=advected,kinetic_luminosity=alpha*w,
                gravitational_luminosity=-alpha/radius,
                sound_speed_ratio=np.sqrt(sound_squared),mach=np.sqrt(w/sound_squared),
                lte_sound_speed_proxy=np.sqrt(lte_sound_squared),lte_mach_proxy=np.sqrt(w/lte_sound_squared),
                velocity_ratio=np.sqrt(w),scaled_radiation_energy=y)


def solve_wind(alpha,energy,beta,base_mach=.8,*,delta=1e-6,rtol=1e-8,method='Radau',samples=1601):
    """Select w_infinity using the stated base Mach closure, with no fitted coefficient."""
    if not 0<base_mach<1:raise ValueError('Published base must be subsonic')
    if beta<10 or energy<=1:raise ValueError('Reference family requires beta>=10 and energy>1')
    if not isinstance(samples,int) or samples<100:raise ValueError('At least 100 profile samples required')
    trials=[]
    def mismatch(w):
        result=integrate_wind(alpha,energy,beta,w,delta=delta,rtol=rtol,method=method)
        if len(result.t_events[0]):
            trials.append(dict(w_infinity=w,status='stalled_before_base'))
            return -base_mach
        if not result.success:raise RuntimeError('Wind integration failed: '+result.message)
        fields=wind_fields(0.,result.y[:,-1],alpha,energy,beta)
        trials.append(dict(w_infinity=w,status='reached_base',base_mach=float(fields['mach'])))
        return float(fields['mach']-base_mach)
    lo=max(1e-4,2*((energy-1)*delta+(alpha*beta*delta/.01)**2))
    hi=energy/alpha*(1-1e-8)
    if lo>=hi:raise ValueError('No positive asymptotic shooting bracket')
    w=brentq(mismatch,lo,hi,xtol=max(1e-12,rtol*.01),rtol=rtol*.01)
    sol=integrate_wind(alpha,energy,beta,w,delta=delta,rtol=rtol,method=method)
    if not sol.success or sol.t[-1]!=0:raise RuntimeError('Selected wind did not reach the base')
    t=np.linspace(0,np.log(1/delta),samples)
    fields=wind_fields(t,sol.sol(t),alpha,energy,beta)
    if np.any(fields['luminosity']<=0):raise ValueError('Selected wind has negative radiative luminosity')
    if max(fields['velocity_ratio'])/beta>.1:
        raise ValueError('Wind exceeds the supported nonrelativistic velocity limit of 0.1c')
    if fields['optical_depth_proxy'][0]<=1 or fields['w'][0]>=1:
        raise ValueError('Selected wind violates published optically thick, bound-base conditions')
    def root_of(key,target):
        v=fields[key]-target;cross=np.flatnonzero(v[:-1]*v[1:]<0)
        if len(cross)!=1:raise ValueError(f'Expected one {key} crossing, found {len(cross)}')
        i=cross[0]
        return float(np.exp(brentq(lambda x:wind_fields(x,sol.sol(x),alpha,energy,beta)[key]-target,t[i],t[i+1])))
    # L=4 pi R_phot^2 sigma T^4 gives y=4 Gamma. This is the paper's
    # photosphere definition, not a solved frequency-dependent atmosphere.
    photo=fields['scaled_radiation_energy']-4*fields['luminosity']
    crossings=np.flatnonzero(photo[:-1]*photo[1:]<0)
    if len(crossings)!=1:raise ValueError('Published photosphere condition has no unique crossing')
    i=crossings[0]
    def photo_residual(x):
        f=wind_fields(x,sol.sol(x),alpha,energy,beta)
        return f['scaled_radiation_energy']-4*f['luminosity']
    photo_t=brentq(photo_residual,t[i],t[i+1]);photo_fields=wind_fields(photo_t,sol.sol(photo_t),alpha,energy,beta)
    sonic=root_of('mach',1.)
    budget=sum(fields[k] for k in ('luminosity','advected_luminosity','kinetic_luminosity','gravitational_luminosity'))
    # Independent integral residuals of momentum and radiation transport.
    momentum=fields['w'][-1]-fields['w'][0]-simpson(np.exp(-t)*(fields['luminosity']-1),x=t)
    radiation=fields['u'][-1]-fields['u'][0]+simpson(2*np.exp(-2*t)*fields['luminosity']*(1+1.5*fields['optical_depth_proxy']),x=t)
    summary=dict(alpha=alpha,energy=energy,beta=beta,base_mach=base_mach,w_infinity=w,
          base_w=float(fields['w'][0]),base_luminosity=float(fields['luminosity'][0]),
          base_optical_depth_proxy=float(fields['optical_depth_proxy'][0]),
          photospheric_radius_ratio=float(np.exp(photo_t)),sonic_radius_ratio=sonic,
          photospheric_optical_depth_proxy=float(photo_fields['optical_depth_proxy']),
          energy_algebraic_closure_relative=float(max(abs(budget-energy))/energy),
          momentum_integral_residual=float(momentum),radiation_integral_residual_relative=float(radiation/fields['u'][0]),
          boundary_delta=delta,rtol=rtol,method=method,
          initial_acceleration_sign=float(np.sign(fields['luminosity'][0]-1)),
          max_velocity_over_c=float(max(fields['velocity_ratio'])/beta),
          scope='Stationary radiation-dominated constant-opacity wind only; no coupled quasi-star evolution',
          trials=trials,v1_validated=False)
    return fields,summary


def run_wind_reference(config,output='results'):
    from . import io
    allowed={'family','alpha','energy','beta','base_mach','delta','rtol','method','samples'}
    if set(config)-allowed:raise ValueError('Unknown wind reference configuration key')
    params={k:v for k,v in config.items() if k!='family'}
    path=io.create_directory(output,'fiacconi2016_wind');io.provenance(path,config)
    try:fields,report=solve_wind(**params)
    except (ValueError,RuntimeError) as exc:
        report=dict(status='reference_failure',reason=str(exc),v1_validated=False)
        report['source_snapshot_verified']=io.finish_provenance(path)
        io.json_write(path/'termination.json',report)
        return path,report
    np.savez_compressed(path/'wind_profile.npz',**fields)
    report['status']='completed_reference'
    report['source_snapshot_verified']=io.finish_provenance(path)
    io.json_write(path/'termination.json',report)
    return path,report
