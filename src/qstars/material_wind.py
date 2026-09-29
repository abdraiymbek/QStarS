"""Local continuum wind equations with an independent material energy equation.

Sources are specified local gas heating per volume and radiative acceleration.
This operator does not compute these sources, transfer, populations, or a
global wind. Its critical condition assumes sources have no explicit dv/dr
dependence; do not substitute it into Sobolev line driving or LTE diffusion.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from .constants import G


@dataclass(frozen=True)
class GasDerivatives:
    pressure: float
    internal_energy: float
    entropy: float
    pressure_dlnrho: float
    pressure_dlnT: float
    energy_dlnrho: float
    energy_dlnT: float

    @classmethod
    def freeeos(cls, raw):
        """Gas-only native FreeEOS values; never add blackbody radiation here."""
        a=np.asarray(raw,float)
        if a.shape != (24,) or not np.all(np.isfinite(a)):
            raise ValueError('One finite 24-field native FreeEOS state required')
        p=float(np.exp(a[0]))
        return cls(p,float(a[3]),float(a[6]),p*a[1],p*a[2],float(a[4]),float(a[5]))


def continuum_coefficients(gas, *, radius, density, velocity, mass,
                           gas_heating_per_volume, radiative_acceleration):
    """Coefficients for y=dlnv/dlnr, t=dlnT/dlnr in a steady flow.

    Material first law: (u_rho-P/rho)*(-2-y)+u_T*t = r Q/(rho v).
    Eliminating t gives (v²-a_ad²)*y = 2*a_ad²-r GM/r²+r*g_rad
                                         - b*r*Q/(rho*v*u_T).
    Q excludes radiative mechanical work and includes only thermal exchange.
    This distinction is necessary for total gas+radiation energy conservation.
    """
    values=[radius,density,velocity,mass,gas_heating_per_volume,radiative_acceleration,
            *gas.__dict__.values()]
    if not np.all(np.isfinite(values)) or min(radius,density,velocity,mass,gas.energy_dlnT,gas.pressure)<=0:
        raise ValueError('Finite state, positive material flow and heat capacity required')
    ct2=gas.pressure_dlnrho/density
    b=gas.pressure_dlnT/density
    U=gas.energy_dlnrho-gas.pressure/density
    cs2=ct2-b*U/gas.energy_dlnT
    if min(ct2,cs2)<=0:
        raise ValueError('Nonpositive gas compressibility or adiabatic sound speed')
    heat=radius*gas_heating_per_volume/(density*velocity)
    numerator=2*cs2-G*mass/radius+radius*radiative_acceleration-b*heat/gas.energy_dlnT
    return dict(numerator=numerator,denominator=velocity**2-cs2,
                adiabatic_sound_speed_squared=cs2,isothermal_sound_speed_squared=ct2,
                temperature_slope_constant=(heat+2*U)/gas.energy_dlnT,
                temperature_slope_velocity_coefficient=U/gas.energy_dlnT,
                material_heating_specific_logradius=heat)


def continuum_rhs(gas, **state):
    q=continuum_coefficients(gas,**state)
    if abs(q['denominator'])<1e-12*q['adiabatic_sound_speed_squared']:
        raise ValueError('Critical point requires regularity and a resolved branch slope')
    slope=q['numerator']/q['denominator']
    return np.array([slope,q['temperature_slope_constant']+q['temperature_slope_velocity_coefficient']*slope])


def adiabatic_eigenvalue_reference(*, mass, radius, density, pressure, gamma=4/3,
                                  outer_ratio=100., samples=501, rtol=1e-9,
                                  critical_offset=1e-6, root_tolerance=1e-12):
    """Accelerating polytropic Parker reference selected from base thermodynamics.

    Calorically ideal material u=P/[(gamma-1)rho], zero radiation and heating,
    a fixed point potential and a vacuum exterior. Restricted to 1<gamma<1.5,
    lambda=GM/(r_base*a_base²)>2. It is NOT primordial gas or an escaping
    radiation wind. No velocity or mdot input. The whole admissible scalar
    sonic-root interval is bracketed, so this reference can distinguish its
    finite hydrostatic branch from a root-finder failure analytically.
    """
    if (not np.all(np.isfinite([mass,radius,density,pressure,gamma,outer_ratio,rtol,critical_offset,root_tolerance]))
            or min(mass,radius,density,pressure)<=0 or not 1<gamma<1.5
            or outer_ratio<=1 or samples<11 or not 0<rtol<1e-3
            or not 0<critical_offset<1e-3 or not 0<root_tolerance<1e-4):
        raise ValueError('Outside declared adiabatic reference domain')
    a2=gamma*pressure/density;lam=G*mass/(radius*a2);h=1/(gamma-1);p=h-1.5
    if lam<=2:
        return dict(classification='model_domain_exit',reason='Base must be below an accelerating sonic point (lambda>2)',
                    lambda_base=lam,physical_quasistar=False)
    if lam>=h:
        surface=radius/(1-h/lam) if lam>h else None
        return dict(classification='hydrostatic_no_wind_reference',mdot=0.,lambda_base=lam,
                    vacuum_surface_radius=surface,finite_surface=surface is not None,
                    proof='Convex sonic residual is nonpositive at both admissible interval endpoints; no interior accelerating root',
                    physical_quasistar=False)
    maximum=(4/lam**2)**(1/p)
    residual=lambda x:lam**4/32*x**(2*p)+h-lam-p*x
    ratio=brentq(residual,0.,maximum,xtol=root_tolerance,rtol=max(root_tolerance,1e-14))
    sonic_radius=G*mass/(2*a2*ratio);sonic_density=density*ratio**h
    sonic_speed=np.sqrt(a2*ratio);mdot=4*np.pi*sonic_radius**2*sonic_density*sonic_speed
    if radius*outer_ratio<=sonic_radius:
        raise ValueError('Reference exterior does not contain selected sonic point')
    slope=(-2*(gamma-1)+np.sqrt(2*(5-3*gamma)))/(gamma+1)
    theta=-(gamma-1)*(2+slope)
    def state(x,z):
        r=sonic_radius*np.exp(x);v=sonic_speed*np.exp(z[0])
        rho=sonic_density*np.exp(-2*x-z[0]);P=rho*sonic_speed**2/gamma*np.exp(z[1])
        gas=GasDerivatives(P,P/((gamma-1)*rho),0.,P,P,0.,P/((gamma-1)*rho))
        return gas,dict(radius=r,density=rho,velocity=v,mass=mass,
                        gas_heating_per_volume=0.,radiative_acceleration=0.)
    def rhs(x,z):
        gas,pars=state(x,z);return continuum_rhs(gas,**pars)
    xbase=np.log(radius/sonic_radius);xend=np.log(radius*outer_ratio/sonic_radius)
    eps=min(critical_offset,-xbase*.01,xend*.01)
    pieces=[]
    for sign,end in [(-1,xbase),(1,xend)]:
        s=solve_ivp(rhs,[sign*eps,end],sign*eps*np.array([slope,theta]),
                    rtol=rtol,atol=rtol*.01,method='DOP853',max_step=.05,dense_output=True)
        if not s.success:
            raise RuntimeError('Numerical reference integration failed; not a no-root claim: '+s.message)
        pieces.append(s)
    x=np.unique(np.r_[np.linspace(xbase,0,samples),np.linspace(0,xend,samples)])
    z=np.zeros((2,len(x)))
    for mask,s in [(x<0,pieces[0]),(x>0,pieces[1])]:
        z[:,mask]=s.sol(x[mask])
    r=sonic_radius*np.exp(x);v=sonic_speed*np.exp(z[0]);rho=sonic_density*np.exp(-2*x-z[0])
    P=rho*sonic_speed**2/gamma*np.exp(z[1]);u=P/((gamma-1)*rho)
    B=v*v/2+u+P/rho-G*mass/r
    return dict(classification='stationary_wind_reference',mdot=mdot,
                sonic_radius=sonic_radius,sonic_density=sonic_density,sonic_speed=sonic_speed,
                base_velocity=mdot/(4*np.pi*radius**2*density),lambda_base=lam,
                terminal_speed=np.sqrt(2*p)*sonic_speed,
                base_density_relative_error=float(abs(rho[0]/density-1)),
                base_pressure_relative_error=float(abs(P[0]/pressure-1)),
                radius=r,velocity=v,density=rho,pressure=P,internal_energy=u,
                entropy_proxy=P/rho**gamma,bernoulli=B,
                mass_flux=4*np.pi*r*r*rho*v,energy_flux=mdot*B,
                physical_quasistar=False,radiation_solved=False,
                sonic_root_residual=float(residual(ratio)))
