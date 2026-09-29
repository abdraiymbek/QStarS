"""Conservative lab-frame fluxes through a moving, spherical matching face.

Newtonian material, arbitrary realizable lab radiation moments. No LTE,
Eddington factor, entropy continuity, wind law, or interface update is assumed.
The potential must be static and use the same energy gauge on both sides.
"""
from dataclasses import dataclass
import numpy as np
from .constants import C


@dataclass(frozen=True)
class InterfaceState:
    radius: float
    enclosed_mass: float
    density: float
    velocity: float
    gas_pressure: float
    gas_internal_energy: float  # specific; includes chemical energy once
    radiation_energy: float
    radiation_flux: float
    radiation_pressure_rr: float
    potential: float
    energy_gauge: str
    gas_temperature: float | None = None
    gas_entropy: float | None = None

    def validate(self):
        values = [self.radius, self.enclosed_mass, self.density, self.velocity,
                  self.gas_pressure, self.gas_internal_energy, self.radiation_energy,
                  self.radiation_flux, self.radiation_pressure_rr, self.potential]
        if not np.all(np.isfinite(values)) or not self.energy_gauge:
            raise ValueError('Finite interface values and an explicit energy gauge required')
        if min(self.radius, self.enclosed_mass, self.density) <= 0 or self.gas_pressure < 0:
            raise ValueError('Positive geometry, mass, density and nonnegative gas pressure required')
        if abs(self.velocity) > .01*C:
            raise ValueError('Newtonian interface is restricted to |v|/c <= 0.01')
        E, F, P = self.radiation_energy, self.radiation_flux/C, self.radiation_pressure_rr
        if E < 0 or P < 0 or P > E*(1+1e-12) or F*F > E*P*(1+1e-12):
            raise ValueError('Radiation moments are not realizable by nonnegative intensities')
        if self.gas_temperature is not None and (not np.isfinite(self.gas_temperature) or self.gas_temperature <= 0):
            raise ValueError('Invalid material temperature')
        if self.gas_entropy is not None and not np.isfinite(self.gas_entropy):
            raise ValueError('Invalid material entropy')


def moving_flux(state, face_velocity=0.):
    """Reynolds flux = stationary lab flux minus face speed times density.

    Radiation enthalpy advection is already contained in the supplied lab F.
    Adding 4 P_rad v would count it twice. Gas work has both relative and
    moving-face terms; a material face still transmits pressure work.
    """
    state.validate()
    if not np.isfinite(face_velocity) or abs(face_velocity) > .01*C:
        raise ValueError('Finite Newtonian face speed |u|/c <= 0.01 required')
    s = state
    area = 4*np.pi*s.radius**2
    mdot = area*s.density*(s.velocity-face_velocity)
    channels = dict(gas_internal_energy=mdot*s.gas_internal_energy,
                    gas_pressure_work=area*s.gas_pressure*(s.velocity-face_velocity),
                    moving_gas_pressure_work=area*s.gas_pressure*face_velocity,
                    kinetic_energy=mdot*s.velocity**2/2,
                    potential_energy=mdot*s.potential,
                    radiation_lab_flux=area*s.radiation_flux,
                    radiation_swept_energy=-area*face_velocity*s.radiation_energy)
    return dict(mass=mdot,
                momentum=mdot*s.velocity+area*(s.gas_pressure+s.radiation_pressure_rr
                                             -face_velocity*s.radiation_flux/C**2),
                total_energy=sum(channels.values()), energy_channels=channels,
                material_entropy=None if s.gas_entropy is None else mdot*s.gas_entropy,
                radiation_entropy=None, face_velocity=face_velocity,
                escaping_radiation_certified=False)


def compare_faces(inner, outer, *, face_velocity=0., relative_tolerance=1e-8):
    """Source-free, zero-thickness jump check; no smoothness condition imposed.

    A resolved finite layer needs its curvature, gravity and storage terms;
    this routine must not be used to compare two different radii directly.
    No statement of a physical wind or branch existence follows from a pass.
    """
    if not np.isfinite(relative_tolerance) or not 0 < relative_tolerance < 1:
        raise ValueError('Invalid interface tolerance')
    for name in ('radius', 'enclosed_mass', 'potential'):
        if not np.isclose(getattr(inner, name), getattr(outer, name), rtol=1e-12, atol=0):
            raise ValueError('Thin matching interface requires identical geometry and potential')
    if inner.energy_gauge != outer.energy_gauge:
        raise ValueError('Gas and gravitational energy gauges must agree')
    a, b = moving_flux(inner, face_velocity), moving_flux(outer, face_velocity)
    defects = {}
    for key in ('mass', 'momentum', 'total_energy'):
        scale = max(abs(a[key]), abs(b[key]))
        if key == 'total_energy':
            scale = max(sum(abs(v) for v in a['energy_channels'].values()),
                        sum(abs(v) for v in b['energy_channels'].values()))
        delta = b[key]-a[key]
        defects[key] = dict(difference=delta, scale=scale,
                            relative=abs(delta)/scale if scale else 0.)
    return dict(conservation_passed=all(v['relative'] <= relative_tolerance for v in defects.values()),
                defects=defects, inner_flux=a, outer_flux=b,
                continuity_imposed=False, physical_solution_certified=False)
