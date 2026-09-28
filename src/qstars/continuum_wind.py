"""A restricted pressure-free grey continuum wind verification family.

No recombination microphysics is replaced: this module cannot be selected by
the stellar evolution driver. The mass supply and base velocity are inputs.
It tests the stationary momentum/escaping-transfer equations, not sonic launch.
Owocki et al. (2017) Eq. 24 is an independent small-drag reference.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from .grey_transfer import solve_moving_grey, TransferConvergenceError
from .constants import G, C, SIGMA


def flux_driven_reference(radius_ratio, gamma, tiring, w_base):
    """Owocki et al. Eq.24 with explicit finite initial kinetic energy.

    L(1)=L0, w=v^2/v_escape^2, m=mdot*GM/(R*L0), x=1-1/r.
    Radiation enthalpy/drag and gas pressure are omitted in this reference only.
    """
    r = np.asarray(radius_ratio, float)
    if (not np.all(np.isfinite(r)) or np.any(r < 1) or not np.isfinite(gamma)
            or gamma <= 1 or not np.isfinite(tiring) or tiring <= 0
            or not np.isfinite(w_base) or w_base <= 0):
        raise ValueError('Positive defined flux-driven reference inputs required')
    x = 1-1/r
    ell = np.exp(-tiring*gamma*x)
    w = w_base-np.expm1(-tiring*gamma*x)/tiring-x
    if np.any(w <= 0):
        raise ValueError('Specified analytic branch stalls in the requested interval')
    return w, ell


def solve_continuum_wind(radius_ratio, *, gamma, tiring, escape_speed_c, w_base=1e-4,
                         tolerance=1e-5, transport_tolerance=1e-6,
                         core_rays=32, max_iterations=60):
    """Iterate Newtonian wind momentum and moving grey angular radiation.

    Fixed point mass, constant comoving isotropic scattering opacity. Gas
    pressure/enthalpy absent. Velocities are restricted to <=.01c; the mechanical
    equation is the Newtonian limit, retaining the comoving radiative force.
    Density is from the rest-mass current rho0*gamma*v; the radiation work and
    Newtonian kinetic-plus-potential flux then have equal opposite sources.
    This is not a relativistic wind solver or a temperature/ionization solver.
    """
    r = np.asarray(radius_ratio, float)
    if r.ndim != 1 or len(r) < 8 or r[0] != 1 or np.any(np.diff(r) <= 0):
        raise ValueError('Increasing radius ratios beginning at exactly one required')
    if (not np.isfinite(escape_speed_c) or not 0 < escape_speed_c < .01
            or not 0 < tolerance < .1 or max_iterations < 1):
        raise ValueError('Invalid low-velocity wind controls')
    x = 1-1/r
    w, _ = flux_driven_reference(r, gamma, tiring, w_base)
    history = []
    for iteration in range(max_iterations):
        beta = escape_speed_c*np.sqrt(w)
        if max(beta) > .01:
            raise ValueError('Pressure-free Newtonian control exceeds 0.01c')
        lorentz = 1/np.sqrt(1-beta*beta)
        # chi0*R = kappa*rho0*R; mass current is invariant, no density clipping.
        chi = tiring*gamma/(lorentz*beta*r*r)
        rt = solve_moving_grey(r, chi, beta, tolerance=transport_tolerance, core_rays=core_rays)
        grad = gamma*r*r*rt['H_comoving']-1
        integral = PchipInterpolator(x, grad).antiderivative()
        candidate = w_base+integral(x)-integral(0)
        if np.any(candidate <= 0):
            raise TransferConvergenceError('An iteration leaves the outward branch; no no-root claim', history)
        error = float(max(abs(candidate/w-1)))
        history.append(dict(iteration=iteration+1, momentum_fixed_point_error=error,
                            radiation=rt['history'][-1]))
        if error < tolerance:
            # Keep the actual state used by the final radiation solve, not a newer
            # unverified velocity iterate. The residual above quantifies its defect.
            ell = r*r*rt['H']
            total = ell+tiring*(w-1/r)
            tau = -PchipInterpolator(r, chi*lorentz*(1-beta)).antiderivative()(r)
            tau -= tau[-1]
            photo = float(np.interp(2/3, tau[::-1], r[::-1])) if tau[0] >= 2/3 else None
            rt.update(w=w, luminosity_ratio=ell, total_energy_ratio=total,
                      optical_depth_radial_to_boundary=tau, photospheric_radius_ratio=photo,
                      radiation_drag_fraction=beta*(rt['J']+rt['K'])/rt['H'],
                      history_wind=history, mass_loss_is_input=True,
                      thermal_structure_validated=False, physical_envelope_match=False,
                      parameters=dict(gamma=gamma, tiring=tiring, escape_speed_c=escape_speed_c,
                                      w_base=w_base, numerical_damping=.5))
            return rt
        # Numerical damping only, fixed and recorded; it changes no equation.
        w = .5*(w+candidate)
    raise TransferConvergenceError('Continuum wind fixed point did not converge', history)


def dimensional_fields(solution, mass, grey_opacity):
    """A dimensional scaling, not an assertion of physical H/He ionization.

    Grey opacity is a declared constant comoving scattering opacity. The
    scattering problem supplies no material T, entropy, enthalpy or population.
    """
    if not np.isfinite(mass) or mass <= 0 or not np.isfinite(grey_opacity) or grey_opacity <= 0:
        raise ValueError('Positive finite mass and constant grey opacity required')
    p = solution['parameters']
    radius_unit = 2*G*mass/(p['escape_speed_c']*C)**2
    luminosity_unit = p['gamma']*4*np.pi*C*G*mass/grey_opacity
    intensity_unit = luminosity_unit/(16*np.pi**2*radius_unit**2)
    mdot = p['tiring']*luminosity_unit*radius_unit/(G*mass)
    radius = solution['radius']*radius_unit
    velocity = solution['beta']*C
    lorentz = 1/np.sqrt(1-solution['beta']**2)
    density_comoving = mdot/(4*np.pi*radius**2*lorentz*velocity)
    photo = solution['photospheric_radius_ratio']
    photo = photo*radius_unit if photo is not None else None
    luminosity = solution['luminosity_ratio']*luminosity_unit
    emergence_angular = 16*np.pi**2*radius[-1]**2*solution['angular_H'][-1]*intensity_unit
    return dict(radius_cm=radius, velocity_cm_s=velocity, density_comoving=density_comoving,
                mass_loss_g_s=mdot, enclosed_mass_g=mass, grey_opacity_cgs=grey_opacity,
                luminosity_erg_s=luminosity,
                radiation_energy_erg_cm3=4*np.pi*solution['J']*intensity_unit/C,
                radiation_pressure_rr=4*np.pi*solution['K']*intensity_unit/C,
                radiative_flux=4*np.pi*solution['H']*intensity_unit,
                radiative_acceleration=grey_opacity*4*np.pi*solution['H_comoving']*intensity_unit/C,
                photospheric_radius_cm=photo,
                effective_temperature_K=(luminosity[-1]/(4*np.pi*SIGMA*photo**2))**.25 if photo else None,
                emergent_luminosity_erg_s=float(luminosity[-1]),
                angular_emergent_luminosity_erg_s=float(emergence_angular),
                gas_temperature_K=None, ionization=None, entropy=None, gas_enthalpy=None,
                sonic_radius_cm=None, temperature_is_spectral_prediction=False)
