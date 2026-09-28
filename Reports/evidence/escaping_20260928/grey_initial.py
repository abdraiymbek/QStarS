"""Static conservative spherical grey transfer with an angular VEF closure.

Hummer & Rybicki (1971), equations 2.1--2.31. This is a bolometric
component, not a Rosseland-opacity wind extension. Extinction must be a
specified grey coefficient. No material temperature or ionization is inferred
for conservative scattering. See docs/escaping_radiation.md.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from .constants import C, SIGMA


class TransferConvergenceError(RuntimeError):
    def __init__(self, message, history):
        super().__init__(message)
        self.history = history


def _linear_weights(dt):
    """Exact formal solution for a source linear in optical path."""
    transmission = np.exp(-dt)
    small = dt < 1e-3
    down = np.empty_like(dt)
    d = dt[small]
    down[small] = d * (.5 + d * (-1/6 + d * (1/24 - d/120)))
    d = dt[~small]
    down[~small] = 1 + np.expm1(-d) / d
    up = -np.expm1(-dt) - down
    return transmission, up, down


def _angular_integral(mu, intensity, power):
    """Integrate mu**power times piecewise-linear I, without a singular p weight."""
    a, b = mu[:-1], mu[1:]
    width = b-a
    # Gauss quadrature exactly integrates the polynomial (degree <= 3).
    mid = (a+b)/2
    d = width/(2*np.sqrt(3))
    t = 1/np.sqrt(3)
    lo = .5*((1+t)*intensity[:-1]+(1-t)*intensity[1:])
    hi = .5*((1-t)*intensity[:-1]+(1+t)*intensity[1:])
    return np.sum(width/2*((mid-d)**power*lo+(mid+d)**power*hi))


class SphericalRays:
    """Tangent rays plus core rays; linear-source formal integration in both directions.

    Extinction between radial nodes is log-linear. Impact parameters include
    every radial node, resolving the increasingly narrow forward radiation cone.
    Arrays use any one consistent length unit; intensities have arbitrary units.
    """
    def __init__(self, radius, extinction, core_rays=24, path_order=8):
        self.r = np.asarray(radius, float)
        chi = np.asarray(extinction, float)
        if (self.r.ndim != 1 or len(self.r) < 8 or chi.shape != self.r.shape
                or not np.all(np.isfinite(self.r)) or not np.all(np.isfinite(chi))
                or self.r[0] <= 0 or np.any(np.diff(self.r) <= 0) or np.any(chi <= 0)):
            raise ValueError('Increasing positive radii and positive finite grey extinction required')
        if core_rays < 8 or path_order < 2:
            raise ValueError('At least eight core rays and second-order path quadrature required')
        self.core_rays = core_rays
        self.p = np.r_[self.r[0]*np.sin(np.linspace(0, np.pi/2, core_rays, endpoint=False)), self.r]
        n, m = len(self.r), len(self.p)
        dt = np.zeros((n-1, m))
        x, w = np.polynomial.legendre.leggauss(path_order)
        for j in range(n-1):
            count = core_rays+j+1
            p = self.p[:count]
            z0 = np.sqrt(np.maximum(self.r[j]**2-p*p, 0))
            z1 = np.sqrt(self.r[j+1]**2-p*p)
            z = (z0[:, None]+z1[:, None])/2+(z1-z0)[:, None]*x/2
            rq = np.sqrt(z*z+p[:, None]**2)
            slope = np.log(chi[j+1]/chi[j])/np.log(self.r[j+1]/self.r[j])
            cq = chi[j]*np.exp(slope*np.log(rq/self.r[j]))
            dt[j, :count] = (z1-z0)/2*(cq @ w)
        self.weights = _linear_weights(dt)

    def moments(self, source, core_intensity):
        """Return J,H,K from independent angular integrals; no flux renormalization."""
        source = np.asarray(source, float)
        core_intensity = np.asarray(core_intensity, float)
        if (source.shape != self.r.shape or core_intensity.shape != (self.core_rays,)
                or not np.all(np.isfinite(source)) or np.any(source < 0)
                or not np.all(np.isfinite(core_intensity)) or np.any(core_intensity < 0)):
            raise ValueError('Finite nonnegative source and outgoing core intensities required')
        n, m = len(self.r), len(self.p)
        outward, inward = np.zeros((n, m)), np.zeros((n, m))
        trans, up, down = self.weights
        for j in range(n-2, -1, -1):
            s = slice(0, self.core_rays+j+1)
            inward[j, s] = trans[j, s]*inward[j+1, s]+up[j, s]*source[j+1]+down[j, s]*source[j]
        outward[0, :self.core_rays] = core_intensity
        outward[0, self.core_rays] = inward[0, self.core_rays]
        for j in range(1, n):
            s = slice(0, self.core_rays+j)
            outward[j, s] = trans[j-1, s]*outward[j-1, s]+up[j-1, s]*source[j-1]+down[j-1, s]*source[j]
            outward[j, self.core_rays+j] = inward[j, self.core_rays+j]
        J, H, K = np.zeros(n), np.zeros(n), np.zeros(n)
        for j in range(n):
            s = slice(self.core_rays+j, None, -1)
            mu = np.sqrt(np.maximum(1-(self.p[:self.core_rays+j+1][::-1]/self.r[j])**2, 0))
            mean = (outward[j, s]+inward[j, s])/2
            diff = (outward[j, s]-inward[j, s])/2
            J[j] = _angular_integral(mu, mean, 0)
            H[j] = _angular_integral(mu, diff, 1)
            K[j] = _angular_integral(mu, mean, 2)
        return J, H, K


@dataclass
class GreySolution:
    radius: np.ndarray
    extinction: np.ndarray
    mean_intensity: np.ndarray
    angular_mean_intensity: np.ndarray
    angular_flux: np.ndarray
    eddington_factor: np.ndarray
    luminosity: float
    optical_depth_to_boundary: np.ndarray
    history: list

    def fields(self, grey_opacity=None):
        """CGS diagnostics. Opacity for acceleration must be the same grey law.

        The tau=2/3 radius is a convention, not a thermalization or color radius.
        No radius/Teff is reported if that surface is outside the supplied domain.
        """
        r = self.radius
        flux = self.luminosity/(4*np.pi*r*r)
        E = 4*np.pi*self.mean_intensity/C
        tau = self.optical_depth_to_boundary
        photo = float(np.interp(2/3, tau[::-1], r[::-1])) if tau[0] >= 2/3 else None
        result = dict(radius=r, extinction=self.extinction, radiation_energy=E,
                      radiative_flux=flux, angular_radiative_flux=self.angular_flux,
                      angular_luminosity=4*np.pi*r*r*self.angular_flux,
                      radiation_pressure_rr=self.eddington_factor*E,
                      eddington_factor=self.eddington_factor, flux_factor=flux/(C*E),
                      optical_depth_to_boundary=tau,
                      photospheric_radius=photo,
                      effective_temperature=(self.luminosity/(4*np.pi*SIGMA*photo**2))**.25 if photo else None,
                      emergent_luminosity=float(4*np.pi*r[-1]**2*self.angular_flux[-1]))
        if grey_opacity is not None:
            kappa = np.broadcast_to(np.asarray(grey_opacity, float), r.shape)
            if not np.all(np.isfinite(kappa)) or np.any(kappa <= 0):
                raise ValueError('Positive finite grey opacity required')
            result['radiative_acceleration'] = kappa*flux/C
        return result


def solve_grey(radius, extinction, luminosity, *, tolerance=1e-7, max_iterations=50,
               core_rays=24, path_order=8, moment_tolerance=1e-9):
    """Solve the static conservative grey Milne problem with a deep diffusion core.

    No incident radiation at the outer sphere; central *net* luminosity specified.
    The inner diffusion asymptotic I+=J+3 mu H must be tested by moving that boundary.
    This API accepts no gas state, Rosseland table or velocity: it cannot silently
    turn a recombining diffusion wind into a validated escaping wind.
    """
    if not np.isfinite(luminosity) or luminosity <= 0 or not 0 < tolerance < .1 or max_iterations < 1:
        raise ValueError('Positive luminosity and valid iteration controls required')
    r_input, chi_input = np.asarray(radius, float), np.asarray(extinction, float)
    # Validate before choosing numerical units.
    rays = SphericalRays(r_input, chi_input, core_rays, path_order)
    unit = r_input[0]
    r = r_input/unit
    chi = chi_input*unit
    t = np.log(r)
    chi_interp = PchipInterpolator(t, np.log(chi))
    f = np.full(len(r), 1/3)
    g = .5
    history = []
    core_mu = np.sqrt(1-(rays.p[:core_rays]/r_input[0])**2)
    previous = None
    for iteration in range(max_iterations):
        f_interp = PchipInterpolator(t, f)
        def rhs(x, y):
            return (1/float(f_interp(x))-1)*y-np.exp(x+chi_interp(x))
        sol = solve_ivp(rhs, (t[-1], t[0]), [f[-1]/g], method='DOP853',
                        t_eval=t[::-1], rtol=moment_tolerance, atol=moment_tolerance*1e-3)
        if not sol.success:
            raise TransferConvergenceError(sol.message, history)
        J = sol.y[0, ::-1]/(f*r*r)
        j, h, k = rays.moments(J, J[0]+3*core_mu)
        new_f, new_g = k/j, h[-1]/j[-1]
        error = max(float(np.max(abs(new_f/f-1))), float(abs(new_g/g-1)))
        change_J = float(np.max(abs(J/previous-1))) if previous is not None else None
        history.append(dict(iteration=iteration+1, closure_change=error,
                            intensity_change=change_J,
                            angular_luminosity_max_error=float(np.max(abs(r*r*h-1))),
                            formal_source_max_error=float(np.max(abs(j/J-1)))))
        if not np.all(np.isfinite(new_f)) or np.any(new_f <= 0) or new_g <= 0:
            raise TransferConvergenceError('Nonpositive angular closure', history)
        if error < tolerance and change_J is not None and change_J < tolerance:
            # Return the moments of the same last formal solution, keeping its flux error.
            factor = luminosity/(16*np.pi**2*unit**2)
            tau = -cumulative_trapezoid(chi_input[::-1], r_input[::-1], initial=0)[::-1]
            return GreySolution(r_input, chi_input, J*factor, j*factor, 4*np.pi*h*factor,
                                new_f, luminosity, tau, history)
        previous = J
        f, g = new_f, new_g
    raise TransferConvergenceError('Variable Eddington factor did not converge', history)
