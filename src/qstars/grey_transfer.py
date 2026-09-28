"""Conservative spherical grey transfer with an angular VEF closure.

Hummer & Rybicki (1971), equations 2.1--2.31. This is a bolometric
component, not a Rosseland-opacity wind extension. Extinction must be a
specified grey coefficient. No material temperature or ionization is inferred
for conservative scattering. See docs/escaping_radiation.md.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.sparse.linalg import LinearOperator, gmres
from .constants import C, SIGMA


class TransferConvergenceError(RuntimeError):
    def __init__(self, message, history):
        super().__init__(message)
        self.history = history


def _moment_propagation(t, f, chi, beta, outer, tolerance):
    """Integrate the linear moment system at every interpolation knot.

    Sixth-order Gauss collocation, with local step doubling of the 2x2
    fundamental matrix. Vectorized stages keep the per-knot method affordable.
    Unlike a single adaptive ODE solve over the entire domain, this cannot skip
    a narrow feature between nearby tabulation knots. No flux correction.
    """
    fi = PchipInterpolator(t, f)
    ci = PchipInterpolator(t, np.log(chi))
    bi = PchipInterpolator(t, beta)
    root = np.sqrt(15.)
    nodes = np.array([.5-root/10, .5, .5+root/10])
    weights = np.array([5/18, 4/9, 5/18])
    rk = np.array([[5/36, 2/9-root/15, 5/36-root/30],
                   [5/36+root/24, 2/9, 5/36-root/24],
                   [5/36+root/30, 2/9+root/15, 5/36]])
    def step(upper, lower):
        h = lower-upper
        x = upper[:, None]+h[:, None]*nodes
        invf, b = 1/fi(x), bi(x)
        interaction = np.exp(x+ci(x))/(1-b*b)**1.5
        matrix = np.empty((*x.shape, 2, 2))
        matrix[..., 0, 0] = invf-1+interaction*b*(1+invf)
        matrix[..., 0, 1] = -interaction*(1+b*b)
        matrix[..., 1, 0] = interaction*b*b*(1+invf)
        matrix[..., 1, 1] = -interaction*b*(1+b*b)
        # Unknown stage states Y_j, each a fundamental 2x2 matrix.
        block = np.eye(6)[None, :, :]-h[:, None, None]*(rk[None, :, None, :, None]*matrix[:, None, :, :, :].transpose(0, 1, 3, 2, 4)).reshape(-1, 6, 6)
        right = np.broadcast_to(np.tile(np.eye(2), (3, 1)), (len(h), 6, 2))
        stages = np.linalg.solve(block, right).reshape(-1, 3, 2, 2)
        return np.eye(2)[None, :, :]+h[:, None, None]*np.sum(weights[None, :, None, None]*(matrix @ stages), axis=1)
    def refined(upper, lower, depth=0):
        middle = (upper+lower)/2
        full = step(upper, lower)
        half = step(middle, lower) @ step(upper, middle)
        error = np.max(abs(half-full), axis=(1, 2))/np.maximum(1., np.max(abs(half), axis=(1, 2)))
        bad = error > tolerance
        if np.any(bad):
            if depth >= 12:
                raise TransferConvergenceError('Moment quadrature refinement exhausted', [])
            half[bad] = refined(middle[bad], lower[bad], depth+1) @ refined(upper[bad], middle[bad], depth+1)
        return half
    propagators = refined(t[1:], t[:-1])
    state = np.empty((2, len(t)))
    state[:, -1] = outer
    for i in range(len(t)-2, -1, -1):
        state[:, i] = propagators[i] @ state[:, i+1]
    return state


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
    def __init__(self, radius, extinction, core_rays=24, path_order=8, beta=None):
        self.r = np.asarray(radius, float)
        chi = np.asarray(extinction, float)
        if (self.r.ndim != 1 or len(self.r) < 8 or chi.shape != self.r.shape
                or not np.all(np.isfinite(self.r)) or not np.all(np.isfinite(chi))
                or self.r[0] <= 0 or np.any(np.diff(self.r) <= 0) or np.any(chi <= 0)):
            raise ValueError('Increasing positive radii and positive finite grey extinction required')
        if core_rays < 8 or path_order < 2:
            raise ValueError('At least eight core rays and second-order path quadrature required')
        self.beta = np.zeros_like(self.r) if beta is None else np.asarray(beta, float)
        if self.beta.shape != self.r.shape or not np.all(np.isfinite(self.beta)) or np.any(abs(self.beta) >= 1):
            raise ValueError('Finite subluminal radial beta on each radial node required')
        self.core_rays = core_rays
        self.p = np.r_[self.r[0]*np.sin(np.linspace(0, np.pi/2, core_rays, endpoint=True)), self.r]
        n, m = len(self.r), len(self.p)
        dt_plus = np.zeros((n-1, m))
        dt_minus = np.zeros_like(dt_plus)
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
            bq = self.beta[j]+(self.beta[j+1]-self.beta[j])*np.log(rq/self.r[j])/np.log(self.r[j+1]/self.r[j])
            gamma = 1/np.sqrt(1-bq*bq)
            dt_plus[j, :count] = (z1-z0)/2*((cq*gamma*(1-bq*z/rq)) @ w)
            dt_minus[j, :count] = (z1-z0)/2*((cq*gamma*(1+bq*z/rq)) @ w)
        self.weights_plus = _linear_weights(dt_plus)
        self.weights_minus = _linear_weights(dt_minus)
        mu = np.sqrt(np.maximum(1-(self.p[None, :]/self.r[:, None])**2, 0))
        b = self.beta[:, None]
        gamma = 1/np.sqrt(1-b*b)
        self.boost_plus = (gamma*(1-b*mu))**-4
        self.boost_minus = (gamma*(1+b*mu))**-4

    def moments(self, source, core_intensity, *, _allow_signed=False):
        """Return J,H,K from independent angular integrals; no flux renormalization."""
        source = np.asarray(source, float)
        core_intensity = np.asarray(core_intensity, float)
        if (source.shape != self.r.shape or core_intensity.shape != (self.core_rays,)
                or not np.all(np.isfinite(source)) or (not _allow_signed and np.any(source < 0))
                or not np.all(np.isfinite(core_intensity)) or np.any(core_intensity < 0)):
            raise ValueError('Finite nonnegative source and outgoing core intensities required')
        n, m = len(self.r), len(self.p)
        outward, inward = np.zeros((n, m)), np.zeros((n, m))
        trans, up, down = self.weights_minus
        sm = source[:, None]*self.boost_minus
        sp = source[:, None]*self.boost_plus
        for j in range(n-2, -1, -1):
            s = slice(0, self.core_rays+j+1)
            inward[j, s] = trans[j, s]*inward[j+1, s]+up[j, s]*sm[j+1, s]+down[j, s]*sm[j, s]
        trans, up, down = self.weights_plus
        outward[0, :self.core_rays] = core_intensity
        outward[0, self.core_rays] = inward[0, self.core_rays]
        for j in range(1, n):
            s = slice(0, self.core_rays+j)
            outward[j, s] = trans[j-1, s]*outward[j-1, s]+up[j-1, s]*sp[j-1, s]+down[j-1, s]*sp[j, s]
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
               core_rays=24, path_order=8, moment_tolerance=1e-11):
    """Solve the static conservative grey Milne problem with a deep diffusion core.

    No incident radiation at the outer sphere; central *net* luminosity specified.
    The inner diffusion asymptotic I+=J+3 mu H must be tested by moving that boundary.
    This API accepts no gas state, Rosseland table or velocity: it cannot silently
    turn a recombining diffusion wind into a validated escaping wind.
    """
    if (not np.isfinite(luminosity) or luminosity <= 0 or not 0 < tolerance < .1
            or not 0 < moment_tolerance < .1 or max_iterations < 1):
        raise ValueError('Positive luminosity and valid iteration controls required')
    r_input, chi_input = np.asarray(radius, float), np.asarray(extinction, float)
    # Validate before choosing numerical units.
    rays = SphericalRays(r_input, chi_input, core_rays, path_order)
    unit = r_input[0]
    r = r_input/unit
    chi = chi_input*unit
    t = np.log(r)
    f = np.full(len(r), 1/3)
    g = .5
    history = []
    core_mu = np.sqrt(1-(rays.p[:core_rays]/r_input[0])**2)
    previous = None
    for iteration in range(max_iterations):
        state = _moment_propagation(t, f, chi, np.zeros_like(r), [f[-1]/g, 1.], moment_tolerance)
        J = state[0]/(f*r*r)
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


def comoving_moments(J, H, K, beta):
    """Exact bolometric Lorentz transformation (Fukue 2018, equations 4--6)."""
    b = np.asarray(beta, float)
    if not np.all(np.isfinite(b)) or np.any(abs(b) >= 1):
        raise ValueError('Subluminal beta required')
    gamma2 = 1/(1-b*b)
    return (gamma2*(J-2*b*H+b*b*K),
            gamma2*((1+b*b)*H-b*(J+K)),
            gamma2*(b*b*J-2*b*H+K))


def solve_scattering(radius, extinction, beta, core_intensity, *, tolerance=1e-9,
                     max_iterations=300, core_rays=32, path_order=8):
    """Steady angular transfer in a specified moving grey scattering medium.

    Isotropic coherent scattering in the comoving frame, no true absorption,
    no thermal or population equation. The formal source is D**4 J_comoving;
    path extinction is chi_comoving/D. The core absorbs inward photons and
    emits the supplied *lab-frame* intensity; exterior incident intensity is zero.
    A matrix-free linear solve replaces slow scattering-source iteration.
    It does not determine the supplied matter velocity or extinction.
    """
    if not 0 < tolerance < .1 or max_iterations < 1:
        raise ValueError('Invalid scattering solver controls')
    rays = SphericalRays(radius, extinction, core_rays, path_order, beta)
    core = np.broadcast_to(np.asarray(core_intensity, float), (core_rays,))
    if not np.all(np.isfinite(core)) or np.any(core < 0) or not np.any(core > 0):
        raise ValueError('A positive core radiation source is required')
    n = len(rays.r)
    zero = np.zeros(n)
    boundary = rays.moments(zero, core)
    rhs = comoving_moments(*boundary, rays.beta)[0]
    # Scale the unknown to remove geometric dilution from the linear solver norm.
    scale = (rays.r/rays.r[0])**2
    history = []
    def operator(s):
        m = rays.moments(s/scale, np.zeros(core_rays), _allow_signed=True)
        return s-comoving_moments(*m, rays.beta)[0]*scale
    op = LinearOperator((n, n), matvec=operator, dtype=float)
    sol, info = gmres(op, rhs*scale, rtol=tolerance, atol=0, restart=60,
                     maxiter=max_iterations, callback=lambda e: history.append(float(e)),
                     callback_type='legacy')
    source = sol/scale
    if info != 0 or np.any(source < 0):
        raise TransferConvergenceError('Scattering source solve failed; no physical nonexistence implied', history)
    J, H, K = rays.moments(source, core)
    J0, H0, K0 = comoving_moments(J, H, K, rays.beta)
    local_error = float(np.max(abs(J0/source-1)))
    if local_error > max(50*tolerance, 1e-10):
        raise TransferConvergenceError('Local scattering source residual exceeds tolerance', history)
    return dict(radius=rays.r, extinction=np.asarray(extinction), beta=rays.beta,
                J=J, H=H, K=K, J_comoving=J0, H_comoving=H0, K_comoving=K0,
                source=source, source_relative_residual=local_error, history=history)


def solve_moving_grey(radius, extinction, beta, *, tolerance=1e-6, max_iterations=80,
                      core_rays=32, path_order=8, moment_tolerance=1e-11):
    """Conservative moments closed by exact-frame grey angular transport.

    Fukue (2018) scattering transfer/moment equations, with net inner H=1 in
    intensity units. The core emits uniform lab intensity and absorbs returning
    photons. This normalization determines the emitted intensity, not its shape.
    No gas temperature, ionization, sonic condition or hydrodynamics is solved.
    The returned radius is in caller units; J,H,K are in normalized intensity units.
    """
    if not 0 < tolerance < .1 or not 0 < moment_tolerance < .1 or max_iterations < 1:
        raise ValueError('Invalid moment solver controls')
    rays = SphericalRays(radius, extinction, core_rays, path_order, beta)
    if np.any(rays.beta < 0) or np.any(rays.beta > .01):
        raise ValueError('Moving VEF validation domain is outward 0 <= v/c <= 0.01')
    r = rays.r/rays.r[0]
    chi = np.asarray(extinction)*rays.r[0]
    b = rays.beta
    t = np.log(r)
    f, g = np.full(len(r), 1/3), .5
    previous = None
    history = []
    for iteration in range(max_iterations):
        state = _moment_propagation(t, f, chi, b, [f[-1]/g, 1.], moment_tolerance)
        if not np.all(np.isfinite(state)) or state[1, 0] <= 0:
            raise TransferConvergenceError('Moving moment integration failed', history)
        y, ell = state/state[1, 0]
        J, H, K = y/(f*r*r), ell/(r*r), y/(r*r)
        source = comoving_moments(J, H, K, b)[0]
        # The first formal sweep obtains the returning core flux. The outgoing
        # intensity is uniform; its contribution to H at the core is I/4.
        _, hin, _ = rays.moments(source, np.zeros(core_rays))
        core_intensity = 4*(1-hin[0])
        aj, ah, ak = rays.moments(source, np.full(core_rays, core_intensity))
        new_f, new_g = ak/aj, ah[-1]/aj[-1]
        error = max(float(np.max(abs(new_f/f-1))), float(abs(new_g/g-1)))
        change_J = float(np.max(abs(J/previous-1))) if previous is not None else None
        history.append(dict(iteration=iteration+1, closure_change=error,
                            intensity_change=change_J,
                            angular_flux_max_error=float(np.max(abs(ah/H-1))),
                            formal_source_max_error=float(np.max(abs(aj/J-1)))))
        if not np.all(np.isfinite(new_f)) or np.any(new_f <= 0) or new_g <= 0:
            raise TransferConvergenceError('Invalid moving angular closure', history)
        if error < tolerance and change_J is not None and change_J < tolerance:
            J0, H0, K0 = comoving_moments(J, H, K, b)
            return dict(radius=rays.r, extinction=np.asarray(extinction), beta=b,
                        J=J, H=H, K=K, J_comoving=J0, H_comoving=H0, K_comoving=K0,
                        angular_J=aj, angular_H=ah, angular_K=ak,
                        eddington_factor=new_f, core_intensity=core_intensity, history=history)
        previous = J
        f, g = new_f, new_g
    raise TransferConvergenceError('Moving VEF solve did not converge', history)
