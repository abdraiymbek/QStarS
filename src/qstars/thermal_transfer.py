"""Static frequency-band absorption, emission and isotropic scattering controls.

The absorption coefficient is constant within each declared frequency band;
Kirchhoff emission uses the *material* temperature. Band boundaries are h nu/k
in kelvin, with endpoints zero and infinity. These are specified-coefficient
verification problems, NOT H/He opacity, population kinetics, or wind models.
No Rosseland opacity is accepted as an absorption or flux mean by this API.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.integrate import simpson
from .constants import C, SIGMA
from .grey_transfer import SphericalRays, TransferConvergenceError

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)


def planck_bands(temperature, photon_temperature_edges):
    """Integrated B_nu and dB_band/dT; 64-point Planck quadrature.

    Planck tail beyond h nu/kT=80 is below 2e-29 of the total. The sum of
    bands is exactly sigma T^4/pi, not an independently renormalized spectrum.
    """
    T = np.atleast_1d(np.asarray(temperature, float))
    e = np.asarray(photon_temperature_edges, float)
    if (T.ndim != 1 or np.any(T <= 0) or not np.all(np.isfinite(T)) or e.ndim != 1
            or len(e) < 2 or e[0] != 0 or e[-1] != np.inf
            or not np.all(np.isfinite(e[:-1])) or np.any(np.diff(e) <= 0)):
        raise ValueError('Positive finite temperatures and ordered 0..infinity bands required')
    x = np.minimum(e[:, None]/T, 80.)
    q = x[..., None]*(_NODES+1)/2
    integral = x/2*np.sum(_WEIGHTS*q**3/np.expm1(np.maximum(q, 1e-100)), axis=-1)
    cumulative = integral*(15/np.pi**4)
    cumulative[0] = 0.; cumulative[-1] = 1.
    fraction = np.diff(cumulative, axis=0)
    # Roundoff in differences of already saturated cumulative integrals.
    if np.min(fraction) < -1e-13:
        raise ValueError('Planck band quadrature lost positivity')
    fraction = np.maximum(fraction, 0.)
    derivative_edge = x**4/np.expm1(np.maximum(x, 1e-100))*(15/np.pi**4)
    derivative_edge[[0, -1]] = 0.
    B = SIGMA/np.pi*T**4*fraction
    dB = SIGMA/np.pi*T**3*(4*fraction-np.diff(derivative_edge, axis=0))
    return B, dB


def opacity_means(absorption, scattering, temperature, edges, J, H):
    """Return four distinct means; a signed band flux is retained in force.

    Inputs are opacities per mass (not extinction per length). A vanishing
    total flux or radiation energy has no defined corresponding mean.
    """
    ka, ks, J, H = np.broadcast_arrays(absorption, scattering, J, H)
    B, dB = planck_bands(temperature, edges)
    if ka.shape != B.shape or np.any(ka < 0) or np.any(ks < 0) or np.any(ka+ks <= 0):
        raise ValueError('Band opacity and moment shapes must agree; positive total opacity required')
    if not all(np.all(np.isfinite(a)) for a in (ka, ks, J, H)) or np.any(J < 0):
        raise ValueError('Finite nonnegative intensity and finite opacities required')
    if np.any(J.sum(axis=0) <= 0) or np.any(abs(H.sum(axis=0)) < 1e-30*np.maximum(J.sum(axis=0), 1e-300)):
        raise ValueError('Zero radiation energy or net flux: corresponding mean is undefined')
    kt = ka+ks
    return dict(planck_absorption=(ka*B).sum(axis=0)/B.sum(axis=0),
                energy_absorption=(ka*J).sum(axis=0)/J.sum(axis=0),
                rosseland_total=dB.sum(axis=0)/(dB/kt).sum(axis=0),
                flux_total=(kt*H).sum(axis=0)/H.sum(axis=0),
                acceleration=4*np.pi/C*(kt*H).sum(axis=0))


class ThermalBands:
    """Frozen coefficients and density on a static spherical grid.

    Prescribed outgoing blackbody core intensity; inward core light absorbed;
    zero incoming exterior radiation. Net core luminosity is an output. Direct
    formal transfer is independently audited for flux error on every result.
    Algebraic convergence alone is never a conservation certification.
    """
    def __init__(self, radius, absorption_extinction, scattering_extinction, edges,
                 core_temperature, *, core_rays=32, path_order=8):
        self.r = np.asarray(radius, float)
        self.ka = np.asarray(absorption_extinction, float)
        self.ks = np.asarray(scattering_extinction, float)
        self.edges = np.asarray(edges, float)
        self.core_temperature = float(core_temperature)
        self.core_rays = core_rays
        self.core_B = planck_bands([core_temperature], edges)[0][:, 0]
        shape = (len(self.core_B), len(self.r))
        if (self.ka.shape != shape or self.ks.shape != shape
                or not np.all(np.isfinite(self.ka)) or not np.all(np.isfinite(self.ks))
                or np.any(self.ka < 0) or np.any(self.ks < 0)
                or np.any(self.ka+self.ks <= 0) or np.any(self.ka.sum(axis=0) <= 0)):
            raise ValueError('Explicit nonnegative absorption/scattering in every band required; thermal coupling cannot vanish')
        self.kt = self.ka+self.ks
        self.rays = [SphericalRays(self.r, chi, core_rays, path_order) for chi in self.kt]

    def fixed_temperature(self, temperature, *, tolerance=1e-10, max_iterations=200):
        if not 0 < tolerance < .01 or max_iterations < 1:
            raise ValueError('Invalid transfer tolerance or iteration limit')
        T = np.asarray(temperature, float)
        B, _ = planck_bands(T, self.edges)
        if B.shape != self.ka.shape:
            raise ValueError('One material temperature per radial node required')
        n = len(self.r); scale = (self.r/self.r[0])**2
        moments = []; histories = []
        for g, rays in enumerate(self.rays):
            albedo = self.ks[g]/self.kt[g]
            core = np.full(self.core_rays, self.core_B[g])
            boundary = rays.moments(np.zeros(n), core)[0]
            thermal = self.ka[g]/self.kt[g]*B[g]
            rhs = (thermal+albedo*boundary)*scale
            def action(value):
                j = rays.moments(value/scale, np.zeros(self.core_rays), _allow_signed=True)[0]
                return value-albedo*j*scale
            history = []
            op = LinearOperator((n, n), matvec=action)
            source, info = gmres(op, rhs, rtol=tolerance, atol=0., restart=min(n, 40),
                                 maxiter=max_iterations, callback=history.append, callback_type='pr_norm')
            if info or np.any(source < 0):
                raise TransferConvergenceError('Band source solve failed; not a physical nonexistence result', history)
            moments.append(rays.moments(source/scale, core)); histories.append(history)
        J, H, K = np.asarray(moments).transpose(1, 0, 2)
        heating = 4*np.pi*(self.ka*(J-B)).sum(axis=0)
        luminosity = 16*np.pi**2*self.r**2*H.sum(axis=0)
        integrated = simpson(4*np.pi*self.r**2*heating, x=self.r)
        defect = (luminosity[-1]-luminosity[0]+integrated)/max(abs(luminosity[0]), abs(luminosity[-1]))
        return dict(radius=self.r.copy(), gas_temperature=T.copy(), J=J, H=H, K=K,
                    radiation_energy=4*np.pi/C*J.sum(axis=0),
                    radiation_temperature=(np.pi/SIGMA*J.sum(axis=0))**.25,
                    radiation_flux=4*np.pi*H.sum(axis=0),
                    radiation_pressure_rr=4*np.pi/C*K.sum(axis=0),
                    luminosity=luminosity, heating_per_volume=heating,
                    energy_integral_relative_defect=float(abs(defect)),
                    histories=histories, population_kinetics_solved=False,
                    wind_eigenvalue_solved=False, coefficients_are_prescribed=True)

    def radiative_equilibrium(self, *, tolerance=1e-7, transfer_tolerance=1e-10,
                              max_iterations=100):
        """Solve sum chi_abs (J-B(Tgas))=0, without imposing Tgas=Trad.

        Static material with no mechanical/chemical heating or advection.
        Coefficients remain frozen: this does not solve temperature-dependent
        H/He populations and must not be applied as an accelerating-wind EOS.
        """
        if not 0 < tolerance < .01 or max_iterations < 1:
            raise ValueError('Invalid equilibrium tolerance or iteration limit')
        dilution = .5/(self.r/self.r[0])**2/(1+np.sqrt(1-(self.r[0]/self.r)**2))
        T = self.core_temperature*dilution**.25
        history = []
        for iteration in range(max_iterations):
            q = self.fixed_temperature(T, tolerance=transfer_tolerance)
            heat = (self.ka*q['J']).sum(axis=0)
            B, _ = planck_bands(T, self.edges)
            error = float(np.max(abs((self.ka*B).sum(axis=0)/heat-1)))
            history.append(dict(iteration=iteration, thermal_balance_error=error,
                                energy_integral_error=q['energy_integral_relative_defect']))
            if error <= tolerance:
                q['equilibrium_history'] = history
                q['thermal_balance_error'] = error
                return q
            # Monotonic local Planck integral: safeguarded log-temperature Newton.
            updated = T.copy()
            for _ in range(60):
                B, dB = planck_bands(updated, self.edges)
                emission = (self.ka*B).sum(axis=0)
                if max(abs(emission/heat-1)) < 1e-12:
                    break
                step = (emission-heat)/(self.ka*dB*updated).sum(axis=0)
                updated *= np.exp(-np.clip(step, -.5, .5))
            else:
                raise TransferConvergenceError('Local material-temperature root did not converge', history)
            T = updated
        raise TransferConvergenceError('Radiative equilibrium iteration exhausted; no no-solution inference', history)
