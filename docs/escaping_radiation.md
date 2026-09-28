# Spherical escaping radiation: supported component and remaining boundary

This component solves **specified grey transfer**, including a moving coherent-scattering field. It is separate from the recombining FreeEOS/MESA diffusion wind. It does not turn a Rosseland table into a flux mean, infer LTE from scattering, or supply a physical envelope mass-loss eigenvalue. No evolution configuration selects this component yet.

## Reproducible targets

The complete [candidate register](../benchmarks/outer_candidates.json) records source availability, equations, EOS, opacity, radiation, sonic and boundary conditions, composition, input scales and available outputs for nine families. Public source availability takes priority, but an available prescribed-flow, line-driven or dust-driven calculation is not the same physical benchmark as a primordial continuum wind. Repository commits and acquired-file hashes are in `Reports/evidence/escaping_20260928/candidate_source_files.json`.

[Hummer & Rybicki (1971)](https://doi.org/10.1093/mnras/152.1.1), Table I, is the independent angular-transfer reference. Its final-iteration values were transcribed from the rendered table, not fitted. The dimensionless test has volume extinction chi=r^-2, outer radius 10, H0=1, and net L=16 pi². The exact central source is at infinite inward optical depth; finite diffusion cores at tau=30,100,300 are numerical truncations to test. Published accuracy is a few tenths of a percent, rising to about 1% at the outer surface. The saved comparison limits are 0.5% internally and 1% at the surface.

The selected **constructed stationary control** uses the precisely specified flux-driven limit in [Owocki, Townsend & Quataert (2017)](https://arxiv.org/abs/1708.07790). With x=1-R/r, w=v²/v_escape², Gamma=L0/L_Edd, and m=mdot GM/(R L0), energy conservation gives ell=L/L0=1-m(w-w_base+x). Neglecting gas pressure and radiation drag gives dw/dx=Gamma ell-1, hence

    ell = exp(-m Gamma x)
    w = w_base + [1-exp(-m Gamma x)]/m - x.

The finite base speed is an explicit boundary condition extending the paper's zero-speed limit, not a coefficient fit. The full numerical control includes moving angular radiation; the analytic solution is therefore an asymptotic comparison, with the omitted drag measured. Exact inputs are in `benchmarks/grey_continuum_wind.json`. The analytic model supplies no material temperature, ionization or sonic eigenvalue. A second parameter set tests the same implementation without retuning. A larger-drag case must not be forced to match this limiting formula.

## Transfer equations and numerical method

Define bolometric moments J=(1/2) integral I dmu, H=(1/2) integral mu I dmu, K=(1/2) integral mu² I dmu. E=4 pi J/c, F=4 pi H and P_rr=4 pi K/c. In the static conservative grey problem, S=J and

    d(r² H)/dr = 0
    dK/dr + (3K-J)/r = -chi H.

We compute f=K/J from angular rays. For y=f r² J and H=H0/r², the second equation becomes

    dy/dln(r) = (1/f-1)y - r chi H0.

No fitted Eddington function or flux limiter enters. A sixth-order Gauss collocation integrator visits each radial interpolation interval and checks local error by step doubling; this avoids skipping narrow tabulated closure features. Straight rays have impact parameters at every radial node plus rays intersecting the opaque core. Inward and outward formal sweeps integrate a source linear in optical path exactly. Extinction is log-linear between radial nodes, with Gaussian quadrature along each chord. Angular moments integrate piecewise-linear intensity with polynomial weights. Two one-sided values at the core limb preserve the radiation discontinuity rather than smearing the stellar disc. The moment equation and the formal angular flux are separately retained: setting H=H0/r² is **not** counted as independent luminosity validation.

For a moving isotropic coherent-scattering medium, beta=v/c, gamma_L=1/sqrt(1-beta²) and D=[gamma_L(1-mu beta)]^-1. The lab path extinction is chi0/D, the lab directional source is D^4 J0, and

    J0 = gamma_L² (J - 2 beta H + beta² K)
    H0 = gamma_L² [(1+beta²)H - beta(J+K)].

Here the subscript 0 means the comoving frame; it differs from the static reference's constant H0. These transformations and the spherical formal solution are given by [Fukue (2018)](https://doi.org/10.1093/mnras/sty358). His relativistic wind figures are not claimed reproduced. Scattering is isotropic in this family; a true polarized Thomson phase function is not implemented.

The conservative moving moments are

    d(r² H)/dr = -r² chi0 gamma_L beta H0
    dK/dr + (3K-J)/r = -chi0 gamma_L H0.

They give equal radiation energy loss and force times velocity. The VEF iteration integrates these equations inward from zero exterior illumination and normalizes the net core flux. The opaque core absorbs returning radiation and emits uniform lab intensity. Its emitted intensity is determined from the specified net luminosity; emitted luminosity and net luminosity are different. Angular transport supplies the Eddington factor and the outer flux factor. Independent integrals using the angular moments test the work budget.

`solve_scattering` also exposes a direct matrix-free scattering-source solve for comparison. A tiny algebraic residual in that method does not guarantee accurate optically thick flux: the archived coarse tau=30 calculation is a deliberate retained counterexample. The moving VEF production candidate is `solve_moving_grey`, explicitly restricted to outward 0 <= v/c <= 0.01. Linear/iteration failure is never a no-wind existence result.

## Stationary continuum control and units

`solve_continuum_wind` iterates the momentum integral with the moving radiation solution. The mass flow is a supplied parameter. The rest-mass current fixes chi0 R=m Gamma/(gamma_L beta rhat²). Newtonian acceleration uses kappa F0/c-GM/r²; mechanical relativity and gas pressure/enthalpy are excluded. With this declared Newtonian mechanical convention, the transported rest-mass current and the radiation work equation give the conserved diagnostic

    L + mdot (v²/2 - GM/r).

Corrections to the mechanical approximation are order beta²; the controls are far below the 0.01c domain limit. This is not a relativistic hydro solver. Fixed numerical damping 0.5 changes iteration speed, not the equations, and is recorded in each run.

`dimensional_fields` supplies E, F, radial radiation pressure, radiative acceleration, density, velocity, luminosity, emergent angular luminosity and a photospheric diagnostic. Radius and luminosity scales follow from a specified M, constant kappa, v_escape/c and Gamma. The 100-Msun, kappa=.34 scaling in the campaign is an illustration of dimensions, **not** a claim that this enormous-radius object remains fully ionized. Through 100,000R it also contains about 3249 central masses of wind, violating neglected wind self-gravity; that scaling cannot be interpreted as a stellar model. Pure scattering has no gas temperature, gas entropy/enthalpy or ionization solution; these fields and sonic radius are returned as `None`.

The moving photospheric convention is outgoing radial optical depth integral chi0 gamma_L(1-beta) dr = 2/3, integrated only to the explicitly stated vacuum boundary. Static solutions use integral chi dr. The effective-temperature diagnostic uses the emergent luminosity and this radius. It is neither a color temperature nor a spectrum or thermalization radius. Boundary-placement tests are necessary before regarding the finite exterior as an approximation to infinity.

## What prevents physical envelope matching

An actual primordial wind needs a material energy/population equation and appropriate absorption/scattering flux means through recombination. The existing Rosseland means and LTE gas EOS remain valuable in their own verified domains, but do not provide these missing equations in a dilute escaping layer. SIROCCO is a possible external non-LTE reference; acquiring its source is not coupling or validating it. No source/equation in this campaign justifies replacing the original recombining physics with constant opacity in a star.

The current continuum control also prescribes mdot and v_base. It cannot determine them from the envelope. Consequently Stage 5 physical matching is not attempted, Stage 7 coupled evolution is not authorized by the science gates, and Stage 8 complete-track replication remains deferred. The existing Mach-one MLT exit remains mandatory. No saturated-convection coefficient, porous opacity, or transition shock has been invented to bypass it.

For a future thin matching layer, conservation supplies the conditions rather than arbitrary equality of every thermodynamic variable. At a sphere moving at u, the Newtonian total mass, radial momentum and energy fluxes per area are, writing Phi for the common gravitational potential,

    j = rho (v-u)
    momentum = j v + Pg + P_rr - u F_lab/c²
    energy = j (ug + v²/2 + Phi) + Pg v + F_lab - u E.

These follow by subtracting the swept conserved density from each lab flux. A source-free layer equates its two fluxes; a layer with stored mass/energy, heating, stresses or shocks needs the corresponding reservoir and jump terms. Entropy need not be equal across a dissipative layer. Gas temperature and radiation energy need not obey E=a Tgas^4 outside thermalization. One must not add radiation enthalpy to F_lab again: it already contains the advected radiation flux. In the isotropic diffusion limit, expanding F_lab=F_diff+4 Pr v and E=3 Pr recovers Report 10's relative enthalpy plus moving pressure-work formula. For a finite spherical layer, curvature, gravity and its inventory must also be integrated. These identities describe a future interface; they do not constitute a solved physical match.

## Reproduction

    python scripts/register_outer_candidates.py
    python scripts/build_registry.py
    python scripts/validate_escaping_radiation.py
    python -m pytest -q tests/test_grey_transfer.py

The saved tolerances precede the formal campaign. Driver, inputs, source snapshot, iteration histories, failed attempts, raw radial arrays and comparisons accompany each run. See Report 11 for measured results and the release gate matrix.
