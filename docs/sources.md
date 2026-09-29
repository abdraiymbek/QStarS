# Source and prescription registry

This register distinguishes implemented laws, derived numerical identities, optional approximations, and unimplemented published models. Primary PDFs are locally cached, excluded from redistribution, and hashed in `sources/acquisition.json`.

## bra2008 — Begelman, Rossi & Armitage (2008)

[Begelman, Rossi & Armitage (2008)](https://arxiv.org/abs/0711.4078) · Sections 2.1-2.4; Eqs 1-14

- Status: `implemented_luminosity_law`
- Assumptions: Radiation dominated convective envelope; low BH/envelope mass; unresolved central accretion flow
- Domain: Hot optically thick envelope; do not apply the analytic photosphere scalings as an EOS
- Variables and units: L erg/s, M g, rho g/cm^3, P dyn/cm^2; alpha dimensionless
- Interpretation: Eq 3 is a selectable luminosity law. Fixed excision geometry is an additional QStarS approximation, not BRA equilibrium replication.
- Implementation: `src/qstars/engine.py`

## ball2011 — Ball et al. (2011)

[Ball et al. (2011)](https://arxiv.org/abs/1102.5098) · Sections 2.1-2.3, 4.1; Eqs 2,5,9-15,24-33; Table 1

- Status: `implemented_reference_and_printed_law`
- Assumptions: BH-only Bondi radius; rho~r^-3/2 cavity; loaded polytrope
- Domain: Reference n=2,3,4 polytropes; complete stellar tracks require different microphysics
- Variables and units: xi and phi dimensionless; accretion g/s; radius cm
- Interpretation: Eq 15 prints 4pi despite Eq 12 plus rB=2GM/cs^2 implying 16pi. Preserve printed mode; do not blend with QUEST. Table 1 stored as unreplicated target.
- Implementation: `src/qstars/reference.py; src/qstars/engine.py`

## dotan2011 — Dotan, Rossi & Shaviv (2011)

[Dotan, Rossi & Shaviv (2011)](https://arxiv.org/abs/1107.3562) · Sections 3.2-3.4; Eqs 19-23,31-32

- Status: `implemented_optional_launch_closure`
- Assumptions: Porous radiation-driven sonic boundary; empirical wind function W; photon-tiring cap
- Domain: Only effective launch prescription, not a solved porous atmosphere/wind
- Variables and units: mdot g/s, L erg/s, gas sound speed cm/s, W dimensionless
- Interpretation: Use gas isothermal sound speed; gamma<=1 gives zero. Wind enthalpy and energy to escape are separately charged.
- Implementation: `src/qstars/winds.py`

## ball_thesis — Ball (2012) thesis

[Ball (2012) thesis](https://arxiv.org/abs/1207.5972) · Chapters 2-5; Bondi and CDAF/ADAF boundaries

- Status: `audited_not_full_evolution`
- Assumptions: Boundary choice changes the evolutionary sequence and its mass limit
- Domain: No single universal quasi-star BH mass limit
- Variables and units: M g, r cm, thermodynamic quantities CGS
- Interpretation: Inspected boundary discussion and contents; thesis has 145 pages and was not audited line by line. Alternatives must remain explicit.
- Implementation: `docs/literature.md`

## fiacconi2016 — Fiacconi & Rossi (2016)

[Fiacconi & Rossi (2016)](https://arxiv.org/abs/1509.07511) · Sections 2.1-2.2; Eqs 12-22,25; Figures 2 and 4

- Status: `implemented_stationary_wind_reference`
- Assumptions: Stationary spherical radiation-dominated wind; radiation enthalpy advection
- Domain: Gas pressure negligible and wind base optically thick
- Variables and units: mdot g/s, L erg/s, specific Bernoulli energy erg/g
- Interpretation: Inward integration retains radiative, advected, kinetic and gravitational energy. Base Mach is explicitly 0.8 +/- 0.05, as in the paper; no gas sonic critical-point selection is claimed. Figure 4 vector paths retain calibration uncertainty and discrepancies; evolving-envelope coupling remains open.
- Implementation: `src/qstars/wind_reference.py; scripts/validate_wind.py; benchmarks/fiacconi2016_figure4.csv`

## coughlin2024 — Coughlin & Begelman (2024)

[Coughlin & Begelman (2024)](https://arxiv.org/abs/2405.00084) · Sections 2-3; Eqs 2-7,11-21; Figs 1,3,4

- Status: `implemented_reference_odes`
- Assumptions: Radiation-dominated saturated inner convection matched to n=3 outer envelope
- Domain: Dimensionless static structure; thin radiative layer approximation
- Variables and units: xi, m, Ki dimensionless
- Interpretation: Direct ODE integration distinguishes corrected Bondi and saturated families. Printed Ki=.46554 and mBH/Mi=.6321 correspond to x_tilde=1, not the stated entropy minimum near 1.015. Retain both and document discrepancy.
- Implementation: `src/qstars/reference.py`

## quest2025 — Campbell, Santarelli & Caplan (2025)

[Campbell, Santarelli & Caplan (2025)](https://arxiv.org/abs/2507.12651) · Eqs 1-5; Fig 1; public code

- Status: `audited`
- Assumptions: MESA EOS/opacity, ML1 convection, cavity, TOV correction
- Domain: Full published run depends on MESA r24.03.1
- Variables and units: CGS; epsilon, eta dimensionless
- Interpretation: Early note Eq 3 lacks (1-epsilon) while using a BH growth convention. Later article and source code have it. No silent normalization correction.
- Implementation: `docs/literature.md`

## santarelli2025 — Santarelli et al. (2025/2026)

[Santarelli et al. (2025/2026)](https://arxiv.org/abs/2510.11772) · Sections V-VI; Eqs 8,11-13; upstream implementation

- Status: `implemented_eq11_law_only`
- Assumptions: Convection-limited accretion with explicit radiative efficiency; MESA full tracks
- Domain: Luminosity law does not reproduce its moving cavity or saturated-convection boundary
- Variables and units: mdot_BH=(1-epsilon)mdot_rest g/s
- Interpretation: Named quest2025_eq11 uses 16pi and distinct rest/BH rates; upstream HEAD inspected at 671eb9b898b28c0d02a6a2455dacc0f5b269e280.
- Implementation: `src/qstars/engine.py`

## hassan2025 — Hassan et al. (2025/2026)

[Hassan et al. (2025/2026)](https://arxiv.org/abs/2510.18301) · Section 2.1, Eqs 5-6; Appendix A cubic

- Status: `implemented_reference_cubic`
- Assumptions: Self-consistent Bondi radius including rho~r^-3/2 cavity
- Domain: Below critical mass, choose small positive radius continuous as MBH tends to zero
- Variables and units: rho g/cm^3, cs cm/s, N dimensionless
- Interpretation: Solve Appendix A cubic directly for stability near zero BH mass; verified equivalent to printed Eq 6. An initial factor-of-two audit suspicion was retracted (Report 4). Full ~0.33 evolutionary limit remains unreplicated.
- Implementation: `src/qstars/reference.py`

## cantiello2025 — Cantiello et al. (2025/2026)

[Cantiello et al. (2025/2026)](https://arxiv.org/abs/2512.17997) · Section 2.2, nonadiabatic GYRE analysis

- Status: `audited_not_implemented`
- Assumptions: GYRE 8.1 nonadiabatic radial modes with fixed inner displacement and MESA microphysics
- Domain: Linear stability results cannot be inferred from nonlinear solver failures
- Variables and units: Eigenfrequency complex; periods time; growth rates inverse time
- Interpretation: QStarS Mach or solver stops are domain/numerical exits, not a pulsational-instability diagnosis.
- Implementation: `src/qstars/validity.py`

## roman2026 — Roman-Garza et al. (2026)

[Roman-Garza et al. (2026)](https://arxiv.org/abs/2603.21714) · Section 2; Appendix A; Table B.1

- Status: `audited_not_replicated`
- Assumptions: Cold net mass gain; MESA 25.10.1; nuclear evolution and GR progenitor; late extrapolation
- Domain: Net-gain sequences do not establish separate feeding/wind thermodynamics
- Variables and units: Mass solar masses, time yr, luminosity solar luminosities in tables
- Interpretation: Do not equate prescribed net gain with a resolved simultaneous two-stream flow; do not mistake analytically extrapolated late evolution for a completed stellar run.
- Implementation: `docs/literature.md`

## begelman2026 — Begelman & Dexter (2026)

[Begelman & Dexter (2026)](https://arxiv.org/abs/2507.09085) · Late-stage quasi-star structural and spectral discussion

- Status: `audited_context`
- Assumptions: Saturated envelope assumptions and late BH-dominated evolution
- Domain: Qualitative context; not a QStarS spectral prediction
- Variables and units: CGS, solar mass where stated
- Interpretation: Spectral interpretation requires atmosphere microphysics absent from this solver.
- Implementation: `docs/literature.md`

## atmosphere2026 — Martins, Schaerer & Marques-Chaves (2026)

[Martins, Schaerer & Marques-Chaves (2026)](https://arxiv.org/abs/2605.12141) · Abstract and atmosphere-model assumptions

- Status: `audited_context`
- Assumptions: Expanding spherical non-LTE CMFGEN atmospheres
- Domain: 5000-12000 K atmosphere models; not interior LTE closure
- Variables and units: T K, L solar luminosities
- Interpretation: Grey flux temperature is not a synthetic spectrum or a validated wind photosphere.
- Implementation: `docs/literature.md`

## mesa_energy — MESA technical documentation, energy equation

[MESA technical documentation, energy equation](https://docs.mesastar.org/en/latest/reference/controls.html#energy-eqn-option) · dedt and eps_grav descriptions; MESA V section 3 Eq 8 cited there

- Status: `audited_numerical_reference`
- Assumptions: Total-energy conservation and thermodynamic consistency are distinct from resolution convergence
- Domain: Stellar evolution numerical practice
- Variables and units: e erg/g, P dyn/cm^2, V_rho cm^3/g
- Interpretation: Use independent energy residuals plus mesh/time convergence. The implementation derives its own spherical discrete work identity, not a transcription of MESA.
- Implementation: `src/qstars/numerics.py; docs/equations.md`

## mesa_install — MESA installation and upstream QUEST

[MESA installation and upstream QUEST](https://docs.mesastar.org/en/latest/installation.html) · Prerequisites, SDK, release pinning

- Status: `audited_backend_decision`
- Assumptions: Validated microphysics requires versioned EOS/opacity/atmosphere infrastructure
- Domain: Existing MESA 24.08.1 and ARM SDK 24.7.1 recovered from earlier runtime; local build changes recorded
- Variables and units: Not applicable
- Interpretation: Direct EOS and opacity calls are implemented for fixed H/He composition. Stellar evolution, atmosphere, and moving cavity production gates remain open. Report 5 corrects the earlier absence claim.
- Implementation: `docs/mesa.md`

## ideal_eos — Classical ideal-gas plus radiation thermodynamics

[Classical ideal-gas plus radiation thermodynamics](https://arxiv.org/abs/1107.3562) · Dotan Eqs 21-23; first-law derivation in docs/equations.md

- Status: `derived_and_tested`
- Assumptions: Nondegenerate fully ionized monatomic mixture in LTE; fixed composition
- Domain: 1e5 <= T <= 1e8 K, electron occupation < 0.1, Newtonian low compactness
- Variables and units: P dyn/cm^2, u erg/g, s erg/g/K
- Interpretation: Heat capacities, adiabatic exponents, and entropy derived consistently from the same EOS; arbitrary entropy zero explicitly fixed.
- Implementation: `src/qstars/eos.py`

## thomson — Electron-scattering opacity

[Electron-scattering opacity](https://arxiv.org/abs/0711.4078) · Section 2.3; Thomson cross section times electrons per gram

- Status: `derived_approximation`
- Assumptions: Fully ionized H/He mixture; scattering dominates; LTE/thermalization separately assumed
- Domain: Hot diffuse envelope; no bound-free/free-free/negative hydrogen/line opacity
- Variables and units: kappa=0.2(1+X) cm^2/g
- Interpretation: No cool extrapolation admitted as validated physics. Optional external opacity tables strictly reject out-of-bounds values and require provenance.
- Implementation: `src/qstars/opacity.py`

## diffusion — Rosseland radiative diffusion

[Rosseland radiative diffusion](https://docs.mesastar.org/en/latest/reference/controls.html) · Temperature-gradient controls; diffusion derivation in docs/equations.md

- Status: `derived`
- Assumptions: Optically thick LTE; F=-ac/(3 kappa rho) grad(T^4)
- Domain: Radiation gradient Knudsen number <=0.1; LTE assumed; no resolved optically thin transfer
- Variables and units: Flux erg/cm^2/s; luminosity erg/s
- Interpretation: Backward-Euler diffusion shares exactly one flux per face; grey extrapolation is an optional outer closure.
- Implementation: `src/qstars/transport.py`

## mlt_lossless — Optically thick lossless mixing-length limit

[Optically thick lossless mixing-length limit](https://arxiv.org/abs/2507.12651) · ML1 reference in Section II; explicit derivation in docs/equations.md

- Status: `optional_approximation`
- Assumptions: Local buoyant blobs, no radiative heat leakage; l=alpha Hp
- Domain: Pe>=100 where convection carries more than 1 percent of transported flux; subsonic blobs; not full ML1 or time-dependent convection
- Variables and units: cp erg/g/K; velocity cm/s; l cm
- Interpretation: Use v^2=g delta l^2 (nabla-nabla_ad)/(8Hp), F=rho cp T v l (nabla-nabla_ad)/(2Hp), capped by P cs. Entropy contrast discretized consistently.
- Implementation: `src/qstars/transport.py`

## compatible_work — Compatible spherical mechanical work

[Compatible spherical mechanical work](https://cnls.lanl.gov/~shashkov/papers/comp.pdf) · Compatible differencing principle; QStarS derivation in docs/equations.md

- Status: `derived`
- Assumptions: Fixed shell masses during mechanical substep; midpoint pressure and velocity; spherical gravity quadrature
- Domain: Smooth subsonic flow; no shock validation
- Variables and units: Total energy erg; nodal velocity cm/s
- Interpretation: Secant area and 1/(r_new r_old) gravity produce exact work cancellation; no energy residual compensation. Article discovered as methodological background; full PDF not audited.
- Implementation: `src/qstars/numerics.py`

## two_stream — Conservative boundary source-cell closure

[Conservative boundary source-cell closure](https://arxiv.org/abs/1509.07511) · First-law open-system derivation; compare wind Eq 3

- Status: `derived_optional_approximation`
- Assumptions: Separate mass sources/sinks at first/last cell; comoving mixing; supplied enthalpy and entropy specified
- Domain: Effective 1D boundary mixing; unresolved geometry, shocks, inflow momentum
- Variables and units: Mass g, enthalpy erg/g, entropy erg/g/K
- Interpretation: Differentiate discrete gravitational energy for boundary potential flux; central-potential work is separately recorded. Wind escape cost debits outward light once.
- Implementation: `src/qstars/mass_transport.py`

## remapping — Conservative overlap remapping

[Conservative overlap remapping](https://cnls.lanl.gov/~shashkov/papers/comp.pdf) · QStarS finite-volume derivation in docs/equations.md

- Status: `derived`
- Assumptions: Piecewise constant volume and total specific energy in mass coordinates
- Domain: Smooth resolved states; remap changes entropy and hydrostatic residual
- Variables and units: Mass g, volume cm^3, energy erg
- Interpretation: New internal energy is remapped total energy minus new discrete gravity and kinetic energies. Decimal overlap endpoints avoid cancellation of large cumulative CGS masses; the sweep is linear in zone count. Reject negative internal energy; log entropy changes.
- Implementation: `src/qstars/mesh.py`

## mesa_microphysics — MESA 24.08.1 EOS and opacity public APIs

[MESA 24.08.1 EOS and opacity public APIs](https://docs.mesastar.org/en/24.08.1/eos/overview.html) · eosDT_get, kap_get; eos_def quantity definitions; source commit 5c166bf174b477823a2ca0031c4f57210c8e67a8

- Status: `implemented_direct_library_interface`
- Assumptions: Fixed H1/He4 composition Z=0; default EOS blends and Type2 opacity with Zbase=0; matched electron state
- Domain: Query rectangle rho=1e-14..1e3 g/cm3 and T=3000..1e8 K is only a bounding box. Reject neutral ideal backstop and first-law defect above 1 percent; no production atmosphere claim.
- Variables and units: P erg/cm3, u erg/g, s erg/g/K, kappa cm2/g; derivatives with respect to natural logarithms
- Interpretation: Versioned C ABI and Python interface cross-checked against separate Fortran process. Wide audit preserves failing blend states. Radiative opacity disables conduction; a separate combined-opacity call measures omitted conduction. Existing binary/build/table provenance must accompany each run.
- Implementation: `backends/mesa/microphysics.f90; src/qstars/mesa_microphysics.py`

## inverse_transport — Inverse radiative/convective flux solve

[Inverse radiative/convective flux solve](https://docs.mesastar.org/en/24.08.1/) · Algebraic derivation in docs/equations.md, using existing diffusion and lossless MLT equations

- Status: `derived_and_tested`
- Assumptions: Same physical closure as forward flux; luminosities retained as nonlinear unknowns; active-piece complex-step derivatives
- Domain: Smooth subsonic hot envelopes; transport validity unchanged by nonlinear method
- Variables and units: L erg/s, logarithmic thermal contrast dimensionless
- Interpretation: Solving thermal contrast from luminosity avoids ill conditioning near adiabatic stratification; no transport coefficient changed. Numerical reference is QStarS derivation, not copied MESA implementation.
- Implementation: `src/qstars/transport.py; src/qstars/nonlinear.py`

## diffusion_domain — Levermore & Pomraning (1981)

[Levermore & Pomraning (1981)](https://adsabs.harvard.edu/pdf/1981ApJ...248..321L) · Section II, Eqs 7 and 18-31; pp321-324

- Status: `implemented_validity_screen`
- Assumptions: Near LTE and slowly varying radiation energy over a mean free path; no flux-limiter implementation added
- Domain: R=|grad E|/(kappa rho E)<=0.1, using E=aT^4. Cell optical depth remains a diagnostic, not a mesh-dependent validity restriction.
- Variables and units: R dimensionless, photon mean free path cm
- Interpretation: Corrects prior tau_cell>=1 rejection, which incorrectly rejected increasingly fine representations of a smooth diffusive field. R<=0.1 is a conservative screening threshold, not a proof of full atmosphere validity.
- Implementation: `src/qstars/validity.py`

## backward_euler_work — Conservative backward-Euler mechanical work

[Conservative backward-Euler mechanical work](https://docs.mesastar.org/en/24.08.1/) · Independent derivation in docs/equations.md; not copied MESA code

- Status: `derived_and_tested`
- Assumptions: Implicit velocity and pressure with secant geometry/gravity; positive numerical kinetic heating retained as internal conversion
- Domain: Smooth secular evolution; first order, requires timestep refinement; not a pulsation solver
- Variables and units: Q_i=dm_i[(Delta v_i)^2+(Delta v_i+1)^2]/4 erg
- Interpretation: Diagnostic heating is derived from momentum work and never from a global energy residual. It is excluded from the external ledger balance. Both midpoint and backward_euler are explicit selectable numerical methods.
- Implementation: `src/qstars/numerics.py; src/qstars/ledger.py`

## mesa_atmosphere_mlt — MESA 24.08.1 atmosphere and native MLT

[MESA 24.08.1 atmosphere and native MLT](https://docs.mesastar.org/en/24.08.1/reference/controls.html) · atm_T_tau_varying.f90; atm_T_tau_uniform.f90; turb set_mlt and mlt.f90

- Status: `implemented_reference_not_evolution_boundary`
- Assumptions: Eddington T(tau); native fixed/varying opacity; Henyey or ML1 with alpha=2 explicit, Henyey nu=8 and y=1/3
- Domain: Native varying atmosphere starts with rho_outer=1e-10; inappropriate as an implicit vacuum boundary at low gravity. MLT subsonic validity checked separately.
- Variables and units: tau dimensionless, Pg dyn/cm2, v cm/s; native automatic derivatives
- Interpretation: At Teff=8000 K, g=3 and kappa=.34, native varying constant-opacity photospheric Pg is 78.96 times the analytic vacuum value. Native MLT availability is not a validated coupled stellar model.
- Implementation: `backends/mesa/microphysics.f90; src/qstars/mesa_microphysics.py`

## freeeos3 — Irwin FreeEOS 3.0.0 (2022)

[Irwin FreeEOS 3.0.0 (2022)](https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/) · mod_free_eos modern API; EOS1 (3,11,-2); mod_free_eos_constants and mod_isotopic_mass_data

- Status: `implemented_audited_reference`
- Assumptions: Free-energy minimization for LTE H1/He4, gas-only evaluation followed by explicit blackbody radiation; no EOS fitting or stitching
- Domain: 153-state component audit is not global EOS certification or non-LTE atmosphere validation. Upstream byte comparisons fail; physical reference fields quantified separately.
- Variables and units: CGS; derivatives in ln rho and ln T; energy zero from upstream
- Interpretation: Source tar SHA256 4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09. ARM portability patch removes duplicate extended/quadruple generic declarations after compiler-kind check; physical formulae unchanged. GPL dependency separately installed.
- Implementation: `backends/freeeos/interface.f90; src/qstars/freeeos_reference.py; scripts/build_freeeos.py`

## grey_mlt_atmosphere_reference — Finite-pressure grey and MLT outer-layer reference

[Finite-pressure grey and MLT outer-layer reference](https://adsabs.harvard.edu/pdf/1965ApJ...142..841H) · Eddington grey relation and hydrostatic first-law derivation in docs/atmosphere.md; native Henyey implementation

- Status: `implemented_candidate_fails_general_use`
- Assumptions: Plane-parallel hydrostatic layer, grey radiative atmosphere to tau=2/3, diffusion plus native MLT below; external Pg explicit
- Domain: No spherical transfer, non-LTE populations, moving mass, wind coupling or physical disruption claim. Geometric thickness and supersonic convection veto promotion.
- Variables and units: dP/dtau=g/kappa, dT/dtau=(T/P)*nabla*g/kappa; column mass g/cm2, energy erg/cm2
- Interpretation: Integrate explicit gas pressure and subtract radiation-pressure gradient; retain depth, mass and energy inventories. Constant opacity has analytic check. Report 10 stops at the first resolved convective Mach-one crossing with a structured validity-domain exit; no supersonic profile is returned. The archived Report 8 continuation remains an invalid diagnostic, not a closure.
- Implementation: `src/qstars/atmosphere_reference.py; scripts/validate_atmosphere_reference.py`

## dilute_gas_screen — Dilute classical-gas limit for H/He

[Dilute classical-gas limit for H/He](https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/) · Ideal-gas limit Pg=rho kT/(mu m_u); Debye-Huckel weak-coupling and nondegeneracy screening

- Status: `implemented_validity_screen`
- Assumptions: H/He only, T>=40000 K, eta<-10 and conservative charge-2 Coulomb Gamma<.01
- Domain: A 1% dilute-gas discrepancy is a rejection threshold, not a correction or global EOS accuracy claim
- Variables and units: Gamma and relative Pg defect dimensionless
- Interpretation: Radiation-dominated total first-law agreement can conceal significant gas-pressure error. The historical MESA table discrepancy is retained and rejected in evolution when the dilute limit applies; its radiation-scaled origin is an inference, not a proven table-generation bug.
- Implementation: `src/qstars/mesa_microphysics.py; src/qstars/validity.py`

## nakauchi2017 — Nakauchi et al. (2017)

[Nakauchi et al. (2017)](https://arxiv.org/abs/1611.09601) · Section 2 equations 1-9, Section 3 Figure 2

- Status: `implemented_diffusion_wind_component`
- Assumptions: Smooth stationary spherical LTE gas plus radiation, negligible wind self-gravity, no wind convection; direct FreeEOS and MESA opacity differ from original microphysics
- Domain: Stop strict calculations at radiation-gradient Knudsen 0.1 or a second sonic approach. Paper photosphere extension is diagnostic only. No terminal escape or coupled evolution certification.
- Variables and units: r cm, rho g/cm3, T K, v cm/s, mdot g/s, L erg/s
- Interpretation: Gas sonic regularity fixes critical luminosity and local accelerating slope. Base matching solves for sonic r,rho,T without supplied mdot. Independent momentum, diffusion and first-law quadratures supplement exact flux algebra. Published profile discrepancies retained without fitting.
- Implementation: `src/qstars/spherical_wind.py; scripts/validate_spherical_wind.py`

## quataert2016 — Quataert et al. (2016)

[Quataert et al. (2016)](https://arxiv.org/abs/1509.06370) · Sections 2-3; steady equations and trapped-radiation wind regime

- Status: `audited_method_alternative`
- Assumptions: Near-surface energy deposition; nearly adiabatic trapped radiation; convection not automatically the dominant transport
- Domain: Not a generic recombining photosphere or microscopic opacity replacement
- Variables and units: Heating erg/s, mass g, speed cm/s
- Interpretation: Provides a physically motivated outflow alternative when hydrostatic convection cannot carry energy; adiabatic approximation needs a photon trapping check. No new heating law added.
- Implementation: `docs/spherical_wind.md`

## owocki2017 — Owocki, Townsend & Quataert (2017)

[Owocki, Townsend & Quataert (2017)](https://arxiv.org/abs/1708.07790) · Equations 5-10; Appendices B-C

- Status: `audited_energy_accounting`
- Assumptions: Steady spherical continuum winds; gas and radiation enthalpy separated from comoving diffusion luminosity
- Domain: Constant opacity limits and simplified sonic sound speed do not solve recombination transfer
- Variables and units: Luminosity erg/s, enthalpy erg/g
- Interpretation: Confirms energy partition and the distinction between observer luminosity and diffusive luminosity. Do not charge wind energy from emergent radiation alone.
- Implementation: `docs/spherical_wind.md; src/qstars/spherical_wind.py`

## owocki2004 — Owocki, Gayley & Shaviv (2004)

[Owocki, Gayley & Shaviv (2004)](https://arxiv.org/abs/astro-ph/0409573) · Porosity-length and power-law porosity models

- Status: `audited_not_implemented`
- Assumptions: Effective opacity from unresolved inhomogeneities; porosity length and distribution are physical closure inputs
- Domain: No unique first-principles porosity parameter for these quasi-star boundaries
- Variables and units: Porosity length cm, optical depths dimensionless
- Interpretation: A defensible alternative family with explicit parameters, not permission to fit an opacity reduction to a failing smooth-wind benchmark. Kept separate from Dotan launch and Fiacconi radiation wind.
- Implementation: `docs/spherical_wind.md`

## jiang2015 — Jiang et al. (2015)

[Jiang et al. (2015)](https://arxiv.org/abs/1509.05417) · 3D radiation hydrodynamics; optically efficient and inefficient convection regimes

- Status: `audited_validity_context`
- Assumptions: Massive-star radiation-dominated iron-opacity layers, multidimensional density and flux correlations
- Domain: Different composition and opacity peak from primordial quasi-stars; cannot transplant a calibrated coefficient without testing
- Variables and units: Optical depth and convective Mach dimensionless
- Interpretation: Supports checking diffusion time versus dynamical time and forbids assuming supersonic local MLT remains valid. Does not prove a universal no-go for 1D quasi-stars.
- Implementation: `docs/spherical_wind.md`

## nakauchi_opacity_audit — MESA historical OPAL / Alexander-Ferguson opacity audit

[MESA historical OPAL / Alexander-Ferguson opacity audit](https://docs.mesastar.org/en/24.08.1/kap/defaults.html) · kap_get_radiative_opacity; kap controls; Nakauchi section 2.1

- Status: `implemented_diagnostic_only`
- Assumptions: Type1 OPAL gs98 or gn93 above 7000 K; MESA lowT_af94_gn93 below; X=.7, Z=0; native interpolation
- Domain: Common table rectangle logR=-7..1, T=3000..1e7 K; Rosseland mean only, no flux-mean or non-LTE inference
- Variables and units: kappa cm2/g; R=rho/(T/1e6)^3
- Interpretation: Separate native handles preserve the default opacity. MESA lowT_af94 header says Alexander & Ferguson (1995); processed tables cannot be certified identical to the authors files. No fitting of sonic inputs.
- Implementation: `backends/mesa/microphysics.f90; scripts/validate_outer_discrepancy.py`

## moving_wind_flux — Reynolds transport for an LTE wind matching sphere

[Reynolds transport for an LTE wind matching sphere](https://arxiv.org/abs/1708.07790) · QStarS derivation from Reynolds transport; stationary energy channels consistent with Owocki et al. equations 5-10

- Status: `derived_and_tested_component`
- Assumptions: Newtonian LTE isotropic radiation pressure, static central potential, lab conserved variables, specified face speed
- Domain: No automatic moving-cavity or evolving atmosphere reconstruction; conservation interface only
- Variables and units: mass g/s, momentum dyn, energy erg/s, face speed cm/s
- Interpretation: Relative enthalpy advection plus P*A*face_speed work equals fixed flux minus swept conserved density. Tested against independently integrated swept volumes and an energy-gauge shift.
- Implementation: `src/qstars/spherical_wind.py; tests/test_wind_moving_boundary.py; scripts/validate_outer_transport.py`

## hummer1971 — Hummer & Rybicki (1971)

[Hummer & Rybicki (1971)](https://doi.org/10.1093/mnras/152.1.1) · Equations 2.1-2.31, Table I final iteration; PDF SHA256 05fdb069368c1b2885358203a6b84dc85eb41eb6e0e4c6026ba17b16f936be20

- Status: `implemented_angular_transfer_reference`
- Assumptions: Static spherical conservative isotropic grey transfer with computed variable Eddington factor
- Domain: Declared grey extinction only; no non-grey H/He flux mean, thermalization or non-LTE inference
- Variables and units: J,H,K bolometric intensities; chi inverse length; net central luminosity
- Interpretation: Independent published moment table, angular flux, diffusion and streaming limits plus boundary/grid refinement. Finite inner diffusion boundary is checked, never silently treated as the exact infinite-depth condition.
- Implementation: `src/qstars/grey_transfer.py; scripts/validate_escaping_radiation.py`

## fukue2018_transfer — Fukue (2018)

[Fukue (2018)](https://doi.org/10.1093/mnras/sty358) · Equations 3-14 and grey scattering moment equations; public indexed primary article

- Status: `implemented_restricted_moving_transfer`
- Assumptions: Isotropic coherent scattering in the comoving frame; frequency-independent specified extinction, no true absorption or material energy equation
- Domain: Outward v/c <= .01 in tested moment solver; no ionization or gas temperature prediction
- Variables and units: D=[gamma(1-mu beta)]^-1; lab source D^4 J0; extinction chi0/D
- Interpretation: Conservative moving moments closed by angular rays; incoming core light absorbed and isotropic emitted intensity determined by net luminosity. Published relativistic figure families are not claimed reproduced.
- Implementation: `src/qstars/grey_transfer.py; docs/escaping_radiation.md`

## owocki2017_analytic — Owocki, Townsend & Quataert (2017) flux-driven limit

[Owocki, Townsend & Quataert (2017) flux-driven limit](https://arxiv.org/abs/1708.07790) · Equation 24 and 26; finite initial kinetic energy generalization independently derived

- Status: `implemented_constructed_stationary_control`
- Assumptions: Constant grey opacity, no gas pressure/enthalpy; analytic reference neglects radiation drag, numerical wind includes moving transfer
- Domain: Mathematical continuum control, supplied mdot and vbase. Not a thermodynamic primordial quasi-star wind or an envelope match.
- Variables and units: w=v^2/vesc^2, x=1-R/r, m=mdot GM/(R L0), Gamma=L0/Ledd
- Interpretation: No coefficient fits. Compare velocity, density, luminosity, force and photospheric diagnostic, energy and momentum; high-drag negative control must not be interpreted through the small-drag analytic approximation.
- Implementation: `src/qstars/continuum_wind.py; benchmarks/grey_continuum_wind.json`

## outer_public_code_audit — Public stationary-wind and transfer candidate audit

[Public stationary-wind and transfer candidate audit](https://github.com/sirocco-rt/sirocco) · Pinned SIROCCO, MPI-AMRVAC and DUSTY source paths, model inputs and hashes in candidate_source_files.json

- Status: `audited_not_coupled`
- Assumptions: Availability ranked separately from physical correspondence; prescribed flow, line driving and dust driving remain different families
- Domain: No installed/coupled non-LTE or dust/line wind claim
- Variables and units: Each candidate records equations, EOS, opacity, radiation, sonic and boundary conditions, composition, M/R/L and available outputs
- Interpretation: Constructed analytic continuum control chosen for exact reproducibility; public transfer packages remain candidates for closing missing gas/radiation microphysics, not substitutes for a solved launch.
- Implementation: `benchmarks/outer_candidates.json; scripts/register_outer_candidates.py`

## sander2017 — Sander et al. (2017) hydrodynamically consistent PoWR

[Sander et al. (2017) hydrodynamically consistent PoWR](https://arxiv.org/abs/1704.08698) · Sections 2-5, critical condition and mass-loss update; PDF hash in Report 12 acquisition

- Status: `audited_not_replicated`
- Assumptions: Coupled non-LTE populations, comoving-frame transfer and hydrodynamic iteration; atmospheric inner optical depth is a matching constraint
- Domain: O/WR-star examples have metals, line driving and explicit turbulence/clumping choices; not primordial quasi-star calibration
- Variables and units: Mass loss g/s, opacity per mass, material temperature K, radiative acceleration cm/s2
- Interpretation: A single critical-point condition is insufficient without consistent inner radiation/thermodynamic constraints and the force response to density. No published PoWR quasi-star run claimed.
- Implementation: `docs/self_selected_outer.md`

## mayer2005 — Mayer & Duschl (2005) primordial mean opacities

[Mayer & Duschl (2005) primordial mean opacities](https://arxiv.org/abs/astro-ph/0411613) · Sections on Rosseland, Planck continuum and line means; conclusion; cached paper

- Status: `audited_not_installed`
- Assumptions: Primordial H/He plus trace D/Li, equilibrium populations; Planck continuum and line sets differ
- Domain: Single-temperature means do not determine the absorption energy mean or flux mean for arbitrary escaping spectra
- Variables and units: rho g/cm3, temperature K, opacity cm2/g
- Interpretation: Paper states full tables available from authors on request. No complete table or non-LTE atom installed by this campaign; trace-Li convention differs from current X=.7 Y=.3 mixture.
- Implementation: `docs/self_selected_outer.md`

## malygin2014 — Malygin et al. (2014) two-temperature Planck means

[Malygin et al. (2014) two-temperature Planck means](https://arxiv.org/abs/1408.3377) · Sections 2-4, equilibrium-temperature degeneracy; cached PDF

- Status: `audited_frequency_requirement`
- Assumptions: LTE chemical populations and high-resolution opacity spectrum; gas and illuminating radiation temperatures distinct
- Domain: Published metallicities -0.3, 0,+0.3 dex are not primordial; multiple thermal roots possible when opacity depends on state
- Variables and units: Tgas, Trad K; Planck/Rosseland cm2/g
- Interpretation: Demonstrates why a Rosseland-only backend is insufficient and why a unique frozen-coefficient thermal root is not evidence for unique real-gas equilibria.
- Implementation: `src/qstars/thermal_transfer.py; docs/self_selected_outer.md`

## thermal_bands_control — Spherical frequency-band thermal transfer verification

[Spherical frequency-band thermal transfer verification](https://arxiv.org/abs/1408.3377) · Standard static transfer and Kirchhoff source; exact chord and thin power-law limits derived in QStarS

- Status: `derived_and_tested_component`
- Assumptions: Specified constant opacity within each band, isotropic coherent scattering, static material energy equilibrium
- Domain: No state-dependent H/He opacity, statistical populations, moving frequency shifts or mechanical energy coupling
- Variables and units: J,H,K band integrated intensities; photon-temperature edges h nu/k in K
- Interpretation: Temperature is solved from absorption/emission balance; Planck, radiation-energy, flux and Rosseland means remain distinct. Independent luminosity/heating integral is never renormalized.
- Implementation: `src/qstars/thermal_transfer.py; scripts/validate_self_selected_outer.py`

## material_energy_critical — Material energy equation and adiabatic Parker eigenvalue control

[Material energy equation and adiabatic Parker eigenvalue control](https://arxiv.org/abs/1207.7348) · Parker polytropic family; independent local first-law elimination in docs/self_selected_outer.md

- Status: `derived_and_tested_component`
- Assumptions: Finite local thermal source without explicit velocity-gradient dependence; caloric ideal gamma between 1 and 1.5 for exact reference; fixed point gravity
- Domain: Reference is not H/He nor radiation escape. Isothermal/diffusive prescribed-gradient critical equations are distinct limits
- Variables and units: mdot g/s selected from base rho,P,M,r; local Q erg/cm3/s
- Interpretation: Gas-only FreeEOS energy derivatives yield the native gas adiabatic sound speed. Reference root interval and hydrostatic alternative have an explicit scoped proof.
- Implementation: `src/qstars/material_wind.py; tests/test_material_wind.py`

## general_radiation_interface — General lab radiation moments at a moving matching face

[General lab radiation moments at a moving matching face](https://arxiv.org/abs/1708.07790) · Reynolds transport derivation in docs/self_selected_outer.md; not an LTE closure

- Status: `derived_and_tested_component`
- Assumptions: Newtonian material and realizable lab radiation moments, same static potential and gas energy gauge
- Domain: Zero-thickness flux jumps only; finite layers require storage, gravity and curvature; no entropy continuity imposed
- Variables and units: Mass g/s, momentum dyn, energy erg/s; E erg/cm3, F erg/cm2/s, Prr dyn/cm2
- Interpretation: Subtract swept gas and radiation energy/momentum. Radiation enthalpy is already in lab flux. Independent normal-shock and energy-gauge controls.
- Implementation: `src/qstars/outer_interface.py; tests/test_outer_thermal.py`

