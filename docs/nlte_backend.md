# Non-LTE backend selection and coupling contract — Report 13

**CMFGEN is selected for component integration. It is not a certified QStarS wind boundary.** The selection follows an audit of source, atomic data, reference inputs and numerical assumptions. QStarS continues to own the stellar structure, evolution, matching, mass-loss eigenvalue, ledger, configuration and validity decisions. A backend's prescribed-flow test is not a selected mass-loss solution.

## Candidate comparison

The audit is pinned in `Reports/evidence/backend_20260928` and the source registry. Availability statements refer to what was actually retrieved, not an assumption that any unavailable item does not exist.

| Approach | Relevant capability | Limitation for this task | Decision |
|---|---|---|---|
| [CMFGEN](https://sites.pitt.edu/~hillier/web/CMFGEN.htm) | Spherical comoving, frequency-dependent transfer; H/He and metal statistical atoms, thermal balance, line/continuum emission and extinction, scattering, moments and force diagnostics. Official source, atomic files and model inputs are downloadable. | Standard tests specify the flow. Cool primordial thermodynamics, material energy/derivatives, radiation work and envelope handoff need independent validation. | Selected; source and diagnostic build executed on Apple Silicon. |
| [PoWR public source](https://github.com/powr-code/PoWR), [hydrodynamic method](https://arxiv.org/abs/1704.08698) | Strong alternative: spherical comoving NLTE, line blanketing, H/He, thermal balance; published hydrodynamically consistent mass-loss method. | Retrieved repository supplies source, not the required atomic/reference input package; supplied build targets Intel Fortran/MKL. Public release cannot be assumed identical to every published hydrodynamic method. | Not rejected physically. More immediate acquisition/porting work than CMFGEN for this experiment. Report 12's source-availability uncertainty is resolved: source is public under GPLv3. |
| [SIROCCO](https://sirocco-rt.readthedocs.io/en/latest/physics/radiative_transfer.html) | Moving-flow Monte Carlo, continuum processes, hybrid macro-atoms, atomic populations and thermal solutions; energy-conserving packet formulation and force estimators. | Its documented line treatment always uses Sobolev. The weak-gradient/subsonic launch and optically thick thermalization region require additional validation; statistical noise complicates an eigenvalue residual. | Valuable independent exterior check, not chosen as the minimum complete launch-region closure. |
| [Lightweaver](https://goobley.github.io/Lightweaver/), [method paper](https://arxiv.org/abs/2107.00475) | Extensible non-LTE multilevel atoms, detailed absorption/emission and scattering, H/He and external atmospheric coupling. | Public formal solvers cover plane-parallel 1D/Cartesian 2D; spherical expanding transfer and a conservative wind energy interface would be new work. | Good atomic/static reference, more new transport code for this objective. |
| [RH](https://github.com/han-uitenbroek/RH) | Mature multilevel transfer and partial redistribution, including a spherical source directory. | Pinned `rhsphere/readatmos.c` explicitly rejects macroscopic spherical flows; `formal.c` rejects the angle-dependent moving-opacity path. | Static reference candidate, not a ready expanding spherical backend. |
| [Cloudy](https://arxiv.org/abs/1705.10877) | Extensive photoionization, chemical and thermal microphysics; useful H/He populations and heating/cooling references. | Escape-probability line treatment is not the required comoving launch-region transfer operator. A detailed atom alone does not determine the wind force response. | Independent atomic/thermal benchmark candidate; not selected as wind transfer solver. |
| Reduced H/He model or group radiation hydrodynamics | Can reduce cost and preserve QStarS ownership; a hot, dilute atomic H/He gas has an explicit ideal translational plus atomic energy. | Report 12's synthetic coefficients do not validate real populations. A group transport library supplies neither a converged atom nor cool-gas material derivatives. Group/level reduction must be checked against resolved calculations. | Use CMFGEN's existing H/He atoms as the first reduction; no sufficiency claim for recombining quasi-stars. |

### Availability, cost and validation evidence

| Candidate | Data/license/reproducibility | Platform/coupling cost | Published evidence and what remains unverified here |
|---|---|---|---|
| CMFGEN | Source 18-Jun-2025 includes GPLv3; atomic 21-Jun-2023, recent models 22-Jun-2023; archive and per-atomic-file hashes recorded. Atomic files retain their own provenance. | Fortran/OpenMP, native arm64 build using gfortran and Accelerate. File-based trial-state coupling plus a read-only spectral export. Local timing measured; Linux build not executed. | [Hillier & Miller (1998)](https://doi.org/10.1086/305350) methodology and supplied stellar models; [Martins et al. (2026)](https://arxiv.org/abs/2605.12141) quasi-star **spectra with prescribed outflows**, not independent hydrodynamic mass-loss selection. No assumption that spectral validation certifies QStarS energy conservation. |
| PoWR | Pinned public commit `b7cb9d27e4970e67de83f9e166476f46b4af7f45`, GPLv3. Atomic/reference inputs not established from the retrieved repository inventory. | Intel `ifort`, 8-byte default integers/reals, MKL configuration and custom model files; Apple port unmeasured. | Sander et al. (2017) solves a hydroconsistent wind. Requires matching source/data/configuration before local replication. |
| SIROCCO | Pinned public GPLv3 source; documented atomic model modes and public examples. | C/GSL, Monte Carlo and parallel workflows; Mac port feasible but not built or timed here. Packet count/seed/error must enter the coupling contract. | Release paper tests against independent photoionization/transfer examples; no locally replicated subsonic quasi-star benchmark. |
| Lightweaver | Pinned MIT source, documented atomic models and Python/C++ interface. | Source build on arm64 plausible; published Intel-Mac wheels do not establish a tested arm64 binary. No local benchmark. | Osborne & Milić (2021) comparisons with established atmospheric transfer; geometry/energy requirements here remain additional work. |
| RH | Source commit `c95bdf66ff6bc1c6cea0ed04dc688100170e8ee6`; repository API did not identify a license. Redistribution terms not resolved. | C source; no local build/timing, spherical moving path explicitly unsupported. | Uitenbroek multilevel transfer methods do not remove the inspected implementation restriction. |
| Cloudy | Public source and atomic documentation; current license text was not successfully retrieved in this audit, so no licensing assertion is made. | C++; no local build/timing. Independent fixed-state reference coupling is conceptually simple. | Ferland et al. (2017) broad atomic/thermal verification; not a validated CMF dynamical launch solver. |
| New reduced/group solver | Would require pinned cross-sections, rate coefficients and reduction errors; no external code's validation automatically transfers. | Small nominal transport cost can be outweighed by atomic/thermodynamic validation work. Not measured. | No new physical validation evidence beyond the executed CMFGEN component. |

CMFGEN is the best **available integration route found in this audit**, not proof that other packages cannot close the problem. Selection does not mean that the reduced model passes, that all CMFGEN options are supported, or that a steady quasi-star solution exists.

## Minimum model actually exercised

The reduced component retains full levels HI=30, HeI=69, HeII=30, with no superlevel compression for those ions, plus fully stripped HII/HeIII. The integer atomic-mass convention is H=1, He=4; number ratio 9.333333333333333:1 gives nominal X=0.7, Y=0.3. Flow radius/velocity and prescribed mass flux come from the supplied Zeta Pup-like reference, with hydro updates, metals, microclumping and X-ray emission disabled. Native level dissolution and the supplied photoionization smoothing settings remain enabled and recorded; consistency of dissolution with material derivatives is an additional bridge requirement. Atomic populations and temperature are solved; adiabatic cooling and population advection are disabled for this initial radiative-equilibrium component.

This is a **hot prescribed-flow atomic test**, not a quasi-star atmosphere, and not a demonstrated minimal atom for cool gas. H−, molecules, missing high levels/continua, nonideal material corrections, relaxation versus flow times, population advection and atom/frequency/radius convergence remain explicit requirements where relevant. No hidden multiplicative opacity, emissivity, heating or force calibration is introduced. Atomic line widths and flow assumptions inherited for a reference calculation are recorded inputs, not fitted quasi-star parameters.

The supplied configuration also has `INC_TWO=F`: two-photon transitions are
omitted in these controls. Their significance must be tested in a cool atom;
the full-level count alone does not imply process completeness. The inherited
photoionization smoothing width is 3000 km/s, and line microturbulence is 20 km/s.
The [physical-control audit](../Reports/evidence/backend_20260928/physical_controls.json)
retains the actual flags. These assumptions have not been justified as a
minimum sufficient quasi-star atom.

The full supplied metal/clumped model is a separate reproduction attempt; changes of numerical controls are recorded. Do not compare its populations or energy directly with the metal-free model as if they represented identical physics.

## Exchange contract

`src/qstars/cmfgen_backend.py` reads completed native material profiles and the opt-in `QSCMF001` spectral stream. `backends/cmfgen/apply_export.py` adds output operations only; original/patched source hashes and the diff are tracked. No change of physical arrays is made. A one-iteration matched-restart comparison of patched/export-enabled and unmodified executables gives identical printed material, population, moment and cooling arrays; its scope is recorded in `export_invariance.json`.

| Quantity | QStarS/interface meaning | Native convention and conversion |
|---|---|---|
| Radius r | cm, strictly increasing outward | Native order outside→inside, units 1e10 cm; reverse and multiply 1e10 |
| Velocity v | cm/s, outward positive | Native km/s; multiply 1e5 |
| Temperature T | K, material temperature | Native 1e4 K; multiply 1e4 |
| Density | g/cm³; mean and inside-clump separately | `DENSITY` is inside-clump; mean = filling factor × density |
| Composition/populations | Species/level number densities cm⁻³ inside clumps; explicit number/mass convention | Full populations from RVTJ, electron density and total nuclei independently checked |
| Frequency | Hz, comoving | Native NU ×1e15. Native FQW already has units Hz; do **not** multiply it again |
| χν, χscat,ν | cm⁻¹, total extinction and coherent scattering inside clumps | Native ×1e−10; stimulated emission is already in net coefficients |
| ην | erg s⁻¹ cm⁻³ sr⁻¹ Hz⁻¹ inside clumps; continuum and line emission excluding coherent χscat J | Native ×1e−10 at the documented export point after line accumulation |
| Jν,Hν,Kν | erg s⁻¹ cm⁻² sr⁻¹ Hz⁻¹, comoving angular moments | Native intensities; H = one-half angular integral μI, outward positive; K uses μ² |
| Boundary conditions | Inner material/grid and incident spectrum or diffusion boundary, outer incoming intensity, luminosity in erg/s | Current preparation is a reference replay. A general QStarS boundary writer is **not** implemented or certified. Native solar luminosity is 3.826e33 erg/s and must not silently use another solar convention |

For coherent isotropic scattering and emissivity,

- E₀ = (4π/c)∫Jν dν, F₀ = 4π∫Hν dν, P₀ = (4π/c)∫Kν dν.
- Q₀ = 4π fvol ∫[(χν−χscat,ν)Jν−ην]dν, positive when **total material energy gains** radiation energy.
- f₀ = (4π/c) fvol ∫χνHνdν, the radial force per mean volume; grad = f₀/ρmean.
- κF = c f₀/(ρmean F₀). Zero flux makes this mean undefined. It is not the Rosseland mean.

These are independently integrated with the native frequency weights. Native printed RVTJ J/H/K are already frequency-integrated cgs quantities. Native mean extinction divided by inside-clump density gives the mass opacity; filling factor must not be applied twice.

The exported scattering split supports **coherent** scattering only. Inelastic/redistributing electron scattering requires its frequency redistribution emissivity and energy exchange to be exported and validated; the adapter rejects that interpretation. Native negative-extinction handling and warning output are retained for audit; QStarS does not repair coefficients by clipping its imported arrays.

With β=v/c and γ²=1/(1−β²), the diagnostic converts bolometric moments:

E = γ²(E₀+2βF₀/c+β²P₀),
F = γ²[(1+β²)F₀+βc(E₀+P₀)],
Prr = γ²(P₀+2βF₀/c+β²E₀).

The lab energy source is γ(Q₀+v f₀). Mechanical work is **not** added a second time to the gas heat. Exact algebraic frame conversion does not make an approximate native transfer equation relativistically exact. The present diagnostic rejects |v|/c≥0.01 and does not certify accuracy below that threshold by itself. Native advective/relativistic transfer options and global energy balance must be tested in the intended regime.

The `GENCOOL` net cooling diagnostic has the opposite sign: cooling minus heating. Its printed percent is 200 net/Σ|individual rates|. It is an **electron thermal** equation, potentially with enabled expansion/nonthermal channels; it is not automatically the total gas atomic-plus-translational source Q₀. A small electron residual does not validate total material/radiation energy conservation.

### Material energy and sonic closure

The diagnostic reconstructs the inside-clump pressure Pgas=(nH+nHe+ne)kT and u=[(3/2)Pgas+Eexcitation+Eionization]/ρclump from the actual full populations and pinned level energies. Neutral H and neutral He ground states define zero atomic energy; successive ionization energies and excited states are explicit. Radiation, molecular, nonideal, degeneracy and turbulence energies are excluded.

FreeEOS uses a different energy reference (including molecular hydrogen), composition masses and equilibrium population response. Its LTE derivatives cannot simply be attached to these non-LTE populations. An energy-zero/composition bridge and derivatives of the same kinetic material model must pass LTE limits, independent heat-capacity tests and finite-rate/flow-time checks before use in the wind equations.

CMF radiation force and heating are nonlocal functionals of the velocity, density, temperature and populations. Report 12's local critical denominator assumes sources without explicit velocity-gradient dependence. Substituting CMFGEN values into that formula would not validate the coupled Jacobian or critical condition. The next physical implementation needs a global boundary-value/eigenvalue formulation (or a justified response reduction), rate/advection and material energy response, with the number of boundary constraints checked. CMFGEN already contains population advection and adiabatic options; their existence is not equivalent to validating this particular coupled model.

### State synchronization and validity

The final spectral loop precedes the final population/temperature correction. Spectral and RVTJ material states must agree within the nonlinear tolerance; an unfinished iterate cannot be exported as a synchronized physical closure. Native errors can return process exit code zero. Completion requires fresh full outputs, exact spectral length, finite grids, retained native diagnostics and independently measured residuals.

Configuration/test preparation is available; production dispatch deliberately raises an error. No physical mass flux enters QStarS through this adapter, no supersonic MLT is re-enabled, and no successful steady-envelope branch is inferred from component convergence. Failed nonlinear iterations are numerical failures. An unsupported regime is a domain exit. A physical absence of a steady solution or multiplicity requires separate continuation/branch evidence.

## Validation sequence and release boundary

Thresholds were recorded in `benchmarks/cmfgen_validation_plan.json` before completed native results: 1e−6 populations/moment normalization, 0.1% nonlinear corrections and radiation-energy budget, 1% electron balance and internal profile refinement. Reference comparisons allow specified percent-level cross-release differences. These are initial component gates; they do not replace eventual wind/evolution error budgets. Reference residuals, including failures in supplied files, must be reported rather than adopted as new tolerances.

Only after independent atomic, force, energy, nonlinear, radius/frequency/atom-size and export-invariance gates pass may the closure be attached to a valid envelope and asked to select mdot. A resulting steady wind must then pass sonic, photospheric, escape and independent ledger checks before evolution/restart/remesh/time studies. The convective handoff and published evolutionary tracks remain downstream gates. QSRD remains deferred.

The executed spatial/angular control reaches 665 radii and 680 rays. The upstream
`RAY(500)` array lacked a bounds check; that run initially failed with a corrupted
DTAU allocation descriptor. A separate source patch raises this capacity to 2048
and guards access before indexing. A hash-matched build capability record is
required to exceed 500 rays; the ordinary 84-point control is unchanged by the
patch. Dust and gamma-ray capacities are not extended. The repaired 665-point
run has a 0.09140% integrated radiation-energy defect and differs from the
333-point radiation/force profiles by at most 0.2019%. These pass the stated
operator targets; the incoming populations remain out of equilibrium, so this
does not close thermal, atom-size or wind convergence. Failures and exact
timings remain in [Report 13](../Reports/Report_13.md).
