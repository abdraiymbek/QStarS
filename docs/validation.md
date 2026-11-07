# Validation report — development release, not QStarS v1

**The general-purpose v1 criteria are not satisfied.** This campaign closes
important numerical checks for a specified hot-envelope problem and adds direct
MESA microphysics and a stationary advective-wind reference. It does not
reproduce complete published quasi-star evolutionary tracks.

The latest recorded full automated suite (2026-10-22) has **341 passing tests and
37 skips**; compiled MESA/FreeEOS and native-compiler tests skip unless their
environment is configured. Report 25 records 357 passing tests after the project
moved to the external SSD. (Report 15 recorded 122 passing tests with MESA and
FreeEOS enabled.) The last recorded `qstars validate --require-v1` campaign
(2026-09-30) has **29 passing checks and two retained printed-value
discrepancies** and correctly exits 2. Passing software tests is distinct from
passing the scientific release gates below.

The [Report 8 atmosphere campaign](../Reports/Report_8.md) now supplies a specific
physical blocker rather than only a list of missing modules. Native MESA's
implicit starting pressure dominates the low-gravity grey boundary. Replacing
that with explicit pressure and audited direct FreeEOS fixes that particular
problem, but the 6000 K, g=1 candidate still demands Mach 2.3367 convection at
tau≈20 under both Henyey and ML1. The 8000 K, g=3 candidate has a depth below
the photosphere of 6.68% of radius, incompatible with one-percent
plane-parallel accuracy. These are converged component-model validity failures,
not physical destruction or completed evolution. See [the equations, sensitivity
tests and resolution criteria](atmosphere.md). General v1 remains blocked.

## Current status (Reports 17–25, the dense-window run and the v1 plan)

**No gate has changed status since Report 23, and QStarS is not v1.** The chain
that v1 now waits on is the cool non-LTE closure at 167 radii, followed by
Gates 1–5 of [Report 23](../Reports/Report_23.md) §8. Gate 6, the Ball et al. 2011
replication, is independent of that closure.

- [Report 17](../Reports/Report_17.md) replaced the Report 16 moment propagation with a
  nonnegative formal-solution transfer. The released iteration still did not
  converge. [Report 18](../Reports/Report_18.md) measured the available
  linearizations against the exact residual. [Report 19](../Reports/Report_19.md)
  passed an analytic transfer tangent, a component result.
- [Report 20](../Reports/Report_20.md) built the exact Jacobian action (RE row conditional),
  and [Report 21](../Reports/Report_21.md) completed it with level dissolution and
  INC_AD/INC_ADV. [Report 22](../Reports/Report_22.md) made the solve faster without
  changing the model. The released closure did not converge in any of them.
- [Report 23](../Reports/Report_23.md) traced the non-convergence to numerical causes and fixed
  them in builds 32–35. The causes were a fixed-J measure blind to trapped
  transitions, the `OLD_DEFAULT` two-photon approximation in single-iteration
  passes, and a wrong free-free Gaunt-factor T derivative. It withdrew every
  earlier "converged" cool state. It also delivered the Gate 1–3 closing tooling,
  the Gate 4 wind module and the Gate 5 node interface (implemented and tested,
  not validated; see below), and the registered Gate 6 execution: **FAIL**.
- [Report 24](../Reports/Report_24.md) found that the Report 23 solves broke the exact
  number and charge conservation rows, and that full trial steps produced
  supercritical masers. `--conserve` and build 36 fix those defects.
- [Report 25](../Reports/Report_25.md): stage 1 (fixed T) lowered the row-scaled residual
  from 2.125 to 0.051 with conservation to 1e-15. It then stopped at a solver
  limit along physically soft front modes. Stage 2 (T released) started from
  there.
- The stage-2 dense hydrogen-front-window run (2026-10-31 to 11-03; recorded in
  `results/nk167_trap/free_r25/README.txt` and `driver.log`; no Report 26
  was written) used about 1,200 exact columns per iterate, with 13 builds
  totalling about 70 h. It lowered the frozen reference norm by only about 3 %
  (0.02930 → 0.02837). At iterate 16 every trial was rejected. The last valid base
  is `results/nk167_trap/free_r25/it16_trial_l0.5`; `history.json` names a
  rejected iterate-17 trial as accepted.
- The [v1 execution plan](../Reports/V1_execution_plan_2026-11-06.md) (2026-11-06) keeps
  the scientific v1 definition. It stops repeats of the dense-front run and orders
  the work: P0 baseline, P1–P2 step construction, P3 component gates, P4 outer
  solution and coupling, P5 Ball repair, and P6 the requirements beyond Gates 1–6.
  It changed no registration or gate.

| Gate (Report 23 §8) | Recorded status (Report 25 §9) | Implemented, not validated |
|---|---|---|
| 167-radius T-released closure converged | FAIL (not converged) | Builds 32–36, driver options in `scripts/cool_newton_krylov.py` and `scripts/cnk_trap.py` |
| 1. Post-convergence cleanup and audits | NOT EVALUATED | `scripts/post_convergence.sh`, `scripts/trace_cleanup.py`, `scripts/native_correction_check.py`; closing rule in `benchmarks/cmfgen_cool_validation_plan.json` |
| 2. 333 radii | NOT REACHED | Nested restart; base pass must be rebuilt (stale traces) |
| 3. Doubled frequency | NOT REACHED | `benchmarks/cmfgen_cool_frequency_plan.json`, `scripts/compare_resolution.py --profiles` |
| 4. Self-selected CMFGEN wind | NOT REACHED | `src/qstars/cmfgen_wind.py` (21 tests), `benchmarks/cmfgen_selected_wind_plan.json`, [contract](cmfgen_selected_wind.md) |
| 5. Node-coupled conservative evolution | NOT REACHED | `src/qstars/wind_node.py` (11 tests, synthetic node; runs report `validated: false`), [contract](wind_node.md) |
| 6. Ball et al. 2011 benchmark | FAIL | Smoke gate numerical blocker at M_BH = 46 M☉; Q2 fails at row 0; campaign not run |

The per-criterion inventory, with artifacts, registrations, gaps and plan phases,
is in [gate_inventory.md](../Reports/evidence/v1_exec_20261106/gate_inventory.md).

## Earlier outer-closure campaigns (Reports 13–16)

[Report 15](../Reports/Report_15.md) has outcome **NUMERICAL BLOCKER**. It repairs
H-minus FF initialization, last-depth opacity and missing continuum/EHB
derivatives, and retains every level in convergence checks. Fixed-J/T cool
population initializations settle on 84/167 radii, with 5/11 recombined cells;
this is not thermal/radiation closure. Released solves and two matched
angular/frequency recovery diagnostics fail. Captured frequency-coupling
operators reproduce runaway J under independent linear algebra. The relaxed
response, global mdot/sonic solve, physical envelope handoff and evolution
remain withheld. See the report's campaign matrix and retained evidence.

[Report 14](../Reports/Report_14.md) has outcome **BACKEND/PHYSICS BLOCKER**.
Cool H−/H/He controls expose native population-output and Rayleigh-indexing
defects, now repaired and tested. The ideal material differential explicitly
includes level energy and population changes. The corrected 84-zone T6 control
contains nine cool partially recombined cells but fails nonlinear convergence,
electron balance (150.59%), global radiation energy (8.6847%) and cumulative
interior radiation energy (14.8590%). T9 develops nonfinite luminosities; the
requested Newton continuation develops a radiation-moment instability. No
nonlocal response from those failed states is fed to a wind solver. Refinement,
all failed inputs and complete evidence are documented in the report. No
physical envelope, evolving wind or published full evolutionary track is
validated; the Mach-one MLT guard remains active.

The preceding [Report 13](../Reports/Report_13.md) records the hot component.
CMFGEN supplies a native spherical non-LTE component with explicit atomic data,
spectral moments, force and energy exchange. Its adapter normalization and
read-only export pass independent checks. Fresh 84/167/333/665-point transfer
controls reduce the radiation-energy defect from 6.312% to 1.532%, 0.3921% and
0.09140%; the last meets the predeclared 0.1% target. The finest common-radius
radiation/force differences are 0.2019%, below 1%. These are fixed incoming-state
operator tests, with populations still out of equilibrium, not converged
quasi-star structure or selected mass loss. The final hot advective continuation
passes 12 per-run checks, with 0.296% electron imbalance and 0.0208% global
radiation-energy defect. Maximum cumulative interior defects are 0.148% in that
case and 0.230% in the finest frozen-state control; these are explicitly
separate from the passing endpoint targets. No cellwise energy or thermal-mesh
certificate is claimed. Cool-gas atomic sufficiency,
consistent material derivatives, the global CMF response and envelope handoff
remain unvalidated; evolution and QSRD were not attempted.

The [backend contract](nlte_backend.md) defines units, energy reference, frames,
signs, source provenance and the boundary between diagnostics and production.
All prior v1 failures remain in the [gate matrix](../Reports/evidence/completion_gates.json).

[Report 11](../Reports/Report_11.md) classified its campaign as **OUTER CLOSURE FAIL**.
The published spherical grey moment table, moving coherent-scattering radiation
work, and constructed small-drag wind comparisons pass their stated component
limits. The 801/1601-node grey profiles differ by at most 0.0329%; wind profiles
at 401/801 nodes differ by at most 0.0471%. Extending the wind exterior from
30,000 to 100,000 base radii changes the photosphere by 0.0695%, below the 0.2%
limit. Shorter exterior failures remain archived. These are fixed-grey controls,
with no material temperature/ionization or mass-loss eigenvalue; they do not
validate a physical quasi-star outer boundary. No coupled evolution was run.

The campaigns retain **190 passing and ten failing numerical checks**, including
coarse-grid, insufficient-depth, finite-exterior and approximation counterexamples.
The report separates resolved finer controls from these retained failures.
See [the method and validity domain](escaping_radiation.md),
[the candidate audit](../benchmarks/outer_candidates.json), and
[the current gate matrix](../Reports/evidence/completion_gates.json).

[Report 10](../Reports/Report_10.md) remains the evidence for the unreproduced
Nakauchi branch, moving-surface identities and mandatory Mach-one MLT exit.
Its attribution campaign was not repeated.

## Spherical outer component campaign

[Report 9](../Reports/Report_9.md) and [the spherical method documentation](spherical_wind.md)
record the new gas-sonic wind, realistic recombination, separate pressure-work
fluxes, independent first-law/diffusion/momentum checks, analytical spherical
limits, known-solution matching recovery and a digitized published profile
comparison. Numerical component convergence passes. The published escaping
branch, complete escaping-radiation transport and a validated Report 8 base
match remain unresolved. No coupled evolution is enabled on that evidence.
This is further component development, not completion of the requested outer
campaign or a general 1D no-go result.

## Controlled spatial and nonlinear convergence

The benchmark has a 10,000-solar-mass envelope, 100-solar-mass BH, initial radii
1e12–1e14 cm, external pressure 1e8 dyn/cm², prescribed boundary luminosity
1e42 erg/s, rest accretion 2e-4 solar masses/year, feeding 0.01 and wind removal
0.008 solar masses/year. It evolves for ten years. These are transparent
controlled boundary conditions, not a full published stellar model.

`configs/hot_resolved_validation.yaml` selects the sparse flux solve,
surface-refined log-radius mesh and conservative backward-Euler mechanics.
The latter's kinetic dissipation is deposited as internal heat, separately
recorded, and excluded from the external energy balance. The historical
midpoint method remains selectable. Its high-resolution acoustic ringing and
convection-domain failure are preserved, not relabeled as physical destruction.

The 256/512/1024/2048-zone backward-Euler runs all complete. Initial density
mass-weighted logarithmic RMS errors against a separate continuum boundary-value
solution are 1.553e-4, 3.881e-5, 9.703e-6 and 2.426e-6: successive factors of
approximately four. The continuum reference is not the evolution solver.

Final interior comparisons use extrema of piecewise cubic interpolants on a
common enclosed-mass domain. Their squared differences are integrated by exact
polynomial quadrature. A separate whole-mass, piecewise-constant overlap norm
is also retained. Interpolation does not turn either representation into an
exact physical solution.

| Relative error bound from logarithmic maximum | 512 vs 2048 | 1024 vs 2048 |
|---|---:|---:|
| Density | 0.248% | 0.0597% |
| Temperature | 0.0221% | 0.0201% |
| Pressure | 0.0880% | 0.0803% |
| Specific internal energy | 0.250% | 0.0207% |
| Radius within the envelope | 0.0220% | 0.0200% |

The common-center comparisons exclude 2.13e-6 and 8.40e-7 of envelope mass,
respectively. Whole-mass overlap log-RMS errors at 1024 vs 2048 are 0.290% in
density, 0.0962% in temperature, 0.385% in pressure, 0.0984% in specific energy
and 0.112% in radius. Entropy difference divided by local Cv is 4.34e-5.
The 256-zone model still has about 1.24% local specific-energy error; it does
not pass the one-percent local-profile target.

All cases use the same maximum timestep, but boundary replacement and state
change limits cause smaller accepted steps at 2048 zones. This is therefore a
combined spatial/adaptive-time campaign. The separate timestep campaign below
bounds its time-discretization contribution. Small global radius or luminosity
differences are not used as substitutes for these interior measurements.

A complete fresh 256/512/1024/2048-zone repeat passes unchanged-source
verification. All saved state arrays, time/BH/step-controller scalars and
ledgers are bitwise identical to the earlier campaign. The
[verified spatial campaign](../benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/campaign.json)
and [repeat comparison](../benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/previous_campaign_comparison.json)
are the primary current-source mesh evidence.

## Time integration, retries, remapping and restart

For the 512-zone backward-Euler benchmark, reducing maximum timestep from
0.02 to 0.005 years changes the largest measured profile quantity, pressure,
by 0.0757%; density changes by 0.0563%. The 0.01-to-0.005 comparison is smaller
(pressure 0.0252%), consistent with the method's first-order time accuracy.
Solver tolerances 1e-8/1e-10/1e-12 agree to approximately 2e-13 in these profile
norms. Changing boundary-remap thresholds 0.4/0.2/0.1 changes density/specific
energy by at most 0.00370%. Restart continuation is bitwise identical, including
volumes, velocities, solved fluxes, radius compensation and the whole ledger.

A deliberately independent time-method comparison matters: at ten years,
midpoint and backward Euler with maximum dt=0.02 years differ by about 13% in
local density/specific energy. Reducing midpoint's maximum step to 0.005 years
reduces the matched-method density discrepancy below **0.801%** and specific
energy below **0.768%**, with pressure below 0.041%. A one-year campaign at
0.02/0.005/0.00125 years also keeps inter-method profile discrepancies below
0.44%, but density differences are not strictly monotone. Large-step midpoint
is not certified by small ledger residuals. These checks do not validate a
pulsation or shock calculation.

The 2048-zone final energy residual is 1.32e-12 of energy throughput; the
512-zone residual is 3.15e-13. Numerical kinetic heating in the 512-zone run
is 2.82e41 erg, about 4.07e-10 of throughput. Its small size alone is not an
accuracy argument: density can be sensitive even when an energy diagnostic is
tiny. Analytical regressions establish that this heating is positive, conserves
total energy and decreases with timestep.

Earlier high-resolution failures exposed dense initialization cost, a Newton
iteration ceiling, an invalid-geometry penalty root, cancellation in cell
volumes and cumulative mass coordinates, depleted boundary cells, and
unresolved acoustic oscillations. Each repair has an independent diagnostic or
regression. Forty 2048-cell remaps preserve mass, volume, energy and geometry.
The old failed and stopped runs remain in the evidence archive.

## MESA microphysics and stellar infrastructure

The optional interface directly calls the public MESA 24.08.1 EOS and radiative
opacity APIs for fixed H1/He4 composition. It includes ionization. It does not
replace these libraries with fitted tables or an invented ionization EOS.
Radiative and conductive opacities are distinguished; omitted conduction is
screened independently. The build records the MESA commit, local changes,
bridge source and library hashes. An additional inventory records 8,135
microphysics data files totaling 15.54 GB.

The 153-state audit exactly matches a separate Fortran process. **Five states
fail the one-percent first-law screen**; the maximum normalized defect is 0.315
near rho=1e-4 g/cm³ and T=6812 K in an OPAL/SCVH–FreeEOS blend. These states
remain excluded. A matching library call is not proof of thermodynamic
consistency or of the correctness of a whole stellar model.

A direct-MESA 128-zone, ten-year run completes with independent energy residual
3.45e-11 of throughput. A timestep/tolerance/restart campaign also completes:
density changes by at most 0.431% for dt=0.02 vs 0.005, and restarts are bitwise
identical. Its 0.1–0.4 remap thresholds do not trigger remapping, so those
identical cases are **not** evidence of MESA remap robustness. General production
microphysics remains conditional on domain checks and missing atmosphere,
convection and composition/stellar-evolution coupling.

Cool-envelope trials terminate at the radial Mach boundary. Their polytropic
initialization is not a relaxed cool quasi-star. These exits are not physical
disruptions or published-track replicas. An ordinary MESA smoke run reaches
100 models; two pinned upstream QUEST seed experiments, including an explicitly
logged startup CGS-conversion experiment, stall during inner relaxation within
the wall budget. Neither constitutes a completed quasi-star track.

## Published results and disagreements

The existing static benchmarks retain Ball n=2/3/4 maximum BH fractions
0.126753/0.101515/0.088521, corrected Bondi N=1/2 limits 0.017040/0.005756,
and the saturated composite printed-cutoff limit 0.617828. They are independent
structure comparisons, not evolutionary validation. CB24's printed Ki and
BH/inner-mass values correspond to scaled radius one rather than the stated
entropy minimum near 1.015; both interpretations remain explicit.

The new Fiacconi–Rossi solver integrates stationary radiation-dominated wind
ODEs with radiative, advected, kinetic and gravitational fluxes. Base Mach
0.8 and the paper's range 0.75–0.85 are explicit; no parameter is fitted.
Changing integrator (Radau/DOP853), tolerance and outer boundary changes the
sampled profiles by less than 5e-6 in log magnitude in the tested cases.
Independent momentum/radiation integrals supplement algebraic energy closure.

Of 461 digitized Figure 4 vertices, 349 fall within the propagated base-Mach
band plus the stated line-width uncertainty under the literal diagnostics.
This is **not a full passing replication**. Thin-region sound speed differs
between continuation of p=U/3 and p=g(tau)U; both proxies are saved, with no
thermodynamic claim in the thin limit. Digitized luminosity sums are roughly
4.82–4.87 versus caption total 5, at a precision comparable to the plotted
line widths. Caption rounding, plotting conventions or unavailable original
data may explain part of the discrepancy; the code does not choose a fitted
energy to force agreement. The evolving envelope is not yet coupled to this wind.

**No complete published quasi-star evolutionary track has been reproduced.**
Ball's tabulated targets, QUEST attempts and the pending independent wind/evolution
families are listed in `benchmarks/manifest.yaml`. No cross-paper evolutionary
prediction without retuning has passed.

## Release gate matrix

PASS means the stated scope has evidence. CONDITIONAL does not certify general
quasi-star use. FAIL marks a required capability or validation result that is
missing; none of the failures below is waived as impossible in spherical 1D.

| User's priority gate | Status | Evidence and remaining boundary |
|---|---|---|
| Spatial/nonlinear convergence | CONDITIONAL | Hot 512–2048-zone profiles, timestep and inter-method comparisons pass the one-percent target; production atmosphere/cavity/wind models untested |
| Production microphysics | FAIL | Direct EOS/opacity and native atmosphere/MLT references work; the explicit-pressure candidate fails geometry/subsonic validity; coupled stellar boundary remains unvalidated |
| Inner cavity and model families | FAIL | Fixed-excision ledger works in the thermal engine; a moving cavity exists only on the quasi-static Ball path (Report 23, verdict FAIL); a conservative moving physical cavity and alternative evolving inner flows remain unvalidated |
| Wind physics | FAIL | Stationary advective wind component exists; the CMFGEN-selected wind module and node interface (Report 23) are implemented and tested, not validated; envelope/atmosphere coupling and complete validation absent |
| Multiple full published evolutions | FAIL | Static results and wind-component comparisons are insufficient; the registered Ball et al. 2011 replication is FAIL (Report 23); zero complete evolutionary replicas |
| Cross-validation without retuning | FAIL | No multi-paper evolutionary predictive validation |
| Numerical robustness | CONDITIONAL | Hot-model campaigns pass; latest full suite 341 passed / 37 skipped (357 passed in Report 25); full production-model campaigns remain |
| Apple Silicon performance | CONDITIONAL | Useful local runs measured; full stellar-track costs unestablished |
| QSRD readiness after independent v1 | FAIL | Deliberately deferred until the preceding framework gates pass |

[Machine-readable gates](../Reports/evidence/completion_gates.json) map these
priorities to the original twelve completion criteria. Source traceability,
finite thermal evolution, separate stream accounting and execution safeguards
have working implementations; they do not cancel the failed general-use gates.

## Reproducibility and evidence

New runs freeze source bytes at import and verify that package files remain
unchanged through completion. Fresh final validation and benchmark runs pass
that check. Earlier exploratory campaigns predate it and require their first
source snapshot and development chronology; later files in a long-running
edited process must not be assumed to be its executed source. A fresh MESA
verification agrees with the earlier ten-year density within 4e-11.

The new evidence archive is `benchmarks/results/2026-09-28-closure/`. It retains
full histories, ledgers, checkpoints, source snapshots, first/final profiles,
figure calibration, comparison reports and selected failures, with hashes in
its manifest. Intermediate profiles and upstream binaries/source are excluded
from that compact archive; all original local result directories remain intact.
Historical evidence remains under `benchmarks/results/2026-09-28/`.
The fresh verified spatial repeat is separately preserved under
`benchmarks/results/2026-09-28-verified-spatial/`. Byte-length/SHA-256 checks
pass for all 1,165 files in the first archive and 52 files in the fresh one.
Both archives preserve their original-byte manifests; compressed logs are
verified after decompression.

See [model families](model_families.md), [validity limits](validity.md),
[equations](equations.md), [MESA](mesa.md), [source registry](sources.md) and
[MacBook performance](performance.md). These remaining requirements cannot be
narrowed away while still calling this a general-purpose validated v1.


## Report 12 self-selected outer prerequisites

The requested outcome is **FAIL**: no physical quasi-star currently selects a
validated escaping thermodynamic wind. `scripts/validate_self_selected_outer.py`
and `scripts/validate_self_selected_refinement.py` provide two frozen campaigns;
`analyze_self_selected_outer.py` compares their archived data. The 139 checks
retain 23 failures, while refined thermal profiles/energy, smooth frequency
limits, adiabatic eigenvalue profiles and native rejected-trial recovery pass
scoped targets. Full native regression: 97 tests pass. All 89 archived files are
SHA256-verified. The user-required physical-envelope, coupled-evolution and
complete published-track stages were not attempted because the joint physical
closure gate remains open. [Report 12](../Reports/Report_12.md) and
[its equations/validity domain](self_selected_outer.md) distinguish these controls
from a production stellar mass-loss model.

Report 16 retains **NUMERICAL BLOCKER**. Exact scaled recurrence attribution and
native fixed-factor eta/chi response pass their bounded checks; shared-face
and unclipped N/H candidates fail. No joint cool energy state, two converged
resolutions or selected wind exists. See [Report 16](../Reports/Report_16.md)
and its [gate matrix](../Reports/evidence/moment_closure_20260930/gate_matrix.json).
