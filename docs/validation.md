# Validation status — development build, not QStarS v1

**The user's v1 completion criteria are not satisfied.** This repository contains
working simulation software and auditable evidence, including negative results.
It must not be described as a generally validated quasi-star evolution code.

## Analytical and software checks

The automated suite currently has **29 passing tests**. It covers EOS first-law
derivatives and inversions, initial hydrostatic residuals, zero spurious convection
in an n=3 polytrope, an exact linear-T⁴ diffusion-flux case, strict opacity-table
bounds, wind energy limits, separate equal feeding/wind, BH rest-energy conversion,
immutable retries, mass/volume/total-energy remap, restart equality, refinement
followed by restart, control-breakpoint mass integration, input rejection,
EOS-domain exits, configurable static reference runs, and bounded-run checkpoints.

Restart continuation is array-exact in the tested environment. This does not
establish bitwise reproducibility across BLAS, platforms, or future versions.

## Evolution and conservation

Six configured examples completed. In the ten-year two-stream pulse example:

- Delivered mass: 0.16 solar masses.
- Escaping wind: 0.08 solar masses.
- BH-consumed rest mass: 0.002 solar masses.
- BH gravitational growth: 0.0018 solar masses.
- Final envelope: 10000.078 solar masses.
- Final mass residual: zero at reported floating-point precision.
- Final energy residual: 4.30e-15 of initial internal energy, 1.16e-13 of
  initial binding energy, and 1.93e-12 of transported energy throughput.

The finite-pressure, concentrated 64-zone example has a larger energy error:
about 8.5e-6 of throughput. A small internal-energy-normalized residual must
not hide this distinction. These budgets are for the defined Newtonian domain
and unresolved engine exchange; they are not a relativistic whole-system mass.

## Resolution and tolerance

The retained final campaign uses mesh counts 16/32/64/128, timesteps from
4e-4 down to 5e-5 years, and nonlinear tolerances 1e-7/1e-9/1e-11.
Zero-pressure models all reach the common 0.01-year comparison time.

- Relative luminosity differences to the smallest timestep decrease from
  3.48e-8 to 1.25e-8 to 3.82e-9.
- Tightening the nonlinear tolerance reduces the energy residual from
  about 8.6e-9 to the few-e-15 floor in this test.
- Mesh internal-energy differences to 128 zones are about **11.1%, 6.50%,
  and 2.80%**. This is monotone improvement, not adequate research convergence.
  Fixed outer radius and a nearly Eddington luminosity can appear converged
  while the internal structure remains inaccurate.
- Concentrating zones at the core and applying finite external pressure improves
  initial structural convergence, but the 128-zone time-dependent case fails
  the nonlinear solve. Its unequal final time is excluded from comparisons.

The trust-region fallback investigation is preserved as a rejected approach.
It made slow progress with repeated timestep reductions, and was interrupted
after 275 seconds at only 2.68e-5 years. Its original interruption reason was
stale; `stop_context.json` records the actual stop without overwriting the
failed artifact. The release retains the earlier solver, and now handles future
interruptions and explicit wall budgets correctly.

## Published results

| Published result | Computed | Reference | Interpretation |
|---|---:|---:|---|
| Ball n=2 maximum BH fraction | 0.126753 | 0.127 | Within printed rounding |
| Ball n=3 maximum BH fraction | 0.101515 | 0.102 | Within printed rounding |
| Ball n=4 maximum BH fraction | 0.088521 | 0.089 | Within printed rounding |
| Corrected Bondi, N=1 | 0.017040 | ~0.017 | Reproduced |
| Corrected Bondi, N=2 | 0.005756 | ~0.0058 | Reproduced |
| Saturated composite, printed cutoff | 0.617828 | ~0.618 | Reproduced static structure |
| Saturated entropy-minimum location | 1.014835 | ~1.015 | Reproduced |
| Ki at entropy minimum | 0.468357 | 0.46554 | Discrepancy retained |
| BH/inner mass at entropy minimum | 0.625102 | 0.6321 | Discrepancy retained |

The last two printed numbers instead match evaluation at scaled radius one.
See [the source interpretation](literature.md). The code preserves both
interpretations; it does not tune a parameter to conceal the mismatch.

These are independent structure/equation benchmarks. **No complete published
quasi-star evolutionary track, radius/temperature history, or advective wind
solution has been reproduced.** Ball's tabulated stellar targets are stored,
explicitly pending. No figure digitization or fitted coefficients were used.

## v1 gate matrix

| Required gate | Status |
|---|---|
| Sourced equations and boundaries | Restricted model documented; production closures incomplete |
| Mass and energy conservation | Demonstrated for tested Newtonian cases, not relativistic total mass |
| Retry/remesh/restart/timestep | Controlled regressions pass; adaptive research mesh incomplete |
| Representative convergence | **Open:** coarse spatial error and concentrated-mesh solve failure |
| Several published quasi-star models | Static results pass; **full evolution/winds open** |
| Explicit selectable prescriptions | Implemented laws and static families are named and traceable |
| Finite thermal evolution | Implemented for restricted hot envelopes |
| Separate feeding and winds | Implemented effective boundary closure; geometry unresolved |
| Validity detection | Temperature, ionization, degeneracy, optical depth, flux, compactness, Mach, conservation and solver checks; not all convection domains covered |
| Reproducible runs | Configurations, snapshots, dependencies, ledger and restart data saved |
| Local performance | Useful restricted examples measured; high-resolution exception documented |
| QSRD as external experiment | **Deferred by the requested prerequisite:** framework not yet fully validated |

The QSRD recovery demonstration and a science campaign have deliberately not
been run. Completing v1 requires resolving spatial/nonlinear convergence,
validated ionization/opacity/atmosphere microphysics, conservative moving
cavities and alternative inner flows, transonic wind integration, and multiple
full published evolution comparisons. MESA integration is the preferred path
for the microphysics/stellar-evolution portion, not an implemented feature.

## Evidence

The final campaign, complete saved example runs, source snapshots, checkpoints,
and selected failed attempts are tracked under `benchmarks/results/2026-09-28/`.
All original local runs are also retained under `results/`. Sequential reports
describe repairs and rejected approaches. The archived campaign JSON includes
every configuration's termination status; incomplete runs are never treated as
common-time convergence data.

- [Final campaign](../Reports/evidence/final_validation.md)
- [Machine-readable evidence](../Reports/evidence/final_validation.json)
- [Measured performance](performance.md)
- [Source register](sources.md)
