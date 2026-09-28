# Report 7 — Numerical validation and remaining v1 release gates

Assessment date: 2026-09-28 UTC. Software version: **0.1.0 development**.
**QStarS v1 has not been achieved.** None of the missing cavity, atmosphere,
wind-coupling or published-evolution requirements is waived as impossible in 1D.

## What this development campaign establishes

The concentrated-envelope failure is no longer a generic inability to solve
the discretized equations. Sparse initialization, an inverse flux formulation,
conservative fine-mesh arithmetic, boundary redistribution and a selectable
implicit mechanical method allow the controlled hot problem to evolve for ten
years through 2048 zones. Independent internal-profile comparisons are now
available, rather than relying on nearly fixed global radius or luminosity.

The initial density error against a separate continuum boundary-value solution
decreases by approximately four with each doubling of resolution. At ten years,
1024 versus 2048 zones differ by at most 0.081% across the compared density,
temperature, pressure, specific energy and radius fields. At 512 versus 2048,
the maximum is about 0.250%. The common-center domain excludes less than
2.2e-6 of the mass; separate whole-mass overlap norms retain the boundary
regions. The 256-zone specific-energy error exceeds 1%, so that resolution does
not pass the adopted local-profile target.

The maximum timestep is shared across meshes, but accepted timesteps are not
identical because boundary-cell transport limits become stricter. A separate
512-zone timestep campaign changes pressure by at most 0.0757% when the maximum
step decreases from 0.02 to 0.005 years. Solver and remapping sensitivity are
smaller. Restart continuation is bitwise identical, including ledger and
compensated radius arithmetic.

A cross-method test exposes a limitation that energy closure alone misses.
Midpoint and backward Euler disagree by about 13% in local density/specific
energy at a 0.02-year maximum step. Reducing midpoint's step fourfold reduces
their matched-step difference below 0.81%. The coarse midpoint result is not
certified. A separate 2048-zone midpoint run exits the lossless-convection
domain; this is not a demonstrated physical instability or envelope destruction.

The backward-Euler option deposits the positive kinetic-energy loss derived
from its discrete momentum equation as internal heat. This conversion is
reported separately and is not counted as external heating. It is not an
energy-residual correction. The implementation and regression tests retain
the historical midpoint method for comparison.

## Production microphysics progress and limits

The pre-existing MESA installation is now used directly for EOS, ionization
and radiative opacity. Report 5 corrects the earlier erroneous statement that
MESA was unavailable. The bridge records release, commit, local changes,
compiler output, source and binary hashes. ABI 3 distinguishes radiative from
radiation-plus-conduction opacity. The independent dataset inventory hashes
8,135 files totaling 15,541,148,301 bytes.

All 153 audit states match a separate Fortran process exactly. Five fail the
one-percent first-law consistency screen, with a maximum normalized defect
of 0.315 near 6812 K in an EOS blend. They remain rejected; exact API agreement
does not establish thermodynamic consistency. No blend coefficients were
adjusted. The present bridge supports fixed H1/He4 composition, not a general
chemical network or a production stellar-evolution backend.

A direct-MESA 128-zone ten-year run completes with an energy residual of
3.45e-11 of throughput. Timestep refinement changes density by at most 0.431%
in the tested campaign, and restart is bitwise identical. Its remapping
threshold variants never trigger remapping, so they supply no MESA remapping
validation. Cool initial trials encounter a radial-Mach domain exit and are
not relaxed published quasi-stars.

An ordinary upstream MESA smoke model reaches model 100. Original and explicitly
modified QUEST startup experiments remain incomplete during inner relaxation.
Their histories and disabled upstream energy-convergence controls are recorded.
Wrapping those executions would not supply QStarS's missing independent whole
stellar budget.

## Published winds and evolutionary validation

The Fiacconi–Rossi stationary reference includes radiation enthalpy advection,
kinetic energy and gravity. Integration tolerance, outer boundary and two
independent ODE integrators agree closely. Momentum and radiation quadratures
are checked independently of algebraic energy closure. Base Mach 0.8 and the
paper's 0.75–0.85 range are explicit inputs, not fitted parameters.

Only 349 of 461 digitized Figure 4 vertices lie within the stated model band
and line-width uncertainty. Thin-region sound-speed conventions and several
percent luminosity differences are preserved. The comparison is partial,
not a full passing reproduction. This stationary component is not coupled
to the evolving envelope or a solved atmosphere.

**Zero complete published quasi-star evolutionary tracks have been reproduced.**
The benchmark manifest includes explicit pending Ball, QUEST, Fiacconi–Rossi,
Dotan and Hassan cases. Static loaded-polytrope agreement and a ten-year
controlled hot calculation do not establish million-year stellar evolution,
lifetimes, wind histories or cross-paper predictive validation.

## Regression, provenance and performance evidence

The finalized test suite passes **52 tests with MESA enabled** in 6.95 seconds.
Four tests require the optional MESA library. The CLI scientific campaign has
29 passing checks and two retained CB24 printed-value discrepancies;
`qstars validate --require-v1` correctly exits 2. The automatic test count
is not used as a proxy for scientific release readiness.

Fresh runs freeze source bytes at process import and check for subsequent
changes. Older exploratory runs predate this protection: their first snapshot
and chronology are needed, and a later on-disk snapshot alone cannot identify
already imported code. A fresh direct-MESA rerun matches the earlier ten-year
density within 4e-11. A fresh complete 256/512/1024/2048-zone spatial campaign
passes unchanged-source verification at every resolution. Its state arrays,
time, BH mass, step controller and whole ledgers are **bitwise identical** to
the earlier campaign. That comparison now supplies the missing current-source
confirmation of the internal-profile convergence results.

Sequential fresh-process MacBook measurements use an Apple M5 with 16 GiB RAM.
The 512-zone resolved hot run takes 24.215 seconds and 88.9 MiB peak RSS;
the 128-zone MESA run takes 31.293 seconds and 1045.4 MiB. Both evolve ten years.
These are useful individual local runs, not measured full published stellar
tracks. Independent CPU grids can benefit from cloud resources. The Linux MESA
bridge and million-year production-track cost remain unverified.

The archive preserves failed approaches, stopped runs, complete histories,
ledgers, checkpoints, first/final profiles and source snapshots. All 1,165 files
in the first new archive, plus 52 in the fresh spatial archive, pass byte-length
and SHA-256 verification, including
decompression where appropriate. Archive integrity does not prove scientific
validity or execution provenance. Full original local result directories remain
intact; upstream binaries and whole upstream source trees are not redistributed.

## Release assessment

| Priority | Gate | Status |
|---|---|---|
| 1 | Spatial and nonlinear convergence | CONDITIONAL: controlled hot problem passes |
| 2 | Production microphysics | FAIL: atmosphere/full stellar coupling absent |
| 3 | Moving cavity and evolving inner families | FAIL |
| 4 | Conserved coupled wind and atmosphere | FAIL |
| 5 | Several complete published evolutions | FAIL |
| 6 | Cross-paper predictions without retuning | FAIL |
| 7 | Numerical robustness | CONDITIONAL: bounded tested models |
| 8 | MacBook performance | CONDITIONAL: full-track cost unknown |
| 9 | QSRD readiness after independent v1 | FAIL: intentionally deferred |

The physical limitation that can properly be narrowed is simultaneous
multidimensional inflow/outflow geometry: 1D can provide independent averaged
streams, not resolve that geometry. This does not excuse any of the other
failed gates. QSRD remains deferred. A general-purpose validated working v1
and a full published-evolution benchmark suite are not delivered by this
campaign.

The detailed [validation report](../docs/validation.md),
[machine-readable gate matrix](evidence/completion_gates.json),
[model-family documentation](../docs/model_families.md),
[validity limits](../docs/validity.md),
[performance report](../docs/performance.md) and
[source registry](../src/qstars/data/sources.json) specify the evidence and scope.
The [fresh spatial comparison](../benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/previous_campaign_comparison.json)
and [spatial campaign](../benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/campaign.json)
retain the exact rerun evidence.
