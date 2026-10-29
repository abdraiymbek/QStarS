# QStarS — Quasi-Star Simulation

A working, modular 1D simulation development framework with finite thermal
evolution, self-gravity, gas and radiation pressure, radiation and convection,
BH growth, separate feeding/wind streams, and conservative energy accounting.

**Status: 0.1 development software. The requested scientifically validated
general-purpose QStarS v1 is NOT complete.** The executable evolves restricted envelopes with either an analytic hot EOS
or direct, audited MESA EOS and opacity calls. Published static structure benchmarks and numerical
conservation tests are implemented. Full published stellar evolution, cool
photospheres, moving cavities, and coupled transonic winds remain open validation gates. A separate
Fiacconi–Rossi stationary wind solver now includes radiation enthalpy advection.
These examples do not establish quasi-star lifetimes or QSRD recovery results.

[Report 24](Reports/Report_24.md) re-examines the stalled stage-1 (fixed-T) closure. **Still not converged; no gate
changes; not v1.** The solves of Report 23 broke exact conservation: the per-component step limiter clips populations
separately, so hydrogen, helium and charge conservation drifted (He up to 8×, charge up to 34×) and the converged deep
interior was damaged while the hydrogen fronts improved. The least squares used deep-layer directions through their
far-field action at the hydrogen front, which the nearly singular local-ALI blocks amplify, and moved T in the fixed-T
stage. Full trial steps failed natively because they made supercritical infrared masers (Λ spectral radius 1.06), not
because of the transfer solver. New: `--conserve` (trial states projected onto the constraints), `--pc-svd-floor`,
`--freeze-t`, and CMFGEN build 36 (Newton solve of the scattering fixed point where the local ALI stalls).
With them the stage-1 solve, restarted from its first accepted state, keeps conservation to 1e-15 and accepted three
full Newton steps (fixed-measure residual 224 → 128). At the fourth, full steps again made supercritical masers and the
solve stalled near 128. A second, row-scaled measure is not monotone, and each step degrades the deep interior. The narrowest remaining
blocker is the ill-conditioned hydrogen ionization front, through which the least squares fits the front with
deep-layer directions.

[Report 23](Reports/Report_23.md) is the requested v1 closure run. **Released cool closure NOT CONVERGED; Gates 1–5 not
reached; Gate 6 (Ball et al. 2011 track) FAIL; not v1.** The closure's non-convergence is numerical, not a missing
equilibrium. The solver measured progress with fixed-J Jacobian blocks, which hide the error of every optically thick
transition (H I n = 1–9, He I and He II resonance levels), and every single-iteration pass ran CMFGEN's `OLD_DEFAULT`
two-photon approximation, three quarters of the radiative-equilibrium row at the H I front. The temperature continuation
of Reports 21–22 therefore heated the cool layer by 3–6 kK while its ionization stayed ten orders of magnitude from
Saha. Build 32 restores the intended two-photon physics; build 33 makes a local-ALI preconditioner for every column
affordable (17.8× fewer instructions, identical output); with them a fixed-T solve reduces the statistical-equilibrium
measure 3× in two iterates. Build 34 corrects the temperature column of the exact action (free-free Gaunt factor);
on it the solves accept full Newton steps and halve the residual, then stall where the local-ALI preconditioner is
nearly singular at the ionization fronts. No closure has converged. The registered Ball et al. 2011 replication evolves
conservatively for 1 255 steps to M_BH = 46 M☉, where its smoke gate is a recorded numerical blocker in the
super-Eddington outer hydrogen opacity peak; its required check Q2 fails at row 0.

[Report 22](Reports/Report_22.md) profiles and speeds up the 167-radius cool-closure solve
without changing the model, the products (≤ 1e-12 per row) or the convergence test. CMFGEN builds
29–31 skip the observer-frame formal solution in solver passes, share one explicit scattering
operator per frequency across tangent directions, and batch the direction sweeps. Base passes are
2.5× faster, Jv products 1.9×, column passes 3.4×. Measured on the same state, the stalled Krylov
residual is a weighting conflict with ultra-trace He levels, not a missing T subspace. With the
user's approval, the resumed solve uses ultra-trace weight 1e-10 (26 instead of 74 products per
iterate). **Released cool closure NOT CONVERGED; not v1.**

[Report 21](Reports/Report_21.md) records **RE-ROW JACOBIAN AND INC_AD/INC_ADV PASS;
released cool closure NOT CONVERGED; not v1**. Level-dissolution derivatives complete the
radiative-equilibrium row of the exact Jacobian action, which now matches the T-released oracle
at its consistency in every direction. Adiabatic cooling and advection are restored and
validated. CMFGEN build 28 fixes a nested-OpenMP deadlock in the threaded tangent. T continuation
reaches the full residual at 84 and 167 radii. There, free Newton drives the H I front layer
toward T → 0 on both grids, as the measured coupled response predicts. Pseudo-transient
continuation heats it, but the iterates are still far from statistical equilibrium and have not
converged. No audit,
resolution pair, wind or evolution result is claimed. **111 tests pass, 23 optional integration
tests skip.**

[Report 20](Reports/Report_20.md) records **JACOBIAN PASS (RE row conditional);
released cool closure NOT CONVERGED**. The analytic transfer tangent now drives the
exact Jacobian action of CMFGEN's cool population/rate/RE residual. It matches
exact-residual oracles at their noise floor on the real 84-radius state, after
completing the native fixed-J blocks, which omit level-dissolution dependences. A
Newton–Krylov solve reduces the fixed-temperature closure to a 1.4e-3 frozen-J
correction. With temperature truly released (a single-iteration CMFGEN restart had
silently kept T fixed), the 84-radius solve does not converge. The evidence is
consistent with an unresolved cool-layer ionization front between grid points. No audit, resolution pair, wind or
evolution result is claimed. **109 tests pass, 23 optional integration tests skip**
(compiled MESA/FreeEOS libraries not configured; native tangent tests run).

[Report 19](Reports/Report_19.md) records **TRANSFER COMPONENT PASS; cool closure
still blocked**. An analytic formal-transfer tangent carries radiation changes
across depth and frequency and passes the 1e-6 fixture target, including coherent
scattering and face flux. Its native population/rate adapter remains open; no
released cool state, physical resolution pair or selected wind is certified.
**101 tests pass, 23 optional integration tests skip** with the native tangent
compiler initialized. [Report 18](Reports/Report_18.md) retains the measured
failure of the existing population Jacobians and the costly exact-residual
Newton–Krylov experiment.

[Report 16](Reports/Report_16.md) records **NUMERICAL BLOCKER**: signed previous-frequency
H/K/N forcing is reproduced, and experimental shared-face and unclipped N/H
repairs are rejected. Actual native fixed-factor nonlocal response passes a
bounded fixture; cool equilibrium and physical selected wind remain invalid.
**123 tests pass**, zero skips. QStarS remains development software, not v1.

[Report 15](Reports/Report_15.md) records **NUMERICAL BLOCKER**. H-minus free-free
entry, opacity and derivative defects are repaired and tested. Fixed-field
cool populations settle on two meshes, but released radiation/rate closure
still fails near the He I continuum edge. No validated self-selected wind or
physical evolutionary track is claimed. **122 tests pass** with native MESA
and FreeEOS enabled. The complete failed-run and gate evidence is retained.

[Report 14](Reports/Report_14.md) records **BACKEND/PHYSICS BLOCKER** in the cool
H/He recombination regime. H− output and Rayleigh-selection defects are repaired;
the resolved-atom gas energy now has tested population and frozen-state
derivatives. Cool nonlinear/thermal transfer still fails validation, so no
physical envelope selects a production wind. The MLT guard remains active.
**117 tests pass** with MESA and FreeEOS enabled.
[Cool material and coupling contract](docs/nlte_cool_closure.md).

[Report 13](Reports/Report_13.md) records the earlier **BACKEND/PHYSICS BLOCKER**. CMFGEN was
selected after comparing existing methods, built natively on Apple Silicon,
and integrated through an audited spectral/population interface. The 665-point
fixed-state transfer control meets its global radiation-energy target, and a
hot H/He continuation passes 12 component checks. Cool material-energy
derivatives and a global wind eigenvalue remained open.
**108 tests passed** with MESA and FreeEOS enabled. This is component validation,
not a physical quasi-star wind. [Backend contract](docs/nlte_backend.md).

[Report 12](Reports/Report_12.md) records **FAIL** for the self-selected outer
boundary. New thermal-band transfer, conservative lab-radiation interface,
material-energy critical equations and base-selected wind controls are verified;
**97 tests pass**. The coupled H/He population/energy/force closure and convective
handoff remain unresolved. No physical-envelope or coupled-evolution gate was
bypassed. [Equations and evidence](docs/self_selected_outer.md).

[Report 11](Reports/Report_11.md) records **OUTER CLOSURE FAIL**. Spherical grey
angular transfer and a constructed stationary continuum control now pass scoped
validation, including radiation work and outer-boundary refinement. **85 tests
pass.** Recombining gas/radiation closure and physical envelope mass-loss matching
remain unresolved; no coupled evolution or QSRD run was performed.

[Report 10](Reports/Report_10.md) retains the earlier Nakauchi discrepancy,
conservative moving-surface accounting and Mach-one MLT exit. Its attribution
campaign was not repeated in Report 11.

[Report 9](Reports/Report_9.md) adds a spherical gas-sonic wind component with
FreeEOS/MESA recombination, full enthalpy/work fluxes and independent integral
checks. It does not yet pass published escaping-wind replication or supply a
validated coupled outer boundary. [Report 8](Reports/Report_8.md) established
the original release blocker: the converged recombining candidate requires supersonic convection,
and the warmer candidate violates thin-atmosphere geometry. Direct FreeEOS,
native MESA atmosphere/MLT references and new gas-pressure validity checks
are available; they are not certified production atmosphere coupling.

## Install and run

Python 3.11+; tested on Python 3.12 on an Apple M5 Mac with 16 GiB RAM.
No GPU or MESA installation is needed for the development solver.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
qstars run configs/hot_accreting.yaml
qstars plot results/<run-directory>
qstars restart results/<run-directory>/checkpoint.npz --until-years 20
qstars run configs/two_stream_pulse.yaml
qstars run configs/hot_resolved_validation.yaml
qstars sweep configs/hot_accreting.yaml configs/sweep.yaml --workers 2
python -m pytest -q
qstars validate
qstars validate --require-v1
qstars reference benchmarks/coughlin2024_saturated.yaml
qstars reference benchmarks/fiacconi2016_advective_wind.yaml
```

Each run prints its unique result directory. `requirements-lock.txt` records
the measured environment, including optional audit/testing packages.
Archiving finished runs, moving to another machine and the disk budget of long
runs: [docs/operations.md](docs/operations.md).
`validate` produces Markdown/JSON evidence and currently exits **2** for two
retained published cutoff discrepancies. The previously failing concentrated
128-zone case now completes; full v1 gates remain open. `--require-v1` also exits 2
while the full scientific gates remain open. A completed example is not v1 validation.

## Examples

| Configuration | Experiment and scope |
|---|---|
| `hot_accreting.yaml` | Ten-year thermal evolution with prescribed BH rest accretion |
| `two_stream_pulse.yaml` | Feed pulse and independent escaping wind; incoming thermodynamics explicit |
| `closed_relaxation.yaml` | Closed energy budget and small temperature perturbation |
| `bra_luminosity.yaml` | BRA08 luminosity law at a fixed excision |
| `dotan_launch.yaml` | Dotan launch/tiring closure at an assumed sonic boundary |
| `pressure_confined.yaml` | Truncated hot envelope for controlled spatial convergence |
| `hot_resolved_validation.yaml` | 512-zone, ten-year controlled hot-envelope problem with documented internal-profile convergence |

These are labeled experiments, **not replicas of complete published stellar
models**. `benchmarks/manifest.yaml` separates executed structure benchmarks
from pending evolution comparisons.

The solver initializes discrete hydrostatic structure and retains radial inertia
during evolution. A coupled nonlinear solve advances momentum and the first law
with finite heat storage. Runs stop outside their subsonic validity domain.
Engine choices are `bra2008`, `ball2011_printed`, `quest2025_eq11` (also named
`ball2012_eq3_15`), `prescribed`, and `off`. This evolution supports **fixed
excision only**. The moving Bondi cavity exists only on the separate
quasi-static path (`numerics.mechanics: quasi_static`, see
[usage](docs/usage.md)). Its registered Ball et al. 2011 replication (Report 23)
relaxes a unique row 0 at three resolutions under four registered amendments, then
its smoke gate stops at M_BH = 46 M☉ (a recorded numerical blocker in the outer
hydrogen opacity peak); no complete track or campaign exists. Corrected Bondi and
saturated-convection models are separate static reference solvers.

The `exploratory`, `research`, and `validation` presets change numerical
defaults in one engine. Their names do not certify scientific validity.

- [Equations and conservation derivation](docs/equations.md)
- [Literature audit](docs/literature.md)
- [Source registry](docs/sources.md) and [JSON registry](src/qstars/data/sources.json)
- [Configuration and outputs](docs/usage.md)
- [MESA backend decision](docs/mesa.md)
- [Atmosphere physics, counterexamples and reproduction](docs/atmosphere.md)
- [Validation and v1 gates](docs/validation.md)
- [Measured performance](docs/performance.md)
- [Supported model families](docs/model_families.md)
- [Validity domain and limitations](docs/validity.md)

## Reproducibility

Runs never overwrite another run. Each directory stores resolved configuration,
commit and source snapshot, dependency versions, source registry, event log
including rejected steps, history CSV, radial NPZ profiles, an atomic checkpoint,
validity assumptions, and termination reason. No checkpoint uses Python pickle.

`completed_requested_interval`, `control_stop`, `model_domain_exit`, and
`numerical_failure` are distinct. None means physical destruction. The unresolved
photospheric radius is `NaN`; effective temperature is a grey flux diagnostic.
Sequential [Reports](Reports/Report_1.md) and compact evidence are tracked in Git;
complete local runs remain under ignored `results/`.

The earlier [hot-envelope validation report](Reports/Report_7.md) records the
resolved hot-envelope evidence and the remaining failed v1 gates. The
[machine-readable gate matrix](Reports/evidence/completion_gates.json) includes
the original twelve completion criteria and the nine current priorities.

`src/qstars/` separates thermodynamics, opacity, structure, transport, engine,
winds, mass transport, ledger, mesh, integration, validity, I/O, reference models,
validation, and plotting. CPU sweeps use independent processes.

MIT-licensed original implementation. Cached papers and audited MESA-QUEST source
retain their own licenses and are not redistributed.
