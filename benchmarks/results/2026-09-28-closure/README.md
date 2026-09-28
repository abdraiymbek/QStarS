# Preserved development validation evidence — 2026-09-28 UTC

**This is not a v1 release archive.** See the repository's
[validation report](../../../docs/validation.md) and
[gate matrix](../../../Reports/evidence/completion_gates.json).

`archive_manifest.json` maps each archived file to its original local path,
original byte length and SHA-256. Files marked `compressed` are gzip encoded;
their hashes refer to the decompressed original bytes. The manifest also
records omitted files and a separately generated plot. `archive_verification.json`
records a complete hash/length check. All original `results/` directories remain
intact locally. Histories, source snapshots, ledgers and checkpoints are retained;
intermediate profiles and upstream binaries/full source are omitted. First and
final QStarS profiles are retained. MESA/QUEST source commits and input hashes
identify their separately installed dependencies.

## Main evidence

| Directory prefix | Evidence and scope |
|---|---|
| `030330.332506Z_spatial_validation` | Ten-year, 256–2048-zone hot-envelope mesh campaign; internal-profile norms and independent continuum initialization |
| `030639.520967Z_robustness_validation` | Final 512-zone timestep/tolerance/remap/restart campaign, unchanged-source checks pass |
| `031253.242612Z_time_method_validation` | One-year midpoint/backward-Euler comparison at three timesteps |
| `031651.519611Z_midpoint_512_refined_time` | Ten-year midpoint refinement and matched-method profile comparison |
| `021532.646459Z_mesa_microphysics_validation` | 153-state direct-library audit; five first-law-screen failures retained |
| `025423.972541Z_robustness_validation` | MESA hot timestep/tolerance/restart evidence; remap thresholds do not trigger remapping |
| `030928.629959Z_mesa_hot_current_verification` | Fresh final-source MESA ten-year confirmation |
| `030927.417120Z_mesa_data_inventory` | 8,135-file, 15.54 GB microphysics dependency hash inventory |
| `024813.279478Z_fiacconi_wind_validation` | Stationary wind convergence and partial published-profile comparison |
| `031903.825105Z_validation` | 29 passing CLI checks, two retained printed-value discrepancies; no v1 claim |
| `031938.396825Z_performance_validation` | Eight fresh sequential MacBook benchmark processes |
| `build_abi3`, `build_verified` | MESA bridge build manifests/logs, without compiled binaries |

The actual directory names include the `20260928T` prefix. JSON report paths
inside original results retain their original `results/...` locations for
provenance. For an archived child, resolve the same run-directory suffix under
this archive instead. Do not silently treat a missing original path as an
absent archived result.

## Preserved counterexamples and earlier approaches

`015005`, `020239`, `020447`, `021050`, `021350`, `023422`, `023549`, `023703`,
`024037` and `024324` preserve the stiff solve, remapping, geometry, iteration
limit and mesh investigations. `025107` retains the midpoint mesh campaign,
including the 2048-zone lossless-convection domain exit. A deliberate stop is
identified in its own saved stop context rather than relabeled as a solver pass.

`021202` is an ordinary MESA smoke run. The two `022006` directories contain
the original and explicit startup-CGS-repair QUEST experiments, neither a
complete quasi-star track. The three `02283...` directories contain cool MESA
initial-model domain exits. Their initial polytropes are not relaxed stellar
models, so these are not physical disruption results.

## Source provenance boundary

Earlier campaigns predate the unchanged-source guard and were run during active
development. They must be read with their first snapshot and chronological
reports. A later snapshot from an already running process does not by itself
prove which module bytes executed. Fresh robustness, method-comparison, MESA
confirmation and performance runs explicitly verify unchanged package sources.
The later verified spatial archive provides a fresh confirmation of the main
mesh campaign. Retaining old files preserves evidence; it does not upgrade an
old run's provenance claim.
