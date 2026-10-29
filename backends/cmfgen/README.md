# CMFGEN external backend

Selected in Report 13 for component integration, **not yet production certified**.
QStarS retains its structure, evolution, conservation ledger, wind eigenvalue,
matching and validity decisions. CMFGEN supplies spherical frequency-dependent
transfer and statistical populations for trial material states. A specified-flow
reference run does not select a wind and is never accepted as a production mdot.

Upstream: <https://sites.pitt.edu/~hillier/web/CMFGEN.htm>. Use the official
18 June 2025 source, 21 June 2023 atomic archive and supplied reference models.
Download hashes and provenance are in
`Reports/evidence/backend_20260928/cmfgen_acquisition.json`.

Source and atomic archives are kept in ignored `sources/cache`; local builds in
ignored `backends/cmfgen/build*`. The source archive includes a GPLv3 LICENSE. Atomic files retain their
individual authorship and provenance. QStarS does not vendor or relicense the
archive. Public availability, successful compilation, numerical convergence and
scientific validation are separate gates.

## Reproduce the component experiment

Use the project's `.venv` (install its `audit` extras for acquisition). The build
below assumes the existing MESA SDK gfortran installation, but does not link the
stellar MESA solver into CMFGEN. It builds an external executable.

```sh
.venv/bin/python backends/cmfgen/acquire.py source atomic models --manifest acquisition.json
mkdir -p backends/cmfgen/build_local sources/cache/atomic_local sources/cache/models_local
tar -xzf sources/cache/cur_cmf_18jun25.tar.gz -C backends/cmfgen/build_local
tar -xzf sources/cache/atomic_data_21jun23.tar.gz -C sources/cache/atomic_local
tar -xzf sources/cache/three_recent_models.tar.gz -C sources/cache/models_local
cp backends/cmfgen/Makefile_definitions.apple backends/cmfgen/build_local/cur_cmf/Makefile_definitions
.venv/bin/python backends/cmfgen/apply_export.py backends/cmfgen/build_local/cur_cmf
.venv/bin/python backends/cmfgen/apply_ray_capacity.py backends/cmfgen/build_local/cur_cmf
.venv/bin/python backends/cmfgen/apply_population_export.py backends/cmfgen/build_local/cur_cmf
.venv/bin/python backends/cmfgen/apply_hminus_rayleigh.py backends/cmfgen/build_local/cur_cmf
.venv/bin/python backends/cmfgen/apply_hminus_ff_repair.py backends/cmfgen/build_local/cur_cmf backends/cmfgen
```

Initialize the installed SDK as documented in `docs/mesa.md`; then set the
absolute extracted source path and build:

```sh
export QSTARS_CMFGEN_ROOT="$PWD/backends/cmfgen/build_local/cur_cmf"
bash backends/cmfgen/build.sh
.venv/bin/python backends/cmfgen/record_build.py "$QSTARS_CMFGEN_ROOT"
.venv/bin/python backends/cmfgen/prepare_hhe.py \
  sources/cache/models_local/sample_grid/zpup_like_model results/new_hhe_component \
  --atomic sources/cache/atomic_local/atomic --iterations 25
touch results/new_hhe_component/QSTARS_EXPORT
.venv/bin/python backends/cmfgen/run_component.py results/new_hhe_component \
  --executable "$QSTARS_CMFGEN_ROOT/exe/cmfgen_dev.exe" --timeout 2400
.venv/bin/python scripts/validate_cmfgen_component.py \
  results/new_hhe_component results/new_hhe_audit --hhe
```

Preparation **must finish before execution starts**. Every destination must be
new. Acquisition verifies pinned SHA256 as well as byte ranges and fails if the
upstream release has changed. The wrapper raises the stack limit to the host
hard limit, records native `/usr/bin/time`, caps runtime, checks fresh complete
outputs, and preserves errors even when Fortran returns zero. A completed run
is still explicitly `completed_component_not_certified`.

The `continue_component.py` command copies a completed trial into a new directory
and requests full Jacobian iterations without changing physical parameters.
Native acceleration and final-iteration matrix-reuse behavior remain active.
Retain its source state. A changed iteration path is not a tolerance relaxation.

`Makefile_definitions.apple` is tested only on arm64/macOS with Accelerate.
Linux must use a checked gfortran/LAPACK configuration and repeat validation;
no Linux execution or scaling measurement is claimed. The spectral stream
requires little-endian 32-bit integers and IEEE 64-bit reals. An incompatible
build is rejected by header/length checks.

The read-only source patch and the modified CMFGEN executable retain upstream
GPLv3 terms; see `NOTICE` and `COPYING.CMFGEN`. General QStarS boundary input
construction, non-LTE material derivatives and production wind dispatch remain
gated. Consult [the scientific contract](../../docs/nlte_backend.md), not just
these commands, before interpreting output.

The upstream relativistic-ray module has a fixed 500-entry array. A 680-ray
request was observed to corrupt the adjacent allocation descriptor and abort.
The separate capacity patch increases that array to 2048 and rejects larger
requests before indexing. The wrapper caps unrecognized builds at 500 rays;
the build record must match the executable hash to permit the guarded capacity.
This extension was tested on the stated non-dust H/He controls; gamma-ray and
dust capacities are unchanged and are outside the extended domain.

## Cool H/He campaign (Report 14)

Use all four patches above before building. The population patch repairs fixed
ion indices in RVTJ when H0 is enabled. The Rayleigh patch corrects stale HMI
naming, the HI array stride and the derivative equation index. Unlike the
read-only export repair, it corrects the affected physical opacity/linearization.
`record_build.py` records those source hashes and sets `hminus_rayleigh_fixed`
only for the repaired source and a newer executable. The runner rejects
H0 plus `INC_RAY=T` without a matching capability record.

```sh
.venv/bin/python scripts/check_cmfgen_rayleigh.py "$QSTARS_CMFGEN_ROOT" results/new_rayleigh_probe
.venv/bin/python backends/cmfgen/prepare_cool.py \
  sources/cache/models_local/sample_grid/zpup_like_model results/new_cool \
  --atomic sources/cache/atomic_local/atomic --case T9 --iterations 60
.venv/bin/python backends/cmfgen/run_component.py results/new_cool \
  --executable "$QSTARS_CMFGEN_ROOT/exe/cmfgen_dev.exe" --timeout 2400
# Only after native completion; this records failed checks without certification:
.venv/bin/python scripts/audit_cmfgen_cool.py results/new_cool results/new_cool_audit
```

`--case T6` selects the other declared control; `--factor 2` geometrically bisects
every original radial interval with the same analytic velocity law. The fresh
`--initialization full` research trial encountered a missing native BA matrix;
its failed inputs are preserved, and it is not a recommended restart recipe.
Completed bootstrap controls can use `continue_component.py`; inspect the
actual native iteration history because native safeguards can still choose
Lambda iterations. These commands reproduce experiments, not a proven cool
convergence path.

The [cool closure contract](../../docs/nlte_cool_closure.md) states the material
energy/derivative constraint and why none of these prescribed-flow controls
selects mdot or validates a convective handoff. See [Report 14](../../Reports/Report_14.md)
for results and retained failures.

## Report 15: H-minus free-free repair and cool numerical recovery

The FF patch preserves the supplied table, interpolation and coefficients. It
initializes the table from either opacity or cooling, fixes last-depth opacity
selection, supplies the missing continuum derivatives, and makes the EHB path
use the same neutral-H ground population. It also retains all levels in native
MAXCH stopping. The native probe compiles the actual source routines; bound-free
terms are zeroed **only in that isolated FF regression**, never in a scientific
run. With the installed compiler environment:

```sh
.venv/bin/python scripts/check_cmfgen_hminus_ff.py "$QSTARS_CMFGEN_ROOT" \
  results/new_ff_probe sources/cache/atomic_local/atomic/HYD/MI/13apr15/hmi_ff
.venv/bin/python backends/cmfgen/recover_cool.py \
  sources/cache/models_local/sample_grid/zpup_like_model results/new_cool_recovery \
  sources/cache/atomic_local/atomic "$QSTARS_CMFGEN_ROOT/exe/cmfgen_dev.exe" --factor 2
```

The second command preserves separate bootstrap, fixed-J/T initialization and
released trials. A stage failure stops the chain with its recorded result;
none is a selected wind. `continue_component.py` can also copy a failed trial's
restart state for explicitly labelled diagnostics; having MOD_SUM does not
certify that checkpoint. No failure directory is overwritten.

Optional `apply_moment_trace.py ROOT` and
`apply_growth_trace.py ROOT OUTPUT` add read-only failure captures. Apply them
before rebuilding and recording executable capabilities. They do not change
operators, clipping, physical coefficients or native stopping. The latter
requires the former. Build 6 has moment tracing, build 7 adds the FF/stopping
repair, and build 8 adds growth tracing. Each executable and modified source
is hashed in Report 15's evidence. Do not assign a newer runner/hash to an
older run. Original upstream GPLv3 and atomic provenance remain applicable.

`validate_cmfgen_component.py` now checks the largest raw and returned
corrections. A native MAXCH can ignore ten levels in an unpatched build.
The runner rejects explicit moment/intensity failures even if native STOP
returns zero. Continuations retain inherited OUTGEN content; new runs save
its initial copy, and the campaign report excludes only a hash-verified prefix.
[Report 15](../../Reports/Report_15.md) documents the remaining numerical
blocker. A native component completion is never production certification.

## Report 16: rejected shared-face moment experiment

`apply_node_moment.py ROOT OUTPUT` requires `apply_moment_trace.py` first.
It reconstructs the K/N frequency contribution using distinct left/right node
weights at one shared spatial face and updates the same native VAR operator.
The default project physics is not changed. This is a rejected candidate;
Report 16's full-cycle request with CHK_J/H/XM=NONE still produces negative J.
Do not use it as a production wind backend or accept the fixed-factor fixture
as a relaxed non-LTE response.

With the compiler environment initialized:

```sh
.venv/bin/python backends/cmfgen/apply_moment_trace.py "$QSTARS_CMFGEN_ROOT"
.venv/bin/python backends/cmfgen/apply_node_moment.py "$QSTARS_CMFGEN_ROOT" results/new_node_patch
bash backends/cmfgen/build.sh
.venv/bin/python backends/cmfgen/record_build.py "$QSTARS_CMFGEN_ROOT"
.venv/bin/python scripts/check_cmfgen_moment_response.py "$QSTARS_CMFGEN_ROOT" results/new_response_probe
```

The probe links actual MOM/VAR/native variation sources. It tests nonlocal
eta/chi derivatives over two frequencies with coherent scattering, relativity
and radiation advection. Formal angular factors are held solely for that
bounded linearization fixture. A missing compiler or failed native fixture
leaves its directory intact. Linux linkage is a replay option, not tested
Linux agreement. The exact read-only seed/edge and rejected N/H patches are in
Report 16's evidence; their original and changed source hashes are explicit.
See [Report 16](../../Reports/Report_16.md) for the remaining coupled angular
realizability/relaxed-response requirement. Native diagnostics and stale
outputs still fail closed; no changed guard is retroactively assigned to a run.

## Report 19: analytic formal-transfer tangent component

`formal_tangent.f90` differentiates the build-17 `QS_FRM_ALI` map on its current
native geometry, including optical-depth/source weights, limited log
interpolation, thick outer extrapolation, coherent scattering and blue-side
ray history. `tangent_map` returns J/H/K/N and their directional changes;
`tangent_scattering` solves the implicit scattering tangent, with either a
dense analytic reference or matrix-free GMRES; `tangent_rsqh` differentiates
the native midpoint r²H interpolation. This CMFGEN adaptation carries GPLv3
provenance; the license text is in `formal_tangent.LICENSE`.

Call after the primal formal/scattering solve at the current frequency. Supply
thermal emissivity to `tangent_scattering`, total extinction, scattering,
their directional changes, and independent B/DBB changes. Include a scattering
change in total dCHI when appropriate. The `prev_p`/`prev_m` arrays hold the
previous frequency's derivatives and stay fixed through the current solve.
They have shape `(maximum native ray length, NP)`; initialize them to zero
at the first frequency and advance them to `out_p`/`out_m` after closure.
Native geometry and coefficient shapes are checked. All geometry and velocity
are held fixed; DIFFUSION and ZERO_FLUX are supported, HOLLOW is rejected.
The component does not modify the native primal ray buffers.

With the compiler environment initialized and recorded formal-ALI libraries:

```sh
.venv/bin/python scripts/audit_formal_tangent.py \
  --cmf-root backends/cmfgen/build_17/cur_cmf --nd 20 --output results/new_tangent_nd20
.venv/bin/python scripts/audit_formal_tangent.py \
  --cmf-root backends/cmfgen/build_17/cur_cmf --nd 39 --output results/new_tangent_nd39
```

The oracle reruns native transfer, rather than finite-differencing the tangent
implementation. The audit retains three standard step sizes and additional
face-limiter probes in the strong-absorption case. Its declared gate-step table
resolves weak signals and stays within the near-zero-H interpolation branch;
it does not select the lowest measured error. The complete larger-step and
roundoff failures are retained. [Report 19](../../Reports/Report_19.md) explains
the limits and remaining native integration work.

This is a component library/probe, not build 18 or a replacement for
`VAR_JREL_V5`. TX/TVX population and continuum/line-history integration remains
open. Neither these fixture grids nor their fast radiation-only products
certify a cool release, a full-model Jacobian, a wind response or an mdot.

## Report 20: exact Jacobian action of the cool residual (builds 18–24)

Builds 18–24 layer read-only hooks on build 17 (patches `jv.patch`,
`jv_replay.patch`, `jv_env.patch`, `jv_fast.patch`, `jv_rays_omp.patch`,
`jv_warm.patch`, `jv_dirs_omp.patch`; each has a manifest recorded by
`record_build.py`). With all of them, one native pass
can return, at a state x, the exact residual R(x), the fixed-J blocks and the
radiation part of J.(x*s) for requested fractional directions s:

| Variable | Effect |
|---|---|
| `QSTARS_JV_DIRS` | int32 `nv,nt,nd` + float64 `s(nt,nd,nv)` (column-major). Writes `QSTARS_JV_RAD.bin`: the exact radiation part of the rate/RE derivative (analytic nonlocal formal tangent with blue-ray history, native dSTEQ/dJ coefficients) |
| `QSTARS_BA_DUMP` | every assembled full block after equation selection, FIXPOP and POPS column scaling (`QSTARS_BA_TRACE.bin`) |
| `QSTARS_REPLACE_FORCE` / trace | fixes the ground-state/ionization equation selection, so the residual definition does not change with the BA |
| `QSTARS_RJ_STORE` / `QSTARS_JV_REPLAY` | store the converged J per frequency; a replay pass rereads it and does one formal solve per frequency (no ALI). The replay residual equals the base residual exactly |
| `QSTARS_JV_STOP` | stop after the solve dump |
| `QSTARS_ALI_WARM` (build 23) | start the formal-ALI from a stored converged J of the same frequency grid; changes only the iteration count (converged J differs by ~1e-9 relative) |
| `QSTARS_JV_EXPORT`, `QSTARS_JV_DUAL`, `QSTARS_JV_CHECKLIN` | per-frequency dchi/deta/dscat/dJ export; dual reference tangent; per-frequency linearized-vs-dual check |

Run with `QS_ALI_LIN=T`, `QS_LIN_FAC=0` so the dumped blocks contain no radiation
response. **The native analytic fixed-J blocks are incomplete**: they omit the
dependence of level dissolution and occupation probabilities on n_e, the total
ion density and T. The complete direct Jacobian replaces the n_e and T columns
by frozen-J central differences (replay passes that reread the base J) and adds
to every ion column the common total-ion-density term measured on the H II
column. `scripts/cool_newton_krylov.py` implements this; `scripts/audit_exact_jv.py`
validates the complete action against exact-residual oracles.

Build 24 (`jv_dirs_omp.patch`) runs the direction loop in parallel only when there are
at least as many directions as threads, so that build 22's ray-parallel `lin_apply` is
actually used for single-direction products (the always-open direction region had
disabled it). Products are identical to build 23 (1.4e-16), but the pass is not faster:
the per-call ray regions are too fine-grained (279 s system time on a loaded machine;
`results/jv/b24_check`). The Newton–Krylov runs use build 23.

**Builds 25–27 (Report 21): level-dissolution derivatives.** The native `VCHI/VETA` hold
the occupation probabilities W and the dissolved fraction T2 constant. Build 25 adds
d(W) in the LTE populations of every bound-free term and d(T2) below each edge. Build 26
moves those terms from `CONT_FREQ` to the actual frequency, as `COMP_VAR_OPAC` does for
`VCHI/VETA`. Build 27 adds the line-opacity ratio `L_STAR_RATIO = W(up)/W(low)`. All
chain through `B/A/X_LEV_DIS(n_e, POPION, T)` as written in `COMP_LEV_DIS_BLK`; W from
the formula equals the native W to 2.2e-16. `QSTARS_JV_NODIS` disables the terms and
reproduces build 23 exactly. With them, the RE row of the complete action matches
T-released exact-residual oracles at the oracle floor (n_e direction: 0.25 → 1.0e-5).

**Build 28 (Report 21): no nested OpenMP team in the tangent.** With fewer directions than
threads, build 22's direction loop ran as an inactive parallel region (`IF` false), so
`lin_apply`'s ray-parallel region became a nested team. On macOS libgomp this intermittently
deadlocked (all threads in team barriers; stack sample in
`Reports/evidence/report21/omp_deadlock_sample.txt`) or crashed with SIGSEGV. `jv_nonest.patch`
moves the loop body into an internal routine and calls it outside any OpenMP construct in that
case. The arithmetic is unchanged: products equal build 27 to 1.2e-16, with much less system time.

**Builds 29–31 (Report 22): solver speed, same arithmetic.**
- *Build 29* (`jv_noobs.patch`): with `QSTARS_NO_OBS_FORMAL` set, the observer-frame formal
  solution (`CMF_FORM_SOL_V2` in `COMP_J_BLANK`, both branches) is skipped. Its only product is
  the observed spectrum (OBSFLUX). The residual, BA, stored J and tangent never read it. It was
  78 % of a column pass. The Newton–Krylov driver sets the variable for its passes;
  verification/final runs do not.
- *Build 30* (`jv_asig.patch`): the per-frequency scattering operator A (ND×ND, the J response
  to a unit scattering emissivity at each node) is built once per frequency with batched sweeps
  (`lin_em_build`, `lin_apply_multi`) and shared by every direction's GMRES, which then
  multiplies by A instead of making a ray sweep per Krylov step. Same map, algorithm,
  preconditioner and tolerance. Active when there are at least `QSTARS_JV_ASIG_MIN` (default 4)
  directions; `QSTARS_JV_NOASIG` disables it.
- *Build 31* (`jv_batch.patch`): with A available, the directions run in one block per thread.
  Each block makes ONE forcing sweep and ONE final sweep for all its directions
  (`lin_scattering_multi`). `lin_setup`'s ray loop runs in parallel; the primal moments are still
  summed serially in ray order. Operator column blocks are sized to fill twice the thread count.
  `QSTARS_JV_NOBATCH` / `QSTARS_JV_NOSETUP_OMP` restore build 30.

Products agree with build 28 per equation row (each row on its own scale) to ≤ 1e-12, residuals
to ≤ 3e-14. See Report 22 and `Reports/evidence/report22/`.

**Build 32 (Report 23): consistent two-photon rate in single-iteration passes** (`two_phot_j.patch`).
Every QStarS pass is a single CMFGEN iteration with `COMP_F=T`, so there is no EDDFACTOR from a
previous iteration. `COMP_OPAC` then silently replaces the intended two-photon method
`TWO_METH=USE_RAD` by `OLD_DEFAULT`, whose partner-photon occupation is exp(−hν_B/kT_e) instead of
J(ν_B)/(2hν_B³/c²) (OUTGEN: "Using OLD_DEFAULT method for two-photon decay"). At the H I front of
the 167-radius cool model this overestimates 1s→2s two-photon absorption 33–101 times, is not
energy-consistent between the RE row and the rate rows, and supplied about 75 % of the RE row.

| Variable | Effect |
|---|---|
| `QSTARS_TWO_PHOT_J` | path of a QStarS stored-J file (`QSTARS_RJ.bin`, written with `QSTARS_RJ_STORE`) on this model's radius and frequency grid. `USE_RAD` is then kept in `COMP_F=T` passes and `TWO_PHOT_OPAC_V3` takes the conjugate-frequency J from that file (`qs_two_phot_get_j`) instead of EDDFACTOR, with the statements of `GET_J_FOR_TWO_PHOT`: linear interpolation in frequency between the bracketing records, the first record above the grid, Rayleigh–Jeans (ν²) below it |

- The file is **fixed (lagged) data**, the equivalent of the previous-iteration EDDFACTOR: the occupation
  `PHOT_OC_TWO` does not depend on the state of the pass. As in native CMFGEN it is evaluated at
  ν_0 − `CONT_FREQ` when the continuum cross-sections are recomputed and reused for the frequencies of that band.
  With the variable set, the stored file is used whenever the method is `USE_RAD` (also if `COMP_F=F`).
- It is loaded and verified once, before any opacity is evaluated (`qs_two_phot_j_init`): `QSRJ0001` header,
  ND, NCF (from the file size), every frequency against the model grid to 1e-14 (as the replay reader), J > 0.
  A missing or mismatched file, `INC_TWO=F`, or a `TWO_METH` other than `USE_RAD` stops the run with a
  `QSTARS_TWO_PHOT_J_FAILURE` file and an OUTGEN error. OUTGEN names the method and the file
  ("Using USE_RAD method for two-photon decay" / "USE_RAD conjugate-frequency J: QStarS stored-J file").
- **Tangent.** With the occupation fixed, the exact Jacobian-vector tangent needs no new term: its two-photon
  section follows the module `TWO_METHOD` of the residual, so the non-`OLD_DEFAULT` branch applies (rate rows:
  dJ(ν)·(1+Y_B) down and dJ(ν)·Y_B up; opacity/emissivity through the native `VCHI/VETA` of
  `TWO_PHOT_VAR_OPAC`; RE row through `CHI_NOSCAT`). `USE_RAD` has no explicit temperature dependence, unlike
  `OLD_DEFAULT`, whose exp(−hν_B/kT_e) is absent from the native `VCHI` and from the tangent's dχ (the frozen-J
  T column of the driver covers only the fixed-J part). `QSTARS_JV_STATS` records the method and the number of
  conjugate-J evaluations. The dependence on J(ν_B) of a self-consistent (same-pass) two-photon J is not
  linearized: it couples each frequency to a lower one that the blue-to-red tangent sweep has not reached.
- Unset, the build is build 31: residual, BA, solution and tangent products are byte-identical (one thread).
- Driver: `Runner.run(..., two_phot_j=PATH)` in `scripts/cool_newton_krylov.py` sets the variable and
  `TWO_METH USE_RAD` for that pass; every pass that shares one residual definition (base, frozen-J columns,
  products, trials) needs the same file. It refuses a missing file, a runner `TWO_METH` other than `USE_RAD`,
  and an executable whose capability record lacks `qstars_two_phot_stored_j` (builds ≤ 31 ignore the variable
  and would silently run `OLD_DEFAULT`). A `QSTARS_TWO_PHOT_J` entry of `extra_env` is the same request: it
  gets the same checks and VADAT setting, and a value that differs from `two_phot_j` is an error (it cannot
  replace the checked file). A pass stopped by the native check raises with the text of
  `QSTARS_TWO_PHOT_J_FAILURE`. `scripts/replay_with_vadat.py --two-phot-j` does the same for a diagnostic
  replay. `record_build.py` flag: `qstars_two_phot_stored_j`.

Validation outputs: `results/two_photon_test/b32/` (`a_result.json`, `b_result.json`, `c_result.json`,
`fail_result.json`; review round 2: `fail_r2_result.json`, `rerun_r2_result.json`, `d_result.json`,
`e_result.json`, `e_bands_result.json`, `e_onoff_result.json`; with the scripts that made them and a README).

Known limit of the tangent, older than build 32 and exposed by its T-direction check (rows 5–6, H I n = 4, 5: 1e-2
at the 84-radius x0 state with the feature on, 2e-4–3e-4 off): the temperature derivative of the infrared continuum
emissivity taken from native `VETA` is 7–9 % too small at 4 µm (26 % below 0.03e15 Hz, exact above 0.36e15 Hz).
It matters where Brα is inverted and the `SRCE_CHK` floor `CHI = 0.1*ETA*CHI_NOSCAT/ETA_CONT` applies, and more with
`USE_RAD`, whose two-photon absorption raises `CHI_NOSCAT` in the infrared 30–400 times at the temperature minimum.
The two-photon rows (2, 3, NT) were affected too, less visibly: build 34 reduced their T-direction errors at depths 14–22
(row 2: 6.1e-5 → 1.5e-5 against the Richardson oracle; rows 3 and NT within that oracle's scatter, and closer against a
third exact-residual oracle, see build 34). See items d and e of that README. Cause and fix: build 34 below (the ion
free-free Gaunt factor's temperature dependence was missing from the tangent).

**Build 33 (Report 23): fast local-ALI linearization pass** (`fast_ali.patch`, capability flag
`qstars_fast_local_ali`). A speed-up of an existing code path; no number changes.

The trapped-level preconditioner needs one frozen-J replay with `QS_ALI_LIN=T`, `QS_LIN_FAC=1` (the local
formal-solution ALI response of `VAR_JREL_V5`: `TX(i,i,2)=Psi`, `TX(i,i,1)=-Psi*S`). Up to build 32 that pass did the
work of a `QS_LIN_FAC=0` pass plus two sweeps over the full `TX`/`TVX(ND,ND,NM)` arrays every frequency
(2 × 135 MB at ND = 167, NM = 606). Profile of a one-thread build-32 pass (167 radii, 43,314 frequencies; seven 10-s
`sample` windows, `results/fast_ali_b33/profile`): 70.6 % (66.7–77.2 % per window) in the QS_FRM_ALI non-finite check
`ANY(.NOT.(ABS(TX).LE.HUGE))`, `ANY(... TVX ...)` at the end of `VAR_JREL_V5`, 22.9 % (17.2–27.2 %) zeroing `TX`/`TVX`
in its QS_ALI_LIN block, 6.5 % everything else (formal solution 2.3 %, `DO_VAR_CONT` loops 1.5 %, BA updates 1.1 %,
`EDD_JREL_VAR_V3` 0.3 %, moment set-up and THOMAS < 0.1 %).

What changed (`plane/var_jrel_v5.f`, `new_main/mod_subs/do_var_cont.f`). Active with `QS_ALI_LIN`, `QS_FRM_ALI`,
`QS_LIN_FAC` ≠ 0, `CONT_VEL`, not `ACCURATE`, `USE_J_REL`, spherical, `SECTION='CONTINUUM'`; otherwise the full path.
- The invariant used: on this branch `VAR_JREL_V5` zeroes `TX` and `TVX` every frequency and sets only `TX(i,i,1:2)`,
  and every later statement of `DO_VAR_CONT` combines elements with equal depth indices. `TX` is therefore zero off its
  diagonal and `TVX` zero everywhere, at every frequency (the hash diagnostic below finds no non-zero word off the
  diagonal of `TX` and none in `TVX` on the full path).
- `VAR_JREL_V5` (flag `QS_ALI_DIAG` in `COMMON/QS_ALI_FAST/`, set by `DO_VAR_CONT` for this call only) does not
  assemble or solve the moment operator whose results the local response discards (tridiagonal set-up, THOMAS,
  advective re-solve, `EDD_JREL_VAR_V3`), zeroes only the diagonals of the matrices selected by `DO_THIS_TX_MATRIX`, and
  checks for non-finite values only the entries that can be non-zero. Executed as before: the INIT block (allocation,
  shared module state, unit 171), J = the stored formal solution with its blow-up and `J_CHK` checks, `GET_RSQH_REL`,
  and the statements that define the response and the `QSTARS_ALI_LIN_Q` statistics.
- `DO_VAR_CONT` skips the rescaling of the previous frequency's matrices (overwritten by the local response), applies
  the clumping, source-only, 3–5, line (in `SIM_INDX` order) and `dSLdT` statements for K = J only (serial loops), does
  not touch `TVX`, and copies `dJ_LOC` from the local band only. `VJ`, the diffusion terms and all BA updates are
  unchanged code. The optimized GIMPLE of every moved statement has the same operations, operand order and fused
  multiply-adds as the full path; every other routine of the executable is instruction-identical to build 32 (apart
  from relocated data offsets).
- A slot freed by `SET_LINE_OPAC`'s forced removal (all `MAX_SIM` line slots taken; it decrements the use count without
  zeroing) keeps its last diagonal until it is reassigned and zeroed by `INIT_LINE_OPAC_VAR_V2`; nothing reads an unused
  slot (the full path zeroes it at the next frequency). Every matrix in use, `TVX`, `dJ_LOC` and `VJ` are the full
  path's numbers at every frequency. In a pass where such a removal happens, the `TX` hash of `QSTARS_ALI_LIN_HASH`
  (and the `TX` count of `QSTARS_VAR_NONFINITE`), which cover unused slots too, can differ between the paths while every
  output agrees; compare outputs there. (The comment above the fast zeroing in `var_jrel_v5.f`, "zeroed when it was
  released", does not hold for this removal; left unchanged to keep the build.)
- Non-finite states: if a factor of the `DO_VAR_CONT` statements is not finite (`CLUMP_FAC`, `-ETA_CLUMP/CHI_CLUMP`
  with `CHI_CLUMP=0`, `LINE_PROF_SIM`, `L_STAR_RATIO`, `LINE_OPAC_CON`, ...), the full path turns the zero off-diagonal
  `TX` entries into NaN (0·Inf) and spreads them through `dJ_LOC`/`VJ` into the BA blocks of neighbouring depths; the
  fast path never touches those entries, so only the local block carries the NaN. The `VAR_JREL_V5` non-finite check
  runs before these statements and catches neither. All finite states are bitwise identical.

| Variable | Effect |
|---|---|
| `QSTARS_ALI_LIN_SLOW` | non-empty: the full path of builds ≤ 32 (in-build reference) |
| `QSTARS_ALI_LIN_HASH=n` | at the first and every n-th frequency, order-sensitive 64-bit hashes of the bit patterns of `TX`, `TVX`, `dJ_LOC`, `VJ` with their non-zero counts (and that of `TX` off its diagonal) to `QSTARS_ALI_LIN_HASH.bin`; diagnostic, either path, one sweep of the arrays per hashed frequency |

A `QS_LIN_FAC` ≠ 0 pass writes `QSTARS_ALI_LIN_PATH` (the path that ran). `QS_LIN_FAC=0` passes take the
build-19/21 branch and execute no new arithmetic (only branch tests and a one-time environment read).

Executable sha256 `6d72dbadb10e3e785751b9c3686e2db681a38ba46136fe7a95583cf2ca7977d6` holds only for the executable
linked in `backends/cmfgen/build_33`: ld embeds the absolute path and timestamp of each `libsubs.a` member (18 `N_OSO`
debug-map entries), an `LC_UUID` and a code signature. For a rebuild elsewhere compare the path-independent fingerprints
in the `reproducibility` block of `fast_ali_patch_manifest.json` (sha256 of every Mach-O section and of `var_jrel_v5.o`,
`do_var_cont.o`; `macho_sections.py`). A rebuild of build 32 + `fast_ali.patch` in a scratch directory had another file
sha256 and all section and object hashes equal (`results/fast_ali_b33/out/r2_fingerprints.json`).

Validation (`results/fast_ali_b33`: README, `out/*.json`; frozen-J replays with the BA dump, one thread):

| check | result |
|---|---|
| `QS_LIN_FAC=1`, 167 radii (pin state it05): build 32 full path vs build 33 fast path | residual, BA, solution, `QSTARS_ALI_LIN_Q`, `TWO_PHOT_SUM` byte-identical |
| same, 84 radii (x0): build 32 full, build 33 full (`QSTARS_ALI_LIN_SLOW`), build 33 fast | all byte-identical; `TX`, `TVX`, `dJ_LOC`, `VJ` hashes of the build-33 passes identical at every 10th frequency (4,332 records) |
| same with the stored two-photon J (`USE_RAD`); with `DO_CL`+`INCL_dSLdT`; with `SRCE_ONLY`+`INCL_dSLdT` (84 radii, full vs fast) | byte-identical, hashes identical |
| `QS_LIN_FAC=0`, 167 radii: build 32 (`two_photon_test/b32/a_r2/a167_b32`) vs build 33 | byte-identical (identical inputs) |
| `QS_LIN_FAC=0` product pass, 84 radii, 3 directions: build 32 (`b32/a/a84_b32`) vs build 33 | residual, BA, solution, `QSTARS_JV_RAD` byte-identical |
| two-photon stored J, 167 radii, fast | runs; residual differs from the `OLD_DEFAULT` pass only in the two-photon level rows (2, 3, 33, 35, 103) and RE, and matches the build-32 two-photon pass to 3e-14 |

With two or more threads no pass is bitwise reproducible, on either path: the formal solution sums J over rays in a
`SCHEDULE(DYNAMIC)` reduction (`cmf_formal_rel_v4.f`), and the local response uses that J. At two threads the same
fast pass run twice differs by 3.5e-13 (residual, of the row maximum) and 2.5e-14 (BA, of the row maximum within a
block); full vs fast differ by the same amounts (3.7e-13, 2.9e-14; 1.4e-18 of each block maximum).

Timing, 167 radii, one thread, back to back on the shared machine (`/usr/bin/time`):

| pass | real s | user s | instructions | cycles |
|---|---|---|---|---|
| `QS_LIN_FAC=1`, build 32 (full path) | 2583 | 2382 | 12.33e12 | 4.20e12 |
| `QS_LIN_FAC=1`, build 33 (fast path) | 97.6 | 91.8 | 0.694e12 | 0.208e12 |
| `QS_LIN_FAC=0`, build 33 | 94.5 | 80.9 | 0.560e12 | 0.183e12 |

The fast path does about 18–20× less CPU work: 17.8× fewer instructions, 20.2× fewer cycles. The 26× wall-time ratio
is a single contended observation and not like-for-like. The slow pass ran for 43 min while the load average moved
between 9 and 25, it was also the profiled pass (seven 10-s `sample` windows), and it ran at a lower effective clock
(1.76 GHz against 2.26 GHz). On one thread the fast pass took 74–182 s of wall time depending on load, at a constant
0.693–0.696e12 instructions. On two threads at low load it took 73.5 s (`QS_LIN_FAC=0`: 78.3 s). The gain grows with
ND: at 84 radii it is 5.7× in instructions. The fast pass still costs 11 s of CPU more than a `QS_LIN_FAC=0` pass. That
cost is the unchanged `VJ` loop of `DO_VAR_CONT` (15 % of the fast pass), two thirds of which adds zeros on the bands
that `dJ_LOC` leaves at zero. Restricting the loop would give the same `VJ`; this was not done. Memory use is unchanged
(`TX`/`TVX` are still allocated).

Driver: `Runner.run(..., lin_fac=1.0)` with this executable, or from the command line (the task's
`--vadat QS_LIN_FAC 1` form works too):

    nice -n 10 .venv/bin/python scripts/replay_with_vadat.py BASE OUT --x0 X0 --exe build_33/.../cmfgen_dev.exe \
        --lin-fac 1 [--two-phot-j base] --depths 33 37 --threads 1

For the pin state `nk167_pin/theta1.00/it05_trial_l1`, both runs (with and without `--two-phot-j base`) are
byte-identical to the validated build-33 passes (`results/fast_ali_b33/cli`). The run without it is therefore also
byte-identical to the build-32 slow pass. The two-photon fast path was checked against the slow path at 84 radii.

**Build 34: Gaunt-factor temperature derivative in the tangent** (`ff_gaunt_dT.patch`, capability flag
`qstars_ff_gaunt_dT`; `new_main/qs_jv.f90` only).

Cause of the infrared deficit of the T column (build 32 known limit above). The residual (`GENOPAETA_V10`, from
`COMP_OPAC`) evaluates the ion free-free Gaunt factor g = `GFF_VEC(CONT_FREQ, T, Z)` at the current temperature. The native
variation routine `VAR_OP_V10` (`VCHI/VETA`, from `COMP_VAR_OPAC`), which the tangent uses, differentiates the free-free
opacity and emissivity of H II, He II and He III with g held fixed. g is bilinear in log10(hν/kT) and log10(Z²Ry/kT). At
the temperature minimum of the 84-radius x0, dln g/dlnT = 0.11–0.22 below 0.05e15 Hz, where ion free-free makes 35–57 % of
the emissivity. The missing term is CHIFF Z² n_e N_ion (dg/dT)/√T ×
{(1 − e^(−hν/kT))/ν³, 2h/c² e^(−hν/kT)}. A split by process (`results/ir_dT_b34/split_processes.py`) shows that this is
the whole deficit. H-minus free-free (table slope in `DO_H0_FF_VAR` since build 7), bound-free with dissolution,
two-photon (`USE_RAD`, fixed stored-J occupation), stimulated emission, the `SRCE_CHK` floor and scattering were already
exact. `COMP_FREE_FREE_V2` has the same fixed-g derivative, but only the electron heat balance calls it (`USE_EHB`/`COMP_EHB`),
and no QStarS pass sets either.

What changed: `qs_gff_dt_vec` uses the table, the statements and the cell of `GFF_VEC`. g equals `GFF_VEC` to the bit
(checked in every pass). It returns the slope within the cell and zero where an argument is clamped. `qs_jv_gff_freq`
evaluates the term at each continuum evaluation frequency for the ions of `COMP_VAR_OPAC`: Z > 0, `PHOT_ID` 1, and
`IMP_VAR(EQION)` or `QSTARS_JV_VALL`. `qs_jv_freq` adds it, multiplied by the direction's dT and moved to FL like the
dissolution terms, to dCHI/dETA, but only where the direction has a T component. Native `VCHI/VETA`, the BA, the residual
and every product of a direction without a T component are unchanged.

| Variable | Effect |
|---|---|
| `QSTARS_JV_NOGFF` | non-empty: no Gaunt term (the build-33 tangent) |

`QSTARS_JV_STATS` has a new last line `ff_gaunt_dT`, with the max |g − `GFF_VEC`|, the number of evaluations and ions,
and the number of free-free resonances (`FF_RES_GAUNT`), which are not differentiated; there are none in H/He.

Executable sha256 `d3e227737611c3bc473eab3510668921951f83a8140cb37272f254d1179829d2` (this directory only). The section
hashes and the `qs_jv.o` hash are in `ff_gaunt_dT_patch_manifest.json`. Every library member and every other `new_main`
object is byte-identical to build 33. Tests: `tests/test_cmfgen_ff_gaunt.py`.

Validation (`results/ir_dT_b34`: README, `out/*.json`):

| check | build 33 | build 34 |
|---|---|---|
| 84 radii, x0, T at depths 16–20: tangent/FD of dη/dlnT, median by band <0.01 / 0.01–0.03 / 0.03–0.07 / 0.07–0.09 / 0.09–0.2 (1e15 Hz) | 0.40 / 0.75 / 0.88 / 0.93 / 0.97 | 1.000000 in every band (p10–p90) |
| same, Br-α floored cells, dχ (depths 17, 18) | 0.930 / 0.910 | 1.000000 |
| 84 radii, complete action vs exact-residual central differences (Richardson), T direction, depths 14–22, rows 5 / 6 (H I n = 4, 5) | 1.06e-2 / 1.32e-2 | 6.9e-8 / 1.1e-7 (oracle 5.2e-8 / 7.9e-8) |
| same, all 131 active rows: median error; rows above 3× their oracle consistency | 1.3e-4; 111 | 6.6e-8 (oracle 7.1e-8); 2 |
| same, rows 3 / NT against a third exact-residual oracle (central, eps 5e-4; review passes `r2/rev_T_c*`) | 2.8e-3 / 2.1e-2 | 7.6e-4 / 1.6e-2 (oracle scatter 1.3e-3 / 2.1e-2) |
| same pair, per-frequency dJ of the tangent against exact-J central differences, depths 16–20, median tangent/FD below 0.01 / 0.01–0.03 / 0.03–0.07 / 0.07–0.09 (1e15 Hz) | 0.07 / −0.13 / 0.27 / 0.73 | 1.000000 in each |
| 167 radii, `nk167_trap/free_ali/it02_trial_l0.5`, T at depths 80–90 (31–39): tangent/FD of dη/dlnT, median below 0.01 / 0.03–0.07 / 0.09–0.2 | 1.97 / 1.54 / 1.31 (0.45 / 0.85 / 0.97) | 1.000000 in every band (p10–p90) |
| same, RE row of the complete action vs exact-residual central differences (2e-3, 1e-3, Richardson), depths 80–90; all 131 rows above 3× their oracle | 3.5e-2; 119 | 1.2e-5 (oracle 1.2e-5); 4 |
| neutrality, one thread, 84 radii, feature on and off | — | residual, BA, solution byte-identical; n=2 / n=1 products bit-identical |

Row 3 (H I n = 2) against the Richardson oracle reads as a regression (4.9e-4 → 2.5e-3, oracle 1.5e-3): that oracle's
FD_a and FD_b differ by factors at depths 21–22. Against the third oracle, build 34 is the closer one (table); its noise
is in the Ly-α band (2.3–3.3e15 Hz), where the exact-J forward and backward differences disagree by a factor of 2.4.

Left as they were: the native BA's `VCHI/VETA` (by design); `DIFFW(NT)` at the diffusion inner boundary, which comes
from the same `VAR_OP_V10` and is used by the tangent only for components at depth ND (added in build 35); and
`COMP_FREE_FREE_V2`. At the nodes of the Gaunt table the residual is not differentiable in T;
the tangent then takes the slope of the selected cell. At depths 80–90 of the 167-radius state, build 33's error in the RE row's T column was
3.5 %. That column is a 15–20× cancellation between the direct part and the radiation part, and the radiation part was
0.19 % off. `results/nk167_trap/free_ali` ran build 33 (`exe_b33`); `free_b34` (started 14 Oct 05:34) runs a copy of
this executable (`exe_b34`). Whether the exact T column removes the oscillating deep T steps had not been tested when
this was written.

**Build 35: the Gaunt term in the diffusion boundary, and a free-free resonance guard** (`ff_gaunt_dT_diffw.patch`,
capability flag `qstars_ff_gaunt_dT_diffw`; `new_main/qs_jv.f90`, three calls added to `new_main/cmfgen_sub.f`). Review of
build 34. Build 34 itself is left unchanged, because `results/nk167_trap/free_b34` runs a copy of it.

- **Diffusion inner boundary.** `DIFFW(NT)` = d(DTDR)/dT(ND) is summed in the DTDR section of `cmfgen_sub.f` from
  T3·`VCHI(NT,ND)`/`CHI(ND)`. Its `VCHI` comes from the same `VAR_OP_V10`, so g is held fixed there too. The tangent
  uses `DIFFW` in dDTDR = `DIFFW`·dPOP(:,ND) for the boundary term of every direction with a component at depth ND.
  `QS_JV_DTDR_GFF` is called three times in that section: at the start, after each `DIFFW` increment (with that
  statement's T3 and `CHI(ND)`, passed as expressions, i.e. copies), and after the scaling by DTDR²/T1. It adds up
  Σ T3·dχ_g(FL, ND)/`CHI(ND)` with build 34's dχ_g, dη_g at depth ND, moved to FL as `COMP_VAR_OPAC` moves `VCHI`, and
  scales the sum as `DIFFW` is scaled. The routine reads the state and writes only `qs_jv_state`. `DIFFW` itself is
  not changed. The tangent adds the term times dPOP(NT, ND) to dDTDR, only where dT(ND) ≠ 0, in both `dir_tangent`
  and `dir_batch`. Ions are those of build 34's term (`IMP_VAR(EQION)`, as `VCHI`; every ion with `QSTARS_JV_VALL`).
  `QSTARS_JV_NOGFF` turns it off along with the build-34 term.
- **Free-free resonances.** `FF_RES_GAUNT` (added to g in `GENOPAETA_V10` for `ION_LEV` = 1) is not differentiated.
  A product pass with a T direction and such a resonance in an ion of the tangent now stops with `QSTARS_JV_FAILURE`
  and an OUTGEN error. Before, it returned a T column that was wrong without any warning. The pass then leaves no
  residual, and the driver raises. `QSTARS_JV_NOGFF` (the build-33 tangent) bypasses the check. There are no
  resonances in the H/He models (counted 0).
- `QSTARS_JV_STATS` has a new last line, `ff_gaunt_dT_diffw`. It records DTDR, the native `DIFFW(NT)`, the Gaunt term
  for `IMP_VAR` ions and for all ions, the number of frequencies, and the number of resonances.

| Variable | Effect |
|---|---|
| `QSTARS_JV_NOGFF` | non-empty: neither Gaunt term (the build-33 tangent); also bypasses the resonance stop |

Executable sha256 `c69836d31fd0425b93ea0cf529e988cdd356ad9a6ea50d0c6949492355a3a6dd` (this directory only). Section hashes
and the `cmfgen_sub.o`/`qs_jv.o` hashes are in `ff_gaunt_dT_diffw_patch_manifest.json`. Every library member and
`cmfgen.o`, `formal_tangent.o`, `formal_tangent_lin.o` are byte-identical to build 34. The disassembly of `cmfgen_sub.o`
has the same ordered sequence of floating-point instructions as build 34 (766), 3 more calls, and a different register
and stack allocation. Tests: `tests/test_cmfgen_ff_gaunt_diffw.py`.

Validation (`results/ir_dT_b34`: README section 8, `out/r2_result.json`, `out/r2_oracle84.json`):

| check | build 34 | build 35 |
|---|---|---|
| d(DTDR)/dT(ND) used by the tangent vs the residual's DTDR (OUTGEN), frozen-J central difference, s_T = 1 at every depth, eps 1e-4: 84 radii (x0) / 167 radii (`free_ali/it02_trial_l0.5`) | +6.97e-4 / +6.97e-4 | +7.3e-6 / +7.3e-6 |
| 167 radii, s_T at every depth (finding 2: depth ~120 and the inner boundary, also in build 34): tangent/FD of dη/dlnT, median and p10–p90 in every band below 7.6e15 Hz, every depth group (1–30, 31–39, 40–79, 80–90, 91–110, 111–130, 131–160, 161–166, 167); dχ/dlnT median below 2.3e15 Hz | same (its dχ, dη exports equal build 35's to 1.6e-15, two threads) | 1.000000 (dη p10 0.999999 at 31–39, 2.3–7.6e15 Hz); 1.000000 |
| same, product change b35 − b34 (the DIFFW term only): largest change / max of the row's complete action, at depth ND; RE row | — | 5.5e-6; 5.9e-8 (falls to 4e-7 at depth 150, 4e-10 at 120) |
| neutrality, 84 radii, one thread, feature on: residual, BA, solution, REPLACE trace, ALI stats, native and tangent spectral exports, products of all three directions (none has a component at ND) | — | byte-identical to build 34 |
| same, feature off (`OLD_DEFAULT` two-photon; the inputs of `fast_ali_b33/r84/prod_b33`) | — | residual, BA, solution, REPLACE trace, ALI stats, products byte-identical to build 34 |
| 84 radii, rows 2 / 3 / 5 / 6 / NT of the T direction (depths 16–20) against the third exact-residual oracle | 1.6e-6 / 7.6e-4 / 1.6e-7 / 1.8e-7 / 1.6e-2 | identical (no component at ND) |

The 7.3e-6 left in d(DTDR)/dT(ND) is not the Gaunt term. Weighting the tangent's own material derivative dχ(ND)/dT (s_T
export) as the DTDR section does reproduces the finite difference to 2.1e-7. So the remainder is a T term that the native
`VCHI(NT,ND)` of the DTDR section lacks and the tangent's continuum loop has. Candidates are the level-dissolution T
terms (build 25) and the line-ratio T terms (the DTDR section's line variation "ignores the effect of a temperature
variation"). It was not identified further: scaled by the product change above, it moves the action by about 6e-8 of a
row's maximum. At the frequencies of Gaunt-table nodes the residual is not differentiable in T, as in build 34.

**Temperature release in single-iteration passes.** `FIX_T=F` alone does not release T
on a restart: `READ_BA_DATA_V3` restores `FIXED_T` from the restart's `BAMATPNT`
("Temperature was held fixed"), and `FIXPOP_IN_BA_V3` then makes the RE row the
identity at all depths. Multi-iteration runs clear it after the first iteration;
single-pass designs (oracles, Newton–Krylov passes) do not. Remove `BAMATPNT`/`BAMAT`
(`scripts/cool_newton_krylov.py:release_t`). Check that row NT of
`QSTARS_RESID_TRACE.bin` is nonzero.

Pitfalls found while building this (all preserved under `results/jv`):
`newunit=` returns negative unit numbers; an environment value longer than the
receiving buffer returns status −1 (present); `ndarray.tofile` always writes C
order — use `qstars.nlte_newton.write_fortran`. NumPy 2 scalar `repr` is not a
Fortran-readable number.

Scope: fixed R, V, geometry and frequency grid; coherent scattering; DIFFUSION
inner boundary. Two-photon: without `QSTARS_TWO_PHOT_J` the passes run CMFGEN's
`OLD_DEFAULT` substitution (builds 18–31 always did; the tangent follows that branch,
without the temperature derivative of its exp(−hν_B/kT_e) factor in dχ). With it
(build 32) the conjugate-frequency occupation is fixed stored-J data and is held
fixed in the tangent.
Velocity/geometry derivatives needed for a wind eigenvalue are not implemented.
