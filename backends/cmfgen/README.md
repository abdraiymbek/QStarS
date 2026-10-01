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
