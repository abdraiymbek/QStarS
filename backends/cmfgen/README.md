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
