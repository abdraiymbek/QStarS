# MESA infrastructure and supported interface

QStarS can call the **MESA 24.08.1 EOS and radiative opacity libraries directly**
for fixed H1/He4 composition (`Z=0`). It keeps its own finite-volume evolution
and independent boundary ledger. This is not yet the MESA stellar-evolution
backend, a production atmosphere, or a validated moving cavity.

Bridge ABI 5 additionally exposes native pressure inversion, Eddington
atmosphere references, full Henyey/ML1 MLT with automatic derivatives, and
radiative opacity using an explicitly supplied external EOS electron state.
ABI 3 and 4 libraries remain readable for their earlier EOS/opacity calls;
ABI 4 also supports the native atmosphere and pressure inversion. Native MLT
requires ABI 5 because exploratory ABI 4 builds had different argument layouts. The
new reference atmosphere has demonstrated validity failures; see
[the boundary audit](atmosphere.md). Evolution still uses the existing
screened envelope transport, not the native MLT reference automatically.

A pre-existing ARM SDK 24.7.1 and MESA 24.08.1 installation were found in the
earlier runtime. The initial development report's absence claim was wrong;
Report 5 records the correction. The existing installation is reused without
modification. Its source commit is
`5c166bf174b477823a2ca0031c4f57210c8e67a8`; the interface build records all local
Git changes. PGSTAR, GYRE, and ADIPLS were disabled in that installation. No
pulsation capability is claimed.

## Build and run

First install the pinned MESA/SDK using the
[official instructions](https://docs.mesastar.org/en/24.08.1/installation.html).
Then, from this repository, with the Python environment activated:

```sh
python scripts/build_mesa_microphysics.py \
  --mesa-dir /absolute/path/to/mesa-24.08.1 \
  --sdk-dir /absolute/path/to/mesasdk \
  --output backends/mesa/build_local
qstars run backends/mesa/build_local/example.yaml
python scripts/validate_mesa_microphysics.py \
  --library backends/mesa/build_local/libqstars_mesa.dylib \
  --mesa-dir /absolute/path/to/mesa-24.08.1
QSTARS_MESA_LIBRARY="$PWD/backends/mesa/build_local/libqstars_mesa.dylib" \
  MESA_DIR=/absolute/path/to/mesa-24.08.1 python -m pytest -q
```

The builder creates a new directory and refuses to overwrite failed builds.
It generates an example with the actual paths. The current build recipe is
validated on Apple Silicon macOS; the ordinary Python solver is portable,
but a Linux shared-library build has not yet been validated. No shell profile
is changed. The macOS SDK build uses deployment target 14.0 unless explicitly
overridden, with that choice recorded. Compiled binaries are local build
artifacts, not distributed under the QStarS license.

## Numerical and physical domain

The interface passes exactly the same composition and electron state to
`eosDT_get` and `kap_get`. It uses MESA's radiation constant and EOS derivatives.
The radiative opacity handle disables electron conduction: the separate
radiation-plus-conduction handle measures the omitted conductive contribution.
Runs reject an omitted contribution exceeding 1% of diffusive heat transport.
This avoids treating conductive opacity as photon optical depth or radiation
force. The initial interface experiment used the combined opacity; those
results are preserved and superseded by ABI 3.

The bounding query rectangle is density `1e-14..1e3 g/cm³` and temperature
`3000..1e8 K`. **This rectangle is not a declaration of validity at every
point.** Every accepted state is checked for a first-law defect below 1% using
symmetric finite differences, and the neutral ideal EOS backstop is rejected.
The first-law test normalizes its residual by `|dU/dlnX| + P/rho` to avoid an
artificial divergence from cancellation. It is a screening limit, not proof of
percent-level evolutionary accuracy. Existing diffusion, optical-depth, Mach,
compactness, geometry, and whole-energy checks remain active.

An additional evolution screen checks gas pressure against the dilute
classical limit where T>=40000 K, eta<-10 and the charge-2 Coulomb bound is
below 0.01. It rejects discrepancies above 1%. The retained test at rho=1e-9,
T=1e6 has an 8.28% gas-pressure discrepancy but a 1.04e-8 total first-law
defect. No table value is adjusted. Low-density neutral fallback and the
table-loader electron floor are separately audited in Report 8.

The 153-point audit matches a separate Fortran process exactly. Five points
exceed the 1% first-law threshold; the maximum is about 0.315 at density
`1e-4 g/cm³`, T≈6812 K, in an OPAL/SCVH–FreeEOS blend. These are retained
counterexamples. MESA documents the difficulty of combining inconsistent
sources in its [EOS overview](https://docs.mesastar.org/en/24.08.1/eos/overview.html).
No coefficients or table entries were adjusted to make the test pass.

The saved dependency inventory hashes 8,135 files totaling 15,541,148,301 bytes
in the installed composition, EOS, ionization and opacity data directories.
This identifies files present in those directories, not a claim that every
file was used by a particular model. Its compressed file manifest and directory
scope are retained in the
[inventory evidence](../benchmarks/results/2026-09-28-closure/20260928T030927.417120Z_mesa_data_inventory/inventory.json).
MESA may use its existing table caches. The inventory and source audit do not
assert that a library execution can never update a cache.

A fresh final-source 128-zone, ten-year QStarS envelope run completes with
energy residual 3.45e-11 of throughput. The separate MacBook benchmark records
31.293 seconds and 1045.4 MiB peak RSS. This is a direct-microphysics envelope
run, not a native MESA stellar-evolution track; see the
[validation scope](validation.md) and [performance report](performance.md).

## Upstream stellar benchmarks

`scripts/run_mesa_reference.py` builds isolated ordinary MESA or upstream
QUEST work directories and records input hashes, compiler logs, explicit
changes, histories, and termination logs. `--saved-seed` uses the upstream
provided stellar model; `--startup-cgs-fix` explicitly corrects the startup
assignment of a solar-mass number to CGS `M_center`. The original and repaired
runs are distinct. A process exit or a short MESA run is not full published
replication, and MESA's energy diagnostic is not the QStarS external ledger.

Remaining production work includes validated atmosphere/convection coupling,
conservative cavity evolution, wind-bearing published tracks, independent
full-system accounting for an actual stellar backend, and cross-validation.
[Upstream QUEST](https://github.com/andysantarelli/MESA-QUEST) is pinned at
`671eb9b898b28c0d02a6a2455dacc0f5b269e280`.
