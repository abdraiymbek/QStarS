# QStarS — Quasi-Star Simulation

A working, modular 1D simulation development framework with finite thermal
evolution, self-gravity, gas and radiation pressure, radiation and convection,
BH growth, separate feeding/wind streams, and conservative energy accounting.

**Status: 0.1 development software. The requested scientifically validated
general-purpose QStarS v1 is NOT complete.** The executable evolves hot,
fully ionized envelopes. Published static structure benchmarks and numerical
conservation tests are implemented. Full published stellar evolution, cool
photospheres, moving cavities, and transonic winds remain open validation gates.
These examples do not establish quasi-star lifetimes or QSRD recovery results.

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
qstars sweep configs/hot_accreting.yaml configs/sweep.yaml --workers 2
python -m pytest -q
qstars validate
qstars validate --require-v1
qstars reference benchmarks/coughlin2024_saturated.yaml
```

Each run prints its unique result directory. `requirements-lock.txt` records
the measured environment, including optional audit/testing packages.
`validate` produces Markdown/JSON evidence and currently exits **2** for two
published cutoff discrepancies and a concentrated-mesh solver failure. `--require-v1` also exits 2
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

These are labeled experiments, **not replicas of complete published stellar
models**. `benchmarks/manifest.yaml` separates executed structure benchmarks
from pending evolution comparisons.

The solver initializes discrete hydrostatic structure and retains radial inertia
during evolution. A coupled nonlinear solve advances momentum and the first law
with finite heat storage. Runs stop outside their subsonic validity domain.
Engine choices are `bra2008`, `ball2011_printed`, `quest2025_eq11`, `prescribed`,
and `off`. Evolution currently supports **fixed excision only**. Corrected Bondi
and saturated-convection models are separate static reference solvers.

The `exploratory`, `research`, and `validation` presets change numerical
defaults in one engine. Their names do not certify scientific validity.

- [Equations and conservation derivation](docs/equations.md)
- [Literature audit](docs/literature.md)
- [Source registry](docs/sources.md) and [JSON registry](src/qstars/data/sources.json)
- [Configuration and outputs](docs/usage.md)
- [MESA backend decision](docs/mesa.md)
- [Validation and v1 gates](docs/validation.md)
- [Measured performance](docs/performance.md)

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

`src/qstars/` separates thermodynamics, opacity, structure, transport, engine,
winds, mass transport, ledger, mesh, integration, validity, I/O, reference models,
validation, and plotting. CPU sweeps use independent processes.

MIT-licensed original implementation. Cached papers and audited MESA-QUEST source
retain their own licenses and are not redistributed.
