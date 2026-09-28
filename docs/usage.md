# Configuration and output reference

Configuration is strict YAML; unknown keys and unsupported physics raise errors.
The safe loader accepts scientific notation and reserves only true/false as
booleans; `off` remains a string. Mass rates use
solar masses/year; other physical controls use CGS and kelvin. Initial mass
means envelope mass excluding the BH. The initializer is hydrostatic n=3, not
thermal equilibrium. A Gaussian temperature perturbation may be specified by
amplitude, enclosed-mass center, and width; the ledger starts after perturbation.

Engine coupling, prescribed accretion, feeding rate/temperature, wind rate, and
surface luminosity accept numbers, linearly interpolated tables, half-open
pulses, or sinusoids. Tables hold endpoint values outside their range.
Timesteps land on table/pulse breakpoints. Negative one-way rates are rejected.
Executable expressions in YAML are unsupported.

```yaml
feeding:
  rate_msun_per_year:
    kind: table
    time_years: [0, 1, 2, 3]
    values: [0.01, 0.01, 0.05, 0.01]
  thermodynamics: specified
  density_cgs: 1.0e-7
  temperature_K: 400000
```

Numerics include mesh count/distribution, solver tolerance, timestep bounds,
fractional-change limit, retries, conservation tolerance, and optional refinement
cadence. The adaptive controller limits state changes; it is not a step-doubling
local-error estimator. Preset names do not establish resolution convergence.

`run --stop-after-steps N` creates a controlled checkpoint stop. `restart`
restores all arrays, BH mass, ledger, next timestep, and step count into a new
directory. Only the final time may change; changing physics requires an
explicitly accounted model transition. Sweeps use Cartesian products of dotted
keys, isolated processes, and separate directories. They do not submit cloud jobs.

## Files

| File | Contents |
|---|---|
| `config.yaml` | Fully resolved controls |
| `metadata.json` | Commit, dirty state, version, dependencies, platform, source hash |
| `source_snapshot.zip` | Exact Python code and registry |
| `working_tree.diff` | Tracked changes relative to the commit |
| `source_prescriptions.json` | Sources, assumptions, implementation status |
| `history.csv` | Global state, mass rates, energy/entropy inventory, complete ledger |
| `profile_*.npz` | Radial CGS arrays with units in keys; face fluxes and velocities |
| `checkpoint.npz` | Atomic state plus configuration and ledger; no pickle |
| `events.jsonl` | Accepted/rejected trials and diagnostics |
| `termination.json` | Status, reason, runtime, process peak RSS, residuals |
| `validity_warnings.json` | Physical limitations |
| `opacity_table.npz`, `external_inputs.json` | Copied external opacity data and hash, when used |

Energies are signed and exclude rest energy except for engine conversion.
`luminosity_erg_s` is gross boundary heat flux; the radiated energy ledger
deducts wind mechanical work. Replacement time uses the sum of flow magnitudes,
so balanced feeding/wind still replaces mass. No-flow replacement time is infinity.
`photospheric_radius_cm=NaN` is deliberately undefined and must not be replaced
by the outer radius. Effective temperature is a flux/area diagnostic.

```python
import numpy as np
history = np.genfromtxt("<run>/history.csv", delimiter=",", names=True)
with np.load("<run>/profile_0000010.npz", allow_pickle=False) as profile:
    radius, temperature = profile["r_cm"], profile["temperature_K"]
```

`qstars plot` reads saved outputs and writes evolution, radial-profile, and
temperature spacetime PNGs. It never executes an evolution run.

Opacity NPZ format: increasing `logrho`, increasing `logT`, rectangular
`logkappa` with shape `(len(logrho),len(logT))`, and nonempty string `source`.
All logs are base 10 in CGS/K. Log-space bilinear interpolation forbids
extrapolation. Table paths resolve relative to the input YAML and are copied
into the output. A provenance string alone is not scientific validation.

## Validation campaigns

```sh
python scripts/validate_spatial.py --mesh log_radius_surface --zones 256 512 1024 2048 --end-time 10 --max-dt .02
python scripts/validate_robustness.py --zones 256 --mesh log_radius_surface
python scripts/digitize_wind_reference.py
python scripts/validate_wind.py
```

Wind digitization requires the audited cached primary PDF and the `audit`
extra. Its CSV and calibration/hash metadata are retained in `benchmarks/`.
The two `fiacconi2016_*_wind.yaml` files run through `qstars reference`; they
do not select an evolutionary wind rate. `docs/mesa.md` gives pinned build,
ionization audit and optional MESA regression commands.

Checkpoint schemas 1–6 remain readable. Current checkpoints include cell
volumes, solved face luminosities, reference mass fractions and compensated
radius updates. A successful restart is tested against an uninterrupted run,
including remaps. Process peak RSS is cumulative within one Python process;
compare memory across fresh processes when making performance claims.

The explicit `mechanics` choice is `midpoint` or `backward_euler`. The latter
records numerical kinetic heating; it is not an unaccounted energy sink.
`configs/hot_resolved_validation.yaml` is the controlled pressure-confined
benchmark, with the accuracy boundary given in the validation report.

Source snapshots freeze the package files seen when a process imports QStarS.
Run metadata reports whether the package remained unchanged through completion.
Editing package files during a campaign invalidates that provenance check;
rerun affected validation cases from a fresh process before claiming exact
reproducibility. This guard also covers lazy imports.
