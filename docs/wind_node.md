# CMFGEN wind node contract and the `selected_nodes` evolution wind

This page defines the file that a converged, self-selected CMFGEN wind (Gate 4)
must deliver. It also defines how the evolution consumes that file (Gate 5).
The software interface is implemented and tested with a **synthetic** node
(`tests/test_node_wind.py`). No CMFGEN-derived node exists yet. Nothing on this
page certifies a physical wind or a coupled quasi-star track.

## `wind_node.json` (schema `qstars_wind_node`, version 1)

All scalars are positive, finite CGS numbers. The field table, the validation
flags and the provenance keys live in one place, `qstars.wind_node`.
`json_schema()` emits them as a draft 2020-12 JSON Schema.
`validate()` additionally enforces the cross-field identities listed below.

| Field | Meaning |
|---|---|
| `M_total_g` | Total mass inside the base. CMFGEN gravity is a point mass. |
| `R_base_cm` | CMFGEN inner boundary, which is the envelope's outer face |
| `L_base_erg_s` | Lab radiative luminosity through the base, 4πR²F |
| `rho_base`, `T_base` | Base density and gas temperature |
| `P_gas_base`, `P_rad_base` | Base gas pressure and lab radiation pressure P_rr (K moment) |
| `mdot_g_s` | Self-selected mass-loss eigenvalue. It is never an input. |
| `v_inf_cm_s` | Terminal speed |
| `L_emergent_erg_s` | Escaping lab luminosity at the outer boundary |
| `R_phot_cm`, `Teff_emergent_K` | Rosseland τ=2/3 radius, and Teff ≡ (L_emergent/(4πσR_phot²))^¼ |
| `energy_gauge` | Gas energy zero, e.g. `neutral_H_He_ground` (`cmfgen_backend.hhe_material`) |
| `base_face` | Every `outer_interface.InterfaceState` field at the base, including `gas_temperature` |
| `validation` | Booleans, one per release-boundary gate ([mapping](#validation-flags-and-the-release-boundary)). They are the producer's claims, not verified evidence. |
| `provenance` | `source_run`, plus sha256 for `cmfgen_executable`, `vadat`, `atomic_data` and `solution` |
| `node_id`, `schema`, `schema_version`, optional `notes` | Identification |

Identities are enforced to a relative 1e−10, because the fields come from one
double-precision export:

- The duplicated base values equal the corresponding `base_face` fields: radius, mass, density, gas pressure, P_rr and temperature.
- The potential is −GM/R.
- 4πR²F = L_base.
- 4πR²ρv = mdot.
- Teff follows from L_emergent and R_phot.
- R_phot ≥ R_base.
- `base_face.energy_gauge` equals `energy_gauge`.

Unknown keys are rejected.

### Validation flags and the release boundary

The flag list is a contract proposal, pending Gate 4 and user approval. It
gives one flag to each gate of the [release boundary](nlte_backend.md#validation-sequence-and-release-boundary)
and to the Gate 1 audits. Two requirements are checked by the evolution run
itself, not by a flag.

| Release-boundary gate | Where it is checked |
|---|---|
| Atomic data | `atomic_data_validated` |
| Radiation force | `radiation_force_validated` |
| Energy: total material plus radiation energy, and the non-LTE energy-zero/composition bridge and material derivatives used in the wind equations | `material_energy_validated` |
| Gate 1 electron, charge and radiation-energy audits | `electron_charge_radiation_audits` |
| Nonlinear convergence | `nonlinear_converged` |
| Radius, frequency and atom size | `radius_converged`, `frequency_converged`, `atom_size_converged` |
| Export invariance | `export_invariant` |
| Sonic, photospheric, escape and independent ledger checks | `sonic_regular`, `photosphere_resolved`, `escapes`, `independent_ledger` |
| Mach-one MLT guard during mdot selection | `mach_one_mlt_guard` |
| Attachment to a valid envelope | No flag. The run applies `validity.check` (Mach-one MLT guard, Pe ≥ 100, Kn, MESA first-law screens) at the start and to every accepted or remapped state. |
| Evolution EOS energy zero against the node gauge | No flag. The start-of-run base-face check refuses a gauge mismatch. The MESA gauge is undeclared. |

Every run reports `termination.json["wind_node"]["validated"]: false`,
whatever the flags say. `node_flags_all_true` reports the flags as claims. A
`validated: true` needs a Gate 4 evidence verifier, which does not exist yet.
It would replace each flag with an evidence path plus sha256 and check on load
that the files exist and the hashes match.

A node set `{schema: qstars_wind_node_set, schema_version: 1, nodes: [...]}`
lists nodes **in track order**. Nodes in a set share one gauge, have unique
ids, and consecutive nodes have distinct base data.

## Evolution use

```yaml
wind:
  model: selected_nodes
  node_path: wind_node.json   # resolved relative to the YAML
  drift_tolerance: 0.01       # required; 0 < delta < 1
```

`rate_msun_per_year` and `terminal_speed_cm_s` must stay at 0, because the node
supplies both. The run copies the node to `wind_node.json` in the run directory
and records its source path and sha256 in `external_inputs.json`, as it does for
opacity tables. A restart re-copies the checkpoint's node.

- **Base data.** These come from the state as (M, R, L, ρ, T): `enclosed[-1]`, `r[-1]`, the reported outer luminosity, and the outer cell's ρ and T. The cell-centre against face offset is a discretisation term and converges with resolution.
- **One node.** The node values are used only while every |Δln x| ≤ δ. Otherwise the run raises `DomainError` (`model_domain_exit`).
- **A node sequence.**
  - The state is projected (Euclidean in ln space, clamped to each segment) onto the polyline of consecutive nodes.
  - Every output is log-linear along the segment, which keeps the Teff identity exact.
  - The state must lie within δ (max norm) of that point.
  - At the end nodes, values are held within δ. Values are never extrapolated, and no functional form is fitted.
- **Explicit lagged coupling.** mdot and v_inf come from the old state, like every other wind mode. `mass_transport.exchange` uses the node v_inf for the escape debit removed·max(0, v_inf²/2 − h − φ − v²/2). The ledger is otherwise unchanged.
- **Accepted states stay inside the neighbourhood.** A candidate step or remap that leaves it is rejected and retried, as for any domain check. The run ends as `model_domain_exit` when the drift comes within 0.01% of δ, which is the same margin rule as the radial Mach exit.
- **Gauge and base-face check at start.** compare_faces(evolution face, node base face) is evaluated at the node's base state. The evolution face uses the evolution EOS gas P and u with LTE radiation E = aT⁴, P = E/3, and the node's flux, mdot and potential.
  - A gauge mismatch is rejected (`model_domain_exit`).
  - The MESA EOS has no declared gauge (FreeEOS uses ground-state H₂), so it is refused until an energy-zero bridge exists. The ideal-gas gauge is `fully_ionized_H_He_at_rest`.
  - Flux defects go to `wind_node_interface.json` and `termination.json["wind_node"]`. They are recorded and not enforced.
- **Per-step diagnostics in `events.jsonl`.** These are not ledger terms.
  - `interface_energy_defect` = (L_out − wind_cost/dt − L_emergent)/L_base.
  - `wind_power_defect` = (wind_cost/dt − (L_base − L_emergent))/L_base.
  - Both depend on the energy gauge and on the discrete outer-cell h and φ.
- **History.** `photospheric_radius_cm` stays NaN in every mode. The node's held or interpolated R_phot belongs to the node's base data, not to the current state. It goes to `wind_node_photospheric_radius_cm`. The node mode also adds the terminal speed, the emergent L and Teff, and the drift.

## Validity boundary and what Gate 5 still needs

- **A validated Gate 4 node.** All `validation` flags must be true and must be backed by archived, hash-checked evidence. Until a verifier exists, every node, synthetic or not, is reported `validated: false`.
- **An admissible QStarS envelope at the node's base data.** The envelope must pass the Mach-one MLT guard, Pe ≥ 100, Kn and the MESA first-law screens. The T9 control scales (M = 10⁶ M☉, R = 2.9×10¹⁵ cm, T_base ~ 10⁴ K) need the MESA EOS, and MESA's energy gauge is undeclared. Coupling therefore also needs an explicit, tested FreeEOS/MESA ↔ `neutral_H_He_ground` energy-zero bridge. Gauge shifts are not applied silently.
- **A boundary condition.** The envelope keeps its own grey or luminosity outer boundary and constant `pressure_cgs`. The node does not impose CMFGEN's base T–L–P relation. The face and energy defects quantify that inconsistency. Thresholds for them must be pre-registered before results.
- **Numerical certification at the node.** This needs mesh, timestep, δ versus δ/2, remap and bitwise-restart studies. The tests here cover only the mechanics: per-step mass and energy residual < 1e−12, transactional retry, bitwise restart, drift and hull exits, remap rejection, self-declared flags never validating, and gauge rejection.
- **Physics scope.** One node certifies a neighbourhood, not a track. Full-track coupling needs further CMFGEN nodes added along the track. Published evolutionary tracks use their own wind laws. No node or interpolation may be adjusted to match them.
