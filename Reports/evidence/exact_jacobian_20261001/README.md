# Report 20 evidence: exact Jacobian action and cool-closure solves (2026-10-01/02)

State x0 = `results/cmfgen_ali_lambda84_long_maxlam1e3` (84 radii, 43 314 frequencies, H/He cool control).

| File | Content |
|---|---|
| `build_18..24_capabilities.json` | executable sha256 and source hashes per build (record_build.py) |
| `lin_vs_dual_fixtures.txt` | linearized vs dual tangent, 84 fixture runs (<= 3e-12) |
| `spectra_audit_jv0.json` | dchi/deta/dscat/dJ vs oracle-exported spectra |
| `residual_audit_native_direct_jv1.json` | native BA + radiation part vs oracle (e = 15.3 before completion) |
| `columns_audit.json/.txt` | frozen-J columns vs native blocks; POPION hypothesis (1.5e-9) |
| `exact_jv_audit.json` | complete action vs Report 18 and T-direction oracles (RE row inactive there) |
| `tfree_jv_audit.json`, `run_tfree.py` | complete action with the RE row ACTIVE (T released) |
| `nodis_jv_audit.json`, `run_nodis.py` | same with level dissolution off (diagnostic) |
| `T_inheritance_README.txt` | FIXED_T inherited from BAMATPNT in single-iteration restarts |
| `*_JV_STATS.txt` | per-pass tangent statistics (closure guard, replay reproduction) |
| `precond_probe_it09_ali.json`, `precond_probe_it10_moment.json`, `mpk_probe_it10.json` | rejected preconditioners / multi-preconditioned directions |
| `build24_product_check.json` | build 24 vs 23 product (1.4e-16) |
| `nk_fixT_history.json` | fixed-T Newton-Krylov history (iterates 0-10) |
| `thermal_curves/*.json` | local frozen-J thermal balance curves (d12 start/iterate 4, with INC_AD, hot side; d17; d22) |
| `superseded_*_README.txt` | preserved failed native-integration runs |

Run directories (not copied; each has README.txt): `results/jv/*`, `results/nk_fixT`, `results/nk_freeT`,
`results/nk_freeT_hot`, `results/superseded_T_inherited`, `results/thermal_curve`, `results/nk167`, `results/nk_dfreq`.
