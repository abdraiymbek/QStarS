# V1 plan execution, 2026-11-06: P0, P1, bounded P2, Ball Q2 audit

Execution record of [V1_execution_plan_2026-11-06.md](../../V1_execution_plan_2026-11-06.md) ("P0 → P1 using existing
products → a bounded P2 comparison; start the independent Ball audit without changing its verdict"). No gate status,
registration or verdict was changed. P2 is in [p2/README.md](p2/README.md); the whole session is written up in
[Report 26](../../Report_26.md). **P2 outcome:** a third defect — the Krylov least squares over nearly collinear raw
directions returned coefficients up to 4e14 (product noise became the model); fixed by `--span-tau`. With all fixes the
released-T solve went 2.837e-2 → 2.185e-2 (−23 %, frozen measure) in 4 accepted iterates (~8.9 h), stopped for review.

## P0: trustworthy baseline (done; exit criterion met)

| Item | Result | Artifact |
|---|---|---|
| Rejected-trial restart hazard | History records now carry `best_trial`, `accepted` (None when no trial lowers the merit) and `status` (`accepted`, `stopped_no_decrease`, `rejected_retry_same_base`, `converged`, ...). Only an accepted trial is written to the restart pointer `OUT/accepted_state.json`. The driver checks the start base **before any native pass** and refuses a trial its run did not accept (`--allow-unaccepted-base` overrides); `--resume` starts from the pointer. Legacy records are judged by `merit_new < merit`. | `scripts/cnk_trap.py` (`trial_decision`, `history_accepted`, `check_restart_base`, pointer I/O), `scripts/cool_newton_krylov.py`; tests `test_trial_decision_*`, `test_restart_check_*`, `test_driver_refuses_to_start_from_the_rejected_iterate_16_trial` (iteration 16's eight trials, verbatim) |
| free_r25 restart pointer | `accepted_state.json` → `it16_trial_l0.5` (history.json left as written; README notes added) | `results/nk167_trap/free_r25/` |
| Frozen reference norm vs local norm | `campaign_norms.py freeze/score`: weights and scales of `it16_from_it16_trial_l0.5` frozen (sha256 `914b7cd4…`); scores by depth region, species group, largest cells, unfloored rows, rows significant locally but unweighted in the frozen measure, conservation and RE net/gross. Reproduces the driver's base merit 2.837128e-2 exactly; it04 → it16: 0.029299 → 0.028371 (−3.2 %). | `scripts/campaign_norms.py`, `frozen_norm_it16.npz`, `frozen_scores.json` |
| Baseline manifest + replay | Source (HEAD, uncommitted diff with sha256, driver hashes), executable/capabilities, 944 input hashes (digest), equation map, two-photon J, VADAT, launch script and command, frozen/local residual. **Replay:** the residual-defining `_tp` pass of the accepted state reproduces **bit for bit** (37 s); a fresh full native pass of the state from x0 reproduces it to 1e-13 relative (frozen difference 2.6e-13). ("numerical_failure" is the normal native status of these one-iteration perturbation passes.) | `scripts/baseline_manifest.py`, `baseline_manifest.json` (+ `_source.diff`), `results/nk167_trap/v1_p0_replay/` |
| README / docs / ledger reconciliation | Chronology brought to Report 25, the dense-window run and this plan; descriptive text only, every status/outcome value byte-identical to before (checked against a saved copy). Unresolved status items listed, not changed. | `README.md`, `docs/{validation,nlte_cool_closure,cmfgen_selected_wind,validity,performance}.md`, `Reports/evidence/completion_gates.json` |
| Gate-to-artifact inventory | 31 criteria (12 original completion criteria, 9 priority gates, closure + Gates 1–6, 3 open items) with recorded status, artifacts, governing plan, gap and plan phase; 5 open registration decisions | `gate_inventory.md`, `gate_inventory.json` |

Unresolved ledger items (need the user's decision; not changed): the active gate "Self-selected wind, coupled evolution,
published evolutionary benchmark" still reads NOT_ATTEMPTED although Gate 6 was run and failed; a stale `outcome`
string of `outer_campaign_completion`; no registered stop rule covers `results/nk167_trap`; no registered atom-size
criterion; `qstars validate` hard-codes its open-gate list (`src/qstars/validation.py`).

## P1: why the useful step disappeared (archived products of iterates 4 and 16, then native)

All offline numbers use the stored exact products and the driver's own row weights (`p1/p1_common.py` rebuilds the
iterate's least squares and reproduces the driver's logged fit 0.8111, μ 6.667e-5 and bounded fit 0.998766 exactly).

**1. Step algebra is numerically unstable.** cond(A) = 4.4e26 (it16) / 2.1e26 (it04); the normal equations formed by
`lm_step_bound` have cond ≈ 1e53. At the driver's own μ the augmented-QR solution of the same penalty problem differs
from the driver's step by ~100 %.

| Bounded step (|s| ≤ 0.5, T 0.03), linear fit ‖Ay − b‖/‖b‖ | it04 | it16 |
|---|---|---|
| unbounded least squares (max |s| ~ 1e3) | 0.832 | 0.811 |
| LM-NE (driver, normal equations) | 0.976 | 0.9988 |
| LM-QR (same family, augmented QR) | 0.938 | 0.969 |
| TSVD (column-scaled, cutoff 1e-2) + LM-QR | 0.924 | 0.950 |
| BOX (explicit component constraints, SOCP) | 0.869 | 0.895 |
| BOX + linearized number/charge rows where active | 0.884 | 0.909 |
| BOX on the projected step (BOXP, n_e response by direct blocks only) | 0.835 | 0.851 |
| **BOX restricted to conservation-consistent steps (BOXE, all products exact)** | **0.896** | **0.912** |

(`p1/p1_methods.py`, `p1/box_conic.py`, `p1/box_proj.py`, `p1/box_exact.py`; logs `methods_*.log`, `box_*.log`.)

**2. The span omits n_e over most of the model, and the projection applies an unscored change.** The driver's activity
mask is `R != 0`. `--conserve` makes the charge residual exactly 0.0 at 143–150 of 167 depths (He number at 5–9), so
the charge equation **and the n_e variable** (same index) drop out of every direction, the bound and the least squares
there, although the native blocks hold a full charge row (34 nonzeros, diagonal −x_ne). Every trial's projection then
sets n_e to the charge sum: for the driver's small steps a 1–2 % n_e change, for a box step up to 48 % — never scored.

**3. Planned vs realized vs native.** With the projection modeled (BOXP), realized and planned agree to 1e-4, but the
native merit was ×4.1 / ×8.0 / ×15.9 at λ = 0.25 / 0.5 / 1 (both bases) — growth linear in λ, i.e. a missing
first-order term, 99.9 % in row 2 (the replaced H ionization balance) at the H front (d22, 23, 30, 38). One exact
product pass of the projection's correction direction (`p1/ne_product.py`) closes it: the model with that exact
product predicts 4.16 / 8.15 / 16.2 against native 4.13 / 8.02 / 15.9; the remainder grows as λ² (0.20, 0.79, 3.16).
The radiation response of the n_e change (‖·‖ = 16× the merit, all in row 2) is what the direct blocks miss.

Native confirmation of the conservation-consistent box step (BOXE; driver trial pipeline exactly — linear map, 0.693
limit, He III trace balance, projection; `p1/p1_native.py`, passes in `results/nk167_trap/v1_p1/`):

| Base | λ | planned | realized | **native** | frozen measure | driver's best trial at this base |
|---|---|---|---|---|---|---|
| it04 | 1 | 0.896 | 0.896 | **0.907** | 2.657e-2 (−9.3 %) | 0.980 (l0.5) |
| it04 | 0.5 | 0.932 | 0.932 | 0.933 | 2.733e-2 | |
| it16 | 1 | 0.912 | 0.912 | 1.130 | 3.205e-2 | all rejected (best 1.00001) |
| it16 | 0.5 | 0.943 | 0.943 | **0.958** | 2.719e-2 (−4.2 %) | |

Conservation after projection: ≤ 1e-12 (the step is consistent); no native pass failed.

**4. Classification.** (a) unstable subspace algebra — yes (normal equations); (b) inadequate feasible span — yes (n_e
and the charge row absent at ~150 depths; directions conservation-inconsistent: the consistency map has rank 140 over
the 1266 directions, 48–53 over the 70 Krylov directions); (c) poor local model — **no**: exact products predict the
native response to O(λ²), and the it04 natives agree with prediction within 1 %. The it16 λ = 1 overshoot is
ordinary nonlinearity at the bound.

**Promotion rule (plan, proposed).** it04 meets it (predicted 10.4 %, native 9.3 %); it16 shows a meaningful native
decrease (4.2 %) but predicted 8.8 % < 10 %. Rule not met on both bases → no campaign launch; bounded P2 test instead.

**Costs and reuse.** The box gain needs the 1196 dense-window columns: on the 70 Krylov directions alone every method
stays at 0.97–0.99 (`p1/nodense_check.py`). Dense columns cannot be recycled across bases: columns of it15 / it12 /
it05 (radiation part from their base, direct part exact) predict 0.91 at it16, but the same step with it16's exact
columns scores 1.47 / 2.39 / 2.55 (`p1/recycle_check.py`). A fresh dense set costs ~5.4 h per base.

**Driver changes from P1** (opt-in; archived iterates replay unchanged): `--struct-active` (number/charge rows and their
variables active wherever their direct-block row is nonzero), `--consistent-dirs` (each direction replaced by its
linearized conservation projection `cnk_trap.project_linear` before its exact product; needs `--struct-active`),
`--step-solver ne|qr|box` (`bounded_step`: augmented-QR LM family or explicit box SOCP; box needs `clarabel`, installed
into `.venv` with uv, not yet in `requirements-lock.txt`). Tests: `tests/test_cnk_trap.py` (27 pass).

## P5: Ball Q2 source audit (started; verdict unchanged)

[ball_q2/README.md](ball_q2/README.md): **mismatch.** The paper (§3.1) states the inversion is apparent in all but the
first profile and appears once T_surf < ~8000 K; that absence at row 0 is visual (log density over nine decades) and the
paper gives no amplitude. Q2's adjacent-cell 1e-3 test is resolution-dependent: the row-0 rise near 29 kK is a
converged ~2 % (0.0085 dex) feature at a fixed location (N = 200/400/800: total rise 1.55/1.89/1.94 %; the adjacent-cell
value 9.1e-3/5.4e-3/2.8e-3 would pass Q2 near N ≈ 2300 with the same physics). The smoke track shows the paper's
recombination inversion reappearing at T_eff ≈ 7.6 kK. A draft source-faithful successor check (Q2-S) is in the audit
README; it is **not registered** and requires the user's explicit approval (thresholds would be chosen after seeing
the data). The registered FAIL stands. The 46 M☉ merge-event reproduction and the 55 M☉ stall were not started.
