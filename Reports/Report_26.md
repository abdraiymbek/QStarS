# Report 26 — executing the v1 plan: a trustworthy baseline, and why the closure solver stopped making progress

**Status: QStarS is NOT v1.** <!-- outcome:start -->The 167-radius cool closure has not converged and no gate changes
status. Executing P0–P2 of the [v1 execution plan](V1_execution_plan_2026-11-06.md) found three defects in the Newton
step construction, each measured on stored exact products and confirmed natively, and with them fixed the released-T
solve moves again: frozen-reference residual 2.837e-2 → 2.185e-2 (−23.0 %) in 4 accepted iterates (~8.9 h wall),
against −3.2 % in the 70 h dense-window campaign before it. The gain per iterate falls (−11.6, −7.4, −3.2, −2.8 %).<!-- outcome:end -->

This session (2026-11-06 18:20 to 2026-11-07, Mac mini, project on the T9 drive) executed the plan's immediate decision:
"P0 → P1 using existing products → a bounded P2 comparison; start the independent Ball audit without changing its
verdict". Evidence: [Reports/evidence/v1_exec_20261106/](evidence/v1_exec_20261106/README.md) (P0, P1, P5) and
[its p2/](evidence/v1_exec_20261106/p2/README.md). The dense-window run of 2026-10-31 to 11-03 that its run README calls
"Report 26" is summarized in the plan and in §1; this report takes the number.

## 0. Outcome

**Three defects, not physics, stopped the solve** (§3–§5):

1. **Normal equations at cond 1e53.** The bounded trial step (`lm_step_bound`) formed AᵀA of a least-squares matrix with
   cond(A) ≈ 3e26; at the driver's own μ the stable augmented-QR solution differs from its step by ~100 %. On the stored
   products of iterates 4 and 16: driver 0.976 / 0.9988 of the merit, stable LM 0.938 / 0.969, explicit box 0.87 / 0.90.
2. **n_e and the charge row dropped from the span at ~150 of 167 depths.** Activity was inferred from `R != 0`;
   `--conserve` makes the charge residual exactly 0.0, so the charge equation and the n_e variable (same index) left
   every direction, bound and fit, and each trial's projection then changed n_e unscored. For a box step that change
   was up to 48 %, and its radiation response in the H ionization-balance row was 16× the merit (native ×4–×16 worse,
   linear in λ, reproduced to O(λ²) by one exact product of the projection's correction).
3. **The Krylov least squares fitted product noise.** Batched pieces of preconditioned residuals are nearly collinear
   (238 directions, 185 independent); the least squares returned coefficients up to 3.8e14 cancelling to |s| = 0.4, so
   the products' own small errors became an O(1) model error (0.47 of the merit). A shrinking-step native directional
   derivative agrees with a fresh exact product of the step to 0.8 % — the derivatives are right, the subspace algebra
   was not. This is the "fit stalls at 0.8–0.9" of Reports 24–26.

**Fixes** (`scripts/cool_newton_krylov.py`, `scripts/cnk_trap.py`, opt-in so archived iterates replay unchanged, each
tested): `--struct-active`, `--consistent-dirs` (directions as the conservation projection applies them, exact products
of those), `--step-solver qr|box` (augmented QR; explicit component-constrained least squares by SOCP, clarabel),
`--span-tau` (least squares in the truncated orthonormal basis of the span). With them, at the base where the old
solver rejected every trial after 70 h, a 7-pass step gives −8.0 % natively and the driver's four iterates −11.6 %,
−7.4 %, −3.2 %, −2.8 % in the frozen measure (§6).

**What did not change:** the fronts are genuinely nonlinear — the driver's full steps at the |s| ≤ 0.5 bound overshot
in iterates 16–19 and λ = 0.5 was accepted (not "always": the corrected 7-pass step at it16 gained 8 % at λ = 1, and
in the council block the re-optimized 0.25 box was accepted at λ = 1, §6b). Dense-window columns cannot be recycled
across bases (§5c). Gates: unchanged (§8).

**P0** (§2) made the baseline trustworthy: the rejected-trial restart hazard is closed (best/accepted/status records, an
explicit restart pointer, refusal before any native pass), a frozen reference norm sits next to the local one, and the
accepted base replays bit for bit. **P5** (§7): the Ball Q2 source audit finds the registered row-0 test stronger than
the paper and resolution-dependent; a successor check needs the user's approval; the FAIL stands.

**Next** (§9): continue the corrected solve with the plan's stop rules and per-hour accounting; then release-T
convergence checks and the registered closing chain (P3). Decisions for the user are listed in §9.

## 1. Starting point

Report 25 left stage 2 (`results/nk167_trap/free_r25`, T released) at 0.0294. The session of 2026-10-31 to 11-03 ran
the dense hydrogen-front window (`--dense-window 19 47`, ~1200 exact single-cell columns per iterate, 13 builds ≈ 70 h)
from iterate 4 to 16: frozen measure 0.02930 → 0.02837 (−3.2 %); at iterate 16 all eight trials, 1/16 fallbacks
included, raised the merit and the driver stopped. Its `history.json` names the rejected best trial `l0.125_d16` as
`accepted` (written before the rejection test). The v1 plan (2026-11-06) reviewed this and set P0 → P1 → bounded P2.

## 2. P0: a trustworthy baseline

- **History and restart.** Records now carry `best_trial`, `accepted` (None unless the merit decreases strictly) and
  `status`; only accepted trials reach `OUT/accepted_state.json`; the driver refuses a start base that its run did not
  accept, before any native pass, and `--resume` starts from the pointer. Regression tests replay iterate 16's eight
  trial merits verbatim. `free_r25/accepted_state.json` → `it16_trial_l0.5`.
- **Frozen vs local norm.** `scripts/campaign_norms.py` freezes one iterate's weights and scales (it16, sha256
  `914b7cd4…`) and reports any pass by region, species group, largest cells, unfloored rows, newly significant rows,
  conservation and RE net/gross. It reproduces the driver's base merit to all digits.
- **Manifest and replay.** `scripts/baseline_manifest.py`: source diff and hashes, executable and 944 input hashes,
  equation map, two-photon J, launch command; the residual-defining `_tp` pass of the accepted state replays **bit for
  bit**, a fresh full pass from x0 to 1e-13.
- **Documents.** README, `docs/validation.md`, `docs/nlte_cool_closure.md`, `docs/cmfgen_selected_wind.md`,
  `docs/validity.md`, `docs/performance.md` and `Reports/evidence/completion_gates.json` brought to Report 25 and the
  dense run (descriptive text only; all 280 status/outcome/boolean fields of the gate file identical to HEAD). A
  31-criterion gate-to-artifact inventory: `evidence/v1_exec_20261106/gate_inventory.md`.

## 3. P1: step algebra on identical products (iterates 4 and 16)

`evidence/v1_exec_20261106/p1/p1_common.py` rebuilds an iterate's least squares from its stored products and reproduces
the driver's fit 0.8111, μ 6.667e-5 and bounded fit 0.998766 exactly. Method comparison, same weights:

| bounded step (|s| ≤ 0.5, T 0.03): linear fit | it04 | it16 |
|---|---|---|
| driver (normal equations) | 0.976 | 0.9988 |
| LM, augmented QR | 0.938 | 0.969 |
| truncated SVD (1e-2) + LM | 0.924 | 0.950 |
| explicit box (SOCP) | 0.869 | 0.895 |
| box, conservation-consistent (all products exact) | **0.896** | **0.912** |

## 4. P1: the n_e hole and the unscored projection

At it16 the charge residual is exactly 0.0 at 150 depths (He number 9), at x0 at 83. The box step on the projected
directions (projection modeled with direct blocks) predicted 0.851 / 0.835 and natively gave ×4.1 / ×8.0 / ×15.9 at
λ = 0.25 / 0.5 / 1 on both bases; 99.9 % of the gap in row 2 (H ionization balance) at d22–38. One exact product of the
projection's correction direction predicts the natives to O(λ²) (4.16 / 8.15 / 16.2 vs 4.13 / 8.02 / 15.9). Restricted
to conservation-consistent steps, the box step natively gave **it04 −9.3 %** (driver at that base: −2.0 %) and **it16
−4.2 %** (driver: all trials rejected); predictions within 1 % at it04.

## 5. P2: deep global Krylov, and the coefficient blow-up

### 5a. Fit by pass count
`results/nk167_trap/v1_p2`, it16 base, P1 fixes, no dense window, 12 passes of 34 directions (12.5 min each). The raw
fit fell 0.992 → 0.858 (5 passes) → 0.661 (12); the 5-pass bounded step predicted 0.868 and gave 0.995 / 0.986 natively.

### 5b. Diagnosis
Native trials at λ = 0.02 and 0.005 of that step: ‖FD + J s‖ = 0.0024 of the merit with a fresh exact product of s,
0.473 with Σ yᵢ J zᵢ from the stored products (‖FD‖ 0.298). Coefficients up to 3.8e14, ‖y·|zᵢ|‖ 5e13 for ‖s‖ 0.38.
In the orthonormal basis of the normalized span truncated at τ (`--span-tau`), the predicted fit equals the fresh-
product fit to four digits for every τ from 1e-1 to 1e-8 (12 candidates verified in one product pass).

### 5c. Dense columns do not transfer
Dense columns of it15 / it12 / it05 recycled at it16 (radiation part from their base) predict 0.91; the same step with
it16's exact columns scores 1.47 / 2.39 / 2.55. Without dense columns the old Krylov span gave nothing (0.97–0.99):
the gain then came from the dense span's near-orthogonality, not its location.

## 6. P2: the corrected solve

| iterate | method | passes | predicted (λ 1) | accepted | local merit | frozen measure | wall |
|---|---|---|---|---|---|---|---|
| 16 (base) | | | | | 2.8371e-2 | 2.8371e-2 | |
| 16 → 17 | raw coefficients, box | 12 | 0.784 | λ 0.5 | 2.5102e-2 (−11.5 %) | 2.5074e-2 (−11.6 %) | 3.0 h |
| 17 → 18 | `--span-tau 1e-6`, box | 8 | 0.851 | λ 0.5 | 2.3544e-2 (−7.5 %) | 2.3224e-2 (−7.4 %) | 2.1 h (shared machine) |
| 18 → 19 | `--span-tau 1e-6`, box | 8 | 0.918 | λ 0.5 | 2.2686e-2 (−3.6 %) | 2.2491e-2 (−3.2 %) | 1.8 h |
| 19 → 20 | `--span-tau 1e-6`, box | 8 | 0.912 | λ 0.5 | 2.1844e-2 (−3.6 %) | 2.1853e-2 (−2.8 %) | 1.8 h |

Frozen measures of the accepted states are those of their `_tp` residual (the next iterate's definition), except the
last (its trial pass). Trials λ 1 / 0.5 / 0.25 / 0.125 relative to the local merit: it17 1.070 / 0.925 / 0.949 / 0.972,
it18 1.127 / 0.964 / 0.972 / 0.984, it19 1.144 / 0.964 / 0.972 / 0.984.

Offline at it16 (driver trial pipeline): 7-pass τ-box step λ 1 → 0.920 (frozen 2.609e-2); 12-pass τ-box λ 0.5 → 0.892.
In iterates 16–19 λ = 1 overshot (×1.07–1.30) at |s| ≤ 0.5. An audit of the stored it18/it19 trials (§6b) shows
the native − linear-model gap growing as λ² (≈ 0.66 λ² of the merit), 65–99 % in row 2 (H ionization balance) at
d19–47: second-order curvature, not a pipeline effect (the λ = 1 projection factors, ≤ 4e-6, are far too small).
Conservation after projection ≤ 4e-6 at λ 1 and ~1e-13 otherwise; RE net/gross at most 1.85e-3 → 3.02e-3 (d33).
Raw-norm audit: ‖R‖ (unweighted) 4.4e15 → 1.5e16 from row 20 (H I, fraction 6e-7) at d151–167, whose imbalance
relative to its gross rates is 1e-5 → 6e-5; not a physical regression.

Residual distribution (frozen measure, share of ‖·‖²), it16 → it20: d1–9 0.120 → 0.065, d19–32 0.014 → 0.027,
d33–47 0.201 → 0.223, d55–61 0.078 → 0.080, d62–100 0.141 → 0.222, d101–167 0.407 → 0.339. The outer He I layers and
the deep interior fall; the inner H front and d62–100 hold their absolute size, so their share grows.

**Review (plan stop rules).** No stop rule fired: every iterate accepted, cumulative −23 % ≫ 1 %. Progress per wall
hour ≈ 2.6 %/h over the four iterates, falling from 3.9 %/h to 1.6 %/h; the predicted bounded fit rose from 0.78 to
0.91 while the unbounded fit stayed at 0.66–0.73. *Corrected after the council review:* that the bound, not the span,
limits the step was not shown — the unbounded fit worsened (0.661 → 0.731), the predicted in-box gain fell (22 % →
9 %) and the significant unbounded step grew (7.4 / 6.3 / 5.9 / 16.3, He II at d21–22). §6b tests the bound directly.
The run was stopped after the planned three corrected iterates for this review; resume with
`results/nk167_trap/launch_v1_p2.sh results/nk167_trap/v1_p2 $(pwd)/results/nk167_trap/v1_p2/it20_trial_l0.5 20 --krylov 8
--eta 0.01 --span-tau 1e-6 --lams 1 0.5 0.25 0.125 --max-outer N` (or `--resume` through the pointer).

### 6b. Council block: a re-optimized smaller bound (2026-11-07)

A three-member council (Claude Opus 5.5, GPT 6.1 Sol, Grok 4.7; three rounds) reviewed this report and agreed one
controlled test, approved by the user: from `it20_trial_l0.5`, only `--lm-bound 0.25 --lm-bound-t 0.015`, machine
exclusive, ≤ 3 accepted iterates, a matched A/B on identical products, a pre-declared radius rule, promotion threshold
and stop rules ([protocol and results](evidence/v1_exec_20261106/p2b/README.md)).

- **Planned vs realized (stored it18/it19 trials):** realized = planned in the model; the native − model gap grows
  as λ² (0.009 / 0.037 / 0.153 / 0.687 at λ 0.125 … 1), 65–99 % in row 2 at d19–47 — curvature, not the pipeline.
- **Matched A/B at it20, same products:** 0.25 box at λ 1 → −7.3 %; old 0.5 box at λ 0.5 → −4.0 % (λ 1 → ×1.33).
- **Three iterates at radius 0.25** (ρ 0.66 / 0.90 / 0.86, radius kept): accepted λ 1 / 0.5 / 0.5; λ 1 lowered the
  merit every time (0.927 / 0.962 / 0.979). Frozen 2.1853 → 2.0721 → 1.9872 → **1.9364e-2** (−11.4 % in 5.46 h,
  **2.09 %/h** against 1.75 / 1.59 %/h for the 0.5 box) — promotion criterion met; the 0.25 box is the new default.
  Since it16: −31.7 %. No stop rule fired. RE net/gross max 3.02 → 3.43e-3 (rising ever more slowly); d33–47 down
  30 %, d62–100 flat.
- **Krylov depth:** the bounded fit keeps gaining through pass 8 (0.963 → 0.888 from pass 3 to 8 at it20) although the
  unbounded fit flattens after pass 3; `--krylov 4` would forfeit most of the step.
- **Alongside (code only, committed):** `--ortho-dirs` (pre-product batch orthonormalization, for the next block),
  the P3 closing-chain plumbing (`closing_chain_check.py`; build 36, stored two-photon J and the solve's VADAT carried
  through `post_convergence.sh`, `trace_cleanup.py`, `verify_cool_state.py`), and the per-iterate ρ diagnostics.
- **Decisions recorded:** Q2 stays FAIL, Q2-S authorized in principle without a threshold; the composite active-gate
  entry split by component (Gate 6 FAIL, Gates 4–5 NOT_REACHED).

The per-iterate gain still decays within a block (5.2 → 2.6 %), as at the 0.5 box: at ~2 %/h the closure is still
far from the |s| < 1e-3 convergence test, so the next block must change the convergence rate, not only the step length.

## 7. P5: the Ball Q2 audit (verdict unchanged)

`evidence/v1_exec_20261106/ball_q2/README.md`. Paper §3.1: the inversion is "apparent in all but the first density
profile … It appears once the photospheric temperature T_surf drops below about 8000 K"; Fig. 1 shows log ρ over nine
decades, with no amplitude. Q2 tests an adjacent-cell outward increase > 1e-3 at row 0. The existing runs show a
converged ~2 % (0.0085 dex) rise near 29 kK at row 0 (N = 200/400/800: 1.55/1.89/1.94 %), whose adjacent-cell value
(9.1e-3/5.4e-3/2.8e-3) would fall below 1e-3 near N ≈ 2300 with unchanged physics; the paper's recombination inversion
reappears in the smoke track at T_eff ≈ 7.6 kK. Verdict: mismatch (the registration is stronger than the paper and
resolution-dependent). A draft successor check (Q2-S) is in the audit; it is not registered. The 46 M☉ merge event and
the 55 M☉ stall were not started (they would have competed with the P2 timing).

## 8. Gates

No gate changes status. Gate 6 (Ball) FAIL stands; Gates 1–5 not reached; the closure (Gate 0) not converged.

## 9. Next, and decisions for the user

1. Continue the corrected stage-2 solve (`launch_v1_p2.sh` + `--resume`), with the plan's stop rules and progress per
   wall hour recorded; consider a smaller bound or an adaptive trust radius, since λ = 0.5 of the 0.5 bound is accepted
   every time. Add orthonormalization of each new batch against the span before its product pass (cheaper than
   discarding dependent directions afterwards).
2. Re-examine Report 25's "physical soft mode" conclusion: the coarse-matrix condition 2.6e19 and the fit stalls were
   measured with the defective algebra; the nonlinearity at the fronts is real, the linear-solve limit was not.
3. Decisions: Ball Q2 successor registration (audit), the unresolved ledger statuses listed in the evidence README,
   `clarabel` as a dependency (installed into `.venv`, not yet in `requirements-lock.txt`), committing this work.

## 10. Reproduction

```
# P0
.venv/bin/python scripts/campaign_norms.py freeze results/nk167_trap/free_r25/it16_from_it16_trial_l0.5 FROZEN.npz
.venv/bin/python scripts/baseline_manifest.py results/nk167_trap/free_r25 --frozen FROZEN.npz --out M.json --replay DIR
# P1 (offline; box steps need clarabel)
.venv/bin/python Reports/evidence/v1_exec_20261106/p1/p1_methods.py it16_from_it16_trial_l0.5
.venv/bin/python Reports/evidence/v1_exec_20261106/p1/box_exact.py it16_from_it16_trial_l0.5
.venv/bin/python Reports/evidence/v1_exec_20261106/p1/p1_native.py it16_from_it16_trial_l0.5 boxe=y_BOXE_it16_from_it16_trial_l0.5.npy
# P2
results/nk167_trap/launch_v1_p2.sh results/nk167_trap/v1_p2 BASE K --krylov 8 --eta 0.01 --span-tau 1e-6 --lams 1 0.5 0.25 0.125 --max-outer N
.venv/bin/python Reports/evidence/v1_exec_20261106/p2/p2_tsvd.py --verify
.venv/bin/python -m pytest -q tests/test_cnk_trap.py
```
