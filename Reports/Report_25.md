# Report 25 — exact columns at the hydrogen front, physical soft modes, and population pseudo-transient continuation

**Status: QStarS is NOT v1.** <!-- outcome:start -->The 167-radius cool closure has not converged and no gate changes status. Stage 1
(statistical equilibrium at fixed T) ran 26 more iterates from Report 24's stall. Its row-scaled residual, the measure
of each equation's imbalance against its own gross rates, fell from 2.13 to 0.051 (about 40×) with conservation exact to
1e-15 and every native pass complete. The solve then reached a solver limit and stopped; stage 2 (T released) was
started from that state and meets the same limit.<!-- outcome:end -->

This session (2026-10-29 21:15 to 2026-10-30, on the always-on Mac mini) did what Report 24 §0 named as next. It put exact
Jacobian columns at the hydrogen ionization front into the preconditioner, then continued the constrained stage-1 solve
(statistical equilibrium at fixed T, 167 radii, `results/nk167_trap/fixedT_b35`). The exact columns were computed and
installed. They did not change the conditioning, because the near-singularity at the front is physical. The solve was
then moved onto a different footing, step by step, each change measured first on exact products at the same base:

1. **The front's soft mode is physical** (§2). With all 1168 exact front columns installed, the front blocks keep the same
   nearly singular mode (relative singular value 2e-6 at depth 23, 1.7e-6 at depth 26, 6e-5 at depth 31, against 5e-6,
   1e-6 and 8e-5 for the local-ALI blocks). It is the slow H I n = 2 + 3 / H II mode: the n = 2 and n = 3 rows are nearly
   parallel through Hα, and the rest is held by trapped Lyα and weak collisions. A coarse space built from these modes
   with their exact products shows that the full exact Jacobian is nearly singular along them too (Galerkin coarse
   matrix condition 2.6e19).
2. **The least squares was fitting the wrong measure** (§3). In every product subspace at the base, the step that lowers
   the preconditioned measure raises the row-scaled measure and the fixed measure in the exact linear model; each
   measure can only be lowered at the cost of another.
3. **Population pseudo-transient continuation** (`--pop-ptc`, §4) bounds the step along the soft modes. All trial
   passes have completed natively since (no supercritical maser), with conservation exact to 1e-15.
4. Four defects of the residual measure and the step showed up on the way and were fixed (§5–§7): frozen trace He III
   in the replaced He II ground rows (`--trace-ion-balance`), equations of gross scale 1e-19..1e-45 dominating the
   row-scaled measure (`--eq-floor 1e-15`), limit cycles of the component limiter (`--adaptive-cap`), and a least
   squares solved in a different measure from the one the trials are accepted on (`--merit rowscale`).

## 0. Outcome

**Where the solve stands.** Stage 1 stopped at iterate 29 (`fixedT_b35/it29_trial_l0.5`): every trial of that iterate,
1/16 fallbacks included, raised the residual. Driver measure 0.0302, row-scaled (iterate-1 row scales, floored) 0.051,
unfloored 0.06; conservation 1e-15; radiative-equilibrium imbalance at the pinned T at most 0.2 % of its gross rates
(depth 33). Stage 2 (`results/nk167_trap/free_r25`, T released, its own step radius 0.03 on T) started there:
0.0302 → 0.0295 in three iterates.

**The narrowest blocker** is the linear solve along the physically soft modes: the inner hydrogen front (depths 33–45,
H I n = 2 and the H ionization balance), the helium front (55–61) and the outer He I layers (1–9). The least squares
in the row-scaled measure, over directions built from the local-ALI preconditioner, now reaches only 0.82–0.98 of the
residual per iterate, and the realized steps gain 0.3–4 %. The preconditioned residual M⁻¹R is wrong exactly along
these modes (its forward model of the Jacobian is good, median error 0.10). None of the following changed the fit:
exact front columns (the requested fix, at iterate 4 and again in stage 2: 0.986 → 0.980), the slow modes of the front
blocks with exact products, 10 or 30 targeted cells, recycled directions of the previous iterate, PTC shifts of
0.03–0.1 in the preconditioner. The fronts' softness is physical (§2–3): in the exact Jacobian the excited and ionized
hydrogen at the 5000 K front is held only by small Lyα escape and collision rates.

**What changed in the method** (all in `scripts/cool_newton_krylov.py` and `scripts/cnk_trap.py`, each with a test):
the least squares and the acceptance in the same row-scaled measure, with the exact number and charge rows weighted
by themselves (§5, §7, §7c); trace He III balanced against the realized update (§5, §7b); equations of negligible
gross scale out of the measure (§5); adaptive component limits against limit cycles (§6); a trial step bounded to
|s| ≤ 0.5 by Levenberg–Marquardt in the subspace, T separately in stage 2 (§7b); targeted cell directions; a PTC shift
in the preconditioner only (§7b). Withdrawn after measurement: the coarse correction, per-depth front pieces, a finer
split, a column-scaled PTC, PTC inside the least squares.

**Next** (each a substantial piece of work; none was started):
1. A dense exact Jacobian on the front sub-domains (full nonlocal columns of the ~40 radiation-responsive rows, one
   depth per direction: about 1200 directions per front, ~6 h at the present product cost), kept for several iterates
   as the sub-domain solver of a two-level preconditioner.
2. A deep GMRES per Newton step (10–20 sequential product passes) instead of 2–3 batched block passes.
3. A nonlinear relaxation of the front populations at fixed radiation field between Newton steps.

The project moved to the external T9 drive during the session (§9b).

## 1. Starting point

Report 24 left stage 1 stalled at fixed measure 128 (`it04_trial_l1`, conservation exact). Its diagnosis: at the
hydrogen front (depths 23–28) the local-ALI preconditioner is so ill-conditioned that a deep direction's far-field
radiative action becomes 97–99.6 % of its preconditioned action, and the least squares fits the front with deep
directions. Resume instructions (run log): trim `history.json` to iterates ≤ 3, move `it05_trial_*` aside, restart
iterate 4 at `it04_trial_l1` with exact front columns in the preconditioner. That was done at 21:28; every later
start, stop and directory move is in `fixedT_b35/README.txt`.

The fixed measure of Reports 23–24 (`scripts/fixed_metric.py`: every state scored by the preconditioner of iterate 1)
stopped being informative in this session: it stays at 125–129 for every trial of every iterate, because the soft
modes it amplifies (front and deep blocks with relative singular values down to 1e-7) dominate it. Progress is
measured here by the row-scaled residual (each equation's imbalance relative to its gross rates, no inverse;
`Reports/evidence/report25/progress.py`), in particular its version with `--eq-floor 1e-15` (§5), column rows*.

## 2. Exact front columns (`--front-probe`)

`scripts/cool_newton_krylov.py --front-probe`, helpers in `scripts/cnk_trap.py` (`front_window`, `front_rows`,
`window_mask`, `block_conditioning`). The window is chosen each iterate as the depths where H II holds less than half
of the hydrogen, widened by two depths (19–47 at `it04_trial_l1`). The probed rows are H⁻, the 30 H I levels, H II and
n_e, plus every column whose local-ALI increment changes a significant equation in the window by more than 1e-4 of its
gross scale. That added the He I n ¹P levels (rows 37, 43, 51, 61, 73, 87), whose 584 Å-series resonance lines are
optically thick at the front although the levels are at 1e-17..1e-23 of helium, and He II ground: 41 rows. One
direction per row and stride-3 colour (123 directions, 7 batched product passes, 3 at a time: 39 minutes) gives the
exact response at each probed depth and its two neighbours, installed in the preconditioner only (1168 columns).
Tests: `tests/test_cnk_trap.py::test_front_rows_take_radiation_responsive_columns_and_front_probes_recover_window_columns`.

Result at `it04_trial_l1`, against the 56 exact products of Report 24's iterate-4 attempt at the same base
(`Reports/evidence/report25/front_check.py`, `front_check.json`):

| preconditioner | q median / p90 / max | front blocks: smallest relative singular value, median / min |
|---|---|---|
| local ALI + floor 1e-2 (Report 24) | 0.36 / 1.58 / 4.4 | 9.9e-3 / 9.8e-3 (floored) |
| + exact front columns, floor | 0.35 / 1.27 / 4.4 | 9.9e-3 / 9.8e-3 (floored) |
| + exact front columns, no floor in the window | 1.95 / 86 / 317 | 7.2e-7 / 1.9e-9 |

The smallest singular values of the front blocks are the same with the local-ALI and with the exact columns
(7.25e-7 → 7.18e-7 median, 1.86e-9 → 1.90e-9 min). Near-null modes (`exact_null` analysis, depths 21–44): the smallest
is a helium scaling mode (the He I ground and helium-number rows near-parallel, ultra-trace He II/He III variables),
handled by the floor; the next is the hydrogen mode above, present in both. Without the floor the exact-column
preconditioner is a much worse approximation of the Jacobian (q p90 86): what stabilizes the slow mode in the full
Jacobian is radiative coupling beyond the neighbour depths, which a block-tridiagonal preconditioner cannot hold.

## 3. The soft mode is physical, and the least squares fits the wrong measure

**Coarse correction** (`cnk_trap.slow_modes`, `CoarseCorrected`; offline `coarse_check.py`): the 52 slow modes of the
front blocks with an H/n_e-dominated left vector (relative singular value < 1e-2) got exact products J Z (3 passes).
A Galerkin coarse correction is singular (coarse matrix condition 2.6e19); a least-squares coarse correction degrades
q to a median of 12–18; with Z in the span the least squares moves the modes by 400–2700. The mode is nearly singular
in the exact Jacobian, nonlocal coupling included: the n = 2 + 3 + H II population at the 5000 K front is held only by
small Lyα escape and collision rates. A preconditioner cannot make that well conditioned; the step along it must be
limited.

**Per-depth front pieces** (`--front-pieces`, one Krylov direction per front depth for the H/n_e part): 51 directions,
rejected (`pieces_check.py`).

**Which step does the subspace contain?** (`objective_check.py`): the least squares solved directly in each measure over
the exact products at the base, each result evaluated by the exact linear model in all three:

| directions | least squares in | row-scaled | fixed | preconditioned |
|---|---|---|---|---|
| Report 24 iterate 4 (56) | preconditioned | 1.29 | 2.12 | 0.17 |
| | row-scaled | 0.68 | 0.92 | 0.86 |
| | fixed | 1.35 | 0.32 | 1.10 |
| all 107 | row-scaled | 0.59 | 4.68 | 2.09 |
| | fixed | 1.36 | 0.29 | 1.14 |

The step the driver took (first row) is predicted to raise both other measures; no subspace contains a step that
lowers all three. That is the stall of Report 24 seen from the linear model.

## 4. Population pseudo-transient continuation (`--pop-ptc`)

`Iterate.add_pop_ptc`: δ × sign(M_ii) × (row l1-norm of the direct blocks) is added to the diagonal of every
statistical-equilibrium row, in the Jacobian (every product, made or reused, is that of the shifted operator) and in
the preconditioner; the H and He number, charge and T rows get none. The residual is unchanged, so a converged state
is a root of the unmodified equations. SER: δ follows the row-scaled residual ratio of the accepted step, divided by
√λ when a shortened step is accepted (§6), and is set to 0 once the step is below `--pop-ptc-switch`. `pairs.npz`
keeps the unshifted products. Test: `test_population_ptc_shifts_jacobian_and_preconditioner_diagonals_except_constraint_rows`.

Offline on the 107 exact products at `it04_trial_l1` (`ptc_check.py`): for δ = 1e-3 … 1e-1 the exact linear model
predicts row-scaled ×0.84–0.88 (undamped ×1.33) and a He-front step component of 1.2–1.7 (undamped 6.7, where the
Report 24 masers formed). With PTC the exact front columns add nothing (×0.83 against ×0.84), so the probe pass
(39 min per iterate) was dropped.

In the run, iterate 4 (δ = 0.1, 20 directions, one split pass) had all four trials complete natively; the full step
lowered the row-scaled residual 2.60 → 1.36 in every depth region. The driver accepted λ = ¼ on its natural measure;
the acceptance rule was changed to the row-scaled residual (`--accept-rowscale`, PTC practice) and iterate 4 amended
to λ = 1 (`history.json`: `accepted_by`, `trials_natural`).

## 5. Replaced ground rows: the measure and trace He III

The forced equation selection (`QSTARS_REPLACE_FORCE`, `results/nk167_adv/transfer/replace_167.txt`) replaces the He II
ground row by the He II / He III ionization balance and the H I ground row by the H I / H II balance. Where He III is
1e-11..1e-43 of helium, that row's variable is He III, whose least-squares weight is 0 (`--ultra-weight 0`) or 0.001–0.03:
no direction moves it, the row cannot be satisfied (linear model 0.337 → 0.341, native 0.344 at depth 20), and the
preconditioner turns its imbalance into He II ground corrections (s 17–61 at depths 42–45).

- **`--trace-ion-balance`** (`cnk_trap.trace_ion_balance`): in every trial, He III where it is trace (≤ 1e-10 of He)
  takes the component that zeroes the linearized He II / He III row. At 75 depths; the unfloored row-scaled residual
  fell 1.13 → 0.46 in one iterate. Test: `test_trace_ion_balance_zeroes_the_linearized_replaced_row_only_where_the_ion_is_trace`.
- **`--eq-floor 1e-15`** (an existing option, now used for the acceptance and SER measure only; the PTC scale keeps each
  row's own gross scale): these rows have gross scales 1e-19..1e-45 of the largest at their depth and made up 62–93 % of
  the row-scaled measure. The same numbers result at floors 1e-15, 1e-12 and 1e-10.

## 6. Limit cycles

He II ground at the inner H-front edge (depths 44–45, 1e-6 of He) took steps +32, −3, +54, −5, −2, +27, +51 in
successive iterates: the factor-2 limiter carried it across its root each time (He II/He 1.5e-6 ↔ 3.6e-6, H II/H at
depth 44 0.16 ↔ 0.28). Its own row is the replaced He II / He III balance (row gross 1e-22..1e-37 of its column's), so
PTC never damped it; it is driven by a 0.2 % imbalance of the He I ground row, whose gross rates are the trapped 584 Å
line.

- **`--adaptive-cap`** (`cnk_trap.adapt_caps`): a limited component whose update reversed direction from the previous
  accepted update gets its limit halved; the others grow back to 0.693. Test:
  `test_adapt_caps_shrinks_only_limit_cycling_components`.
- **λ-aware SER**: δ × ρ / √λ_accepted, so a shortened step raises δ.
- *Withdrawn*: a PTC scale of max(row gross, column gross) for the replaced-row ground states (`--pop-ptc-colscale`,
  removed). Pinning He II ground locally made the He I ground and helium-number rows parallel (block conditioning
  6e-9 → 4e-16); the preconditioner broke at depth 43 (max|M⁻¹R| 1e12) and the run was stopped within a minute.

## 7. The least squares in the acceptance measure (`--merit rowscale`)

With the least squares in the preconditioned measure, the exact linear model predicted, on the iterate-11 products, an
increase of the acceptance measure for every δ from 0.003 to 0.416 (×1.23 to ×9): the fit chased the He II ground
component at depth 44 (+60 to +90). The H ionization balance at depth 18 stayed at −0.06 for every δ, so this was not
PTC stagnation. Solved in the row-scaled measure itself (directions still preconditioned residual pieces): ×0.85 at
δ = 0.1, ×0.93 at 0.416. Adopted at iterate 12; since then every trial has been monotone in λ.

A finer split of the residual pieces (`--split-fine 1 60 5`, iterate 14) lowered the linear residual (0.50 against
0.56) but the full step gave ×1.00: the extra freedom put large amplitudes into components with little row-scaled
action (step max|s| 752 off the significant set, 252 limited variables). Dropped.

## 7b. Past the PTC: what the remaining residual needed

From iterate 16 the floored measure fell only 0–6 % per iterate. Each cause was measured on the exact products of the
iterate in question before the driver was changed:

- **He III balanced against the realized update** (iterate 18). The first `--trace-ion-balance` set He III against the
  *planned* step of the other components; the limiter then truncated those while the He III change went through, so at
  the He front (depths 59–60) the balance pushed row 102 the wrong way each iterate (0.03 → 4.4, He III ×15; the direct
  linear prediction 0.572 against native 0.563 shows the Jacobian itself is right). He III is now balanced inside each
  trial against the realized, limited update of every other component (direct action, unshifted). The unfloored
  residual came back 5.33 → 2.39 → 0.71.
- **A trapping-aware PTC scale** (`--pop-ptc-scale pc`, iterate 19). What held the floored residual (75 %) was the H
  ionization balance at depth 18, which only H I n = 2 can satisfy. With the PTC scale = fixed-J gross rates (untrapped
  Lyα, 3e4 × the local-ALI row scale there) the preconditioner's n = 2 correction at depth 18 was exactly 0, so no
  direction contained it. With the row scale of the preconditioner blocks (local ALI included) M⁻¹R has n = 2
  corrections −0.35..+0.17 and its largest significant component falls 187 → 44.
- **Targeted directions** (`--target-dirs 10`, iterate 20): per pass, the preconditioned correction of each of the ten
  cells with the largest row-scaled residual, alone. Up-weighting the depth-18 cell showed the subspace could fix it only
  by damaging the rest (×2.2); the fix is a coordinated local change (n = 1 and n = 2 +21 % together through the Lyα,
  H-number and charge rows) that region-wide pieces cannot express. First-pass fit 0.43 against 0.74.
- **No PTC inside the least squares** (iterate 20). True-Jacobian predictions on iterate 20's 58 products: ×0.92 at
  δ = 0.156, ×1.04 at 0.05, ×1.83 at 0.02, ×0.405 at δ = 0. In a least squares the shifted operator absorbs residual
  through δ P s instead of damping the step. Iterate 20 was redone as pure Newton on the same products (its PTC step
  discarded, `superseded/it20_ptc_step`): **0.1075 → 0.0490 in one iterate (×0.456)**, every trial native-complete and
  monotone in λ.
- **PTC in the preconditioner only** (`--pc-ptc 0.1`, iterate 21): without any shift the directions are again
  dominated by amplified soft modes (first-pass fit 0.874). The shift is kept in M (the directions) and removed from the
  Jacobian, the products and the least squares. Test: `test_preconditioner_only_ptc_leaves_the_jacobian_exact`.
- **A bounded step** (`--lm-bound 0.5`, iterate 23). Without the Jacobian shift the full steps overshot on the inner
  front's soft mode (H I n = 2 at depths 37–38, s = 3.6; full steps ×5.5 and ×2.4 worse, λ = ¼ and ⅛ accepted). The
  line-search step is now the Levenberg–Marquardt step in the subspace with the smallest damping that keeps every
  |s| ≤ 0.5 (the convergence test keeps the unbounded Newton step). On iterate 22's products: ×0.813 against ×0.778
  unbounded, nothing to clip. In iterate 23 every trial improved, the full step too.

## 7c. The constraint rows of the row-scaled measure

The row-scaled objective weighted each equation with the significance weight of the variable of the same index. For
the three exact rows that is wrong: the helium-number row is the He III row, and He III is ultra-trace (weight 0) in
the neutral-helium layers, so helium conservation had no weight there. In iterate 26 the least squares changed He I
ground by +183 % and −152 % at the adjacent depths 56 and 57; every trial's projection restored the helium totals
(corrections up to 27 %), so the realized update differed from the step the linear model scored (fit 0.758, full step
×1.045). With the three exact rows at weight 1 (`Iterate.row_weights`; test
`test_row_scaled_measure_weighs_the_constraint_rows_whatever_their_variables_weight`) the same products predict ×0.964
with helium held to 1e-3, which is what the native trials of iterates 23–26 delivered (×0.95–0.97). The better fits
since iterate 12 were in part this defect. The projections of iterate 27 are ≤ 0.7 %.

Two cheaper-subspace attempts in the same period: 30 targeted cells instead of 10 (iterate 25, no gain: the residual
is diffuse, the largest 100 cells hold 61 % of its square) and recycling the previous iterate's directions with their
radiation part from their own base and the direct part recomputed (`--recycle-stale`, iterate 26: fit 0.80 with one
fresh pass instead of 0.89, but the stale radiation parts make the prediction optimistic).

## 8. Trajectory

Accepted bases of the constrained stage-1 solve, measured by `Reports/evidence/report25/progress.py` (iterate-1
reference row scales; rows* with `--eq-floor 1e-15`; conservation: largest |factor − 1| of the H, He and charge
projections). Passes of iterates 5–8 were packed by retention after they were measured (`progress_log.txt` keeps the
latest table).

<!-- trajectory:start -->
| base | method (from this iterate on) | rows* | rows (unfloored) | fixed | conservation |
|---|---|---|---|---|---|
| it04_trial_l1 (Report 24) | local-ALI preconditioner, floor | 2.125 | 2.60 | 127.9 | 1e-15 |
| it05_trial_l1 | population PTC δ 0.1, natural acceptance (amended to row-scaled) | 0.809 | 1.36 | 128.7 | 7e-16 |
| it06_trial_l1 | row-scaled acceptance | 0.507 | 1.18 | 128.3 | 1e-15 |
| it07_trial_l1 | | 0.293 | 1.14 | 127.2 | 8e-16 |
| it08_trial_l0.5 | `--eq-floor 1e-15` | 0.241 | 1.13 | 127.6 | 6e-16 |
| it09_trial_l0.125 | `--trace-ion-balance` | 0.235 | 0.46 | 126.2 | 8e-16 |
| it10_trial_l1 | λ-aware SER | 0.192 | 0.27 | 124.7 | 9e-16 |
| it11_trial_l1 | | 0.187 | 0.20 | 129.0 | 1e-15 |
| it12_trial_l0.125_d16 | `--adaptive-cap` | 0.187 | 0.20 | 129.0 | 6e-16 |
| it13_trial_l1 | `--merit rowscale` | 0.170 | 0.18 | 126.2 | 8e-16 |
| it14_trial_l1 | | 0.172 | 0.22 | 127.0 | 1e-15 |
| it15_trial_l0.25 | (`--split-fine`, dropped) | 0.163 | 0.66 | 126.6 | 1e-15 |
| it16_trial_l1 | `--krylov 2` | 0.149 | 1.41 | 129.3 | 1e-15 |
| it17_trial_l0.125 | | 0.150 | 2.38 | 128.4 | 6e-16 |
| it18_trial_l0.5 | | 0.149 | 5.33 | 126.8 | 9e-16 |
| it19_trial_l1 | He III against the realized update | 0.148 | 2.39 | 130.5 | 7e-16 |
| it20_trial_l0.5 | `--pop-ptc-scale pc` | 0.144 | 1.45 | 130.7 | 1e-15 |
| it21_trial_l1 | `--target-dirs 10`, no PTC in the least squares | 0.074 | 0.71 | 132.9 | 2e-15 |
| it22_trial_l0.25 | `--pc-ptc 0.1`, `--krylov 3` | 0.066 | 0.29 | 133.2 | 2e-15 |
| it23_trial_l0.125 | | 0.064 | 0.22 | 133.4 | 1e-15 |
| it24_trial_l0.5 | `--lm-bound 0.5` | 0.058 | 0.20 | 133.8 | 8e-16 |
| it25_trial_l0.5 | (`--target-dirs 30`: no gain) | 0.057 | 0.19 | 133.8 | 7e-16 |
| it26_trial_l0.5 | `--pc-ptc 0.03` | 0.053 | 0.17 | 133.9 | 1e-15 |
| it27_trial_l0.5 | `--recycle-stale 1` | 0.052 | 0.15 | 133.8 | 9e-16 |
| it28_trial_l0.5 | constraint rows weighted (§7c) | 0.052 | 0.09 | 133.5 | 8e-16 |
| it29_trial_l0.5 (stage-1 end) | `--krylov 2`; iterate 29: slow-mode span, no decrease, stop | 0.051 | 0.06 | 133.2 | 1e-15 |

Stage 2 (`free_r25`, T released, driver measure with the RE row): 0.0302 → 0.0295 (iterate 1, λ ¼) → 0.0295
(iterate 2, λ ⅛) → 0.0294 (iterate 3 with exact front columns, λ ⅛); stopped at the start of iterate 4.
<!-- trajectory:end -->

Every trial pass since iterate 4 has completed natively; none met a supercritical maser. The fixed measure (iterate-1
preconditioner) does not follow the solve; its rise from 128 to 134 is the moving front scored by the iterate-1 blocks
(§1).

## 9. Gates

<!-- gates:start -->
No gate changes status; Report 24 §6 stands.

| Gate | Status |
|---|---|
| 167-radius T-released closure converged | **FAIL**: stage 1 stopped at a solver limit (§0), stage 2 started |
| 1. Post-convergence cleanup and audits | NOT EVALUATED |
| 2. 333 radii | NOT REACHED |
| 3. Doubled frequency | NOT REACHED |
| 4. Self-selected CMFGEN wind | NOT REACHED |
| 5. Node-coupled conservative evolution | NOT REACHED |
| 6. Ball et al. 2011 benchmark | **FAIL** (Report 23 §6). It cannot pass under its registration: check Q2 ("no density inversion at row 0") fails at every resolution because the model has a resolved inversion near 28 kK, and the registration forbids relaxing it. Moving it needs the owner's decision on Q2, not more computation |
| v1 | **NO** |
<!-- gates:end -->

## 9b. Operations: power loss and the move to the T9

The Mac mini switched off about 12:00 on 2026-10-30 (no cause recorded; restart after power failure is off, which
needs the owner's `sudo pmset -a autorestart 1`). Iterate 25's second Krylov pass was lost (`superseded/it25_powerloss`);
its first pass was reused.

The internal disk was then down to 26 GB free, so the whole project moved to the external T9 SSD (exFAT, which also
holds unrelated user folders). It lives in two APFS sparse-bundle images, `/Volumes/T9/QStarS/QStarS.sparsebundle`
(repository, `.git`, `.venv`, `results/`) and `qsrd-runtime.sparsebundle` (MESA and its SDK, the executables'
`@rpath`). Each is attached at its original path (`/Volumes/T9/QStarS/attach_qstars.sh`, also `scripts/attach_qstars.sh`,
after every restart), so no absolute path in the code, run records or executables changed. The run records compare
resolved paths when they reuse passes, so a symlink to the drive would have broken same-base reuse. Copies were
made with `ditto` and checked entry by entry with BLAKE2b (repository 626 967 entries, 79.5 GB; runtime 20 685
entries, 21.7 GB). `ditto` had folded 414 AppleDouble `._` files into extended attributes; they were copied back byte
for byte. A native replay pass takes 39.3 s from the T9 against 37.6 s internally (raw exFAT: 45.1 s), with an identical
residual; the test suite passes (357).

## 10. Mistakes and superseded attempts

- A second driver ran for 22 minutes (22:32–22:56): the coarse-space driver was not stopped after its offline rejection
  before the PTC driver was started. Its first Krylov pass was moved to `superseded/` by the new driver; history was not
  touched (iterates ≤ 3 at the time). `results/nk167_trap/driver_running.sh` now guards every launch.
- Removing `--pop-ptc-colscale` also deleted the `--trace-ion-balance` option definition; the 02:02 restart failed at
  argument parsing and 10 minutes were lost. The full launch command is now parsed before every restart.
- Every stopped pass is kept in `fixedT_b35/superseded/` (`README.txt` there).

## 11. Reproduction

```bash
# exact front columns vs the Report 24 preconditioner (56 exact products at it04_trial_l1)
.venv/bin/python Reports/evidence/report25/front_check.py
# the least squares in each measure; per-depth front pieces; coarse space; population PTC
.venv/bin/python Reports/evidence/report25/objective_check.py
.venv/bin/python Reports/evidence/report25/pieces_check.py
.venv/bin/python Reports/evidence/report25/coarse_check.py
.venv/bin/python Reports/evidence/report25/ptc_check.py [--front]
# progress of every accepted base in the row-scaled (rows*), fixed and conservation measures
.venv/bin/python Reports/evidence/report25/progress.py it04_trial_l1 it05_trial_l1 ...
```
`results/nk167_trap/fixedT_b35/README.txt` lists every start, stop and directory move of the run.
