# Report 28 — larger steps for low-population levels (option 1), and what actually limits the step

**Status: QStarS is NOT v1.** The 167-radius cool closure is not converged. Option 1 of Report 27 (let low-population
levels change by large factors per step) was implemented and tested offline; it does not help. The diagnosis it
produced points to targeted exact Jacobian columns, which are being tested natively (§4).

## 1. What was implemented

`scripts/cnk_trap.py`:
- `update_factor(..., knee=...)`: per-component knee of the linear update map (the exact linear range x(1 − s)
  extends to s = knee; the positive C1 tail beyond starts there).
- `lowpop_bounds(frac, active, rad, thresh, factor)`: asymmetric bounds −(F − 1) ≤ s ≤ 1 − 1/F for levels below
  `thresh` of their element (never H⁻, H I n = 2, H II, He III, n_e, T), with the matching knee and cap ln F.

`scripts/cool_newton_krylov.py`: `bounded_step` accepts a pair (lower, upper) of per-component bounds (QR and box
solvers; a symmetric pair reproduces the scalar radius). Tests: `test_lowpop_bounds_allow_factor_f_changes_with_an_exact_linear_update`,
`test_bounded_step_respects_asymmetric_bounds_in_both_solvers` (31 tests pass). Not yet wired into the driver.

## 2. Offline result: no gain (iterate 32 products, `evidence/report28/lowpop_offline.py`)

| low-population bounds | components | predicted bounded fit |
|---|---|---|
| none (run setting: 0.25, T 0.015) | 0 | 0.9904 |
| fraction < 1e-8, F = 10 | 17 635 | 0.9904 |
| < 1e-6, F = 10 | 20 877 | 0.9904 |
| < 1e-6, F = 100 | 20 877 | 0.9904 |
| < 1e-4, F = 10 | 21 038 | 0.9904 |

The low-population components use the room (|s| up to 9 / 25) without any gain. Widening the bound of the major
components as well (`radius_scan.py`): 0.25 → 1 → 3 → 10 → 30 gives 0.9904 → 0.9887 → 0.9877 → 0.9855 → 0.9849. The
unbounded fit (0.9217) is reached only through absurd components — the largest is He III at d14 (fraction 2.5e-19,
|s| 8e4), a trace ion the trace-ion balance sets anyway. **Within any physical step the span offers ~1–1.5 %.**

## 3. The span is the limit, and so is the preconditioner (`split_test.py`, one exact product pass at it32)

| directions | unbounded | bounded |
|---|---|---|
| region × group pieces of W M⁻¹R (as the driver), alone | 0.9991 | 0.9991 |
| the same split into major / low-population parts, alone | 0.9960 | 0.9969 |
| it32 span (8 passes) | 0.9217 | 0.9904 |
| it32 span + split pieces | 0.9144 | 0.9855 |

The preconditioned residual is now almost useless as a direction (0.999): the local-ALI preconditioner no longer
predicts the correction; the span's gains come from repeated exact products. Splitting by population class adds little.

**Where the needed corrections are** (it32 unbounded step, `step_convergence`-style): major populations (> 1e-6 of
their element, 907 variables) median 2.7e-3, 90 % 0.10, max 2.3; low-population (8 904) median 5.9e-3; n_e max 0.08;
T max 0.44 (d40, 33, 37, 8, 17 — the inner H front and the outer layers). The 154 major components with |s| > 0.05 are
the H ionization balance (row 2, 55 cells), He II ground (30), He I ground (27), T (26) and H II (14), at d1–61 and a few
at d120–166. The atmosphere's own ionization/temperature structure at the fronts is not converged; the low-population
levels are not the only issue.

## 4. Test in progress: targeted exact columns

Exact columns gave the largest native steps of this work (Report 26 P1: −9.3 % / −4.2 % per step with 1196 dense
front columns), but the full front window cost 5.4 h per base. Restricted to the seven rows that carry the large
corrections — rows 2, 3, 32, 33, 102, 133, 134 (H ionization balance, H I n = 2, H II, He I and He II ground, n_e, T) — over
d1–61 it is ~430 exact columns (~2 h). Iterate 33 from `it33_trial_l0.5` runs with them plus 8 Krylov passes
(`--dense-window 1 61 --dense-rows 2 3 32 33 102 133 134`, otherwise the Report 27 setting). Result: §5.

## 5. Iterate 33 with targeted exact columns

427 exact columns (117 min) fit the base residual to 0.916 on their own (d1–9 0.048 → 0.007 of the residual²,
d33–47 0.089 → 0.034); with 8 Krylov passes (596 directions) the unbounded fit is 0.836 and the 0.25-box step predicts
0.912 (without the columns at it32: 0.990). Trials λ 1 / 0.5 / 0.25 / 0.125: 1.141 / **0.957** / 0.968 / 0.983 of the local
merit; accepted λ 0.5 (ρ 0.73): local 1.5994e-2 → 1.5308e-2 (−4.3 %), frozen 1.6470e-2 → **1.6117e-2 (−2.1 %)**, 3.82 h
→ 0.56 %/h (it31/it32 without the columns: 0.48 / 0.20 %/h). RE net/gross max 3.48e-3; significant max|s| of the
unbounded step 235 (soft front modes made reachable by the exact columns; the λ 1 overshoot is again curvature).

Continuation: the run goes on with the targeted columns every iterate, stop rule unchanged (< 0.5 %/h twice).

| iterate | dense cols | accepted | trials λ 1 / 0.5 / 0.25 / 0.125 | frozen measure | drop | wall | %/h | sig. max|s| |
|---|---|---|---|---|---|---|---|---|
| 33 | 427 | λ 0.5 | 1.141 / 0.957 / 0.968 / 0.983 | 1.6117e-2 | −2.1 % | 3.82 h | 0.56 | 235 |
| 34 | 427 | λ 0.25 | 1.260 / 1.002 / 0.987 / 0.991 | 1.5928e-2 | −1.2 % | 4.71 h | 0.25 | 738 |
| 35 | 427 | λ 0.25 | 1.215 / 1.001 / 0.989 / 0.992 | 1.5783e-2 | −0.9 % | 4.29 h | 0.21 | 1582 |

## 6. Stop (2026-11-08 23:57) and the state of the closure

Iterates 34 and 35 ran below 0.5 %/h: stop rule. Best state `results/nk167_trap/v1_p2/it36_trial_l0.25` (pointer
`accepted_state.json`), frozen it16 measure **1.5783e-2** — Report 26 start 2.837e-2: **−44 %**; conservation ~1e-16;
RE net/gross max 3.45e-3. Not converged: the unbounded significant step grows (235 → 738 → 1582 over iterates 33–35) and
every full step overshoots (×1.14–1.26) while λ 0.25–0.5 gains ~1–2 %.

**What this means.** With the step algebra repaired (Reports 26–27), every remaining lever of the Newton–Krylov
method has been measured on this problem: bound size and shape (no gain), low-population factor bounds (no gain),
local relaxation (diverges), deeper or split spans (small gain), targeted exact columns (−2.1 %, then ~1 % per 4 h).
The linear model is accurate for small steps (ρ 0.6–0.9), but within any step where it is accurate it offers only
~1–2 %; the large corrections it asks for (T up to 40 % at the inner H front; ionization balances at d1–61; soft
front modes, now unbounded steps of 10³) lie where the equations are strongly nonlinear. The method has reached its
practical limit on this model: further progress needs a globalization that crosses the nonlinearity, not more linear
algebra.

**Options (none started):**
1. **Continuation across the front nonlinearity** — e.g. a homotopy in a physical parameter that softens the fronts
   (a hotter or more strongly irradiated start, or the T-pin/frontpin homotopies of Reports 22–23 re-run with the
   corrected algebra), each stage solved by the now-reliable Newton–Krylov. Days of compute; the most direct attack on
   what limits the solve.
2. **A mature nonlinear solver comparison** (plan route C: PETSc SNES trust region / nonlinear preconditioning on the same
   exact residual and products). A week of adapter work; independent check of whether the limit is method or problem.
3. **The convergence-scope decision** (Report 27 §6.3): whether levels below ~1e-8 of their element belong in the
   Newton convergence test or in the closing chain's cleanup — a scientific decision for the user.
