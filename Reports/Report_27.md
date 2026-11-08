# Report 27 — continuous corrected closure run (from 2026-11-07)

**Status: QStarS is NOT v1.** The 167-radius cool closure (released T) is being solved continuously. Live record;
updated as the run progresses.

## 1. Setup

Continues Report 26 (§6b: the re-optimized 0.25 box adopted). At the user's request the run goes on without
review checkpoints; it stops only on its own stop conditions.

- Run directory `results/nk167_trap/v1_p2`, start `it23_trial_l0.5` (frozen it16 measure 1.9364e-2), iterate 23 on.
- Method: Report 26 flags (`--struct-active --consistent-dirs --step-solver box --span-tau 1e-6`), box 0.25 (T 0.015),
  plus two changes: `--ortho-dirs` (each batch orthonormalized against the span before its product pass) and
  `--krylov 12` (the bounded fit still gained 0.075 from pass 3 to 8 in iterate 20).
- Stop conditions: all trials rejected (driver); progress < 0.5 %/h on two consecutive iterates; RE net/gross max
  > 5e-3; significant max|s| > 16 with frozen gain < 2 %.
- Progress is the frozen it16 measure (`scripts/campaign_norms.py`), per-iterate evaluation by
  `Reports/evidence/v1_exec_20261106/p2b/block_eval.py`.

## 2. Progress

| iterate | passes | accepted | frozen measure | drop | wall | %/h | RE net/gross max | sig. max|s| |
|---|---|---|---|---|---|---|---|---|
| 22 (end of Report 26) | 8 | λ 0.5 | 1.9364e-2 | | | | 3.43e-3 | 7.8 |
| 23 | 12 (ortho) | λ 1 | 1.8900e-2 | −2.4 % | 2.46 h | 0.97 | 3.55e-3 | 12.2 |
| 24 | 12 (ortho) | λ 0.5 | 1.8281e-2 | −3.3 % | 2.45 h | 1.34 | 3.55e-3 | 7.3 |

*Change after iterate 24:* 12 passes gained about what 8 did per iterate (2.4 / 3.3 % against 2.6 % at it22) for
~35 % more wall time, so the run continues from it25 with `--krylov 8` (orthonormalized directions kept).

| iterate | passes | accepted | frozen measure | drop | wall | %/h | RE net/gross max | sig. max|s| |
|---|---|---|---|---|---|---|---|---|
| 25 | 8 (ortho) | λ 0.5 | 1.7939e-2 | −1.9 % | 1.68 h | 1.11 | 3.54e-3 | 9.4 |
| 26 | 8 (ortho) | λ 0.5 | 1.7896e-2 | −0.24 % | 1.69 h | **0.14** | 3.53e-3 | 13.3 |

## 3. Iterate 26: the slowdown, diagnosed (2026-11-07/08)

Iterate 26 gained only 0.24 % (0.14 %/h): its bounded step predicted almost nothing (linear fit 0.984).

- **What the span can still reduce** (`evidence/report27/region_reach.py`, it26 products): each region's residual fitted
  alone — d1–9 → 0.15, d10–18 → 0.06, d19–32 → 0.14, d33–47 → 0.27, d48–54 → 0.33, d55–61 → 0.38, **d62–100 → 0.86**,
  d101–167 → 0.54 (of itself). d62–100 holds 33 % of the residual² and the global fit leaves it at 1.03.
- **Preconditioner test** (`pc_region_test.py`, one exact product pass of region × species pieces of W M⁻¹R): at the He
  front d55–61 the run's preconditioner (local ALI + PTC 0.03 + SVD floor 0.01) reaches 0.37 with the span, the local-ALI
  blocks without floor/PTC **0.10**, fixed-J blocks 0.16. For d62–100 no preconditioner helps (0.995 alone, 0.84–0.85
  with the span).
- **Not ALI noise** (`ali_tol_test.py`): the it26 state recomputed with QS_ALI_TOL 1e-13 (1602 formal solutions max)
  instead of 1e-10 changes the residual by 3.6e-9 in the frozen measure (of 1.79e-2); the d62–100 residual is real.

Change from iterate 28 on: preconditioner without the SVD floor and PTC shift (`--pc-svd-floor 0 --pc-ptc 0`); with
orthonormalized directions the large front modes that the floor tamed no longer distort the least squares.
| 27 | 8 (ortho) | λ 1 | 1.7900e-2 | +0.03 % (local −1.1 %) | 1.78 h | **−0.01** | 3.61e-3 | 10.1 |

Iterate 27 is the second iterate below 0.5 %/h: stop rule. Its diagnosis is §3; the run continues from it28_trial_l1 with the change named there.

| iterate | passes | accepted | frozen measure | drop | wall | %/h | RE net/gross max | sig. max|s| |
|---|---|---|---|---|---|---|---|---|
| 28 | 8 (ortho, no floor/PTC) | λ 1 | 1.7417e-2 | −2.7 % | 1.67 h | 1.62 | 3.56e-3 | 7.9 |
| 29 | 8 (ortho, no floor/PTC) | λ 1 | 1.6844e-2 | −3.3 % | 1.67 h | 1.97 | 3.51e-3 | 23.2 |
| 30 | 8 (ortho, no floor/PTC) | λ 0.5 | 1.6662e-2 | −1.1 % | 1.66 h | 0.65 | 3.49e-3 | 24.2 |

## 4. Where the solution stands (after iterate 30, 2026-11-08 07:00)

- **Stop condition at iterate 30:** significant max|s| 24.2 with a frozen gain of 1.1 % (< 2 %). The large step is the
  He II ground state at d25–26 (row 102), a near-trace ion there (2e-9 – 5e-8 of He, He III 1e-47) governed by the
  replaced He II/He III balance row — the same pathology as Report 25's d44–45. Every other species' unbounded step is
  ≤ 2.3. Not a bulk divergence; the run continues, but this variable will block the final convergence test and needs a
  treatment (e.g. the trace-ion balance extended to a near-trace He II), to be done before convergence.
- **Residual by region (frozen measure, absolute, 1e-4), it16 → it31:** d1–9 0.97 → 0.08, d10–18 0.17 → 0.06, d19–32
  0.11 → 0.06, d33–47 1.62 → 0.55, d48–54 0.15 → 0.06, d55–61 0.62 → 0.30, **d62–100 1.14 → 1.06**, d101–167 3.28 → 0.62.
  d62–100 now holds 38 % of the residual².
- **d62–100 is nearly converged in the convergence criterion** (`evidence/report27/step_convergence.py`): its unbounded
  Newton step is ≤ 0.027 (median 6e-4; 65 % of variables already below 1e-3). Its residual is spread over the highest
  H I levels (rows 23–31, ~3 % each, |R|/gross ~1e-4, alternating sign with depth) where hydrogen is fully ionized, and
  two near-trace He I levels; the row-scaled merit weights these near-LTE equations heavily. Its large share of the merit
  overstates its distance from the solution.
- **The convergence criterion itself** (all significant unbounded steps < 1e-3), corrected solver only: median step
  1.4e-2 (it20) → 8.6e-3 (it25) → 7.1e-3 (it30); share already below 1e-3 26 % → 29 % → 31 %; largest steps 10.5 / 9.4 /
  24 (He II trace above). The median halved in 10 iterates (~17 h); at that rate it reaches 1e-3 in roughly two more
  days, the largest outliers later. (Unbounded steps are least-squares estimates within each iterate's span.)

## 5. Why the gain per iterate falls, and three relaxation attempts that failed (2026-11-08)

At it31 the bounded step predicted 1.5 % against 9 % unbounded. **The bound is not the limit**
(`evidence/report27/radius_test.py`, it31 products): radius 0.25 / 0.35 / 0.5 and T radius 0.015 / 0.05 / 0.15 / 1 all
predict 0.984–0.985. Most components at the bound are near-trace levels; the unbounded gain comes from huge relative
changes of near-trace populations (He II ground at the H front, s ≈ 24; near-trace He I in d62–100).

**Where the residual is** (it32, frozen measure): rows of levels with 1e-10 < fraction < 1e-8 hold 41 % of the
residual², < 1e-6 56 %, and trace levels (≤ 1e-10) another 32 %: ~88 % of the merit is in the equations of
low-population levels; the populations that carry the atmosphere are a small part of it.

**Relaxing the low-population levels directly** (each one native pass from it32_trial_l1, frozen measure 1.653e-2):

| attempt | relaxed | result |
|---|---|---|
| structural trace solve at fixed J, fraction ≤ 1e-6 (`trace_cleanup.py --sig-frac 1e-6`), one sweep | 21 283 variables, max s 5.5e5 | **2.5e1 (×1500)**: trapped levels (H I n=2 ~1e-8 of H at the front) cannot be solved at fixed J |
| same with the local-ALI blocks, factor-2 limit, conservation projection (`ali_relax_test.py`) | 21 283 | **1.44 (×87)**: d19–47 fronts destroyed |
| local-ALI relaxation of d62–167 only | 13 493, max s 2.2 | **5.9e-2 (×3.6)**: d62–100 ×10, d101–167 ×40 worse |

Local relaxation (Route B in its simplest form) does not work here: even in the inner atmosphere the local radiation
response is not accurate enough for these levels; their equations need the nonlocal exact response that only the
Newton–Krylov products carry. The relaxed states were separate test passes; the main run was not touched.
| 31 | 8 (ortho, no floor/PTC) | λ 1 | 1.6529e-2 | −0.8 % | 1.67 h | 0.48 | 3.46e-3 | 10.1 |
| 32 | 8 (ortho, no floor/PTC) | λ 0.5 | 1.6470e-2 | −0.36 % | 1.76 h | 0.20 | 3.45e-3 | 22.8 |

## 6. Stop (2026-11-08 10:07) and where this leaves the closure

Iterates 31 and 32 ran below 0.5 %/h: stop rule. With the diagnoses of §3–§5 there is no further evidence-based change
of the present method to make, so the run is stopped, not restarted. Best state: `results/nk167_trap/v1_p2/
it33_trial_l0.5` (pointer `accepted_state.json`), frozen it16 measure **1.6470e-2** (Report 26 start 2.837e-2: −42 %);
conservation ~1e-16; RE net/gross max 3.45e-3. Convergence criterion: median unbounded significant step 5.5e-3,
28.7 % of significant variables below 1e-3, largest 22.8 (near-trace He II at the H front).

What the run established:
- The corrected Newton–Krylov method (Report 26 fixes, 0.25 box, orthonormalized directions, unfloored local-ALI
  preconditioner) reduces the residual of the populations that carry the atmosphere; the frozen measure fell 42 % and
  every region except d62–100 dropped by 2–13×.
- Its limit is the equations of low-population levels (88 % of the remaining merit): they need large relative changes
  that a bounded linear-update step cannot make, and they cannot be relaxed locally (§5).

Candidate next methods (none started):
1. **Log-variable steps for low-population levels** (`--update-map exp` for levels below ~1e-6 of their element, bound
   in log space): Newton steps can then change a level by large factors; products stay exact. Lowest risk, a code
   change of a day plus a few iterates to judge.
2. **A trace-level Newton solve with exact products** between outer iterates: Krylov over the low-population variables
   only (significant state held), exact nonlocal response, until their equations converge — the nonlinear
   preconditioner of Route B with the exact operator instead of the local one.
3. **Reconsider which levels the convergence test must cover**: levels below ~1e-8 do not feed back on the atmosphere
   (the closing chain's cleanup already treats ≤ 1e-10 separately); a justified, registered threshold would change what
   "converged" means and needs the user's decision, not a solver change.
