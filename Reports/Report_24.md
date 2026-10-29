# Report 24 — under the stall: broken conservation, supercritical trial masers, and a constrained restart

**Status: QStarS is NOT v1.** The 167-radius cool closure has not converged, and no gate changes status. This report is
one night's work (2026-10-21 22:00 to 2026-10-22 06:00 CDT) on the stage-1 solve that Report 23 left stalled. It
replaces Report 23's account of the stall with four measured defects. Three are fixed and tested; the fourth is a
property of the trial states that the fixes are meant to avoid:

1. **The solver broke exact conservation laws.** The per-component step limiter clips each population's change on its
   own. A Newton step that satisfies the linear number- and charge-conservation equations no longer does after
   clipping. Every accepted state of every Report 23 solve violates them: helium number by up to 8×, hydrogen by up
   to 2.3×, charge balance by up to 34× (§4). The deep interior, converged in the starting state, was damaged
   iterate by iterate (fixed-measure residual 0.3 → 72). The "224 → 118" progress of Report 23 was the hydrogen fronts
   improving while this happened. Fixed by projecting every trial state onto the constraints (`--conserve`, §5).
2. **The least squares used distant directions as blunt instruments.** At iterate 1 the piece of the deep layers
   (He II, He III, n_e at depths 120–167) got the coefficient 172. 97 % of its preconditioned action lay at the
   hydrogen front (depths 24–28), where the nearly singular local-ALI blocks amplify a tiny far-field coupling
   (§4.2). That put s = 4.4 at the deep even depths, which the limiter turned into the conservation drift of item 1.
   A singular-value floor on the preconditioner blocks bounds the amplification (`--pc-svd-floor`, §2).
3. **The fixed-T stage did not hold T fixed.** T was a variable held by a weakly weighted pin row. At iterate 4 the
   least squares gave a deep temperature-only direction the coefficient 358, and the full-step trial put T at
   0.50–1.22 of the pinned profile. Fixed by removing T from every direction (`--freeze-t`, §4.3).
4. **Full Newton steps failed natively because the trial states contained supercritical masers**, not because of the
   preconditioner. The steps inverted H I infrared lines at the He front (12→8, 16→9, 15→9 at depths 60–64).
   CMFGEN's negative-opacity floor then left the total opacity below the electron-scattering opacity, and the
   radiation problem at those frequencies has no nonnegative solution: the assembled Λ operator has spectral
   radius 1.063 and 1.058, and its fixed point is negative at depths 57–67 (§3). Build 36 adds a Newton solve of the
   scattering fixed point for frequencies where the iterative one fails. It converges the frequencies that have a
   solution and shows why these do not.

## 0. Outcome

**Where the solve stands.** With the three fixes, stage 1 was restarted from its first accepted state (iterate 1). In
its first three iterates every trial pass completed natively, every accepted state kept hydrogen, helium and charge
conservation to 1e-15, and the full Newton step was accepted each time. The fixed-measure residual went
224 → 179 → 144 → 128. The unconstrained path of Report 23 reached 118, with conservation broken by up to 34×. At the
fourth iterate the full and half steps met supercritical masers again; the quarter step was accepted on the driver's
natural measure but is worse in the fixed measure (140). Stage 1 is stalled at about 128 with exact conservation. The
run was stopped at the end of the session (05:16) and is resumable (§5).

This is not yet a converging solve. A second measure with no inverse, the row-scaled residual, went
2.32 → 1.39 → 1.69 → 2.60 → 2.85, and every accepted step degrades the deep interior (fixed measure at depths
119–167: 1.5 → 27 → 67 → 76 → 79).

**The narrowest remaining blocker** is the hydrogen ionization front (depths 23–28). Its linear problem is so
ill-conditioned in the preconditioner that a deep-layer direction's far-field radiative action becomes 97–99.6 % of
its preconditioned action. A radiation part that is 22 % of the direction's net exact action, because trapping
cancels most of the rest, is enough. The least squares then fits the front with deep-layer directions and pays for it
in the deep layers (§4.2, §5). Next:
1. Exact Jacobian columns for the front depths in the preconditioner, from one batched probe pass.
2. Then a front-local least squares, compared on the same products and adopted only if the true linear residual holds.
3. Then the constrained solve to convergence, T release, and Gates 1–5.

**What was delivered tonight:**
- CMFGEN build 36: a Newton fallback for the scattering fixed point, with a diagnostic dump; bitwise neutral where the
  ALI converges, verified at one thread.
- Driver options `--pc-svd-floor`, `--conserve`, `--freeze-t` and `--trial-workers`.
- Tests for the floor, the projection and the build-36 record and patch (full suite: 341 passed, 37 skipped).
- The diagnosis of all four defects, with the scripts that reproduce every number (`Reports/evidence/report24/`).
- A summary page (`Reports/report24_page.html`).

`--freeze-t` has no unit test of its own: its directions are built inside the driver loop. A runtime check stops the
driver if a step has a temperature component.

## 1. Starting point

Report 23 ended with stage 1 of the staged closure: statistical equilibrium at fixed temperature, 167 radii, build 35,
local-ALI preconditioner, `results/nk167_trap/fixedT_b35`. The fixed-measure residual (`scripts/fixed_metric.py`, the
preconditioner of iterate 1) had fallen 224 → 153 → 130 → 118 and then stalled, with every full and half step failing
inside CMFGEN. Report 23 attributed the stall to the local-ALI preconditioner, nearly singular where hydrogen and
helium change ionization. The plan was to regularize that preconditioner and resume.

The machine was shared all night with four to five other coding agents (load average 20–50 on 10 cores, swap nearly
full at the start). Every process of this work ran at `nice` 15, the trial passes two at a time until 00:29 and four
at a time afterwards. A native pass that took 11 minutes on the idle machine took 20–45 minutes.

## 2. The preconditioner floor

`--pc-svd-floor TAU` (`scripts/cnk_trap.py:floor_blocks`, `Iterate.floor_preconditioner`). Each diagonal block of the
preconditioner is equilibrated as `BlockTridiag.factor` does (rows by their largest entry, then columns). Its singular
values below TAU times the largest are raised to that floor, with the singular vectors kept, and the scaling is
undone. Well-conditioned blocks are unchanged. A nearly singular direction is amplified by at most 1/TAU.

Measured against exact Jacobian products at eight archived states (`Reports/evidence/report24/floor_q.py`,
`floor_q.json`). q = ‖W(M⁻¹Jz − z)‖/‖Wz‖ for each exact pair (z, Jz); 1 means the preconditioner explains none of
the action.

| state | fixed-J q med / p90 / max | local ALI q med / p90 / max | ALI, floor 1e-2: q med / p90 / max | ‖W M⁻¹R‖: ALI → floored |
|---|---|---|---|---|
| fixed-T it1 | 0.69 / 1.16 / 1.6 | 0.60 / 2.37 / 5.0 | 0.34 / 0.93 / 0.96 | 216 → 25 |
| fixed-T it2 | 0.73 / 1.28 / 3.4 | 0.51 / 3.23 / 11.4 | 0.40 / 0.88 / 8.4 | 302 → 42 |
| fixed-T it3 | 0.72 / 1.36 / 2.0 | 0.63 / 4.23 / 19.7 | 0.55 / 0.73 / 0.93 | 576 → 240 |
| fixed-T it4 | 0.66 / 2.38 / 5.8 | 0.70 / 8.83 / 10.3 | 0.45 / 0.80 / 0.96 | 1204 → 270 |
| fixed-T it5 | 0.69 / 3.85 / 34.0 | 1.21 / 19.24 / 76.4 | 0.49 / 0.86 / 0.90 | 3477 → 126 |
| T-released it1 | 0.64 / 1.54 / 6.9 | 0.44 / 3.93 / 54.4 | 0.34 / 3.06 / 6.6 | 216 → 25 |
| T-released it2 | 0.70 / 2.79 / 7.3 | 0.77 / 24.3 / 236.8 | 0.41 / 4.27 / 9.8 | 873 → 30 |
| T-released it3 | 0.67 / 1.17 / 4.3 | 0.41 / 6.86 / 40.9 | 0.36 / 1.60 / 6.2 | 192 → 28 |

TAU = 1e-2 was the largest value tried (1e-6 … 1e-2). It raised 805–854 singular values per state (of 167 blocks) and
gave the best median and worst case at every state. The floor does not make the preconditioner exact: the T-released
p90 stays at 1.6–4.3. It removes the amplification that inflated the measure 10–30× (Report 23 §2.9). Tests:
`tests/test_cnk_trap.py::test_floor_blocks_bounds_amplification_and_keeps_well_conditioned_blocks`.

With the floor, iterate 4 of stage 1 was restarted at its best base (`it04_trial_l0.25`, fixed measure 118.2). Its
46-direction least squares reached a linear residual of 0.051. The native trials then completed (build 35), and all
three raised the fixed measure: λ = 1: 194.7, ½: 158.6, ¼: 137.8. The driver accepted ¼ on its natural measure
(576 < 755). The next two sections explain both results.

## 3. Why full steps failed inside CMFGEN

Every failed trial of Report 23 and of this night ends the same way (`OUTGEN`): `QS_FRM_ALI: scattering not
converged` after 2000 formal solutions at one to three infrared frequencies, then `Mean intensity blowing up` or an
invalid (NaN) formal J.

**What the trial state contains.** CMFGEN's negative-opacity log (`NEG_OPAC`) of the failed trial
`it04_trial_l0.5` lists inverted H I lines at the He front, at the frequencies that fail:

| frequency (10¹⁵ Hz) | transition(s) | depths |
|---|---|---|
| 0.028547 | H I 12→8 (and He II 24→16) | 60–65 |
| 0.027749 | H I 16→9 | 60 |
| 0.025982 | H I 15→9 (and He II 30→18) | 60 |

The inversions sit at the even (inserted) depths 60, 62 and 64. The trial step changed H I n = 12 by +68 to +82 %
across depths 60–65 and n = 8 in a zig-zag (×0.91, 1.46, 1.10, 1.69, 1.03, 1.15), so n_u g_l / (n_l g_u) of 12→8 went
from 0.56–0.81 to 1.03–1.09 at the even depths (`Reports/evidence/report24/sawtooth_check.json`).

**Why the transfer fails there.** With the `SRCE_CHK` option CMFGEN replaces a negative line opacity by
CHI = 0.1·ETA·CHI_NOSCAT/ETA_CONT (`cmfgen_sub.f:2352`). That keeps CHI positive but can leave it below the
electron-scattering opacity CHI_SCAT, so the net absorption is negative. In the last scattering iterate at 0.028547
(`QSTARS_ALI_TRACE`), CHI_SCAT/CHI is 1.07, 1.58 and 2.49 at depths 60, 62 and 64
(`albedo_profile_it04_trial_l0.5_f35609.txt`). The local ALI step then has no contracting form (CHI_SCAT·L* ≥ 1),
and it falls back to plain Λ iteration.

**Build 36: a Newton solve of the scattering fixed point** (`backends/cmfgen/ali_newton.patch`, capability flag
`qstars_ali_newton_fallback`, `new_main/mod_subs/comp_j_blank.f` only). Where the local ALI has not converged after
`QS_ALI_MAX` formal solutions, the fixed point J = Φ(ETA + CHI_SCAT·J) is solved by Newton's method. dΦ/dJ is assembled
from one-sided differences (relative step 1e-6) of the same formal solution, ND + 1 formal solutions per step, and
solved with dense LU. At most 10 steps are taken, each halved (down to 1e-3) to keep J positive. The result is accepted
only if every formal J is finite and nonnegative and the relative residual is below `QS_ALI_TOL`. Otherwise TC and the
ray store are recomputed from the last ALI iterate, which is the build-35 result. Opacities, emissivities, the
negative-opacity option and the fixed point are unchanged. Frequencies where the ALI converges execute no new
arithmetic. The tangent already solves its scattering system by GMRES, which does not need a contracting local
operator (`formal_tangent.f90:scattering_gmres`). `QSTARS_ALI_NEWTON_DUMP` writes the first Newton Jacobian of each
fallback frequency for offline analysis. Tests: `tests/test_cmfgen_ali_newton.py` (manifest chain, patch shape,
reproduction of the source from build 35, capability record).

**Result on the failed trial** (`results/nk167_trap/b36_tests/it04_trial_l0.5_b36`, `QSTARS_ALI_NEWTON`): the fallback
ran at all three frequencies and failed at all three ("iterations" twice, "positivity" once). The pass then ended as
before. The dumped Jacobians (`Reports/evidence/report24/newton_dump_analysis.py`, `newton_dump_it04_trial_l0.5.json`)
say why:

| frequency | spectral radius of dΦ/dJ | depths carrying the Perron vector | depths where the fixed point J* < 0 | min J*/J |
|---|---|---|---|---|
| 0.028547 | 1.063 | 62–66 | 61–67 | −0.0019 |
| 0.027749 | 1.058 | 59–61 | 57–62 | −1.94 |

dΦ/dJ is (up to interpolation terms) a nonnegative operator. A nonnegative fixed point for a positive source exists
only if its spectral radius is below 1. At these frequencies it is above 1: the trial opacity makes the He front a
supercritical maser, and no steady radiation field exists. That is a property of the trial state, which no solver
of the transfer can remove. The remedy is on the step side (§4–5).

**Neutrality.** The same converged pass (`it04_trial_l0.25`) run with builds 35 and 36 at one thread
(`results/nk167_trap/b36_tests/neutral_1thread_b35`, `_b36`): residual and BA traces byte-identical, ALI statistics
identical, no fallback used. At two threads the build-36 re-run differed from the original build-35 pass by 7e-13 of
each row's maximum, with 645 against 814 iterations in the pass's second J computation. That is OpenMP run-to-run
variation of the ALI iteration path, not the build.

## 4. The conservation defect

The model's equations include three exact linear constraints at every depth. The H II row is hydrogen number
conservation, the He III row helium number conservation, and the n_e row charge balance. At iterate 4 their residuals
equal ΣH − H₀, ΣHe − He₀ and (n(H II) + Σn(He II) + 2n(He III) − n(H⁻)) − n_e to every printed digit at eight depths,
with H₀, He₀ the totals of x0. A Newton step of the linearized system keeps them satisfied, because they are
linear. The driver then applies the update
x·f(λs) with a per-component limit |log f| ≤ ln 2 (`--component-log-cap 0.693`). Clipping one component of a
constraint and not the others breaks the constraint. Nothing restored it.

Measured against the totals of x0, which satisfies all three to 1e-9 (`cnk_trap.conservation_totals`):

| run / state | H total / x0 | He total / x0 | (H II + He II + 2 He III − H⁻) / n_e |
|---|---|---|---|
| fixedT_b35 start_base (transferred) | 0.885–1.005 | 0.743–1.005 | 0.985–1.000 |
| every run, it01_trial_l0.5 (first accepted step) | 0.887–1.451 | 0.753–2.396 | 0.828–2.369 |
| fixedT_b35 it04_trial_l0.25 (fixed measure 118) | 0.880–2.350 | 0.128–1.236 | 0.059–33.9 |
| free_b34 it03_trial_l0.25 (Report 23's "87") | 0.884–1.519 | 0.500–1.349 | 0.394–6.82 |
| free_b34 it06 trials | 0.89–1.94 | 0.69–3.43 | 0.24–27.4 |
| free_ali it04 trials | 0.92–1.39 | 0.65–1.95 | 0.35–1.65 |

The worst damage is at the inserted even depths of the deep interior. In fixedT_b35, from iterate 2 on, n_e at depths
160, 162, 164 and 166 halved every iterate (the limiter's factor): 3.96e13 → 1.98e13 → 9.89e12 → 4.95e12 cm⁻³ at
depth 166. The odd depths beside them did not change. He/H at depths 124–166 fell from 0.107 to 0.012–0.061.

The cause is in the steps. At iterate 1 the preconditioned residual at those depths asked for a 2 % correction of
He II 1s and n_e. The least-squares combination of the 18 directions put s = +4.4 there, 200× larger and past
s = 1, where the linear update would make the population negative. The limiter halved He II and n_e and could not
move the rest. Each later step repeated this. By iterate 4, the helium conservation residual at those depths was
6.8 times the helium present, the preconditioner turned it into an 88-fold He III correction, and that pulled n_e
the wrong way through the shared direction piece (He II / He III / n_e group). He III itself is within 0.25–2.4× of
its Saha ratio to He II at those depths; the 88-fold correction is the conservation residual's image, not a
He III deficit.

**What it did to the measure** (fixed measure by depth region, `Reports/evidence/report24/collect_page_data.py`):

| state | 1–22 | 23–45 (H fronts) | 46–75 (He front) | 76–118 | 119–150 | 151–167 | total |
|---|---|---|---|---|---|---|---|
| start_base | 11.3 | 350.0 | 12.9 | 19.6 | 0.2 | 0.1 | 351.0 |
| it01_trial_l0.5 | 8.3 | 212.6 | 37.6 | 59.0 | 1.5 | 0.5 | 224.0 |
| it02_trial_l1 | 5.6 | 137.7 | 27.0 | 48.9 | 30.0 | 16.8 | 152.6 |
| it03_trial_l0.25 | 4.8 | 100.6 | 17.7 | 50.3 | 55.6 | 29.9 | 130.3 |
| it04_trial_l0.25 | 4.3 | 77.7 | 12.4 | 50.1 | 66.0 | 29.9 | 118.2 |

The hydrogen fronts improved 350 → 78. The deep interior, converged at the start, went from 0.3 to 72.

### 4.2 Where the large step components came from

`Reports/evidence/report24/leakage_check.py` (`leakage_check.json`) solves each iterate's least squares again from the
run's exact products and reports, for every direction, its coefficient and the share of its preconditioned action
W M⁻¹ J z that lies outside the depths where z is nonzero.

| subspace | linear residual | largest coefficient (direction: depths, rows) | its action outside its depths | largest action at | max \|s\| at the deep even depths |
|---|---|---|---|---|---|
| iterate 1, raw local ALI (18 products) | 0.240 | 172 (depths 120–167, He II–n_e) | 97 % | depth 24 | 17.7 |
| iterate 1 products, floor 1e-2 | 0.845 | −105 (depths 1–9, He II–n_e) | 20 % | depth 9 | 3.1 |
| iterate 4, floor 1e-2 (46 products) | 0.051 | 358 (depths 120–167, T only) | 17 % | depth 166 | 407 |

At iterate 1 the deep direction was used to fit the hydrogen-front residual (H I n = 2, n_e, He II 1s and H II rows at
depths 24–28) through its far-field action, which only the nearly singular front blocks of the preconditioner make
large. A Levenberg–Marquardt penalty on ‖W s‖ does not remove it: the direction's own weighted norm is small, and
penalties large enough to stop it raise the linear residual from 0.24 to 0.54–0.91. With the floor, the deep
directions of iterate 4 have coefficients 1.0–1.4 and leak 1–7 %. The floored metric needs its own directions: the
iterate-1 products, made for the raw metric, fit it only to 0.845.

### 4.3 Temperature in the fixed-T stage

`--pin-all` replaces the RE row by c(T − T_fix)/T_fix at every depth, and T stays a variable of the step. Its
significance weight is small (`--pc-t-floor 1e-3`), so the least squares can move T to fit other rows. Until iterate 3
it did not (every accepted base has T = T_fix to 1e-4). At iterate 4 the deep temperature-only direction got the
coefficient 358 (table above), the step put |s_T| = 0.70 at depth 166, and the full-step trial had T at 0.50–1.22 of the
base. `--freeze-t` removes T from every direction: residual pieces and previous steps have their T row zeroed, reused
products with a T component are skipped, and the driver stops if a step has a T component.

## 5. The fix, and the constrained restart

`--conserve` (`cnk_trap.project_conservation`): every line-search trial state is projected onto the constraints
before it is run. At each depth the H populations (H I levels and H II) are scaled by one common factor so their total
equals x0's, likewise all He populations, and n_e is set to the charge sum. Population ratios within a species, T and
H⁻ are unchanged. A `--lams` value of 0 tries the projected base itself. The projection changes no equation and no
solution, since the converged state must satisfy the constraints exactly. It undoes what the limiter did to the trial
state, and the merit decides acceptance as before. Tests:
`tests/test_cnk_trap.py::test_project_conservation_restores_totals_and_charge_and_keeps_ratios` (totals and charge
exact, ratios, H⁻ and T unchanged, idempotent, a constrained state unchanged).

**First use, at the broken iterate-4 base** (build 36, 00:29–01:13). Trials λ = 1, ½ and ¼ of the iterate-4 step,
projected, failed natively: a supercritical maser at 0.36 × 10¹⁵ Hz (CHI_SCAT/CHI up to 1.70). Build 36's fallback
converged the slow-ALI frequencies at 3.19 × 10¹⁵ Hz in two Newton steps (residual 1e-12) and correctly failed at the
maser. The projected base itself (λ = 0) completed with fixed measure 171.0 against 118.2 (H fronts 78 → 143). Its
populations had adapted to totals off by up to 8× (He) and 33× (n_e); restoring the totals unbalances the rate
equations. A state that far off the constraints cannot be repaired by projection.

**Restart at iterate 1** (01:20, `it01_trial_l0.5`, the first accepted state: He ≤ 2.4×, charge ≤ 2.4× off, deep
interior 1.5 / 0.5) with the floor, `--conserve` and `--freeze-t`. The iterate-2 to -5 directories of the
unconstrained trajectory are kept in `fixedT_b35/branch_unconstrained_it2plus/`.

Trials of the constrained solve (`Reports/evidence/report24/restart_summary.py`, `restart.json`; fixed measure with the
iterate-1 reference preconditioner as in §4; "deep" is depths 119–167; the conservation column is the largest
|factor − 1| of H, He and charge):

<!-- constrained-table:start -->
| iterate | trial | native pass | driver natural merit | fixed measure | H fronts (23–45) | deep (119–167) | conservation |
|---|---|---|---|---|---|---|---|
| 1 | λ = 1 (accepted) | completed | 18.69 | 179.2 | 163.7 | 27.1 | 2e-15 |
| 1 | λ = 0.5 | completed | 23.48 | 201.5 | 188.6 | 13.7 | 1e-15 |
| 1 | λ = 0.25 | completed | 27.59 | 216.1 | 204.4 | 7.3 | 1e-15 |
| 1 | λ = 0 | completed | 34.56 | 238.9 | 228.4 | 1.4 | 1e-15 |
| 2 | λ = 1 (accepted) | completed | 22.18 | 144.2 | 103.4 | 67.1 | 2e-15 |
| 2 | λ = 0.5 | completed | 25.42 | 142.8 | 111.0 | 53.7 | 1e-15 |
| 2 | λ = 0.25 | completed | 29.04 | 151.8 | 128.5 | 40.4 | 1e-15 |
| 2 | λ = 0 | completed | 34.53 | 177.3 | 161.8 | 27.1 | 3e-16 |
| 3 | λ = 1 (accepted) | completed | 19.26 | 128.1 | 83.4 | 75.9 | 1e-15 |
| 3 | λ = 0.5 | completed | 23.21 | 134.2 | 91.3 | 70.9 | 1e-15 |
| 3 | λ = 0.25 | completed | 25.80 | 138.2 | 96.1 | 68.8 | 1e-15 |
| 3 | λ = 0 | completed | 33.80 | 142.4 | 101.4 | 67.1 | 4e-16 |
| 4 | λ = 1 | failed: maser, CHI_SCAT/CHI 3.11, no solution | — | — | — | — | — |
| 4 | λ = 0.5 | failed: maser, CHI_SCAT/CHI 2.15, no solution | — | — | — | — | — |
| 4 | λ = 0.25 (accepted) | completed | 57.14 | 139.9 | 84.1 | 78.5 | 9e-16 |
| 4 | λ = 0 | completed | 66.72 | 127.9 | 83.1 | 75.9 | 3e-16 |
<!-- constrained-table:end -->

In iterates 1–3 every trial pass completed natively. Every state conserves H, He and charge to roundoff. The steps are far tamer than in the unconstrained solve: largest |s|
9.8, 16.1 and 12.3 (against 17.7 and 407), limiter active on 278, 116 and 155 variables (against 2309 at iterate 4). The
full step was accepted each time; the hydrogen fronts went 213 → 164 → 103 → 83 and the total 224 → 179 → 144 → 128 in
the fixed measure. The least squares fitted less well at iterate 3 (linear residual 0.455 after three passes, against
0.206 and 0.080). Projecting a base without a step (λ = 0) never helped: it is the worst trial of every iterate,
and at iterate 1 it is worse than the unprojected base (239 against 224).

A second measure without any inverse does not show steady progress. The row-scaled residual (each equation's
imbalance relative to its gross rate sensitivities; `rowscale_by_region.py`, `rowscale_by_region.txt`) went
2.32 → 1.39 → 1.69 → 2.60 along the constrained path. After iterate 3 its hydrogen-front part rose 0.64 → 2.19 while the
preconditioned measure there fell 103 → 83. The two measures weight trapped transitions very differently (Report 23
§2.7), and both must reach zero at the solution. So the constrained solve is not yet demonstrably converging; it is a
solve that keeps the constraints and completes full steps. Along the unconstrained path the row-scaled residual went
2.32 → 7.43 while the preconditioned measure fell 224 → 118. That "progress" was progress in a measure that hid what the
solve broke.

**What is still wrong.** The deep interior gets worse with each accepted step (fixed measure 1.5 → 27 → 67 at depths
119–167). The least squares is still using deep directions through their far-field action at the hydrogen front. At
constrained iterate 2 the deep hydrogen piece (depths 120–167) has 99.6 % of its preconditioned action at depth 23
and the coefficient −13.6. The preconditioned residual there asks for a 14 % decrease of H I n = 2; the step gives a
3× increase. The far-field action is partly physical. In the row-scaled exact action of that direction, the direct and
radiation parts nearly cancel inside its depths (0.418 and 0.421 before cancellation, the trapping of Report 23), and
the radiation reaching other depths (0.018) is 22 % of what remains (`floor_sweep_it2.py`, `floor_sweep_it2.txt`). The
preconditioner turns that into 99.6 %. Floors from 3e-2 to 3e-1 do not fix it cleanly: the coefficient stays near −12
at 3e-2, and at 1e-1 the deep components fall to 0.79 while the linear residual rises to 0.50. Restricting each
direction's action to its own depths in the least squares (`local_lsq_test.py`, `local_lsq_test.txt`) cuts the deep components 3.7 → 1.2 but raises
the true linear residual 0.080 → 0.18, and the largest front component to 114. It is not adopted.

**Iterate 4: the stall returns** (04:13–05:15). With 56 directions the least squares reached 0.173. The full and half
steps met supercritical masers at 0.0585 and 0.0342 × 10¹⁵ Hz (CHI_SCAT/CHI 3.11 and 2.15). Build 36's Newton fallback
found no nonnegative solution there, as it should. The quarter step completed and was accepted on the driver's natural
measure (57.1 < 66.7). In the fixed measure it is worse than its base (139.9 against 128.1). The row-scaled residual
also rose (2.60 → 2.85), its hydrogen-front part to 2.47. Stage 1 is stalled again, now with exact conservation, at the
hydrogen front.

The run was stopped at 05:16, at the end of the session. `fixedT_b35/history.json` holds iterates 1–4; the partial
iterate-5 passes are in `superseded/it05_started_0515`. The best base in the fixed measure is `it04_trial_l1` (128.1).
Resume from it after the front fix of §0: trim the history to iterates ≤ 3, move `it05_trial_*` aside and start
iterate 4 there.

## 6. Gates

No gate changes status; Report 23 §8 stands.

| Gate | Status |
|---|---|
| 167-radius T-released closure converged | **FAIL**: not converged. Stage 1 (fixed T) not converged; see §5 |
| 1. Post-convergence cleanup and audits | NOT EVALUATED |
| 2. 333 radii | NOT REACHED |
| 3. Doubled frequency | NOT REACHED |
| 4. Self-selected CMFGEN wind | NOT REACHED |
| 5. Node-coupled conservative evolution | NOT REACHED |
| 6. Ball et al. 2011 benchmark | **FAIL** (Report 23 §6, unchanged) |
| v1 | **NO** |

## 7. Withdrawn and corrected statements

- *Report 23 §0 and §2.8*: "the corrected solves cut the residual in a fixed measure to about half (224 → 118)" and,
  for the T-released solve, "to 87 after the second iterate". The measure did fall. The states behind those numbers
  violate number and charge conservation by large factors (§4), and the deep interior was damaged while the fronts
  improved. They are not progress toward the closure in the sense the report implied.
- *Report 23 §0*: "full Newton steps break CMFGEN's formal solution" because of the nearly singular preconditioner.
  The native failures are supercritical masers in the trial states (§3). The preconditioner's singularity was real
  and is bounded by the floor (§2), but it was not the mechanism of the native failures.
- *Report 23 §2.9*: "regularizing that preconditioner, then the fixed-T and T-released solves" as the remaining path.
  Regularizing it was necessary and not sufficient; the conservation defect had to be fixed first.

**Failed or superseded attempts preserved:**
- build 36's first test pass without the dump (`b36_tests/it04_trial_l0.5_b36_nodump_stopped`) and a superseded
  neutrality pass (`b36_tests/superseded`), stopped when the dump option was added;
- the unprojected iterate-4 trials and the start of iterate 5 (`fixedT_b35/branch_b35_noconserve`, with the full
  history);
- the pre-floor iterates 5–6 (`fixedT_b35/branch_pre_floor`);
- the unconstrained iterates 2–5, the projected iterate-4 trials and their 1/16 fallbacks
  (`fixedT_b35/branch_unconstrained_it2plus`);
- two driver starts stopped minutes in: the constrained iterate 1 before `--freeze-t`
  (`superseded/it01_from_it01_trial_l0.5_killed_0120`) and iterate 2 before `--col-workers 6`
  (`superseded/it02_restart_0230`).

## 8. Reproduction

```bash
# preconditioner floor, offline (exact products of archived states)
.venv/bin/python Reports/evidence/report24/floor_q.py
# conservation and the fixed measure by region, for any passes of fixedT_b35
.venv/bin/python Reports/evidence/report24/collect_page_data.py start_base it01_trial_l0.5 it04_trial_l0.25
# a failed trial re-run with build 36 and the Jacobian dump, then the spectral analysis
.venv/bin/python Reports/evidence/report24/replay_trial.py results/nk167_trap/fixedT_b35/it04_trial_l0.5 OUT --exe results/nk167_trap/exe_b36/cmfgen_dev.exe --env QSTARS_ALI_NEWTON_DUMP=1
.venv/bin/python Reports/evidence/report24/newton_dump_analysis.py OUT
# the iterate-4 least-squares step and its singular values; far-field leakage of the directions (iterates 1 and 4)
.venv/bin/python Reports/evidence/report24/step_svd.py
.venv/bin/python Reports/evidence/report24/leakage_check.py
# constrained iterate 2: floor sweep and exact-action split; region-local least squares
.venv/bin/python Reports/evidence/report24/floor_sweep_it2.py
.venv/bin/python Reports/evidence/report24/local_lsq_test.py
# the constrained restart: trial table (restart.json) and the row-scaled second measure
.venv/bin/python Reports/evidence/report24/restart_summary.py
.venv/bin/python Reports/evidence/report24/rowscale_by_region.py
# summary page
.venv/bin/python Reports/evidence/report24/build_page.py
# the constrained stage-1 solve as run from 2026-10-22 01:20 (build 36), resumed by the same command at its latest base
results/nk167_trap/launch_fixedT_b36.sh it01_trial_l0.5 1 --pc-svd-floor 1e-2 --trial-workers 4 --eta 0.1 --conserve --freeze-t --lams 1 0.5 0.25 0
```
`results/nk167_trap/fixedT_b35/README.txt` lists every start, stop and directory move of the run.
