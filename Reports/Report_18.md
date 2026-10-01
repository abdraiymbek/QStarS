# Report 18 — measuring the Newton linearizations against the exact formal-J residual

**Outcome: NUMERICAL BLOCKER, now measured.**

- **Diagnosis.** I tested every available CMFGEN linearization (native moment response, and the
  local formal-ALI response at strengths 0, 0.5 and 1) for whether it predicts the change of the
  *exact* rate-equation residual along its own Newton step. None does. Three of the four
  point uphill.
- **Mechanism.** In the cool layer, a change in populations is almost cancelled by the
  radiation field's response. With the radiation response left out, the predicted residual
  change is about 16 times too large. The local and moment responses get this cancellation
  wrong.
- **First working method.** A Newton–Krylov step built only from exact-residual evaluations
  (GMRES with 12 finite-difference products) is the first method that lowers the exact
  residual without diverging. The gain is 11–26% per outer step, at about 2 hours per step.
- **Not delivered.** No released cool closure, no energy/electron gates, no resolution pair,
  no validated Jacobian, no self-selected wind and no production mdot.

QStarS remains **0.1 development, not v1**.

Starting point is Report 17. The radiation state is still the nonnegative formal
characteristic solution (build 13 onward). All runs in this report start from the same
state x0, the endpoint of the 40-iteration fixed-T Λ stage
(`results/cmfgen_ali_lambda84_long_maxlam1e3`, 84 radii). Required physics is unchanged.
Every hook added in builds 16–17 is diagnostic only and sits behind an environment
variable.

## 1. Inexact Newton with a scaled local response (`QS_LIN_FAC`) — closes Report 17 §4.5

Build 15 (`103640e8…`), `QS_ALI_LIN=T`.

| Run | Settings | Result |
|---|---|---|
| `cmfgen_alilin_fac05_fixT` | fac 0.5, T fixed | Steps stay mild (scale ≈ 1, 30–60% corrections). Major-species max\|STEQ\| sits at 0.3–0.7 over iterations 56–69 and does not decrease. Stall; killed at iteration 69 |
| `cmfgen_alilin_fac05_release` | fac 0.5, T free | He II trace levels run away at d11–12: max increase 1e5% → 3.7e7% by iterations 57–63. Killed |
| `cmfgen_alilin_fac075_release` | fac 0.75, T free | Same runaway, 1.1e7% at d12 by iteration 63. Killed |

Scaling the local response down removes the first-step blow-up (Report 17 §4.4), but it
does not give convergence. Released temperature reproduces the He II trace runaway seen
with full strength.

## 2. Does each linearization predict the exact residual change? (build 16)

Build 16 sha256 `a155d4270812ed6f608f2994c000b026327b58c2bd67d48ddd99338e8838312c`;
`backends/cmfgen/jac_check.patch`. It adds two diagnostic hooks:

- a dump of the raw assembled residual R just before the block-tridiagonal solve
  (`QSTARS_RESID_TRACE.bin`);
- a one-shot perturbation at the first iteration, POPS ← POPS·(1−ε·s), after which all
  derived state (ion sums, n_e, level dissolution, full-atom arrays) is rebuilt.

**Test.** CMFGEN solves BA·s = R and then updates POPS ← POPS·(1−s). If BA is right along s,
then R(x0·(1−εs)) ≈ (1−ε)·R0. The error measure is
e = ‖ΔR + εR0‖ / ‖εR0‖: e ≪ 1 means accurate along s, e ≥ 1 means useless.

- The steps s are the first BA corrections from x0 (record 41 of each run).
- I picked ε so that max|εs| = 1e-2 and 1e-3.
- All four runs are T fixed. The perturbed state was reproduced to ≤3e-16.

Script: `scripts/jacobian_direction_check.py`. Output: `results/jac_check/jacobian_direction_check.json`.

| Linearization (source run) | max\|s\| | e (ε→1e-2 / 1e-3) | cos(ΔR, −εR0) | ‖R‖ along step (prediction 1−ε) | Nonlinearity |
|---|---|---|---|---|---|
| frozen J, fac 0 (`cmfgen_alilin_fac0`) | 0.98 | 0.939 / 0.939 | **+0.993** | 0.99938 (pred 0.98983) | 1.6e-2 |
| local ALI, fac 0.5 (`cmfgen_alilin_fac05`) | 1.14 | 1.72 / 1.72 | **−0.988** | 1.0063 (pred 0.9912) | 9.6e-4 |
| local ALI, fac 1 (`cmfgen_alilin_fixT84_b`) | 5.4e3 | 50.1 / 50.3 | **−0.979** | 1.00009 | 6.0e-2 |
| native moment (`cmfgen_ali_linfixT84_fromlam`) | 75.9 | 117 / 117 | **−0.767** | 1.0118 | 4.1e-4 |

All four are in the linear regime: e agrees between the two ε levels, and the nonlinearity
measure is ≤ 6e-2. So these are properties of the Jacobians themselves.

- **Frozen J.** This direction goes downhill, but the exact residual falls by only about 6%
  of the predicted amount. The radiation response cancels about 94% of the population effect.
- **Local ALI (0.5 or 1) and native moment.** These directions go **uphill** on the exact
  residual. The fac-1 error is concentrated at d11–15, where the Report 17 blow-up sits. The
  fac-0.5 error peaks at d71–75.

**Correction to Report 17 §4.4.** Report 17 said the local operator made the rate equations
nearly singular on trapped lines, and called that the cause. These measurements support a
narrower statement. Even with no radiation response at all, the linearization is about 16×
too strong along its step (e = 0.94). Every analytic response available here then gets the
near-cancellation wrong in sign or size. The q = 1 resonance-line trapping that Report 17
measured is one part of that cancellation, not a separate cause.

## 3. Newton–Krylov on the exact residual (build 17)

Build 17 sha256 `bff4b237a495ae2c51603bc9f4361bc1be6d53ff7215e1ca32b049cf8830efb2`;
`backends/cmfgen/krylov.patch`. It adds two more hooks:

- with `QSTARS_KRYLOV_R0` set, the right-hand side becomes (R(x0·(1−εv)) − R0)/ε;
- the raw block-tridiagonal solution is dumped (`QSTARS_SOL_TRACE.bin`) and the run stops
  after the solve.

**Method.** One native run gives one product M⁻¹J·v. M is the frozen-J BA (fac 0); J is the
exact Jacobian, applied by finite difference. GMRES is run offline with
`scripts/jfnk_krylov.py`, which I checked on a synthetic operator (it recovered the solution to
2e-5 relative). The equation solved is M⁻¹J(x0·s) = M⁻¹R0 = s_fac0.

**Noise floor.** Two identical build-17 runs give exact residuals that differ by up to
4.6e-5 relative. Most cells differ by 1e-12 to 1e-9; the cancellation-dominated ground-state
ionisation rows amplify this. The solution reproduces to 1.4e-9 absolute. I used the
build-17 residual as R0 for every product. This floor limits *native* finite-difference
Jacobian validation to roughly 1e-4 to 1e-3. The 1e-6 target can only be checked in the
operator fixtures.

**GMRES convergence (preconditioned relative residual).**

| Vectors | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|
| GMRES residual | 0.593 | 0.132 | 0.0793 | 0.0792 | 0.0434 | 0.0416 | 0.0387 | **0.0326** |

- The k = 5 solution has coefficient y₁ = 38.6 on the frozen-J direction, so the exact Newton step
  is tens of times longer than the frozen-J step. This is the same cancellation as in §2.
- The resulting steps are large: max s = 1.8 at k = 4 and 8; at k = 12, max s = 4.8 and
  min s = −9.0. With the linear update (1−λs), λ = 1 makes populations negative; at k = 12
  even λ = 0.25 does. I therefore also tried the always-positive update x0·exp(−λs).

**Exact residual after the step** (`results/jfnk/step/evaluate.json`). Two measures:

- **Raw ‖R‖** is 99.95% H I high-level rows 20–22 at the innermost depths d76–83.
- **The frozen-J correction ‖s_fac0(x)‖ in the cool layer d10–30** is dimensionless and carries
  98% of ‖b‖. This is the more meaningful measure for the cool-closure problem.

| Step | Raw ‖R‖/‖R0‖ | Frozen-J correction, d10–30, relative to x0 |
|---|---|---|
| frozen J, λ = 1 (baseline) | 0.939 | 0.931 |
| Krylov k=4, λ = 0.5, linear update | 0.866 | 0.926 |
| Krylov k=8, λ = 0.5, linear update | 0.880 | **0.890** |
| Krylov k=8, λ = 1, exp update | 0.760 | **0.886** |
| Krylov k=12, λ = 0.5, exp update | **0.736** | 0.966 |
| Krylov k=12, λ = 0.1, linear update | 0.947 | 0.992 |

With the exp update, the residual relative to its own size at x0 (the row-scaled measure)
grows 5–58× in a few cool-layer cells (d18–19). The large trace-level changes in those steps
are therefore not yet trustworthy.

**What this shows.** A Newton step built from exact-residual products works in the linear
sense: GMRES gets to 3e-2 with 12 vectors. But one step reduces the residual by only 11–26%.
The reason is that the exact step changes populations by factors of 2 to e⁹, far beyond the
linear range. Repeating this as an outer loop would cost about 1.5–2 hours per step. It is
feasible as a demonstration, not as a production solver. I did not iterate it, so no
convergence claim is made.

## 4. Gate decision

| Gate | Status |
|---|---|
| Nonnegative, realizable radiation field | PASS (every transfer pass in builds 15–17) |
| Released cool populations + T + radiation converged | FAIL. Fixed-T fac 0.5 stalls at major-species corrections 0.3–0.7; T-free runs diverge (§1) |
| Electron balance / radiation energy gates | NOT EVALUATED (no converged state) |
| Two spatial and frequency resolutions | NOT ATTEMPTED (no converged 84-radius state) |
| Validated nonlocal Jacobian | FAIL. Measured: every analytic linearization has e ≥ 0.94 along its own step, three of four uphill (§2). A finite-difference exact product exists but is noise-limited to about 1e-4 (§3) |
| Global mdot / sonic / self-selected wind | NOT ATTEMPTED. No production mdot |

**Physical existence.** Nothing here indicates that the cool solution does not exist. The
exact residual goes down along the frozen-J and Krylov directions, so x0 is not a local
obstruction. Every divergence traces to a Jacobian that I measured to be wrong.

## 5. Narrowest blocker and next steps

**Blocker.** The Newton system needs an analytic radiation response that reproduces the
near-cancellation between populations and the radiation field to better than the 6% that
survives it. Along the dominant direction, the true Jacobian is about 0.06 times the
frozen-J Jacobian. Neither the moment recursion nor a local or scaled Λ* achieves this. The
same response is what the wind eigenvalue problem needs.

Next, in order:

1. **Tangent-linear formal transfer.** Carry dI/dη and dI/dχ along each ray, in space and
   through the comoving-frame frequency term, and fill TX/TVX with them. Two points are
   specific to this code:
   - The local-ALI branch currently zeroes `TX(:,:,3:NM)` at every frequency. Those slots
     accumulate the continuum and line responses across frequency (`do_var_cont.f`
     L231–233 and L442–459), so line variables lose their response at bluer frequencies.
     A correct operator must propagate these accumulators.
   - Validate in the fixtures against finite differences to 1e-6. In the native model,
     validate with `jacobian_direction_check.py`, targeting e ≪ 1 (the floor is about 1e-4).
2. **Preconditioner.** Use that operator to precondition the Krylov solve of §3. GMRES should
   then need 2–4 vectors instead of 12 or more. Add a damped outer loop in log variables
   with a merit function on the cool-layer frozen-J correction.
3. Only after a converged 84-radius release: run the electron/energy audit
   (`scripts/audit_cmfgen_cool.py`), then 167 radii and a doubled frequency grid, then the
   global mdot/sonic eigenvalue problem with the Mach-one MLT guard.

## 6. Reproduction and housekeeping

**Environment.**

- `source ~/qsrd-runtime/mesasdk/bin/mesasdk_init.sh`, `MACOSX_DEPLOYMENT_TARGET=14.0`, two
  threads per run.
- Builds 15, 16 and 17 sit in `backends/cmfgen/build_1{5,6,7}`. Each layers its patch on the
  previous one, and `record_build.py` reports every capability flag as true.
- Run `make d_new_main` twice; the first pass may only archive without relinking.
- Kill runs by PID; BSD `pkill -f "a\|b"` matches nothing.

**Driver scripts.**

- `results/jfnk/drive_krylov.sh` runs the sequential Krylov loop.
- `results/jfnk/launch_eval.sh` runs the line-search residual evaluations.

**Cost.**

- One 84-radius iteration takes 450–820 s.
- One exact product M⁻¹J·v costs one such iteration.

**State.** Nothing was committed to git. Run `pytest` before committing. No native processes
are running.
