# Report 23 — the v1 closure run: why the 167-radius closure did not converge, and Gates 1–6

**Status: QStarS is NOT v1.** The 167-radius T-released cool closure has not converged. This report shows why: the
cause is numerical, not a missing equilibrium of the H/He model. Three defects of the solver and its passes were
found, measured and fixed (builds 32–35 and new driver options). With them the solves make real progress (the residual
halves in a fixed measure) and then stall: the local-ALI preconditioner is nearly singular at the ionization fronts
(§2.8–2.9). No solve has converged. Gates 1–5 cannot close before the closure does; their tooling,
tests and pre-registered rules are in place. Gate 6 (the Ball et al. 2011 evolutionary benchmark) was implemented,
registered and executed. Its smoke gate is a recorded numerical blocker at M_BH = 46 M☉, and its registered verdict is
FAIL (§6).

## 0. Outcome

**The single narrowest remaining blocker:** the preconditioner at the ionization fronts. Every defect that made the
earlier solves wrong is identified and fixed. The corrected solves cut the residual in a fixed measure to about half
(224 → 118) and then stall, because the local-ALI preconditioner is nearly singular where hydrogen and helium change
ionization (depths 23–24, 42–43 and the He front near 60). There it inflates the measure to 1e3–3e3, the least squares
chases those modes, and full Newton steps break CMFGEN's formal solution. Regularizing that preconditioner, then the
fixed-T and T-released solves, then Gates 1–5, is the remaining path (§2.8–2.9).

**What was wrong**, each measured (§2):
1. **The solver could not see trapped transitions.** It preconditioned with fixed-J Jacobian blocks and measured
   progress as ‖W·M⁻¹·R‖. In an optically thick transition the radiation response cancels the fixed-J rates (H I n = 2
   at the front: exact diagonal 3e-4–6.5e-4 of the fixed-J one). The measure therefore understated the error of every
   such level: H I n = 1–9, He I ground and ¹P series, He II n = 1–3. States it called converged were not. The PTC
   iterates of Reports 21–22 raised the cool layer's temperature by 3–6 kK while its ionization stayed at the cool
   values. Neutral hydrogen there was up to 7e10 times its Saha value, and n₂/n₁ was 6e-9 of Boltzmann.
2. **Every single-iteration pass ran the wrong two-photon physics.** `COMP_F=T` silently replaced
   `TWO_METH=USE_RAD` by `OLD_DEFAULT` (`comp_opac.f:123-124`). That gives 33–101× too much 1s→2s absorption at the
   H I front: about 75 % of the radiative-equilibrium (RE) row there, energy-inconsistent, and the row's only
   temperature dependence. With n = 2 unbalanced on top, the RE row was not a thermal balance. CMFGEN's own electron
   balance said the front cools, while the RE row said it heats. Two predictive tests confirmed the diagnosis to 1 %
   (§2.6).
3. **The exact Jacobian action had a wrong temperature column.** The tangent held the ion free-free Gaunt factor
   fixed in dη/dT and dχ/dT. That put the radiation part of the T column off by up to 2× in the infrared and the RE
   row's T column 3.5 % off deep inside the envelope, where direct and radiation parts cancel 15–20×. It became
   decisive once `USE_RAD` made two-photon absorption the main infrared opacity. Found by the build-32 tangent check,
   fixed in build 34 (now at the oracle floor, 1.2e-5) and at the inner boundary in build 35.

**What fixes them:**
- build 32: `USE_RAD` with the state's own stored J;
- a local-ALI preconditioner for every column, made affordable by build 33 (17.8× fewer instructions, identical
  output);
- H I n = 2 made significant, and a linear update map;
- build 34: the Gaunt-factor derivative.

**Effect.** On build 33 the T-released solve stalled. Its steps were depth-oscillating in T, most trials failed
natively, and the measure at successive bases rose 26.7 → 33.4. On build 34 from the same state, every trial completes
and the full Newton step is accepted. The residual in one fixed preconditioned measure falls 224 → 131 in the first
iterate, with every species group and the RE row decreasing, and to 87 after the second iterate. Later iterates and the
staged fixed-T solve stall near half the starting residual (§2.8).

**Numerical or physical?** Numerical. "No root up to 12.3 kK" came from an unconverged n = 2 plus the `OLD_DEFAULT`
term. Whether the front's finite-temperature equilibrium lies near 7–9 kK, as the corrected electron balance suggests,
needs the converged state.

**Earlier statements withdrawn** (§9):
- the front as a "thermally unstable stretch" with "net heating ∝ T" (Report 21 §4.3–4.4);
- the 10.2 kK local root as evidence of a hot equilibrium (this report's first draft);
- the premise that levels below 1e-10 of their element cannot affect other equations (Reports 20–22);
- every statement that a cool state was converged (Reports 20–22).

**Gates** (§8): none of the six passes.
- Gate 6 is evaluated: FAIL.
- Gates 1–5 are not reached.
- Delivered: the blocker's diagnosis and its fixes (builds 32–35 and driver options, tested in `test_cnk_trap.py`,
  `test_cmfgen_fast_ali.py`, `test_cmfgen_ff_gaunt.py` and `test_nlte_newton.py`);
  Gate 1–3 tooling with five defects fixed and pre-registered closing rules; a Gate 4 wind module with 21 analytic
  tests and a Gate 5 node-coupling interface with 11 tests; the Ball replication in seven reviewed stages plus its
  registered execution.
- Tests: 334 passed, 37 skipped (Report 22: 111 and 23).

**Compute.** The request was made at 2026-10-06 21:38. In the first 167 h the laptop was asleep for 154 h (`pmset`
log). The last day ran with the machine kept awake; the solves, five native builds and the Ball campaign ran then.

## 1. Starting point (Report 22)

Report 22 left the 167-radius T-released cool closure (`results/nk167_hot`) on the fast path: build 31, ≈ 28 min per
iterate, user-approved ultra-trace weight 1e-10. The front was heating under pseudo-transient continuation (PTC) and had
not converged. Gates 1–6 had not started.

## 2. The 167-radius released closure

### 2.1 Least-squares weights: three corrections

The Newton-Krylov step minimizes ‖W M⁻¹(R − J s)‖ over the span of exact products.

**Ultra-trace rows.** Levels at ≤ 1e-20 of their element (He II/He III at the front, fractions down to 1e-34) carry
preconditioned corrections that grew as the front heated: 2e6 at iterate 14, 2e8 at iterates 20–23, 9e9 at iterate 29,
3e10 at iterate 34. With any finite weight they eventually dominate the least squares. Offline on the archived exact
pairs of iterates 20, 22 and 23:

| Weight | Front RE-row residual reduced to | Non-ultra LS residual | Front T step |
|---|---|---|---|
| 1e-10 (approved, Report 22) | 0.83–0.95 | 0.33–0.57 | 0.8–1.8 % |
| 1e-12, 1e-14 or 0 (identical) | 0.011–0.017 | 0.07–0.12 | 4.3 % |

The run went to 1e-12 at iterate 25 and to 0 at iterate 31, after a fixed-T diagnostic run stalled with 90 % of its
residual in those rows. These levels are converged afterwards by the trace cleanup of the closing chain (§3). This is
the same solver-path setting the user approved in Report 22; products and the convergence test are unchanged.
Evidence: `ultra_weight_scan_it20_23.txt`, `residual_structure_it22.txt`, `pin_structure.py`.

**The jump at the significance threshold.** Levels above 1e-10 had weight 1, levels just below had 1e-3. He I levels
whose fraction crosses 1e-10 at depths 91–92 flipped category between iterates and carried 20–98 % of the base merit.
Since iterate 33, `--weight-ramp 2` makes the weight a log-linear ramp from 1e-3 at 1e-12 to 1 at 1e-10
(`significance_weights`, unit-tested). Evidence: `merit_classes.py`.

**Neither correction touched the real problem** (§2.5): H I n = 2 sits at ~1e-13, below both thresholds.

### 2.2 Trajectory under PTC (`results/nk167_hot`, iterates 0–34)

| Base | Max significant step | δ | Front max T (K) | Front mean T (K) | Front max RE net/gross | Outside front |
|---|---|---|---|---|---|---|
| start | 2.40 | 0.039 | 7 900 | 7 900 | 7.9e-3 | 2.8e-3 |
| it08 | 0.81 | 0.032 | 8 919 | 8 311 | 6.6e-3 | 1.2e-3 |
| it16 | 0.33 | 0.031 | 9 057 | 8 623 | 6.9e-3 | 9.0e-4 |
| it24 | 0.11 | 0.030 | 9 644 | 9 108 | 7.5e-3 | 9.6e-4 |
| it28 | 0.12 | 0.030 | 10 335 | 9 617 | 8.6e-3 | 9.7e-4 |
| it34 | 0.09 | 0.030 | 10 489 | 9 925 | 9.0e-3 | 1.0e-3 |

The front heated by 2 600 K while its RE imbalance rose. Each step was accepted at λ = 1 with a 10–50× drop of the
natural measure. The run was paused at iterate 34 (`front_progress_nk167_hot.json`). Its direction was set by the
corrupted row (§2.5), so it is not resumed.

### 2.3 Registered tests, and what became of them

- **Stop rule** (`benchmarks/cmfgen_cool_closure_stop_rule.json`, sha256 `a2711719…`; amendment 1 `cfe15815…` counts
  awake hours only, registered when the machine woke from 31.8 h of hibernation). The run was paused by decision, not
  by the rule; the rule did not fire.
- **Coupled existence test** (`cmfgen_cool_existence_test.json`, `c55b1f13…`; amendment 1 `93bb1c8e…`). T was held by a
  pin row (`--hom-mode frontpin`, new, unit-tested), everything else Newton-converged, and the RE imbalance read at the
  front.
  - *Attempt 1* (T free outside the window) did not converge and is kept.
  - *θ = 1.0* (all T pinned at the iterate-29 profile) reached significant step 0.09 at iterate 5; front RE +7.2e-3 to
    +8.8e-3.
  - *θ = 1.15* did not contract (step 3.8, merit rising); indicative RE +9.4e-3 to +11.6e-3 at 12.1–12.3 kK.
  - **The readout is invalid.** The registered test read the sign of the RE row as the thermal direction at
    "statistical-equilibrium populations". §2.5 shows the populations were not in statistical equilibrium (n = 2), and
    the row is dominated by the two-photon artifact. No existence statement is made from either point. Both runs are
    stopped and kept (`results/nk167_pin`).

### 2.4 The contradiction that exposed the problem

One full native iteration at iterate 28, and again at the fixed-T state, gives CMFGEN's per-process electron balance
(GENCOOL): net **cooling** at depths 32–39, 0.2–2.7e-9 erg cm⁻³ s⁻¹, almost all H I collisional. The RE row of the same
states gives net **heating**. With `COMP_EHB=T` the native electron heating balance `STEQ_T_EHB` equals −GENCOOL to four
digits (depth 35: −1.1734e-9), while the RE row is unchanged. Evidence: `gencool_front_it28.json`,
`gencool_front_pin100_it05.json`, `re_vs_ehb_pin100_it05.json`.

### 2.5 Diagnosis

A three-lens investigation (equations, numbers, physics) with an independent synthesis that re-verified each claim
(`re_row_synthesis.md`, `re_row_numbers_analysis.md`, scripts in `re_row_scripts/`).

**No sign or unit error.** Row 134 of the trace is raw `STEQ_T`. Its continuum term is
`FQW·(CHI_NOSCAT·J − ETA_NOSCAT)` (`cmfgen_sub.f:2427-2429`), positive for net absorption. Lines enter as
`−ETAL·ZNET`, and no collisional term enters. In physical units, erg cm⁻³ s⁻¹ = 4π·1e-10·R_T. A direct integral of the
spectral export reproduces the trace to 0.1 %.

**The energy identity.** With R_n the rate-row residuals, E the level energies, and D the two-photon energy defect:

4π·1e-10·R_T = EHB + h·Σ(E₁ − E_n)·R_n + D

It closes to 3.5e-4 at depth 35 (5.9567e-9 against 5.9588e-9); without D it is 24 % short. The RE row equals the
electron balance only when every level row vanishes and D = 0.

**Numbers at the fixed-T state** (`results/nk167_pin/theta1.00/it05_trial_l1`):

| Depth | T (K) | R_T / l1 | 4π·1e-10·R_T | Native EHB | R₃ (n = 2), cm⁻³ s⁻¹ | hν_Lyα·R₃ share of R_T | Two-photon share of R_T | 2γ up-rate: OLD_DEFAULT / exact |
|---|---|---|---|---|---|---|---|---|
| 33 | 10 279 | 7.2e-3 | +2.51e-9 | −4.34e-10 | 147.7 | 0.963 | 76.0 % | 81.3 / 0.80 |
| 35 | 10 387 | 8.75e-3 | +5.96e-9 | −1.17e-9 | 356.9 | 0.979 | 74.9 % | 186.6 / 3.09 |
| 37 | 10 250 | 6.0e-3 | +1.37e-8 | −2.66e-9 | 806.8 | 0.960 | 73.1 % | 401.7 / 12.2 |

- **The row scale is the Lyα gross scattering rate**, not a heating scale: at depth 35, dR_T/dln n₁ = +273.4 and
  dR_T/dln n₂ = −267.2 of l1 = 542. The whole electron budget is 1.7e-3 of it. "|RE net/gross| ≤ 1e-3" therefore
  certifies nothing about thermal balance near the front.
- **n = 2 is frozen.** n₂(d35) stays within 4.3649e-5–4.3754e-5 cm⁻³ over 60 passes spanning 7.7–12.9 kK. In the
  archived exact products the n = 2 row's fixed-J action is +0.41615 and its radiative response −0.41549: Lyα coupling
  cancels 99.8 %. This is Λ-iteration-type false convergence in the frozen-J metric.
- **The two-photon term.** The reconstructed `OLD_DEFAULT` up-rate matches CMFGEN's `TWO_PHOT_SUM` to four digits
  (1.5599e-6 against 1.5601e-6 s⁻¹ per atom). The driver's frozen-J differences give dR_T/dlnT = +3.408, +7.383,
  +14.997 at depths 33/35/37; the derivative of the `OLD_DEFAULT` term alone gives 3.407, 7.380, 14.985. The native
  analytic entry is −4e-3, because `two_phot_var_opac.f` has no dT term.
- **The rest of the row** (24–29 %) is continuum photons redshifting into the trapped Lyman lines. It is real within
  the model and temperature-independent. In statistical equilibrium it leaves through Balmer photoionization of n = 2,
  so it drives ionization and ~1 eV of photoelectric heating per event, not 10.2 eV of heat.
- **The 10.2 kK local root** (`thermal_curve_d35_it25.json`: +1.47e-3 at 9 469 K, −1.21e-3 at 10 889 K) is where real
  Lyα collisional cooling balances the spurious two-photon defect D. It is not a hot equilibrium.
- **Excluded:** state mismatch, `SCL_LN`, adiabatic/advection terms (5e-13 against 6e-9), non-coherent scattering,
  frequency quadrature.

### 2.6 Two predictive tests

**Test 1: frozen-J replay with `TWO_METH=NOSTIM`** (`scripts/replay_with_vadat.py`; `two_photon_nostim_replay.json`).
Predicted from the reconstruction before the run:

| Depth | R_T now | R_T predicted | R_T measured | R₃ now | R₃ predicted | R₃ measured |
|---|---|---|---|---|---|---|
| 33 | 1.994 | 0.48 | 0.477 | 147.7 | 66 | 66.4 |
| 35 | 4.742 | 1.19 | 1.190 | 356.9 | 170 | 170.3 |
| 37 | 10.925 | 2.94 | 2.935 | 806.8 | 405 | 405.1 |

**Test 2: one native pass with n = 2 multiplied by its first-order balance factor** (10–720× at depths 26–46),
`NOSTIM`, `COMP_EHB=T` (`n2_boost_test.json`):

| Quantity at depths 33 / 35 / 37 | Before | Predicted | Measured |
|---|---|---|---|
| Σ R_n, n ≥ 2 (cm⁻³ s⁻¹) | 145 / 351 / 793 | within ±10 / ±25 / ±60 of 0 | 0.2 / 0.3 / 0.4 |
| Ionization-row imbalance | 1e-3 to 8e-3 | ≈ 64 / 165 / 391 | 63.5 / 163.9 / 390.3 |
| Electron balance (erg cm⁻³ s⁻¹) | −4.3e-10 / −1.17e-9 / −2.66e-9 | −3.5e-10 / −0.95e-9 / −2.0e-9 | −3.2e-10 / −8.7e-10 / −1.9e-9 |

At depth 30 (8.1 kK) the electron balance becomes +6.5e-12, as predicted. The thermal balance of an n = 2-balanced
front lies near 8 kK, before the ionization has relaxed.

### 2.7 The wider cause: the solver's measure is blind to every trapped transition

The n = 2 finding is one case of a general defect. The driver preconditions with the fixed-J ("direct") Jacobian
blocks M and measures progress as ‖W·M⁻¹·R‖. A fixed-J block holds the full, untrapped rates of every line and
continuum. Where a transition is optically thick, the radiation response cancels them. M⁻¹R then understates the
error of the levels involved by the trapping factor, and Krylov directions built from M⁻¹R do not contain their
modes. Report 20 §4.2 recorded this cancellation as the reason the linear solve stalls; its consequence for the
convergence measure, and for the rule that levels below 1e-10 of their element are "trace", was not drawn.

Measured on the 167-radius model (`Reports/evidence/report23/trapped_levels/`):

| Measurement | Result |
|---|---|
| Exact H I n = 2 column (21 single-depth exact products) against the fixed-J column, depths 26–46 | Diagonal ratio 3.1e-4–6.5e-4. Depth 35: −9.95 against −2.06e4; nearly all of the exact diagonal is Hα absorption into n = 3 (+9.2) |
| Native local-ALI linearization (`QS_LIN_FAC 1`), one pass | Its n = 2 column equals the exact local column to 0.1–5 % |
| Levels whose local-ALI diagonal is below half the fixed-J one | H I n = 1–9 (ground: down to 1.1e-5), He I ground (1.1e-6) and the ¹P series, He II n = 1–3 (n = 2: 1.2e-6) |
| Largest M⁻¹R over significant variables at the state the old measure called nearly converged (pin state) | 7.7e-3 with fixed-J blocks; 1.63 with the exact n = 2 column (H II and n_e at the front) |
| Cool start state, H I ground at depth 18 | 1.7 with fixed-J blocks; 366 with exact n = 1 and n = 2 columns |

**What the earlier iterates really were.** Departure coefficients b₁ = n₁/n₁(Saha) and the ratio n₂/n₁ against
Boltzmann, from the archived states (`departure_coefficients.txt`):

| State | Depths 20–47 (the layer PTC heated) | Depths 59–130 (interior) |
|---|---|---|
| Cool start (Λ-iterated; T minimum 3.8 kK) | b₁ = 0.16–7.5; n₂/n₁ within 7 % of Boltzmann | b₁ = 1.00–1.27; n₂/n₁ Boltzmann to 1e-3 |
| `nk167_hot` iterate 34 (T = 7.9–10.5 kK there) | b₁ = 3.5e5–6.9e10; n₂/n₁ = 5.7e-3 down to 5.7e-9 of Boltzmann | b₁ = 0.39–0.89; n₂/n₁ = 1.1–2.4 |

The T-homotopy of Report 21 and the PTC of Reports 21–23 raised the temperature of the cool layer by 3–6 kK while
its ionization and excitation stayed at their cool values, ten orders of magnitude from Saha–Boltzmann at the new
temperature, and the measure reported convergence to a few 1e-3. The same holds, less visibly, for every state that
measure called converged: the fixed-T closures of Report 20 and the "full residual" states of Report 21 are converged
in the fixed-J measure only. The cool start itself is near LTE in the interior because it remembers the LTE start of
its Λ iteration, not because statistical equilibrium was solved there.

### 2.8 Solver changes made and tried

All are options of `scripts/cool_newton_krylov.py` (helpers in `scripts/cnk_trap.py`, tests in
`tests/test_cnk_trap.py`); with none of them given the driver is unchanged.

1. **Build 32: the intended two-photon physics in single-iteration passes** (`--two-phot-j`). `USE_RAD` is kept, with
   the conjugate-frequency J read from a stored-J file (`QSTARS_TWO_PHOT_J`). With the variable unset, build 32 is byte-identical to build 31 (residual, blocks, products). With it set, the
   two-photon up-rate at depth 35 is 3.081 cm⁻³ s⁻¹ (186.64 with `OLD_DEFAULT`; 3.081 from a band-exact replica of
   the native algorithm), the energy identity closes without D to 2e-5–9e-4 over depths 30–41, and the exact action
   matches central differences to 2e-9–5e-6 in population directions, where the `OLD_DEFAULT` Jacobian is off by 26 %
   on the n = 2 row. `NOSTIM` is not adopted: the
   true two-photon absorption it drops is 1 % of the row at depth 35 but comparable to collisional excitation at
   depths 39–41. In the driver the residual and blocks of each base come from a frozen-J replay of that base with its
   own stored J as the two-photon J (43 s), and its columns, products and trials use the same file; a converged state
   therefore satisfies the equations with its own conjugate-frequency J. The Jacobian action holds that J fixed.
2. **Exact preconditioner columns from exact products** (`--probe-rows`, `--probe-stride`). One batched pass perturbs
   a level on every k-th depth; the response at the perturbed depth and its two neighbours replaces that level's
   block-tridiagonal column in M only. Contamination by the other perturbed depths: 4.6 % at k = 4, 2.5 % at k = 6
   (median, n = 2).
3. **Local-ALI increment for every column** (`--ali-precond`): M + (BA[`QS_LIN_FAC 1`] − BA[`QS_LIN_FAC 0`]) from one
   frozen-J replay per base. **Build 33** makes that pass cheap: 70 % of it was a non-finite scan of the
   full ND×ND×NM response arrays and 23 % their zeroing, every frequency, although only the diagonal is used. The
   fast path does 17.8× fewer instructions (74–117 s instead of 43 min at 167 radii, one thread) with byte-identical
   output, and `QS_LIN_FAC 0` passes are byte-identical to build 32 (implementation, adversarial review with its own
   native passes, fix; `results/fast_ali_b33/README.txt`).
4. **Explicitly significant rows** (`--sig-rows`, here H I n = 2) and a **linear update map** (`--update-map linear`:
   x(1 − s) for increases; exp(−s) overshoots a ×35 increase to the cap).
5. **Build 34: the temperature column of the exact action.** The tangent took dη/dT and dχ/dT from the native
   `VAR_OP_V10`, which differentiates the ion free-free terms at a fixed Gaunt factor, while the residual evaluates
   g(ν, T) at the current temperature. A per-process split shows this one term is the whole infrared deficit. The
   tangent now adds dg/dT from the same table cell and clamps as `GFF_VEC`. Tangent/finite-difference ratios of
   dη/dlnT go from 0.40–0.97 (84 radii, front) and 1.31–1.97 (167 radii, depths 80–90) to 1.000000 in every band. The
   complete action's RE-row T column at depths 80–90 goes from 3.5e-2 to 1.19e-5 (oracle consistency 1.24e-5), and
   at 84 radii the rows outside 3× their oracle consistency fall from 111 to 2. Residuals and non-T products are
   byte-identical to build 33. The reviewer reproduced every figure with independent harnesses. **Build 35** adds the
   same term at the diffusion inner boundary (d(DTDR)/dT error 7e-4 → 7e-6; largest product change 5.5e-6 of a
   row's maximum, at depth 167) and stops on free-free resonances, which H/He lacks
   (`results/ir_dT_b34/README.txt`).
6. **Native stops.** `STOP_IF_MW F` keeps CMFGEN's flux-accuracy and inner-optical-depth checks as warnings, so
   that a trial pass is judged by the merit; the base states of the corrected solve sit at the 0.5 edge of the
   flux check. A converged state must pass both checks.

Preconditioners measured against exact products at the pin state, q = ‖W(M⁻¹Jz − z)‖/‖Wz‖ over 56 archived Krylov
pairs (`precond_q.txt`; 0 is perfect, 1 is useless):

| Preconditioner | Median q | Range | n = 2 unit probes (median) |
|---|---|---|---|
| Fixed-J blocks (Reports 20–22) | 0.88 | 0.23–1.41 | 2.2 |
| Fixed-J + exact n = 2 column | 0.88 | 0.23–1.41 | 0.02 |
| Fixed-J + local-ALI increment, all columns | 0.40 | 0.00–4.09 | 0.05 |
| The same, population columns only / H columns only | 0.86 / 0.93 | up to 13.9 / 6.1 | 0.05 |

The local-ALI increment is the best measured but not good: nine of the 56 pairs are amplified 2–4×, all through the
same nearly singular blocks at the two ionization fronts (depths 23–24 and 42–43, where the second singular value of
the block falls from 1.5e-3–5.7e-3 to 4e-5–8e-5). Report 20 rejected it for that reason on a random direction. A
partial increment is worse than none.

**Attempts on the fixed-T problem** (T pinned at every depth, populations and n_e solved; each kept):

| Run | Configuration | Outcome |
|---|---|---|
| `results/nk167_trap/fixedT`, first start | Exact n = 1, n = 2 columns; "coupling significance" (a level is significant if a unit error changes a significant equation by > 1 % of its gross scale) | Stopped before the first product. The rule promoted 5 100 levels, among them He II n = 2 and He III at fraction 1e-33 with M⁻¹R ≈ 1e6: an ionization equation's row belongs to the neighbouring ion, not to the variable of the same index. Rule withdrawn |
| same, restart | Exact n = 1, n = 2 columns; n = 2 significant; probe pairs in the least-squares span | Linear residual 4e-4, but all three trial states refused by CMFGEN's inner-luminosity guard. The stride-6 probe directions move one level by a common factor on every sixth depth of the whole model; the fit used them, and the trial luminosity has a period-6 sawtooth in depth (5.8e9–8.4e9 L☉ between neighbours). Probe pairs removed from the span. Second defect: with fixed-J columns left for n ≥ 3, their untrapped Lyβ, Lyγ… entries (3.9e7 at depth 59) sit in a ground-state row whose exact diagonal is 3.1e4 and produce corrections of order 100 in the ionized interior |
| `results/nk167_trap/fixedT_ali` | Local-ALI increment for every column, n = 2 significant, linear update capped at ×3, start from the transferred 84-radius fixed-T state | Natural measure 140 → 97 → 46 and Newton step 59.5 → 14.2 in two iterates. Then every trial was refused by CMFGEN's inner-luminosity stop: L(r)/L(ND) fell to 0.04–0.34. At fixed T the relaxed populations cannot carry the luminosity, so the populations cannot relax further without T |
| `results/nk167_trap/free_ali`, `free_b34` | T released from that state; the same preconditioner, n = 2 significant, two-photon `USE_RAD` | On build 33: depth-oscillating T steps from the wrong T column; most trials failed natively and the measure at successive bases rose 26.7 → 33.4 before a local-ALI pass stopped the run. On build 34 (same start): every trial completes and the fixed-measure residual falls 224 → 131 → 87 → 77, then stalls at 81. From the first T step the step cools the deep interior (depths 90–167, 15–21 kK) by 5–10 % per iterate with a sign-changing profile; a T dip forms at depth 119 and the radiative flux there turns inward. Cause: residual pieces below 1e-3 of the total, the interior T pieces among them, were merged into one remainder direction that the least squares scaled by ~110 |
| `results/nk167_trap/free_b35`, `fixedT_b35` | Build 35, every depth-region × species piece a separate direction (`--split-tol 0`); then the staged route: statistical equilibrium at fixed T first | T released: linear residual 0.080, no deep-interior cooling, full steps accepted (natural measure 28.9 → 13.6 → 6.5), then the T request at depths 100–130 grows (|s_T| 0.75 → 1.27); clipping T, a T-only Levenberg–Marquardt penalty (population components 3× larger) and a 10× PTC weight (no effect) were tried; whole-step scaling gives 5 % per iterate. Fixed T: fixed-measure residual 224 → 153 → 130 → 118, then 119, 120: helium stalls at the He front and every full or half step fails natively (invalid formal-solution J); stopped, resumable |

### 2.9 What closes the blocker

1. A preconditioner and measure that carry the trapped response of every column and stay bounded at the ionization
   fronts. The local-ALI increment is the candidate; its front blocks need regularizing or the exact neighbour
   coupling.
2. With it, the fixed-T problem (do the trapped levels equilibrate and the ionization relax?), then the T-released
   closure from the cool start, with the electron balance (`COMP_EHB`) as an independent check on every state.
3. A convergence test stated in that measure. Convergence in the fixed-J measure must not be reported as convergence.

All QStarS cool passes since Report 13 used `OLD_DEFAULT`. The Jacobian audits of Reports 20–21 remain valid: they
compare the action with the residual actually computed. The physical conclusions drawn from front temperatures, and
every statement that a cool state was converged, do not.

**Open: the tangent's infrared temperature derivative** (found by the build-32 tangent check; older than build 32, not
fixed). In the T direction at the 84-radius x0 state the exact action reproduces only 92 % of the nonlocal response of
H I n = 4 and n = 5 (rows 5–6: e = 1.1e-2 and 1.3e-2 at depths 14–22, oracle consistency 5e-8). The n = 1 and n = 2
directions agree on these rows to 4e-6 and 1e-7, and rows 2, 3 and NT are not affected. Per-frequency exports
isolate the cause (`results/two_photon_test/b32`, items d and e of its README):
- The temperature derivative of the continuum emissivity that the tangent takes from native `VETA` is too small in the
  infrared: tangent / frozen-J finite difference is 0.93 and 0.91 at 4 µm (depths 17–18), 0.74 below 0.03e15 Hz, and
  1.0000 above 0.36e15 Hz. Candidate, untested: `COMP_FREE_FREE_V2` differentiates free–free at fixed Gaunt factor.
- It reaches the residual where Brα is inverted (depths 17–18) and CMFGEN floors the opacity to
  0.1·η·χ_abs/η_cont. The tangent's floor derivative is exact and inherits the factor.
- With `OLD_DEFAULT` the same direction gives 2e-4–3e-4 on these rows, because the response is 30–60 times smaller.
  Once `USE_RAD` is honoured, two-photon absorption (1s + infrared + ultraviolet photon → 2s) is the infrared
  absorption opacity at the temperature minimum: 32–66 times the `OLD_DEFAULT` value at 4 µm, 130–400 times at 10 µm.
  The native multi-iteration restart model, which ran `USE_RAD`, agrees with the build-32 pass within 9 %. Every
  `OLD_DEFAULT` single-iteration pass therefore also had far less infrared absorption there than the model it
  restarted from.

## 3. Gates 1–3: closing chain made ready (cannot pass before the closure)

A read-only survey of the six gates, then an implement → adversarial review → fix workflow, produced:

- **`scripts/post_convergence.sh`** fixed and hardened:
  - It exited under `set -e` when no cleanup sweep was needed (reproduced).
  - It hard-coded build 28. It now defaults to build 31 with an override.
  - Destinations are now per state, so reruns no longer collide.
  - It refuses to verify unless the cleanup history shows convergence.
- **`scripts/trace_cleanup.py`**: a view-mutation bug made row selection order-dependent (8 844 against 8 905 rows at
  iterate 24). A new opt-in `--pairing structural` follows CMFGEN's own equation assembly, so every trace level uses
  its own equilibrium row and a conservation row is never used to solve for a trace variable. At iterate 24 all
  12 393 trace variables pair; backward error 2.3e-13. Diagonal-only solutions differ by up to 70 % there, so the
  earlier "negligible advection" assumption held only for the rows the old rule kept.
- **`scripts/native_correction_check.py`** (new): the audit's "last native correction" per class, computed from a base
  pass without a 20-minute verify run. Validated on iterate 23.
- **Pre-registered closing rule** (`benchmarks/cmfgen_cool_validation_plan.json`, `closing_phase`). The audit's last
  correction is CMFGEN's native one, which amplifies the exact Newton T step 16× (median) to 8e5×. So passing the
  Newton-Krylov test at 1e-3 does not pass the audit. Rule: continue pure Newton until the predicted post-cleanup
  correction is ≤ 1e-3 on every class the audit counts.
  - The first declaration of the floating-point floor was wrong by orders of magnitude; the reviewer caught it. The
    floor is now *measured* from repeated identical-state passes: at iterate 14, T 1.5e-12–9.5e-12, significant
    ≤ 1.9e-10, all variables ≤ 1.7e-4.
  - Floor-limited levels are reported as FAIL-at-floor, never exempted; the exemption question is marked
    `pending_user_approval`.
- **333 radii**: the 333 restart carried five stale bootstrap trace files, so every 333 pass would have read the wrong
  record, and the "ready" base pass of Report 21 was defective. Fixed: traces moved aside, copy-ignore added,
  single-record guards added in the driver and tools, the stale base moved to `superseded/`. `split_regions` is now
  nested for ND = 333. The base pass still has to be rebuilt (15–20 min).
- **Frequency**: `benchmarks/cmfgen_cool_frequency_plan.json` (all nine recipe keys, 85 359 frequencies; force, heating
  and luminosity criteria with explicit normalizations, two of them `pending_user_approval`).
  `scripts/compare_resolution.py --profiles` is new; it reproduces the archived hot frequency pair exactly and refuses
  to report a pass unless both states passed their own audits.
- **ALI visibility**: the driver now logs unconverged ALI frequencies per pass (one accepted base had one at 1.56e-9
  against 1e-10).

**Consequence of §2 for Gate 1.** Its tooling is ready, but the n = 2 finding means a converged state must be judged by
the all-level audit *and* by the native electron balance. `COMP_EHB` should be added to the verify pass.

## 4. Gate 4: self-selected CMFGEN wind (module and tests; no node yet)

`src/qstars/cmfgen_wind.py`, `tests/test_cmfgen_wind.py` (21 tests), `benchmarks/cmfgen_selected_wind_plan.json`
(sha256 `6c63a583…`, registered before any result), `docs/cmfgen_selected_wind.md`.

- **Equations.** Continuity and momentum with g_rad from the non-LTE field, isothermal critical-point coefficients, the
  regular sonic slope, and a hydro step solving for ṁ in two modes (optical depth, the primary; radius, the
  cross-check). Everything is carried in logarithms; tables are never extrapolated or clipped.
- **Classification**, fixed in advance: `stationary_transonic_wind`, `hydrostatic_no_wind` with an explicit bound,
  `photon_tiring_exit`, `mach_one_mlt_exit`, `reduction_invalid`, `numerical_failure` (never relabelled as physics).
- **Review** found four major defects, all fixed:
  - an ill-conditioned parameterization of the primary mode (the unknowns are now (ln ṁ, ln(τ_c − τ_o)));
  - an uncaught integration stall;
  - a fail-open convergence check;
  - screens read on the wrong grid.
- **Input-chain audit.** On eight archived states the spectral g_rad agrees with CMFGEN's own HYDRO value once one
  definitional difference is corrected (`grad_input_audit.json`).
- **Open.** Every real input must come from a converged released closure. Each outer hydro iteration then needs a new
  closure, and resolution checks must be repeated on the self-selected structure. On the iterate snapshots
  Γ ≈ 0.1–0.6 in the dense subsonic layers, so `hydrostatic_no_wind` is a live possible outcome.

## 5. Gate 5: node-coupled conservative evolution (interface and tests; no node yet)

`src/qstars/wind_node.py`, `selected_nodes` wind mode, `tests/test_node_wind.py` (11 tests), `docs/wind_node.md`.

- **Node contract.** Base face, ṁ, v∞, emergent L, photosphere, energy gauge, provenance hashes, validation flags;
  cross-field identities checked to 1e-10.
- **Coupling.** Single-node drift guard; multi-node log-linear interpolation with no extrapolation (leaving the hull is
  a `model_domain_exit`). Per-step ledger residuals stay < 1e-12; restart is bitwise.
- **Review fix.** A node could declare itself validated. Evolution now always reports `validated: false` until an
  independent Gate 4 evidence verifier exists.
- **Latent defect fixed.** Control mappings for `wind.terminal_speed_cm_s` and `boundary.pressure_cgs` raised an
  uncaught `TypeError` at runtime.
- **Open.** A validated Gate 4 node; an admissible QStarS envelope at the CMFGEN base data; the energy-gauge bridge
  between MESA's EOS zero and CMFGEN's.

## 6. Gate 6: the Ball et al. 2011 evolutionary benchmark

**Why Ball 2011.** It is the strongest published quasi-star *evolution* with a printed table (Table 1; thesis
Table 3.1). It uses no wind and an LTE grey photosphere, so it does not depend on the CMFGEN closure.

**Pre-registration.** `benchmarks/ball2011_track_plan.json`, sha256 `a10b4166…`, recorded in
`ball2011_track_plan.registration.json` before any Ball run. The registration is by sha256 because the user has not
asked for commits; the driver refuses to run if the plan or either target CSV differs. Registered: the model, a
declared replication mode, 44 comparison cells with per-column tolerances derived from stated allowances, the endpoint
classification protocol, conservation and bitwise-restart requirements, an eight-run convergence campaign, and verdict
and stop rules.

- **Accretion law.** The thesis prints 16π in its eq. 3.15, the paper 4π in its eq. 15. The printed table contradicts
  the 4π form by a factor of 3.94–4.03, so the 16π form is registered (unit-tested).
- **Tolerances** (|ln ratio| at matched M_BH): 0.05 for t, Ṁ, L, T_eff, T₀; 0.08 for r_B; 0.10 for R*; 0.15 for ρ₀ and
  M_cav; plus the paper-versus-thesis spread and printed rounding. A strict 1 % tier is reported per cell.
- **Replication mode** `published_assumption`. QStarS's own validity screens (convective Mach > 1, Knudsen and
  free-streaming at τ < 10, opacity-table clamp) are recorded, not enforced. Ball's surface is super-Eddington at
  4.5–5 kK, where they will fire. A PASS in this mode means "published-assumption replica", not physics within the
  QStarS domain.

**Implementation**, seven stages, each implemented → adversarially reviewed → fixed:

| Stage | Content | Notable review findings (fixed) |
|---|---|---|
| 0 | Plan, registration, thesis table CSV, U1 table-consistency tests | Record-keeping only |
| 1 | MESA bridge ABI 7 with the STARS opacity switch and clamp flags; domain sweep (3 628 states inside the domain) | Minor |
| 2 | ML1 convection (`src/qstars/mlt.py`), matching MESA's own routine on 1 000 states | Velocity not checked on all states |
| 3 | Cavity mass in gravity and ledger, remap/split/merge, coloured Jacobian | Split cells lost energy; cut-cell remainders survived |
| 4 | Cavity, photosphere, quasi-static step, screen registry | Controller metric pins dt; re-split jumps at the base |
| 5 | Driver, relaxation, fold classification, `scripts/validate_ball_track.py` | 11 issues: mislabelled outcomes and unenforced registered checks |
| 6 | Config, checkpoint schema 7, CLI | None above minor |

Decisions taken within the plan's amendment rule (numerical method only, never tolerances or physics):
- **Merges conserve total energy exactly**, as registered (≤ 1e-12 per event). "Merge entropy ≥ 0" cannot hold at the
  same time, because a merge moves discrete gravitational energy into heat. The entropy change is reported per event.
- **Endpoint momentum form** instead of the spec's secant form, which left every accepted state 2·Δr/r off static
  balance (a test asserts the difference).
- **Relaxation mesh cycles**, without which two seeds cannot agree to the registered 1e-6.

Default behaviour is byte-identical: 128 files from seven example runs and 164 existing checkpoints were compared
against the pre-change source.

**Stage 7, the registered execution** (2026-10-14; command log `Reports/evidence/report23/ball/stage7_commands.log`, run
directories under `benchmarks/results/2026-10-14-ball-track/`). Each amendment was registered by sha256 only after a
preserved registered run showed its defect, in the order the final runbook prescribed:

| Step | Command (amendments in force) | Outcome | Consequence |
|---|---|---|---|
| 0 | `relax` C1, C0, C2 (none) | `relaxation_blocker` at N = 200, 400, 800: seeds 2000 and 1000 R☉ have Γ_ph ≥ 1, 12 000 R☉ is below the 3000 K MESA floor, 6000 R☉ never takes a step | Amendment 3 (grey seed surface; seed radii by a rule using only registered radii) |
| 1 | `relax` C1, C0, C2 (3) | `relaxation_blocker`: the 3464 R☉ seed fails on a false "Γ_ph ≥ 1" of the registered photosphere bracket | Amendment 4 (bracket scan) |
| 2 | `relax` C1, C0, C2 (3, 4) | `not_unique`: seed spreads 1.5e-2, 3.1e-3, 6.5e-4 | Amendment 2 (relaxation mesh cycles) |
| 3 | `relax` C1 (3, 4, 2) | `converged_unique`, spread 3.9e-8 | — |
| 4 | `smoke`, and `smoke --fallback F1` (3, 4, 2) | Both: 0 accepted steps; four rejections on the registered Lagrangian controller metric, then dt < t_dyn | Amendment 1 (fixed-mass base region and mapped metric; F1 recorded as tried, it cannot act on a controller rejection) |
| 5 | `relax` C1, C0, C2, C5, C6, C7 (3, 4, 2, 1) | All `converged_unique` (C1: 5.5e-11), no enforced-screen error at row 0 | — |
| 6 | `smoke` (3, 4, 2, 1) | **Failed.** 1 255 accepted steps (median dt 139 yr, 0.6 % rejected) to M_BH = 46.0 M☉, then eight consecutive rejections whose controller change (1.5–1.9 % in T) does not shrink with dt, and dt < t_dyn | **Numerical blocker** (plan stop rule; F1, F2 and F3 used up) |

**Why the smoke gate stopped** (`smoke_jump/diagnosis_and_verification.md`; investigator plus an adversarial
verifier, both reproducing from the preserved checkpoint):
- The registered split/merge rule merges two cells 0.05 M☉ below the surface, inside the outer hydrogen opacity peak.
  There the local Eddington factor reaches 3.4, density rises outward ninefold, and κ changes by a factor 1.8–2.4
  between neighbouring cells. The merge changes the discrete solution of the next cell inward by +1.9 % at every
  dt, which is the change of the N = 200 discretization error there (refinement converges at second order to 2.6 %
  below the 174-zone value). The same step without the merge is smooth and proportional to dt. `fold.classify` finds
  no fold (σ_min flat, det J of constant sign over 44.8–49.2 M☉). Class: numerical method.
- The proposed Amendment 5 (measure an event step on the pre-event mesh) was **rejected** and not registered. It
  removes mesh-event jumps from error control with no bound: its dry run accepted jumps up to 11.7 × the registered
  δ_target, which amounts to a tolerance relaxation the amendment rule forbids
  (`benchmarks/ball2011_track_plan_amendment_5.rejected_draft.json`).
- No merge fix reaches M_BH = 100. A diagnostic dry run with merges switched off and the registered controller stops
  at 55.0 M☉: Newton iterations grow from 4 to 40 between dt = 0.1 and 2 yr in the same front (Γ up to 4.4), and 14 of
  its 20 rejections are Newton stalls just above the registered 1e-12 tolerance. That stop is unclassified.

**Known in advance, from the dry runs, and not acted on:** the registered verdict would be FAIL whatever the tracks
did. Required check Q2 ("density inversion absent at row 0") fails at row 0 at every resolution (largest outward
density increase 9.1e-3, 5.4e-3, 2.8e-3 at N = 200, 400, 800 against 1e-3; a resolved inversion at ≈ 28 kK). The
row-0 Richardson error of R* is 0.0099 against the 0.01 limit.

**Gate 6 outcome:** FAIL. The registered smoke gate is a recorded numerical blocker at M_BH = 46 M☉ in the
super-Eddington outer hydrogen opacity peak; the convergence campaign C0–C7, the restarts and the 44-cell comparison
were not run; and Q2 fails at row 0 by the registered rule. No published evolutionary benchmark is reproduced. What
is established: a moving-cavity quasi-static quasi-star evolves conservatively on MESA microphysics for 1 255 steps
(41 M☉ of accretion): per-step energy residual changes ≤ 1.2e-14 of the initial internal energy and mass ≤ 3.6e-16
(registered budgets 5e-13 and 1e-14), cumulative energy residual ≤ 6.3e-14, from a unique relaxed row 0 at three
resolutions.


## 7. How the work was checked

Seven multi-agent workflows ran. Each implementation was followed by an adversarial review that re-ran the tests and
re-derived the mathematics; every confirmed finding was fixed before the next stage.

| Workflow | Agents | Outcome |
|---|---|---|
| Chain readiness survey | 7 | Per-gate status; five real defects; "v1 cannot be declared on the ledger's own terms" |
| Gates 1–5 development | 9 | Tooling, wind module, node coupling; one critical and six major findings fixed |
| Ball design (judge panel) | 4 | Three independent designs scored, one spec synthesized |
| Ball stages 0–3 | 12 | Registered plan, ABI 7, ML1, core extensions |
| Ball stages 4–6 | 10 | Step, driver, integration |
| RE row against thermal balance | 4 | The diagnosis of §2.5, re-verified by the synthesizer |
| Build 32 | 3 | Accepted by the reviewer with no critical or major finding; two driver guards added; a pre-existing infrared temperature-derivative defect of the tangent isolated (§2.9) |
| Ball Stage 7 unblock | 4 | Registered seed shown unbuildable on MESA; grey-seed and photosphere-bracket amendments with failing-before tests; uniqueness and row-0 screens enforced; records always written; F1 implemented (rescues nothing) |
| Build 33 (fast local-ALI pass) | 3 | 17.8× fewer instructions, byte-identical output; accepted with minor caveats (path-dependent executable hash, wording) |
| Ball smoke-gate stop | 2 | Mesh-merge jump in the outer hydrogen opacity peak diagnosed; the proposed controller amendment rejected by the verifier |

Test suite: **334 passed, 37 skipped** (Report 22: 111 passed, 23 skipped). The skips need compiled MESA/FreeEOS libraries or the native compiler environment.

## 8. Gate table

| Gate | Status | Evidence |
|---|---|---|
| 167-radius T-released closure converged | **FAIL**: not converged. Causes found and fixed (trapped-transition measure, `OLD_DEFAULT` two-photon, Gaunt-factor T derivative); the corrected solves halve the fixed-measure residual and then stall at the ionization fronts (local-ALI preconditioner nearly singular) | §2 |
| Numerical or physical? | **Numerical** (trapped transitions invisible to the solver; `OLD_DEFAULT` two-photon artifact; wrong T column of the action); the front equilibrium is not yet demonstrated | §2.5–2.8 |
| 1. Post-convergence cleanup and audits | **NOT EVALUATED** (no converged state); tooling fixed, closing rule registered | §3 |
| 2. 333 radii | **NOT REACHED**; stale-trace defect fixed, base pass to rebuild | §3 |
| 3. Doubled frequency and resolution comparison | **NOT REACHED**; plan registered, comparison tool extended | §3 |
| 4. Self-selected CMFGEN wind | **NOT REACHED**; module, 21 tests, registered plan | §4 |
| 5. Node-coupled conservative evolution | **NOT REACHED**; interface, 11 tests | §5 |
| 6. Published evolutionary benchmark | **FAIL**: registered smoke gate a numerical blocker at M_BH = 46 M☉ (outer hydrogen opacity peak, Γ up to 3.4); campaign not run; Q2 fails at row 0 by the registered rule | §6 |
| Other ledger criteria (multiple replicas, cross-paper prediction, H⁻ sensitivity, atom size, moving-cavity conservation) | **OPEN**, outside this chain | `completion_gates.json` |
| v1 | **NO** | — |

## 9. Withdrawn statements, failed and superseded attempts

**Withdrawn.**
- *Report 21 §4.3–4.4 and Report 22 §5*: "on the cool branch the front's net heating is ∝ T", "thermally unstable
  stretch", "PTC moves the front in the physical (heating) direction". The quantity was the RE row at frozen
  populations. It is net radiative absorption into an unbalanced n = 2, three quarters of it the two-photon artifact.
  The electron balance was cooling throughout.
- *Reports 20–22*: levels below 1e-10 of their element are "numerically irrelevant to the other equations". False for
  H I n = 2.
- *This report's first draft*: the 10.2 kK local root as evidence of a finite-temperature equilibrium, and "a strong
  coolant is available, so the equilibrium is near 10 kK". The root is artifact-set. The coolant is real, and it puts
  the equilibrium *cooler*.
- *Report 22 §4.2*: the offline gain from ultra weight 1e-10 was real for the least squares it measured. The larger
  front T steps it allowed went in the wrong direction.
- *Reports 20–22*: every statement that a cool state was converged, or "converged to a few 1e-3" (the fixed-T
  84-radius closure of Report 20, the full-residual states of Report 21). They were converged in the fixed-J measure,
  which is blind to every trapped transition (§2.7).
- *This report*: the "coupling significance" rule (§2.8). It promotes ultra-trace levels through ionization rows.

**Failed or superseded, all kept on disk.**
- PTC iterates 15–34 of `results/nk167_hot` (direction set by the corrupted row); paused.
- Existence test attempt 1 (free T outside the window): did not converge.
- Existence test θ = 1.0 and θ = 1.15: readout invalid (§2.3); θ = 1.15 also did not contract, and its pinned T
  overshot because steps are applied as x·exp(−s).
- Ultra weights 1e-10 and 1e-12: each stopped decoupling the ultra-trace rows as their corrections grew.
- The first pre-registered floating-point floor (wrong units and model); replaced by a measured floor before use.
- `results/nk167_trap/fixedT` (two configurations, §2.8) and `fixedT_ali` (stopped after two iterates when the fixed-T
  populations could no longer carry the luminosity).
- `results/nk167_trap/free_ali`: a Levenberg–Marquardt trust region on T (shrank every component; no predicted
  progress), and two PTC restarts whose trials were all ended by CMFGEN's optional flux-accuracy stop before
  `STOP_IF_MW F`. Kept in `superseded/`.
- Ball Amendment 5 (rejected, §6) and the merge-free and dt-continuation dry runs that diagnosed the smoke-gate stop.

**Process slips, each logged in the run READMEs.**
- A restart without `--max-outer` ended the run at the default 30 iterates (40 min lost).
- A switch armed from zsh passed its option list as one word, and the restart failed at argument parsing (12 min).
- Two waiting commands matched the `pgrep` pattern of the queue they watched.
- A driver edit (front-pin weights) broke broadcasting over product directions; caught at the first product pass,
  fixed, and the pass reused.
- The probe directions of the exact-column preconditioner were first added to the least-squares span (period-6
  sawtooth in the trial luminosity, §2.8).
- A Stage-7 command line again passed zsh array words as one argument (refused at argument parsing, no run made), and
  a bash script used an associative array, which macOS /bin/bash 3.2 lacks (stopped during the C1 relax, which
  completed on its own).
- A PTC restart at the same base did not reuse the exact products of that base (the reuse looks only in other
  iterate directories); the pass it started was stopped and the products were moved so that they are reused.
- The laptop was asleep for 154 of 167 h before this session's final day.

## 10. Reproduction

Registered plans (sha256 prefix): closure stop rule `a2711719` + amendment `cfe15815`; existence test `c55b1f13` +
amendment `93bb1c8e`; Ball track `a10b4166`; selected wind `6c63a583`; cool frequency `7326bbce`.

```bash
# RE row against the electron balance, and the two predictive tests (§2.4-2.6)
.venv/bin/python scripts/verify_cool_state.py STATE OUT --x0 results/cmfgen_ali_lambda167 \
  --exe backends/cmfgen/build_31/cur_cmf/exe/cmfgen_dev.exe --replace results/nk167_adv/transfer/replace_167.txt \
  --free-t --vadat INC_AD T --vadat INC_ADV T --vadat COMP_EHB T [--vadat TWO_METH NOSTIM]
.venv/bin/python scripts/replay_with_vadat.py BASE OUT --x0 results/cmfgen_ali_lambda167 \
  --exe backends/cmfgen/build_31/cur_cmf/exe/cmfgen_dev.exe --vadat TWO_METH NOSTIM
.venv/bin/python scripts/front_progress.py results/nk167_hot
.venv/bin/python scripts/thermal_curve.py BASE DEPTH OUT --x0 ... --exe ... --replace ...

# closing chain, once a state converges
EXE=backends/cmfgen/build_31/cur_cmf/exe/cmfgen_dev.exe scripts/post_convergence.sh RUN BASE X0 REPLACE
.venv/bin/python scripts/native_correction_check.py BASE [--noise-pair P1 P2]

# Ball 2011 (refuses to run if the plan or target sha256 differ from the registration)
.venv/bin/python scripts/validate_ball_track.py --plan benchmarks/ball2011_track_plan.json relax|smoke|track|campaign|compare
```

```bash
# trapped-level solver (builds 32-33); launch scripts record every restart in the run README
results/nk167_trap/launch_fixedT_ali.sh            # fixed T, local-ALI preconditioner
results/nk167_trap/launch_free_ali.sh BASE K --t-homotopy --hom-mode ptc --ptc-delta0 0 --ptc-target-dT 0.05 \
  --t-log-cap 0.1 --component-log-cap 0.693 --vadat STOP_IF_MW F
.venv/bin/python scripts/monitor_closure.py results/nk167_trap/free_ali
.venv/bin/python scripts/replay_with_vadat.py BASE OUT --x0 results/cmfgen_ali_lambda167 \
  --exe backends/cmfgen/build_33/cur_cmf/exe/cmfgen_dev.exe --lin-fac 1 --two-phot-j base

# Ball Stage 7: register an amendment after its preserved failing run, then one registered command per call
.venv/bin/python Reports/evidence/report23/ball/register_amendment.py N RUN_DIR [RUN_DIR ...]
Reports/evidence/report23/ball/stage7_run.sh --amendment benchmarks/ball2011_track_plan_amendment_3.json ... relax --run C1
```

Driver options added in this report: `--hom-mode frontpin --pin-depths D0 D1 --pin-factor θ [--pin-all]`,
`--weight-ramp`, `--ultra-weight 0`, `--two-phot-j`, `--ali-precond`, `--probe-rows/--probe-stride`, `--sig-rows`,
`--update-map`, `--sig-coupling` (withdrawn); `scripts/switch_at_base.sh` restarts a run at its next base (call it
from bash). Evidence: `Reports/evidence/report23/`. Nothing is committed.
