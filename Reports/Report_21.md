# Report 21 — closing the cool-closure gates

**Status: QStarS is NOT v1.** Two requested gates are closed. The released cool closure has not
converged and remains the single narrowest blocker. Every later gate depends on it.

## 0. Outcome

**Closed this run**
1. **Level-dissolution derivatives in the Jacobian/RE response** (builds 25–27). The complete
   action, RE row included, now agrees with the exact T-released residual at the oracle's own
   consistency in every direction: RE row 5.9e-6 (T), 6.8e-5 (n_e), 2.3e-5 (populations), against
   0.25 in n_e before. Spectral opacity-change mismatches fell from 841 314 cells to 160 (§2).
2. **`INC_AD` and `INC_ADV` restored and validated.** They are on in every pass. The
   block-tridiagonal direct part is validated by an odd-depth test where only the depth coupling
   responds: 7.8e-6, against 1.0 without the off-diagonal blocks (§3).

Also delivered:
- Build 28 removes a nested-OpenMP defect (deadlock and intermittent SIGSEGV). Products are
  unchanged to summation order (1.2e-16), with less threading overhead: the same pass took 360 s
  instead of 810 s under similar load, and system time fell from 981 s to 361 s.
- A 333-radius nested grid with its base pass, for the spatial check.
- Pseudo-transient continuation on T and split Krylov directions in the driver.
- An exact measurement of the front's coupled thermal response.

**Not closed (narrowest blocker): the T-released cool closure.** T continuation (residual shift,
λ: 0 → 1) carries the corrected RE row to λ = 1 at 84 and 167 radii without a retreat. At λ = 1,
free Newton steps then drive the H I front layer down the cool branch, toward T → 0, on both grids
(§4.1–4.2). Refinement does not remove this. The exact coupled response explains it. On the cool
branch, at frozen populations, the front's net heating is proportional to T: coupled slope ≈ R_RE,
5–7e-3 of the RE row's gross Lyman rate (§4.4). Newton's only target there is T = 0, a root that
does not exist. Pseudo-transient continuation moves the front in the physical (heating)
direction, from 5.5 kK and from 7.9 kK and 10.5 kK hot starts, with accurate split-direction
solves. The iterates are still far from statistical equilibrium, however: CMFGEN's thermal and
radiative energy balances disagree there by orders of magnitude. The continuation moves the
front only 3–5 % in T per hour-long iterate (§4.5). **No converged released state exists yet**, so the
electron/energy audits, the spatial (167/333) and frequency resolution checks, the self-selected
wind, coupled evolution and the evolutionary benchmark were not reached. **v1 is not declared.**

**What closes the blocker.** The 167-radius PTC runs `results/nk167_hot` (7.9 kK start) and
`results/nk167_hot2` (10.5 kK start) were left running at the end of this session. Their natural
measure falls 2–15× per iterate, but the front populations are not yet in statistical
equilibrium. If either converges, the chain is ready to run:
1. `scripts/post_convergence.sh` (trace cleanup, verification, audits);
2. transfer to 333 radii (`scripts/transfer_state_nested.py`) and the same solve there;
3. the doubled-frequency solve.

If instead the coupled front proves to have no finite-temperature root, the question becomes
physical rather than numerical. In the H/He-only cool model, the dense neutral layer is cooled
only by H I collisions, H⁻ and free-free (GENCOOL, §4.5): there are no metal coolants. Adding
physics would be a separate, explicit model decision, not a convergence aid.

## 1. Starting point (Report 20)

Report 20 delivered an exact Jacobian action for CMFGEN's cool residual. Its RE row was conditional:
it missed level-dissolution derivatives and had 16–25 % net errors in the thick interior. The 84-radius
T-released solve did not converge from either branch, with the d11/d12 front unresolved. `INC_AD` and
`INC_ADV` were still off. No audits, resolution pair, wind or evolution existed.

## 2. Level-dissolution derivatives completed (builds 25–27)

| Build | sha256 (executable) | Adds |
|---|---|---|
| 25 | `c9270eac6886…1d08` | d(occupation probability W) in the LTE populations of every bound-free term and d(dissolved fraction T2) below each edge, chained through `B/A/X_LEV_DIS(n_e, POPION, T)` as written; evaluated at `CONT_FREQ` and cached like the native continuum |
| 26 | `e54dd83b52a3…62e9` | dissolution terms moved from `CONT_FREQ` to the actual frequency exactly as `COMP_VAR_OPAC` moves `VCHI/VETA`; optional `VCHI_ALL` (shown identical: every variable is "important" in this model) |
| 27 | `558395f7336f…0bb8` | dissolution derivative of the line opacity: `L_STAR_RATIO = W(up)/W(low)` (`SET_LINE_OPAC`, one-to-one super levels) |
| 28 | `107289a7ab29…ed67` | threading only (§9): no nested OpenMP team in the tangent. Products equal build 27 to 1.2e-16 (`build28_vs_build27_product.json`); the same 3-direction pass took 360 s instead of 810 s under similar load (system time 361 s vs 981 s) |

The occupation probability W evaluated from the formula equals the native W to 2.2e-16 at every
level and depth. `QSTARS_JV_NODIS` reproduces build 23's products to 8e-17.

**Spectral layer** (`results/jv/dis_tests`; n_e direction, replay passes at x0 with n_e perturbed
±1e-3 exporting CMFGEN's own opacity/emissivity, against the tangent's exported changes;
evidence `spectral_ne_direction.json`):

| Build | opacity-change cells off by > 1e-3 (of 3 638 376) | emissivity change, rel. L2 |
|---|---|---|
| 23/25 without the terms | 841 314 | 7.5e-2 |
| 25 (bound-free W and T2) | 89 586 | 1.5e-4 |
| 26 (moved to FL) | 89 586 | 6.7e-5 |
| 27 (line W ratio) | **160** | 6.7e-5 |

The remaining 160 cells are not discontinuities: none shows a jump or sign change across the ±1e-3
perturbation. They lie between 1.7e13 and 8e14 Hz and their relative error reaches 13 %. Their
origin is not identified. Their weight in the residual is below the oracle floor (next table).

**Complete action vs the T-released oracle of Report 20** (`results/jv/tfree`,
`scripts/audit_action_oracle.py`; RE row = radiative-equilibrium row; evidence
`action_audit_b23.json`, `action_audit_b27.json`):

| Direction | all rows (ε_a / ε_b) | RE row, build 23 | RE row, build 27 (ε_a / ε_b) | oracle consistency (all / RE) |
|---|---|---|---|---|
| T | 4.8e-5 / 4.8e-6 | 1.1e-4 | 5.9e-6 / 1.6e-6 | 4.3e-5 / 5.4e-6 |
| n_e | 3.0e-4 / 3.0e-5 | 0.25 | 6.8e-5 / 1.0e-5 | 2.7e-4 / 5.9e-5 |
| random populations | 1.8e-5 / 1.8e-6 | 4.5e-3 | 2.3e-5 / 2.3e-6 | 1.7e-5 / 2.1e-5 |
| T + n_e + populations | 2.4e-5 / 2.4e-6 | 3.7e-3 | 6.1e-6 / 1.8e-6 | 2.1e-5 / 5.6e-6 |

Every row of the complete action, the RE row included, now agrees with the exact residual at the
oracle's own consistency, and the error falls tenfold with the oracle step (the oracle's O(ε)
nonlinearity). Build 25 alone made the RE row *worse* (n_e 0.25 → 0.32), because it added the
bound-free terms without the line ratio; that intermediate result is kept in
`build25_first_audit.json`. The direct (fixed-J) part keeps its frozen-J completion
(Report 20 §2.2), which is exact and already validated.

## 3. Adiabatic cooling and advection restored

`INC_AD=T` (adiabatic cooling in the RE equation: expansion work and advection of internal
energy, `EVAL_ADIABATIC_V3`) and `INC_ADV=T` (advection of every level population in the rate
equations, `STEQ_ADVEC_V4`, relaxation parameter 1) are now on in every Report 21 pass. Both are
linear and appear analytically in the native blocks as off-diagonal depth bands (22 302 nonzero
upper-band entries at x0). `cool_newton_krylov.Iterate` therefore keeps the native lower/upper
blocks and completes only the diagonal blocks. The all-depth frozen-J columns (n_e, T, H II) are
corrected by subtracting the native off-diagonal contribution of the same column before the
n_e/T replacement and the POPION term.

Validation (`results/jv/tfree_adv`, `results/jv/adv_offdiag`; evidence `action_audit_adv.json`,
`adv_offdiag_test.json`):
- complete action vs a T-released exact-residual oracle with `INC_AD`/`INC_ADV` on: identical
  agreement to the table above (advection is small in the global norm: it changes the residual by
  3e-8 relative);
- **sharp test of the off-diagonal blocks**: the T + n_e + populations direction perturbed at odd
  depths only, with frozen J. The response at even depths comes only from the depth coupling. The
  direct action agrees with the finite difference there to **7.8e-6** (1.0 without the off-diagonal
  blocks); RE row at even depths 8e-9.

`INC_AD` and `INC_ADV` are therefore restored in the residual and validated in the Jacobian. Their
physical effect is judged on the converged state (§5).

## 4. T-released solve with T continuation

**Driver** (`scripts/cool_newton_krylov.py`): Newton–Krylov on the full population + n_e + T
residual with the native radiation field. Each Krylov direction's product J·z is the build 27
complete action, which comes from one tangent replay pass (the exact radiation part) plus the
block-tridiagonal direct part of §3. The preconditioner is the factored direct part M. The merit
is the significance-weighted ‖W M⁻¹ R‖, with weights 1 / 1e-3 / 1e-8 for significant, trace and
ultra-trace rows (trace rows are fixed from the start state, not re-tuned). Steps are accepted by
the natural monotonicity test (Deuflhard), with a Levenberg–Marquardt trust region in the Krylov
subspace and parallel line-search trials, each a full base pass. Nothing is clipped. The
convergence criterion is native CMFGEN's own: the maximum relative Newton correction over
significant levels must be below 1e-3. A trial that does not lower the merit retreats λ and is never
accepted.

**T continuation (residual shift)**: `R_T(x; λ) = R_RE(x) − (1 − λ) R_RE(x_s)`, with x_s the start
state (frozen in `homotopy_shift.npz`). At λ = 0 the start state satisfies the T row exactly, and
the populations, n_e and radiation are solved at the start T. Each stage adds Δλ = 0.25
(predictor–corrector: at most 2 corrector iterates per stage, stage tolerance 1e-2). The step
doubles after a fast stage and halves after a stage that fails to decrease the merit. At λ = 1 the
residual is CMFGEN's own, untouched. The Jacobian is the same at every λ (the shift is constant
in x), so the validated action applies without modification.

### 4.1 84 radii (`results/nk_adv84`, residual shift from iterate 5)

Iterates 0–4 used the pinning homotopy (§9). The start of the shift is the converged λ = 0 pin
state `it05_trial_l1` (T = Report 20's best fixed-T profile, populations converged at that T).
Natural measure ‖W s̄‖ is shown at the start of each iterate and after the accepted trial;
"linear" is the relative residual of the augmented minimal-residual solve after 4 products.

| it | λ | natural measure | max significant step | linear | accepted | after |
|---|---|---|---|---|---|---|
| 5 | 0 | 8.47 | 2.50 | 0.81 | trust 0.5 | 7.52 |
| 6 | 0 | 26.1 | 7.93 | 0.65 | trust 0.25 (radius 1 refused by CMFGEN) | 18.5 |
| 7 | 0.25 | 2.87 | 0.90 | 0.91 | trust 0.25 | 2.34 |
| 8 | 0.25 | 0.664 | 0.233 | 0.87 | full | 0.340 |
| 9 | 0.5 | 0.491 | 0.273 | 0.89 | full | 0.152 |
| 10 | 0.5 | 1.93 | 0.575 | 0.99 | trust 0.25 | 0.084 |
| 11 | 0.75 | 2.41 | 0.811 | 0.64 | full | 0.807 |
| 12 | 0.75 | 0.570 | 0.249 | 0.92 | full | 0.209 |
| 13 | **1** | 1.05 | 0.463 | 0.18 | full | 0.163 |
| 14 | 1 | 1.42 | 0.692 | 0.82 | full | 0.390 |
| 15 | 1 | 1.03 | 0.384 | 0.90 | full | 0.123 |
| 16 | 1 | 1.11 | 0.481 | 0.79 | trust 0.25 | 0.659 |
| 17 | 1 | 3.20 | 2.87 | 0.28 | trust 0.25 | |
| 18 | 1 | | 0.546 | | trust 0.5 | |
| 19 | 1 | 0.633 | 0.226 | 0.88 | full / half | 0.334 |
| 20 | 1 | 0.435 | 0.210 | | full | 0.294 |
| 21 | 1 | | 0.119 | | stopped | |

The continuation carried the corrected RE row from the start state to λ = 1 in eight iterates
without a single retreat. On the way, the coupled exact step heated the cool layer
(d16–22: 3.8–5.1 kK → 5.2–6.6 kK) and cooled the outer three depths (12.4 → 9–10 kK). This is
the direction of the frozen-J local balance of Report 20 §4.5, which Report 20's incomplete RE row
had contradicted. At λ = 1 the iteration stagnates. Each accepted step lowers the natural measure
3–8-fold, but the fresh linearization of the next iterate returns it to ~1. The front cells
d18–20 settle into a grid-scale zig-zag (T ≈ 5.4 / 6.5 / 5.1 kK). At iterate 17 the step turns back
toward the cool branch (s_T = +0.41…+0.48 at d17–19, i.e. ×0.62). The dominant significant
component throughout is the He II ground state at d16–20 (1e-7–1e-5 of helium), whose
population follows T as ≈ e^(36 Δln T). **The 84-radius released solve has not converged.** The
evidence (zig-zag at the front, a step that alternates between branches) is again consistent
with a front narrower than the 84-grid spacing. It is not a demonstration that no solution
exists on this grid. From iterate 17 the iteration followed the cool branch. The front cells
d17–19 cooled monotonically (3.9/3.3/3.9 → 2.4/2.0/2.5 → 2.2/1.8/2.1 kK) while the measures fell,
and every new step asked for another 10–15 % cooling there. This is Report 20's cold-branch
runaway (T → 0 at the front), now with the complete RE row and `INC_AD`/`INC_ADV`. It is not
convergence toward a finite-temperature state. The run was stopped at iterate 21 so the
167-radius solve has the machine.

### 4.2 167 radii (`results/nk167_adv`)

Grid `results/cmfgen_ali_lambda167` nests the 84 grid (84-depth k = 167-depth 2k−1). The start
state is Report 20's best fixed-T state, transferred (shared radii copied, midpoints
log-interpolated). Iterates 0–2 used the pinning homotopy (λ = 0, 0.25 pinned). The residual shift
starts at iterate 3 from the λ = 0.25 pin state. Iterate 4 restarted on build 28 after the
OpenMP deadlock (§9).

| it | λ | natural measure | max significant step | linear | accepted | after |
|---|---|---|---|---|---|---|
| 0 | pin 0 | 2.15 | 0.404 | 0.21 | full | 0.093 |
| 1 | pin 0 | 1.82 | 0.922 | 0.23 | full | 0.030 |
| 2 | pin 0.25 | 0.894 | 0.681 | 0.30 | full | 0.024 |
| 3 | 0 | 11.1 | 7.49 | 0.65 | trust 0.5 | 10.2 |
| 4 | 0 | 4.08 | 2.18 | 0.84 | trust 0.25 | 3.22 |
| 5 | 0.25 | 1.73 | 0.757 | 0.81 | trust 0.25 | 0.404 |
| 6 | 0.25 | 0.430 | 0.236 | 1.00 | full | 0.033 |
| 7 | 0.5 | 1.95 | 1.27 | 0.85 | trust 0.5 | 0.033 |
| 8 | 0.5 | 0.957 | 0.704 | 0.95 | full | 0.037 |
| 9 | 0.75 | 2.05 | 1.20 | 0.60 | full | 0.264 |
| 10 | 0.75 | 0.244 | 0.149 | 0.93 | full | 0.0031 |
| 11 | **1** | 1.22 | 0.675 | 0.14 | full | 0.021 |
| 12 | 1 | 1.34 | 0.764 | 0.70 | full | 0.190 |
| 13 | 1 | 0.818 | 0.460 | 0.90 | full | 0.025 |
| 14 | 1 | 0.525 | 0.186 | 0.93 | full | 0.026 |

Each stage at 167 radii ends one to two orders of magnitude lower than the 84-radius stages. The
coupled step is smooth across the cool layer (§4, first shift step: s_T heats d30–42 with no
sign alternation between neighbours). At λ = 0.75 the dip d31–39 is at 5.0–6.1 kK (3.8 kK at the
start), with a shallower minimum at d40 (4.7 kK). The outer layers are unchanged (12.5 kK at d1).

At λ = 1 the 167-radius iteration does what the 84-radius one did, only smoothly, without a zig-zag.
Each step cools the warm plateau of the front: d33–38 went 5.5 → 4.3 → 3.7 kK, and the
iterate-14 step asks for another 14 %. The natural measure falls slowly (1.22, 1.34, 0.82, 0.53 at
the start of iterates 11–14). The run was stopped after iterate 14. **Grid refinement does not
remove the cold-branch slide**, so the blocker is not only under-resolution of the front.

### 4.3 Why free Newton runs away on the cool branch, and pseudo-transient continuation

Report 20 §4.5 found that the frozen-J local balance of the cool layer is *net heating at every
cool temperature, rising slowly with T*. The RE function there has a positive T slope and no cool
root. A Newton step reduces a positive net heating by *lowering* T, toward the zero of its linear
extrapolation. On such a branch it chases a root that does not exist and drives the layer to
T → 0. That is exactly what 84-radius iterates 17–21 did. Physical relaxation goes the other way:
net heating raises T, toward the stable root on the hot side (≈ 7.9 kK at d12 in Report 20). The
residual-shift path did heat the layer at λ < 1, and the free Newton iteration at λ = 1 turned
back to the cool branch.

**Pseudo-transient continuation (PTC) on T** (`--hom-mode ptc`): the RE row is replaced, *per
iterate*, by `G_T = R_RE − δ c (T − T_k)/T_k`, where T_k is the current base temperature and c the
RE row's gross scale. This is one backward-Euler pseudo-time step of dT/dτ ∝ net heating. The
Jacobian gains −δc on the T diagonal. The radiation products and every other row are unchanged.
At any fixed point T = T_k, so the term vanishes and the converged state is a root of CMFGEN's
unmodified residual. δ starts at 0.1 and follows switched evolution relaxation (δ scales with the
RE net/gross imbalance). δ is set to 0 when the regularized significant step falls below 1e-2 or
δ below 1e-4. Convergence is judged on the pure Newton step at δ = 0 with the native 1e-3
criterion. Unit test: `test_ptc_residual_jacobian_and_heating_direction` (finite-difference
Jacobian, the term vanishes at the base, net heating raises T).

Practical settings, each found in the 84-radius test and recorded in `results/nk84_ptc/README.txt`:
- **δ scale.** δ0 = 0.1 pinned T (steps ≤ 1e-3), because the RE net/gross imbalance is only ~1e-4
  away from the front. The automatic δ0 = max|R_RE/c| / 0.05 is dominated by the front (≈ 0.13).
  It still moved T by only 0.2 % with the unsplit Krylov solve.
- **Split directions.** Every unsplit direction is dominated by population components, with T
  ~1–2 % of each. With `--split` (T its own group per depth region, all products exact) the linear
  residual falls to 0.11 (62 directions, 3 passes), against 0.6–0.9 before.
- **Component cap.** The accurate split step asked He I excited levels at d16–18 (1e-16–1e-20 of
  helium) for e^29. With the default cap of e^9.2 every trial raised the measure, so the update
  cap was tightened to ×10 per iterate (`--component-log-cap ln 10`). It binds only on
  non-feeding-back trace components, and convergence is still judged on the unlimited step.
- **SER direction.** Heating the 84-radius front from 5.8 to 5.95 kK *raised* its net heating
  (RE/c 6.6e-3 → 7.0e-3): an unstable stretch, with the stable root hotter. SER was therefore made
  decrease-only (δ held while the imbalance grows, ×4 on a rejected trial), and δ was reset to
  ≈ 0.03 so a pseudo-step can move T by ~20 %, within the 0.25 T trust region.

### 4.4 Coupled local T sensitivity at the front (measurement)

From the exact radiation response of the T columns (`probe_columns.py`, 84 radii, iterate-5
state, populations frozen, J responding), the RE row's coupled slope at the front is positive
and almost equal to the net heating itself:

| depth | T [kK] | direct ∂R_RE/∂ln T / c | radiation / c | coupled / c | R_RE / c |
|---|---|---|---|---|---|
| 16 | 3.98 | +5.3e-4 | +9.5e-4 | +1.5e-3 | +1.4e-3 |
| 17 | 3.80 | +1.9e-3 | +3.5e-3 | +5.4e-3 | +5.1e-3 |
| 18 | 3.80 | +2.4e-3 | +4.4e-3 | +6.8e-3 | +6.5e-3 |
| 19 | 3.97 | +1.9e-3 | +3.9e-3 | +5.8e-3 | +4.7e-3 |
| 20 | 4.27 | +7.4e-4 | +1.1e-3 | +1.8e-3 | +1.7e-3 |

With populations frozen the front's net heating is roughly proportional to T. Its only root is
T = 0, which is where free Newton steps go. A finite root needs the populations to change with T:
hydrogen ionizing near 8–10 kK, where Report 20's locally re-equilibrated frozen-J balance turns
negative (≈ 7.9 kK at d12). The gross scale c of the front's RE row comes entirely from H I
population columns: thick Lyman terms that nearly cancel. With PTC at δ ≈ 0.03–0.1, the
split-direction solves therefore reduce the imbalance through H I population/J adjustments rather
than through T. Over three PTC iterates (84 radii: 2, 167 radii: 1) the natural measure fell 2–4×
per iterate, but the RE net/gross norm stayed at 1.06e-2 (84) and 1.14e-2 (167) with T nearly
unchanged. Crossing to the hot root at that pace would take many hour-long iterates.

### 4.5 Pseudo-transient continuation runs

| Run | start | iterates | natural measure (accepted) | RE net/gross | front T |
|---|---|---|---|---|---|
| `nk84_ptc` | 84 r, λ = 1 state `it14` (front 5.2–6.6 kK) | 2 | 2.23 → 0.66; 1.25 → 0.47 | 1.02e-2 → 1.07e-2 → 1.06e-2 | d17–18 5.80/5.71 → 5.95/5.87 kK |
| `nk167_ptc` | 167 r, `nk167_adv/it12` (plateau 5.4–5.5 kK) | 2 | 2.08 → 0.50; 1.16 → 0.18 | 1.14e-2 → 1.14e-2 → 1.18e-2 | heated 1 %, then 2.5–3.8 % |
| `nk167_hot` | 167 r, T ≥ 7.9 kK at d19–47 (populations re-equilibrated) | 5 accepted, running (iterate 4 starts at 2.77: not monotone) | 4.45 → 4.16; 4.22 → 2.57; 2.36 → 0.16; 1.47 → 0.26; 2.77 → 2.63 (trust 0.031) | 1.63e-2 → 1.52e-2 → 1.53e-2 → 1.56e-2 → 1.57e-2 → 1.51e-2 | d35 7.90 → 7.74 → 8.15 → 8.87 kK (it04) |
| `nk167_hot2` | 167 r, T ≥ 10.5 kK at d30–40 only | 1 accepted, running | 1.75 → 0.46 | 2.30e-2 → 2.34e-2 | step still heats d33–37 (~1 %) |

With PTC, every accepted step moves the front in the physical direction (heating). The linear
solves with split directions are accurate (0.13–0.43 after one pass). The front's net heating,
however, stays at 5–8e-3 of its gross Lyman rate. In the hot start, hydrogen at d35 is still 99 %
neutral at 8.15 kK. Its RE/c is insensitive to T between 7.7 and 8.2 kK (7.3–7.9e-3). The stable root
of the unmodified residual, if one exists, therefore lies where hydrogen at the front begins to
ionize, at least several hundred kelvin and possibly kilokelvins hotter. The PTC approaches it at
3–5 % in T per 35–60 min iterate.

The frozen-J local thermal curve at 167-radius d35 (`scripts/thermal_curve.py` on the hot-start
state, populations re-equilibrated at each T; `thermal_curve_167_d35_hot.json`) locates that root:

| T [K] | 7900 | 9085 | 10 270 | 11 455 |
|---|---|---|---|---|
| RE/gross | +2.0e-3 | +2.2e-3 | +1.4e-3 | −4.7e-3 |
| neutral H | 0.991 | 0.991 | 0.989 | 0.982 |

A stable local root lies at ≈ 10.5 kK (net heating below it, net cooling above). The 12.6 kK point
failed numerically and is not used. Frozen J is the thin limit of the local response, so this
estimates the target and does not establish the coupled solution. A second hot start with
T ≥ 10.5 kK on the net-heated front only (d30–40, `results/nk167_hot2`) tests it directly.

**Caution on the sign of the front balance (GENCOOL check).** At the 10.5 kK hot start, with J
recomputed by the formal-ALI, the RE row reads +0.9…+1.1e-2 of gross at d33–37 (net heating). One
complete native iteration at the same state (`scripts/verify_cool_state.py`, diagnostic only)
writes CMFGEN's thermal-balance breakdown (`scripts/gencool_breakdown.py`,
`gencool_hot2_start_front.json`). It shows strong net *cooling* at the same depths
(3.5e-10–3.0e-8 erg cm⁻³ s⁻¹, dominated by H I collisional terms). The radiative form
Σ(χJ − η) and the thermal form of the energy balance agree only when the populations are in
statistical equilibrium with J. Their gross disagreement shows that this state, whose populations
were equilibrated to the old J, is far from it. The same holds for the original restart's
GENCOOL (3.8 kK front: net heating ≈ 135 % of total from H I and H⁻ bound-free;
`gencool_x0_restart_front.json`). Values of the RE row at the start states and early iterates
therefore describe unconverged states, not the thermal balance of a solution. The measured
statement that stands is the coupled T-slope at frozen populations (§4.4). Whether a
finite-temperature root of the coupled problem exists at the front is decided only by the
coupled solve, which has not converged.
The same check at the 7.9 kK hot run's iterate-4 state (front 8.0–8.9 kK, after four PTC iterates;
`gencool_hot_it04_front.json`) gives the same picture:

| d30–40 | value |
|---|---|
| GENCOOL net cooling | 8e-12 to 8e-10 erg cm⁻³ s⁻¹ (≈ 190 % of total, H I collisional) |
| RE row / gross | +4.7e-4 to +7.9e-3 (net heating) |

Radiative heating and collisional cooling both flow into the H I level populations, so the front
populations are still not in statistical equilibrium.

## 5. Audits

*(Pending a converged released state.)* The closing chain is ready:
`scripts/post_convergence.sh RUN BASE X0 REPLACE` runs `trace_cleanup.py` (all-level
correction), then `verify_cool_state.py` (one complete native iteration with spectral export) and
`audit_cmfgen_cool.py`. The gates are unchanged from Report 20: electron balance ≤ 1 %, global
and cumulative radiation energy ≤ 0.1 %, last native correction ≤ 0.1 % on all levels,
populations/charge ≤ 1e-6, spectral material temperature ≤ 0.1 %. Build 28 carries every
capability flag the audit checks (H-minus repairs included).

## 6. Resolution

*(Pending.)* Prepared:
- **Spatial.** The 84-radius released solve has no converged state (§4.1), so the 84/167 pair
  cannot be formed. A third nested grid was built: `results/cmfgen_ali_lambda333`
  (`scripts/prepare_mesh333.py`; same plan, source model and atomic data, mesh factor 4, numerical
  controls copied from the 167 restart). Radius and velocity are identical to the 167 and 84 grids
  at every shared point (`mesh333_control.json`). The 333 equation table and the T-released
  reference base pass (`results/nk333/it00_base`, 1162 s on 4 threads) are done.
  `scripts/transfer_state_nested.py` carries a 167 state onto it.
- **Frequency.** The doubled-frequency settings of Report 20 apply at either grid; the
  doubled-frequency solve starts from the converged state of the same grid.
- **Comparison.** `scripts/compare_resolution.py --nested` (T, n_e and dominant-ion fractions at
  shared radii; threshold `profile_refinement_relative` = 1 %).

## 7. Wind, evolution, benchmark, v1

Not reached: the self-selected wind is gated on a converged, audited, resolution-checked released
state. A survey of the code shows that the pieces exist but are not connected:
- `cmfgen_backend.SpectralSnapshot.moments` returns the CMFGEN radiative acceleration;
- `material_wind.continuum_coefficients` takes a local g_rad;
- `spherical_wind.match_base` solves the mass flux as an eigenvalue, but for an LTE diffusion
  wind only.

No code selects ṁ and the sonic point with CMFGEN's g_rad. Coupled evolution uses `winds.rate`
(off / prescribed / Dotan 2011). No published evolutionary benchmark has been run
(`fiacconi_evolution`, `dotan_evolution` pending). **v1 is not declared.**

## 8. Gate table

| Gate (this request) | Status | Evidence |
|---|---|---|
| Level-dissolution derivatives in the Jacobian/RE response | **PASS** | §2; `action_audit_b27.json`, `spectral_ne_direction.json`, `build_27_capabilities.json` |
| `INC_AD`/`INC_ADV` restored and validated | **PASS** | §3; `action_audit_adv.json`, `adv_offdiag_test.json` |
| Threaded tangent robustness (nested-OpenMP deadlock/SIGSEGV) | **FIXED** (build 28, products equal to 1.2e-16) | §9; `omp_deadlock_sample.txt`, `build28_vs_build27_product.json` |
| T continuation to λ = 1 (84 and 167 radii) | **PASS** (continuation), no retreat | §4.1–4.2; `nk_adv84_history.json`, `nk167_adv_history.json` |
| 167-radius T-released cool closure solved | **FAIL** (not converged; cold-branch slide at λ = 1, PTC in progress) | §4.2–4.5 |
| 84-radius T-released cool closure solved | **FAIL** (cold-branch runaway; stopped) | §4.1, §4.5 |
| Electron balance, global/cumulative radiation energy, last correction, charge, spectral material temperature | **NOT EVALUATED** (no converged state; chain ready) | §5; `scripts/post_convergence.sh` |
| Spatial resolution (167 vs 333) | **NOT REACHED** (333 grid, table and base pass ready) | §6; `mesh333_control.json`, `nk333_README.txt` |
| Doubled-frequency resolution | **NOT REACHED** | §6 |
| Self-selected ṁ / sonic wind | **NOT ATTEMPTED** (gated) | §7 |
| Coupled evolution and published evolutionary benchmark | **NOT ATTEMPTED** (gated) | §7 |
| v1 | **NO** | — |

## 9. Failed and superseded attempts

Every attempt below is kept on disk with its driver log; partial passes interrupted by a restart
are under `superseded/` with a README line giving the reason.

- **Build 25 alone** (bound-free dissolution terms without the line ratio): the RE row of the action
  got worse in the n_e direction (0.25 → 0.32). Kept in `build25_first_audit.json`; fixed by build 27.
- **Pinning T homotopy** `R_T = λ R_RE + (1−λ) c (T − T_fix)/T_fix` (`--hom-mode pin`,
  `results/nk_adv84` iterates 0–4, `results/nk167_adv` iterates 0–2): the λ = 0 stage converged quickly
  (natural-measure reductions ×15, ×440, ×22 at 84 radii; ×23, ×62 at 167). At λ = 0.25 the
  coupled step moved T by at most 3e-4. Cause: the pin scale c is the gross scale of the RE row (the
  l1 norm of its native blocks), which exceeds the net RE imbalance by orders of magnitude. T
  therefore stays pinned until λ ≈ 1, and the continuation adds nothing. Superseded by the residual
  shift.
- **Residual shift without a preconditioner floor**: at λ = 0 the T row of the residual is zero by
  construction. The preconditioned merit was still 889, because the direct-part preconditioner's
  RE-row T diagonal is near singular at many depths (|M_TT| / ‖row‖₁ between 7e-8 and 2e-2), so
  M⁻¹ turns population residuals into spurious T corrections. Restarted with `--pc-t-floor 1e-3`,
  which adds 1e-3 of the row scale to M_TT in the **preconditioner only**. The Jacobian, the products
  and the residual are unchanged, so the converged state does not depend on the floor. Merit at the
  same state: 889 (floor 0), 29.2 (1e-4), 5.79 (1e-3), 0.907 (1e-2).
- **Exact radiation columns in the preconditioner** (`scripts/probe_columns.py`,
  `results/nk_adv84/probe_it05_trial_l1`): one batched tangent pass (10 directions, depth period 5)
  measured the exact radiation part of the T and n_e columns. The response beyond ±1 depth is small
  (T 0.1 %, n_e 6.7 % of the column norm). Added to the block-tridiagonal M, it made the
  preconditioner *worse* on the six exact products of iterate 5:
  q = ‖W(P⁻¹Jz − z)‖/‖Wz‖ went from 0.24–1.15 (M without floor) to 0.51–1.35. Partial radiation
  response over-cancels, as the local ALI did in Report 20. The same measurement shows that the
  1e-3 T floor also degrades q (0.99–1.92). The floor is kept only because it removes the merit
  artifact; its cost is recorded here. Not adopted.
  As *augmentation* directions in the least squares, the 10 exact probe products lower the
  iterate-5 linear residual only from 0.81 to 0.75 (`scripts/probe_augment_test.py`,
  `probe_augment_test.txt`). The stall therefore does not come from the T/n_e subspace. It is the
  distributed thick-continuum cancellation of Report 20 §4.2 (H I, He I rows).
- **Intermittent native faults, then a deadlock → build 28.** Two 3-direction product passes of
  the 84-radius run died with SIGSEGV (iterates 6 and 8, 732 s and 448 s into the frequency
  loop). Each identical rerun completed (`Runner.run` keeps a failed pass as `<name>_segv` and
  reruns it once, logged). A repeat of the iterate-8 pass with identical inputs reproduces its
  radiation products to 8.5e-17 relative (`scripts/repeat_product_pass.py`,
  `product_pass_repeat.json`), so completed products were not corrupted. A 167-radius product pass
  then hung at 0 % CPU. A stack sample (`omp_deadlock_sample.txt`) shows every thread waiting in
  libgomp team barriers inside the ray-parallel region of `lin_apply`. That region was nested in
  `qs_jv_freq`'s direction loop, which ran as an *inactive* parallel region (`IF` clause false)
  whenever there were fewer directions than threads. The nested team is the defect: libgomp on
  macOS deadlocked or crashed in it intermittently. **Build 28** (`jv_nonest.patch`, sha256
  `107289a7ab29…ed67`, threading only) runs the direction loop outside any OpenMP construct in that
  case. The pass and its sample are kept (`superseded/it04_jv02_omp_deadlock`).
- **Process slips**: one relaunch script passed its option list as a single argument under zsh and
  failed at argument parsing before running anything. Relaunched with the arguments written out.

## 10. Reproduction

Native builds 25–27 are layered on build 23/24, with the same toolchain as Report 20:

```bash
export MESASDK_ROOT=$HOME/qsrd-runtime/mesasdk; source $MESASDK_ROOT/bin/mesasdk_init.sh
export MACOSX_DEPLOYMENT_TARGET=14.0
export QSTARS_CMFGEN_ROOT=$PWD/backends/cmfgen/build_27/cur_cmf
make -C $QSTARS_CMFGEN_ROOT d_new_main
python backends/cmfgen/record_build.py --help     # capability record (build_27_capabilities.json)
```

Patches and manifests are `backends/cmfgen/jv_dis.patch` (build 25: `qs_jv.f90` dissolution terms,
`CONT_FREQ` passed from `cmfgen_sub.f`), `jv_vall.patch` (build 26: FL rescaling and the optional
`VCHI_ALL`) and `jv_dis_lines.patch` (build 27: line `L_STAR_RATIO`), each with its
`*_patch_manifest.json`. `QSTARS_JV_NODIS=1` disables the dissolution terms (it reproduces build 23
to 8e-17).

Validation of §2–3:

```bash
.venv/bin/python scripts/audit_action_oracle.py results/jv/tfree --out Reports/evidence/report21/action_audit_b27.json
.venv/bin/python scripts/audit_action_oracle.py results/jv/tfree_adv --out Reports/evidence/report21/action_audit_adv.json
```

T-released solves (§4), 84 and 167 radii; the 167 run adds
`--x0 results/cmfgen_ali_lambda167 --replace results/nk167_adv/transfer/replace_167.txt --col-workers 3`:

```bash
.venv/bin/python scripts/cool_newton_krylov.py --out results/nk_adv84 \
  --exe backends/cmfgen/build_27/cur_cmf/exe/cmfgen_dev.exe --threads 6 --krylov 4 --eta 0.25 --augment 2 \
  --warm --free-t --t-homotopy --hom-mode shift --lam0 0 --dlam 0.25 --lam-tol 1e-2 --stage-max-iters 2 \
  --pc-t-floor 1e-3 --lm-trust 0.5 --ultra-frac 1e-20 --ultra-weight 1e-8 --natural \
  --vadat INC_AD T --vadat INC_ADV T --start-base results/nk_adv84/<base> --start-k <k>
```

Every restart is logged with its reason in each run's `README.txt`.

Regression suite (2026-10-04, MESA SDK compiler initialized, `MACOSX_DEPLOYMENT_TARGET=14.0`):
**111 passed, 23 skipped** (the skips need compiled MESA/FreeEOS libraries). New this report:
`test_split_pieces_partition`, `test_lm_step_newton_limit_and_trust_radius`,
`test_release_t_removes_inherited_ba_pointer`, `test_t_homotopy_residual_and_jacobian_term_consistent`,
`test_ptc_residual_jacobian_and_heating_direction`. Closing chain on a converged
state: `scripts/post_convergence.sh RUN BASE X0 REPLACE`, which runs `trace_cleanup.py`, then
`verify_cool_state.py` and `audit_cmfgen_cool.py`. Resolution comparison:
`scripts/compare_resolution.py` (`--nested` for 84 vs 167).
