# Report 20 — the exact Jacobian action of the cool closure, and what stops the released solve

**Outcome: JACOBIAN PASS (RE row CONDITIONAL); released cool closure NOT CONVERGED —
narrowest blocker: the cool-layer thermal/ionization front at 84 radii. QStarS remains 0.1
development software, not v1.**

## 0. Outcome in one paragraph

The Report 19 formal-transfer tangent is now inside CMFGEN's population/rate/RE
residual (builds 18–23). Every pass can return the exact directional derivative of
the cool residual: populations, n_e and T; extinction, emissivity, scattering, all
lines, the diffusion boundary and blue-frequency history. It agrees with
exact-residual finite differences on the real 84-radius cool state at the oracles' own
consistency floor (1e-5 → 1e-6 as ε shrinks; Report 18's directions, where every
native linearization failed, at their floor). Getting there exposed and fixed the
native fixed-J blocks, which omit every level-dissolution dependence on n_e, the ion
density and T. A dissolution-off diagnostic confirmed that this omission is also the
one significant gap left in the radiative-equilibrium row of the action. With that
action a Newton–Krylov solve took the fixed-T closure from merit 22 to 0.011 (largest
significant frozen-J correction 0.63 → 1.4e-3), not yet to the predeclared criterion.
Releasing T exposed a silent CMFGEN behaviour: a single-iteration restart inherits
`FIXED_T` from `BAMATPNT`, so every earlier "FIX_T=F" pass had held T. With T truly
released, the solve does not converge from either a cold or a hot start. The coupled
Newton direction points to a cool neutral layer whose d11/d12 ionization front falls
between two 84-grid points. The evidence is consistent with that front cell having
no consistent RE temperature on this grid: the cold branch drives it toward T → 0,
and its frozen-J local balance has no cool-side root. This is an interpretation,
not a demonstrated non-existence or instability. No
audits, resolution pair, wind, evolution or benchmark were run, because each needs a
converged released state. The next step is a 167-radius released solve (prepared)
with T continuation, plus the omitted wind energy terms (`INC_AD`, `INC_ADV`).

## 1. What Reports 17–19 left open

Report 19 validated an analytic tangent of the native `QS_FRM_ALI` formal-transfer
map (depth- and frequency-nonlocal, with blue-ray history) on bounded fixtures.
It was not connected to the rate equations. Report 18 had measured every native
linearization against the exact residual of the 84-radius cool state x0 and found
none usable (e ≥ 0.94 along its own step; three of four uphill), and showed that
exact-residual Krylov products were noise-limited (~1e-4) and cost a full native
iteration each.

## 2. Native integration (builds 18–24)

All builds layer read-only hooks on build 17; every patch has a manifest and
`record_build.py` capability flag. The primal state, residual and native step
are unchanged (residual reproduction below).

| Build | sha256 (executable) | Adds |
|---|---|---|
| 18 | `de772b0c495c…ffaa` | `QSTARS_JV_DIRS`: exact radiation part of J·(x∘s) per direction inside the native frequency loop; full-block dump; forced/traced equation selection; stop-after-solve; per-frequency δχ/δη/δσ/δJ export |
| 19 | `bc18b46d9a93…3d79` | linearized tangent (one primal setup per frequency), OpenMP over directions, stored converged J and replay passes (one formal solve per frequency, no ALI), moment linearization skipped when its output is identically zero |
| 20 | `94e8126c3962…4abe` | environment presence checks (a truncated value is present) |
| 21 | `53bfbfe1efaa…a916` | fixed-J fast path: ~49 % of a pass had been `memset` of identically-zero TX/TVX |
| 22 | `151f79054898…2898` | ray-parallel linearized tangent; correct (≤1.5e-12), but inactive in the model because the direction loop always opened a parallel region |
| 23 | `ea0251cfa7ab…c7c1` | `QSTARS_ALI_WARM`: formal-ALI warm start from a stored J (iteration count only; converged J differs by ~1e-9) |
| 24 | `f876a520b581…ffd2` | direction loop parallel only when it fills the threads: products identical (1.4e-16), **no speed-up** (fine-grained ray regions, 279 s system time). Recorded, not used |

### 2.1 Radiation part (analytic)

At every frequency, for each direction, the adapter forms the exact directional
changes of total extinction, thermal emissivity and coherent scattering
(continuum from the native VCHI/VETA, electron scattering, Rayleigh, **all**
lines in resonance — weak lines included — and the native negative-opacity
floor differentiated as written), the inner-boundary B and dB/dR (through
DIFFW), and propagates them with the Report 19 tangent, carrying the ray
derivatives to the next (redder) frequency. The resulting δJ_ν enters the
residual through exactly the coefficients the native accumulation uses:
band-integrated bound-free JREC/JPHOT with QFV_R/QFV_P, the low-temperature
recombination branch, the line ZNET sums closed at the end of each resonance
zone (with the native line-cooling scaling), the radiative-equilibrium row and
the two-photon rates. The two-photon conjugate occupation that the native
residual reads back from the J file during the sweep is held fixed (small,
documented omission).

### 2.2 Direct (fixed-J) part — a new finding

With `QS_ALI_LIN=T`, `QS_LIN_FAC=0` the native blocks contain no radiation
response, and they are exactly block-diagonal in depth. **They are nevertheless
not the fixed-J Jacobian.** Frozen-J central differences (replay passes at
x(1±h e) that reread the base J — no transfer noise) show

| Column (all depths) | ‖FD − native‖/‖FD‖ |
|---|---|
| electron density | 1.65 |
| H II | 0.40 |
| He II ground | 1.00 |
| temperature | 0.063 |

The omitted dependence per unit absolute population is identical for H II and
He II ground to **1.5×10⁻⁹**: it enters only through the total ion density
(POPION), i.e. through level dissolution and occupation probabilities
(`COMP_LEV_DIS_BLK`, `OCCUPATION_PROB`, `LTEPOP_WLD_V2`), which the native
linearization treats as constants. The complete direct Jacobian therefore
replaces the n_e and T columns by frozen-J differences and adds the common
POPION term to every ion column (H II, all He II levels, He III): 4 (fixed T)
or 6 frozen-J passes per state.

## 3. Validation chain on the real cool state

Everything below is on the 84-radius Report 18 state x0
(`results/cmfgen_ali_lambda84_long_maxlam1e3`, 43 314 frequencies), against
archived exact-residual finite differences that rerun the full converged
formal-ALI transfer (Report 18 `results/jac_check`, ε = 1e-2 and 1e-3 of the
largest component; new T-direction runs `results/jv/oracle_t`). Evidence:
`Reports/evidence/exact_jacobian_20261001/`.

**Primal reproduction.** Build 18 reproduces the build-16 residual to 1.2e-11
(L2); build 19 reproduces build 18 to 8.0e-12; a replay pass reproduces its base
residual exactly (difference 0.0). Skipping the identically-zero moment
linearization changes the fixed-J blocks by at most 2.9e-15.

**Tangent implementation.** The linearized tangent equals the dual-number
reference to ≤3.0e-12 on all 84 fixture runs (7 families × 6 directions × 20/39
radii) and to 3.9e-13 per frequency in the full model; the radiation parts of the
residual derivative agree to 4e-16–8e-16. The ray-parallel version agrees to
≤1.5e-12 with 4 threads.

**Spectral layer** (all frequencies and depths, analytic vs finite differences of
the oracle's own exported spectra; relative L2 of the response):

| Direction | δη | δσ | δχ | δJ (ε=1e-2) |
|---|---|---|---|---|
| frozen-J Newton step (fac0) | 4.3e-8 | 2.2e-7 | 4.3e-4 | 3.8e-5 |
| local-ALI fac 0.5 step | 4.4e-8 | 1.1e-7 | 5.2e-4 | 2.6e-5 |
| local-ALI fac 1 step (max s 5400) | 1.5e-6 | 1.2e-8 | 7.8e-4 | 3.9e-3 |
| native moment step | 3.9e-7 | 1.5e-9 | 9.5e-5 | 6.1e-4 |

The δχ residual is the level-dissolution derivative that the native opacity
derivatives omit (largest just below the H I Lyman edge at depths 16–20); it
affects δJ at the 1e-4 level.

**Complete Jacobian action** (direct part completed by frozen-J columns; error
e = ‖J·s − FD‖/‖FD‖ over active rows):

| Direction | e (larger ε) | e (smaller ε) | oracle self-consistency |
|---|---|---|---|
| T profile (T column; RE row inactive) | 4.8e-5 | 4.8e-6 | 4.3e-5 |
| n_e profile | 3.0e-4 | 3.0e-5 | 2.7e-4 |
| all populations, random | 1.8e-5 | 1.8e-6 | 1.7e-5 |
| T + n_e + populations | 2.4e-5 | 2.4e-6 | 2.1e-5 |
| Report 18 fac0 step | 1.3e-3 (ε=1e-2) | 1.6e-2 (ε=1e-3) | 1.6e-2 |
| Report 18 fac0.5 step | 1.4e-4 | 1.0e-3 | 9.6e-4 |
| Report 18 fac1 step | 7.5e-3 | 6.2e-2 | 6.0e-2 |
| Report 18 moment step | 6.2e-5 | 4.3e-4 | 4.1e-4 |

**Correction (found during the solve, §4.3).** The four "T released" oracle runs
requested `FIX_T=F`, but a single-iteration restart inherits `FIXED_T=.TRUE.` from
the restart's `BAMATPNT`, so the radiative-equilibrium row was still the identity
with zero residual. These four rows therefore validate the **T column** (every
other row's response to T, n_e and populations) and not the radiative-equilibrium
row. The RE row is validated separately with T actually released (§4.3).

In these four directions the discrepancy falls tenfold with the oracle step:
it is the oracle's O(ε) nonlinearity; the action itself is accurate to ~1e-6. In
the Report 18 directions the agreement equals the oracle's own floor. For
comparison, Report 18 measured e = 0.94 (frozen J), 1.72 and 50 (local ALI at
0.5 and 1) and 117 (native moment) along the same directions. With the native
fixed-J blocks instead of the completed ones, e = 15.3 (fac0) — the direct gap,
not the radiation response, dominated that error.

## 4. Newton–Krylov solve

### 4.1 Method (`scripts/cool_newton_krylov.py`)

The unknowns are the populations, n_e (and T when released) at all 84 radii,
updated multiplicatively, x ← x·exp(−λs), so positivity holds by construction
and nothing is clipped. Each outer iterate k does the following at the accepted
state x_k:

1. **Base pass**: one complete native iteration at x_k (reached from the
   archived restart through the read-only perturbation hook), stopped after the
   linear solve. It gives the exact residual R(x_k) with the radiation field
   converged by the formal-ALI to 1e-10, the native fixed-J blocks, and the
   converged J stored for replay.
2. **Direct Jacobian**: frozen-J central differences in the n_e and H II columns
   (and T when it is released), h = 1e-4. Replay passes reread the base J, so
   these contain no transfer noise. The result completes the native blocks as in §2.2.
   M = the complete fixed-J Jacobian (block-diagonal in depth, no radiation
   response), factored with per-pivot equilibration.
3. **Linear solve** for the exact Newton step J s = R. Products J·z = J_dir z +
   (radiation part from one tangent replay pass) are exact actions of the
   complete Jacobian, validated in §3. The objective is left preconditioned and
   significance weighted: ‖W M⁻¹(R − J s)‖. W = 1 on rows of levels holding
   more than 1e-10 of their element (and on n_e, T) and 1e-3 on trace levels
   (the T-released runs add a third tier, 1e-8 at or below 1e-20; §4.4). The raw residual is unusable as an objective because its rows span
   ~50 decades. Iterates 0–5 used GMRES. From iterate 7 an augmented
   minimal-residual (GCR-type) solve batches the preconditioned residual with
   the two previous Newton steps in one pass. Later passes each add the weighted
   preconditioned residual of the current least-squares solution.
4. **Line search**: three trials λ ∈ {1, ½, ¼} run in parallel, each a complete
   native base pass (ALI warm-started from the base J, build 23; the warm start
   changes the converged J by 5.9e-10 relative). The trial with the smallest
   merit ‖W M_k⁻¹ R‖ is accepted if it decreases the merit. A per-component
   log-step cap ln(10⁴) applies only to the update of ultra-trace levels, which
   ask for e⁵⁷ changes at d14. It was never active after iterate 0.
5. **Convergence criterion** (fixed before any results): the exact Newton step
   satisfies max|s| < 1e-3 over all levels holding more than 1e-10 of their
   element, plus n_e. The all-level maximum is reported separately, and so are the
   native frozen-J correction M⁻¹R (which also vanishes at a solution) and the
   unweighted residual.

The merit changes with M_k between iterates, so cross-iterate progress is also
measured with the iterate-0 blocks fixed (`scripts/nk_metrics.py`, ‖M₀⁻¹R‖ by
region).

### 4.2 Fixed-temperature solve (`results/nk_fixT`)

The temperature is held at x0's profile; the populations, n_e and the radiation
field are solved together. Every accepted state is a complete native pass with
the formal-ALI converged to 1e-10.

| it | merit ‖W M_k⁻¹R‖ | linear residual (products) | max\|s\| significant / all | λ | merit after | ‖M₀⁻¹R‖ / max_sig\|M₀⁻¹R\| at accepted state |
|---|---|---|---|---|---|---|
| 0 | 22.0* | 0.049 (2) | — / 11.2 | ½ | 5.63* | 22.0 / 0.634 (x0) |
| 1 | 29.8* | 0.058 (4) | — / 1.94 | 1 | 12.3* | 5.67 / 0.483 |
| 2 | 70.3* | 0.076 (1) | — / 4.64 | ½ | 28.0* | 2.61 / 0.235 |
| 3 | 12.7* | 0.073 (4) | — / 1.23 | 1 | 4.41* | 2.28 / 0.212 |
| 4 | 11.0* | 0.096 (5) | 0.90 / 1.49 | 1 | 4.12* | 1.27 / 0.237 |
| 5 | 10.4* | 0.068 (5, recycled) | 3.90 / 6.22 | ¼ | 8.19* | 0.94 / 0.124 |
| 6 | 0.343 | 0.272 (4) | 1.92 / 1.92 | ½ | 0.191 | 1.08 / 0.235 |
| 7 | 0.197 | 0.535 (4, augmented) | 0.174 / 0.612 | 1 | 0.111 | 2.03 / 0.136 |
| 8 | 0.113 | 0.213 (4) | 0.201 / 0.364 | 1 | 0.0256 | 1.23 / 0.0259 |
| 9 | 0.0267 | 0.540 (4) | 0.569 / 0.569 | 1 | 0.0164 | 1.24 / 0.0204 |
| 10 | 0.0168 | 0.617 (3 passes, 56 split directions) | 0.115 / 0.229 | 1 | 0.0105 | 1.48 / 0.0059 |
| 11 | 0.0114 | | | | | 1.34 / 0.00139 |

\* unweighted merit (iterates 0–5); from iterate 6 the significance-weighted merit.
The merit is measured with each iterate's own M_k. The last column uses the
fixed x0 blocks M₀. "max_sig" covers levels holding more than 1e-6 of their
element (`scripts/nk_metrics.py`).

**Why the linear solve stalls at 0.2–0.6.** At the outer He I ground row (d6) the
direct and radiation parts of the exact action are −15.007 and +15.003 for a
random direction (§3 oracle; exact sum −3.28e-3, oracle −3.28e-3). In optically
thick continua the population is in detailed balance with a J that follows it, so
the fixed-J preconditioner overstates these rows by up to ~1e4. In the random
direction about 200 significant entries cancel by more than a factor of 100
(H I at d4–12 and d22–52, He I at d4–9 and d24–53, He II at d68–83). This
cancellation is not diagonal: per-entry ratios change sign between neighbouring
depths. Three alternatives were measured on exact products at iterates 9–10 and
rejected:
- native local-ALI blocks (factor 0.5 or 1) as preconditioner;
- the native moment linearization as preconditioner (q from 0.001 to 3.1;
  near-singular modes);
- region × species split directions (0.74 against 0.75 for the same pass).

Multi-preconditioned directions (moment/ALI) added nothing (0.735 → 0.733).
Recycling earlier iterates' directions helps modestly (0.73 → 0.42 with 18
stale products). The solve therefore converges linearly, with the merit falling
by ~0.6 per iterate. It reached the significant-level 1e-3 scale at iterate 11
(fixed-M₀ max 1.4e-3), with the remaining merit two-thirds in levels near the
1e-10 significance edge (H I n=3–4 at d24–27) and in trace He II. The solve was
stopped at iterate 11 to release T, **not converged** to its criterion: the last
complete Newton step (iterate 10) was 0.115 over significant levels and 0.229 over
all levels.

### 4.3 The temperature was never released in a single-iteration restart

Before starting the T-released solve, the radiative-equilibrium row of every
"`FIX_T=F`" pass turned out to be identically zero in the residual and the identity in
the BA. `READ_BA_DATA_V3` restores `FIXED_T` from the restart's `BAMATPNT`
("Temperature was held fixed", written by x0's last lambda iteration), and
`FIXPOP_IN_BA_V3` then fixes T at all depths. A multi-iteration native run
clears this after its first iteration (`FIXED_T=RD_FIX_T`), so earlier reports'
multi-iteration runs are not affected. Every single-pass design here
(oracles, NK passes) is. The fix (`release_t`) removes `BAMATPNT`/`BAMAT` from a
T-released pass. Every pass computes its own BA, and `FIXED_T` keeps its
initial `.FALSE.`. Affected passes are preserved in
`results/superseded_T_inherited/` (README).

**Validation with the RE row active** (`results/jv/tfree`, `scripts/audit_tfree_jv.py`,
evidence `tfree_jv_audit.json`). One T-released base pass at x0 carries the
products of the four §3 directions. Six frozen-J column passes (n_e, H II, T)
and eight exact-residual passes follow, all with T actually released.

| Direction | e all rows (ε_a / ε_b) | e RE row | e RE row, direct part only | oracle consistency (all / RE) |
|---|---|---|---|---|
| T | 4.8e-5 / 4.8e-6 | 1.1e-4 | 1.0 | 4.3e-5 / 5.4e-6 |
| n_e | 3.0e-4 / 3.0e-5 | 0.25 | 35 | 2.7e-4 / 5.9e-5 |
| random populations | 1.8e-5 / 1.8e-6 | 4.5e-3 | 26 | 1.7e-5 / 2.1e-5 |
| T + n_e + populations | 2.4e-5 / 2.4e-6 | 3.7e-3 | 4.3 | 2.1e-5 / 5.6e-6 |

The population, n_e and charge rows are at the oracle floor. The RE row is a
near-cancellation: its direct part alone is wrong by 1–35 because the radiation
part cancels it. The radiation part of the RE row is accurate to ~0.1 % of its
own size for T and ~0.25–0.75 % for n_e. The error grows inward through the
thick interior (d55–84), and in the cool layer (d10–25) it is ~0.4 % of the net
response. Because of the 30–300× cancellation, the net RE-row response at d31–50
(T direction) and in the deep interior (n_e direction) is wrong by up to 16 % and
25 % respectively. The most likely source, consistent with every measurement here,
is that the native opacity/emissivity derivatives (`VCHI`/`VETA`) feeding the
tangent omit the level-dissolution dependence on n_e, the ion density and T, the
same omission found in the fixed-J blocks (§2.2). The evidence is that using the
exact (frozen-J) H II column in the RE row makes the RE error *larger*
(3.2e-4 → 4.5e-3 for the random direction), so the native H II entry partly
compensates an omission in the radiation part. **Confirmed by a dissolution-off diagnostic**
(`results/jv/tfree_nodis`, VADAT `DO_DIS=F`, same x0 and directions, evidence
`nodis_jv_audit.json`). Without level dissolution:
- the complete action agrees with the oracle to 4.0e-9 (random populations) and 4.0e-5 (n_e);
- the RE row agrees to 2.3e-5 (random) and 6.8e-5 (n_e);
- the native fixed-J blocks need no completion at all (results are identical with or without the POPION term).

A small T-direction RE-row discrepancy (1–2e-4) remains with dissolution on and
off. Level dissolution is therefore the one significant gap in both the native
direct blocks (§2.2) and the radiation part. Completing it means differentiating
the occupation-probability opacities/emissivities (`COMP_LEV_DIS_BLK` inputs n_e,
ion density, T) inside the tangent; that is the identified next Jacobian
improvement. This is a measured limitation of the action, not a fitted correction. For Newton it means at most linear
convergence of the deep and mid-depth temperatures, and convergence is always
judged on the exact residual.

With T actually released at x0 (`results/jv/tfree/base`), the RE residual is
large. The native frozen-J T correction is −1.1 (fractional) at d1, −4.8 at
d13, −5.2e2 at d17, −26 at d21, −0.32 at d25 and +4 % at d33–45. That is the same
stiff signature as the populations: with J fixed, the thermal balance of the
cool layer barely responds to T. x0's temperature came from a grey/lambda stage
and is not in radiative equilibrium in the cool layer.

### 4.4 Temperature-released solve (`results/nk_freeT`)

The T-released solve starts from the best fixed-T state (`nk_fixT/it11_trial_l1`), with the
RE row active. Three globalization changes were adopted (two further ones failed and are listed in §9).
All are recorded in the run README, and none touches the convergence criterion.

- **Trust region instead of clipping.** The first exact Newton step asked for T changes up
  to e^5.95: heating the outer wind (d1–4) toward 4×10⁴ K and cooling the cool layer
  (d13–22) toward 10³ K. Clipping only the T components (|s_T| ≤ 0.2) left the population
  components sized for the unclipped T change, and the full step raised the merit 84×.
  The driver now takes Levenberg–Marquardt steps in the Krylov subspace already
  computed: min ‖A y − b‖² + μ‖W S y‖², with μ chosen by bisection so that max|s_T| ≤ Δ.
  No extra products are needed. The line-search trials sit at Δ, Δ/2, Δ/4 on this curved
  path. Δ doubles after a full step with actual/predicted reduction > 0.5 and shrinks to
  the accepted radius otherwise.
- **Three-tier significance weights.** At iterate 1 the radius-0.5 trial had merit 289,
  92 % of it from He II/He III levels at d11 holding 10⁻⁴³–10⁻⁵⁴ of helium. These
  levels cannot affect any other equation. Levels at or below 1e-20 of their element now
  weigh 1e-8 (1e-3 between 1e-20 and 1e-10, 1 above). Re-scored this way, that trial
  improves the merit (5.97 → 5.28). The final all-level native audit still covers these
  levels; they must be converged in a closing fixed-J phase (§5), which was not reached.
- **Natural monotonicity test** (from iterate 4, `--natural`). At iterate 3, 86 % of the
  preconditioned merit was one entry: the RE row at d12 (T = 2870 K), where the fixed-J
  T derivative nearly vanishes (M⁻¹R = +15.9). Row-scaling instead hides T entirely,
  because that row is balanced to ≤ 0.5 % of its gross rates. Trials are therefore judged
  by Deuflhard's test: the simplified Newton correction ‖W J_k⁻¹R(trial)‖, estimated in
  the iterate's Krylov subspace from the existing exact products, must decrease.

**History** (`results/nk_freeT/history.json`; T in K at selected native depths, 1 = outermost):

| it | step | accepted | max\|s\| significant (Newton) | T d1 | T d12 | T d17 | T d22 | T d25 |
|---|---|---|---|---|---|---|---|---|
| start (fixed-T state) | — | — | — | 12398 | 4974 | 3798 | 5085 | 6447 |
| 0 | T clip 0.2 | λ = ¼ | 5.95 | 13034 | 4732 | 3613 | 4837 | 6777 |
| 1 | LM 0.5/0.25/0.125 | r 0.25 | 1.29 | 15733 | 3685 | 3182 | 3788 | 7082 |
| 2 | LM 0.25/0.125/0.0625 | r 0.25 | 4.33 | 16817 | 2870 | 3019 | 3424 | 7183 |
| 3 | LM 0.5/0.25/0.125 | r 0.25 | 6.33 | 16933 | 2235 | 3000 | 3381 | 7190 |
| 4 | LM + natural test | none useful | **89.5** (T at d12) | — | — | — | — | — |

Temperatures are those of the state accepted at the end of each iterate.

Every Newton step points the same way: heat the outer wind, cool the cool layer.
At iterate 4 the exact step at d12 asks for s_T = +89.5 (linear target T < 0)
together with full recombination (s ≈ +30 for n_e and H II). Every trial, at any
radius, increased the natural measure, except steps cut 16-fold (which changed it
by 1e-4 relative). The solve was stopped there for the diagnostic below.

### 4.5 Local thermal balance of the cool layer (diagnostic)

`scripts/thermal_curve.py`: at one depth, with the stored converged J and every
other depth fixed, T is set to f·T, and the local populations and n_e are
Newton-solved with frozen-J replay passes. The RE residual is reported relative to
the RE row's gross scale (row l1-norm). Frozen J is the optically thin limit of the
local response, so this is a diagnostic of the local balance, not a solution.
Evidence: `Reports/evidence/exact_jacobian_20261001/thermal_curves/`. All points listed
are locally converged (local frozen-J step ≤ 2e-6). Unconverged points are omitted:
at d17 and d22 every point below 6 kK kept a local step of 0.8–1.0 (one level being
driven toward zero), and at d17 the points above 10.5 kK did not converge.

| Depth (state) | T [K] → RE/gross |
|---|---|
| d12 (start state, before any T step) | 2487 → +7.9e-8; 4974 → +9.5e-8; 7462 → +9.2e-8; **9949 → −3.2e-6** (local step 4e-3) |
| d12 (iterate 4) | 1118 → +7.4e-8; 2235 → +7.0e-8; 4470 → +8.2e-8; 7152 → +1.1e-7; 7823 → +3.8e-8; **8940 → −6.4e-7** |
| d12 with adiabatic cooling (`INC_AD=T`) | 1118 → +7.4e-8; 2794 → +6.9e-8; 7152 → +1.1e-7; **10058 → −3.5e-6** |
| d17 (iterate 4) | 6000 → +1.2e-3; 9000 → +1.9e-3; 10500 → +7.7e-4 |
| d22 (iterate 4) | 6762 → +2.4e-5; **10144 → −7.2e-5** |

(RE > 0 is net heating, STEQ_T = Σ FQW(χJ − η).) With J held fixed, d12 is
**net-heated at every temperature from 1.1 to 7.8 kK**: the balance is almost flat
and rises slowly with T. d17 and d22 are net-heated at every converged point below
10.5 and 6.8 kK respectively. The only local root is on the hot side: ≈7.9 kK at d12 (stable:
heating minus cooling falls through zero) and between 6.8 and 10 kK at d22. At d17 the
balance is still heating at 10.5 kK; hotter points did not converge locally. Adiabatic
cooling, which the closure omits (`INC_AD=F`), changes the d12 balance by ≲30 % below
7 kK and does not create a cool root. The T-released Newton iteration was therefore
moving the cool layer *away* from the only local equilibrium: its coupled Jacobian,
dominated by the extremely thick Lyman continuum and lines (τ ~ 10⁵ at n_H ~ 10⁸ cm⁻³),
had the opposite T sensitivity. Frozen J misrepresents exactly those thick frequencies,
so this diagnostic does not establish the coupled solution. It does show that the cool
branch the solve was following has no local equilibrium above ~1.1 kK at d12.

### 4.6 Temperature-released solve from the hot branch (`results/nk_freeT_hot`)

Initial guess, not a fitted parameter: the cold-branch iterate-4 state with T raised
to ≥ 7900 K (the d12 local root) at d10–24. Populations and n_e at all depths were
re-equilibrated with J frozen and T held (`scripts/equilibrate_at_T.py`, 5 sweeps,
final local step 5e-5). The same solver as §4.4 was used.

| it | natural measure | trials (radius → measure) | accepted |
|---|---|---|---|
| 0 | 83.3 | 0.5 → 867; **0.25 → 21.1**; 0.125 → 74.5 | r 0.25 |
| 1 | 55.0 | 0.25 → 82.3; 0.125 → 56.6; **0.0625 → 54.75** | r 0.0625 |
| 2 | 56.1 | 0.0625 → 56.4; **0.031 → 55.8**; 0.016 → 55.9 | r 0.031 (stopped) |

The first step removed most of the inconsistency of the guess. The coupled exact
Newton step then points back toward a cool layer: s_T = +0.74…+1.18 at d15–21
(×e⁻¹, about 2.5–3.7 kK), opposite to the frozen-J local balance of §4.5. The step
is dominated by the He II ground state at d15–19 (1e-6–1e-5 of helium, s = 10–14).
Near 7.9 kK its population scales as roughly e^(36 Δln T) (24.6 eV against
kT ≈ 0.68 eV), so any useful T step is strongly nonlinear for it. He III at d1
(1.6e-10, s = −45) is a smaller contributor.

**Reading of §4.4–4.6 (interpretation, not a demonstrated result).** From both
starting branches, the coupled exact Newton step, with the very thick Lyman radiation
responding, points to a cool neutral layer (~2.5–4 kK). Its ionization front
(d11/d12) falls between two 84-grid points, and the evidence is consistent with that
front cell having no consistent RE temperature on this grid: the cold-branch
iteration drives it toward T → 0, and its frozen-J local balance has no root on
the cool side. Progress per iterate is limited to Δln T ≈ 0.03–0.25 by the
exponential T sensitivity of the ionization balance. No converged released state
was obtained, and none is claimed. Nothing here shows that a cool solution does not
exist, or that the layer is physically unstable. Deciding that needs the
167-radius solve and the omitted wind energy terms.

## 5. Audits

The electron-balance, global and cumulative radiation-energy, last-correction and
charge gates (`scripts/verify_cool_state.py` → `scripts/audit_cmfgen_cool.py`)
require a converged released state. None exists, so **no audit was run**. The
gates and thresholds are unchanged: electron balance ≤ 1 %, global and cumulative
radiation energy ≤ 0.1 %, last native correction ≤ 0.1 % on all levels, charge
≤ 1e-6. Running them on an unconverged state would only certify a mixed-iteration
state. `verify_cool_state.py` now releases T correctly (`release_t`) and accepts
the doubled-frequency VADAT settings, so it is ready for a converged state.

`scripts/trace_cleanup.py` is the closing phase for the all-level correction
gate. It solves the fixed-J trace-level subsystem at the converged significant
state, excluding trace rows driven by significant variables (CMFGEN's
species-conservation rows, §9). On the cold-branch iterate-4 state the filtered
solve keeps 4888 of 6244 trace rows, with largest step 1.5 instead of 1e40
unfiltered. It has not been run, because there is no converged state to clean.

## 6. Resolution

Not reached: there is no converged released 84-radius state to compare.
Prepared and preserved:
- `results/nk167`: T-released 167-radius base pass at the nested restart
  (`results/cmfgen_ali_lambda167`; shared radii, density and velocity identical
  to the 84 grid), 1.2 GB, 682 s with 2 threads;
- `results/nk167/replace_167.txt` (equation-selection map);
- `scripts/transfer_state_167.py` (84 → 167 state transfer);
- `results/nk_dfreq`: T-released doubled-frequency base pass at x0
  (`FRAC_SP`, `AMP_FAC` square-rooted; `MAX_BF`, `dV_LEV`, `FRAC_DOP`,
  `dV_CMF_PROF`, `dV_CMF_WING` halved; 58 MB stored J, about twice the
  frequencies), 587 s with 2 threads.

The driver splits depth regions consistently on the nested grid (tested). The
d11/d12 thermal front sits between two 84-grid points (6466 K and 2235 K at
iterate 4 of the cold branch). That makes the 167-radius solve important for
the cool layer, not a formality.

## 7. Wind, evolution, v1

Not attempted. Each depends on a validated released cool closure, which this
campaign did not reach (§4.4–4.5). The Mach-one MLT guard, the prescribed-flow
scope of the CMFGEN controls and the absence of a production mass-loss rate are
unchanged. The closure still omits adiabatic cooling and advection (`INC_AD=F`,
`INC_ADV=F`); the validation plan requires both before any wind use. No
published evolutionary benchmark was run. **QStarS is not v1.**

## 8. Gate table

| Gate (task order) | Status | Evidence |
|---|---|---|
| Report 19 tangent integrated into the full rate/RE Jacobian action (populations, T, n_e; emissivity, extinction, scattering, lines and continuum, diffusion boundary, blue history, nonlocal response, direct terms) | **PASS** | §2; builds 18–23; `qs_jv.f90` |
| Complete action vs exact-residual oracle on the real cool state | **PASS** (population, n_e, charge rows; T column): e at the oracle floor (1e-5 → 1e-6 with ε); Report 18 directions at their floor | §3; `exact_jv_audit.json`, `tfree_jv_audit.json` |
| RE row of the action (T released) | **CONDITIONAL**: cool layer ~0.4 %, but up to 16–25 % net error in the thick interior. The cause is missing level-dissolution derivatives in the native opacity derivatives, confirmed with dissolution off (2e-5–7e-5) | §4.3; `nodis_jv_audit.json` |
| Native fixed-J blocks | **CORRECTED**: n_e/T columns and the POPION term completed (1.5e-9 hypothesis test); not needed without dissolution | §2.2, §4.3 |
| Fixed-T Newton–Krylov (populations + radiation) | **NOT CONVERGED** to the fixed criterion (|s| < 1e-3). Fixed-x0-block significant correction 0.634 → 1.4e-3; merit 22 → 0.011; stopped to release T | §4.2; `nk_fixT_history.json` |
| Released populations + T + radiation converge together | **FAIL**. Cold branch: Newton drives the d12 front to T → 0; the frozen-J local balance there has no root from 1.1 to 7.8 kK. Hot branch: see §4.6 | §4.4–4.6 |
| Electron balance, global/cumulative radiation energy, last correction, charge | **NOT EVALUATED** (no converged state) | §5 |
| No negative / non-finite radiation; no stale-state convergence | **PASS for every accepted pass** (complete native formal-ALI passes to 1e-10, strict reuse checks, no clipping of J or populations) | §4 |
| 84 / 167 radii and doubled-frequency resolution | **NOT REACHED** (T-released base passes prepared) | §6 |
| Self-selected wind, coupled evolution, published evolutionary benchmark | **NOT ATTEMPTED** | §7 |
| v1 | **NOT DECLARED** | §7 |

**Narrowest remaining blocker.** The temperature of the cool recombination layer
when T is released, most plausibly the unresolved d11/d12 thermal and ionization
front at 84 radii (§4.4–4.6).
The populations-plus-radiation problem at fixed T is solved to about 1e-3, and the
Jacobian is exact apart from the documented dissolution gap in the RE row.

**Next, in order:**
1. 167-radius T-released solve from the transferred cold-branch state
   (`transfer_state_167.py`, `results/nk167` base and table ready). The front
   gets two more cells, and the doubled-frequency run (`results/nk_dfreq`) follows
   the same path.
2. T continuation instead of free trust-region steps: converge populations and
   radiation at fixed T (the solve that works, §4.2), then update T on the reduced
   RE system and repeat. Each inner solve reuses the exact action.
3. Restore the wind energy terms the validation plan requires before wind use:
   adiabatic cooling (`INC_AD`, local in T but couples d to d+1 through dT/dR)
   and advection (`INC_ADV`, block-tridiagonal). This needs the direct-part
   completion generalized from block-diagonal to block-tridiagonal.
4. Differentiate the level-dissolution factors (`COMP_LEV_DIS_BLK` inputs) in the
   tangent's opacities/emissivities and the direct blocks. That removes the
   POPION completion and the RE-row gap (§4.3).
5. Only after a converged released state: the electron/energy audits, the
   resolution pair, then the wind eigenvalue problem with the Mach-one MLT guard.

## 9. Failed and superseded attempts (preserved)

Nothing below was used for an accepted result. Each item is archived with a README.

**T-released solve** (`results/nk_freeT/superseded/README.txt`, `results/superseded_T_inherited/`)
- *Inherited fixed T*: every pass with `FIX_T=F` before the `release_t` fix
  (§4.3), including `results/jv/oracle_t` (still valid as T-column evidence).
- *T-only clip* (iterate 0): clipping |s_T| ≤ 0.2 left the populations sized for the
  full T change; the full step raised the merit 84×. Replaced by the LM trust region.
- *Two-tier weights* (iterate 1): a trial's merit of 289 was 92 % He II/He III
  levels at 10⁻⁴³–10⁻⁵⁴ of helium. Third tier added.
- *Row-scaled merit without floor* (iterate 4): dominated by He II↔He III
  balance rows with gross scale 10⁻¹³–10⁻⁴⁴. Floor added. With the floor, the
  merit no longer measured the ill-conditioned T distance; replaced by the
  natural monotonicity test.
- *Interleaved trace sweep* (iterate 4, `it04_trace_INVALID`): solving all
  trace rows with trace unknowns only is wrong. CMFGEN stores the species
  number-conservation equation in the last ion's row (He III, H II), so He III
  absorbed the whole He conservation error (negative at d4, ×1e40 at d16; trace
  blocks cond ≤ 1.6e50; significant imbalance rose 3.18 → 5.16). The option is
  marked DO NOT USE; `trace_cleanup.py` now excludes such rows.
- *Cold-branch continuation* (iterates 0–4 of `nk_freeT`): drove the d12 front
  toward T → 0 (§4.4–4.5). Preserved as the evidence for the thermal diagnosis.

**Native integration** (`results/jv/`)
- `failed_b19base_unitbug`: build 19's first link. `newunit` returned a negative
  unit number that the J store treated as unset, so the store stopped after
  frequency 1. Fixed with explicit open flags; the final build 19 executable
  replaces it.
- `frozenj_b19_replay_inactive`: with long paths, `get_environment_variable`
  returned status −1 (truncated). Build 19 read that as "unset", so the
  "frozen-J" passes silently reran the full transfer. Build 20 treats a
  truncated value as present. Every column was rerun.
- `frozenj_failed_eps`: a `np.float64` repr in the perturbation-ε environment
  string. Fixed by casting to `float`.
- `superseded_transposed_dirs`: direction files written in C order
  (`ndarray.tofile`) were read as Fortran column-major. Fixed with
  `write_fortran`, plus a layout round-trip test.

**Newton–Krylov** (`results/nk_fixT/superseded/README.txt`)
- *Uniform log-step cap.* Ultra-trace He II levels at d14 asked for e⁵⁷ changes,
  so the cap shrank every step to λ≈0.1–0.4. It was replaced by a per-component
  cap that leaves the exact direction unchanged elsewhere.
- *Stale reuse.* The first driver reused column and product passes made at a
  different base. Caught before acceptance, and the trials were killed. Reuse
  now requires byte-identical state, replay source and directions.
- *Right-preconditioned GMRES* minimized the raw residual, whose rows span ~50
  decades. Replaced by left preconditioning.
- *Recycled low-rank preconditioner* (iterate 5): exact (z, Jz) pairs from the
  current and previous iterates in a Woodbury-corrected preconditioner. The
  merit rose (λ = 0.25 accepted, step 3.9). Disabled.
- *Unweighted merit/objective*: ~100 ultra-trace He II/He III cells at d18 held
  M⁻¹R at 1.0 and dominated both. Replaced by significance weighting
  (fixed rule: 1e-3 below 1e-10 of the element).
- *Local-ALI preconditioners* (iterate 9, `precond_probe`): adding the native
  local-ALI increment to the completed fixed-J blocks does not help.
  q = ‖W(M⁻¹Jz − z)‖/‖Wz‖ over the six exact products at that state is
  0.82–1.01 for fixed J, 0.75–1.00 for ALI factor 0.5, and 35–250 for factor 1
  (nearly singular). The local operator over-cancels thick continua and
  under-cancels lines, matching its Report 18 failure as a Jacobian.

## 10. Reproduction

Native builds (each layered on the previous build tree; manifests list every
patched file and its sha256):

```bash
source ~/qsrd-runtime/mesasdk/bin/mesasdk_init.sh
export MACOSX_DEPLOYMENT_TARGET=14.0
python backends/cmfgen/record_build.py --help
```

Builds 18–23 apply `backends/cmfgen/{jv,jv_replay,jv_env,jv_fast,jv_rays_omp,jv_warm}.patch`
and add `qs_jv.f90`, `formal_tangent.f90` and `formal_tangent_lin.f90` to `new_main`
(see each `*_patch_manifest.json`).

Regression suite (2026-10-02, MESA SDK compiler initialized): **109 passed, 23 skipped**
(the skips need compiled MESA/FreeEOS libraries). New tests cover the block solver,
the Fortran layout, the native linearized-vs-dual tangent (1 and 4 threads), the
nested region split, the LM trust-region step and `release_t`. `tests/conftest.py`
pins `MACOSX_DEPLOYMENT_TARGET=14.0` as `build.sh` does. Without it, the SDK
gfortran rejects this macOS version and the native tangent tests fail to compile.
`qstars validate --require-v1` still exits 2 (29/31 CLI checks; the two CB24
printed-value discrepancies are unchanged).

Validation chain (§3):

```bash
.venv/bin/python -m pytest tests/test_nlte_newton.py -q
.venv/bin/python scripts/audit_native_jv.py --help
.venv/bin/python scripts/audit_direct_columns.py --help
.venv/bin/python scripts/audit_exact_jv.py --help
```

Newton–Krylov solve (§4). The exact command lines are in each run's
`stdout.log`/`driver.log`, and every native pass directory holds its
`run_manifest.json` (executable sha256, environment, VADAT):

```bash
.venv/bin/python scripts/cool_newton_krylov.py --out results/nk_fixT --exe backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe --threads 6 --krylov 4 --eta 0.25 --augment 2 --warm
```

T-released validation (§4.3) and solves (§4.4–4.6):

```bash
.venv/bin/python results/jv/tfree/run_tfree.py
.venv/bin/python scripts/audit_tfree_jv.py
.venv/bin/python results/jv/tfree_nodis/run_nodis.py
.venv/bin/python scripts/audit_tfree_jv.py --dir results/jv/tfree_nodis
.venv/bin/python scripts/cool_newton_krylov.py --out results/nk_freeT --exe backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe --threads 6 --krylov 4 --eta 0.25 --augment 2 --warm --free-t --lm-trust 0.5 --ultra-frac 1e-20 --ultra-weight 1e-8 --natural --start-step results/nk_fixT/it11_trial_l1/QSTARS_STEP.bin
.venv/bin/python scripts/thermal_curve.py results/nk_freeT/it04_trial_r0.25 12 results/thermal_curve/d12_it04
.venv/bin/python scripts/equilibrate_at_T.py results/nk_freeT/it04_trial_r0.25 results/nk_freeT_hot/prep --t-floor 7900 --depths 10 24
```

The exact sequence of restarts and options for each run is in its `README.txt`
(`results/nk_fixT`, `results/nk_freeT`, `results/nk_freeT_hot`).

Progress and fixed-reference metrics:

```bash
.venv/bin/python scripts/nk_progress.py results/nk_fixT
.venv/bin/python scripts/nk_metrics.py results/nk_fixT it00_base it08_trial_l1
```
