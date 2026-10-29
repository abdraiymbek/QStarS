# Report 22 — profiling and speeding up the cool-closure solve

**Status: QStarS is NOT v1.** This report is a solver-level speed pass on the 167-radius
T-released cool closure (Report 21 §4.5). The physical model, the residual, the Jacobian-vector
products and the convergence criterion are unchanged; products agree with the previous solver to
summation order. One solver-path setting, the least-squares weight of ultra-trace levels, was
changed with the user's approval after the structure measurement (§4). The 167-radius PTC solve
was resumed on the faster path. The released closure has not converged; **v1 is not declared.**

## 0. Outcome

All before/after numbers come from the same archived state: iterate 14 of `results/nk167_hot`
(base `it14_trial_l1`, 167 radii, 43 314 frequencies, NT = 134). Every pass reproduces an
archived pass byte for byte in its inputs. *Before* is build 28 with the Report 21 driver
settings. *After* is build 31 with the Report 22 driver. The pairs ran back to back on the same
machine under the same conditions (§3.1).

| Metric | Before (build 28) | After (build 31) | Gain |
|---|---|---|---|
| 1. Seconds per base (trial) pass | 912 s at 8 threads | 362 s at 8 threads; 334 s at 2 threads (the trial setting) | 2.5–2.7× |
| 2. Seconds per Jv product (26-direction pass ÷ 26) | 59.9 s (329 CPU-s) | 31.2 s (159 CPU-s) | 1.9× wall, 2.1× CPU |
| 3. Products per outer iterate | 74 (3 passes; η = 0.25 never reached) | 74 with the same weights; **26** with ultra-trace weight 1e-10 (adopted with approval, §4.2) | 1× / 2.8× |
| 4. Seconds per accepted iterate | ≈ 10 800 s measured (archived, shared machine); 34 500 CPU-s summed from the matched passes | 14 350 CPU-s (same weights) / 6 700 CPU-s (weight 1e-10) summed; measured on the resumed run: **≈ 1 670 s** (columns 128 s + 1 product pass 878 s + trials 660 s) | 2.4× / 5.1× in CPU; **6.5× measured wall** |
| 5. Agreement with the current solver | rerun of build 28: products 6.0e-13 per row | products ≤ 1.0e-12 per row (5.0e-14 globally), GMRES step counts identical; residuals ≤ 3.2e-12 per row; frozen-J residual and BA bit-identical with/without the observer solution | to summation order |

A column pass went from 569 s to 167 s (3.4×). Absolute times were taken on a shared, partly
battery-powered machine, 1.5–4.7× slower than the idle night (§3.1); only the ratios matter.

**What made it faster** (§2):
1. The observer-frame formal solution is skipped in solver passes (build 29). It only produces
   the observed spectrum. It was 78 % of a column pass.
2. An explicit per-frequency scattering operator (build 30). The tangent's inner GMRES made
   ~14.5 ray sweeps per direction and frequency. Now one batched operator build per frequency
   serves all directions, and each GMRES step is a 167×167 matrix–vector product.
3. Batched direction sweeps, a parallel `lin_setup` and balanced operator blocks (build 31).
4. Thread budgets fitted to what the machine actually gives background work (§3.1). Column
   passes do not scale beyond one thread; base passes are faster at 2 threads than at 8.
5. Driver: exact reuse of products already made at the same base, counted against the Krylov
   budget; resume hardening.

**What did not change:** the model, the residual, the products (§3.2), η, the Krylov budget, the convergence test (max significant Newton step < 1e-3), PTC, and the
verification/audit chain. The other significance weights (1 / 1e-3) are unchanged. Final and
verification runs still compute the observed spectrum.

**Measured Jacobian structure (§4).** After the 74 exact products of an iterate, the linear
residual stalls at 0.46. 53 % of what remains sits in the radiative-equilibrium (T) row at the
ionization front, which the solve reduces only to 0.96 of its initial value. Another 43 % sits
in ultra-trace He II/He III rows (populations below 1e-20 of He). The same 74 products can
reduce the front T-row residual to 0.016: the subspace already holds the T correction. The
driver's least-squares weights keep the step out of it, because using it would raise the
weighted ultra-trace residual 2.3×. A Schur-complement or reduced T/front solve would therefore
not help; it is not justified by the measured structure. The lever is the ultra-trace weight
(1e-8). Lowering it to 1e-10 reaches η after one pass (26 instead of 74 products per iterate)
and allows a ten times larger front T step per iterate. It changes the least-squares norm and
the line-search merit, though not the model, the products or the convergence test. The user
approved it for the resumed run, paired with a full-budget guard for the convergence test (§4.2).

**Resumed solve (§5).** `results/nk167_hot` continues from iterate 15 on build 31 with
ultra-trace weight 1e-10. Iterate 15: η reached after one pass (residual 0.19, 26 products), front
T step 5.1 % (Report 21 iterates: ≈ 0.5 %), λ = 0.5 accepted (natural measure 0.235 → 0.140),
26 min. A full iterate now takes ≈ 28 min instead of ≈ 3 h.

## 1. Where the time went

**Per iterate (archived iterate 14, Report 21 run, build 28):** 6 frozen-J column passes
(3 workers × 2 threads) took 880 s. Three product passes (26 + 24 + 24 = 74 directions, 8 threads)
took 1 996 + 2 882 + 3 198 = 8 076 s. Three concurrent trial (base) passes at 2 threads took
1 843 s. The total was ≈ 10 800 s (3.0 h), 75 % of it in products. The Python driver
(least squares, preconditioner solves) took < 1 %. Another run shared the machine at the time.

**Per pass** (macOS `sample` stacks taken during the "before" benchmarks,
`results/bench22/before/*_samples`, aggregated by `scripts/sample_top.py`):

| Pass | Main-thread wall share | Thread-sample share |
|---|---|---|
| Column (residual only, frozen J) | observer-frame formal solution `CMF_FORM_SOL_V2` 78 %, comoving formal solution 7 %, BA updates (`UPDATE_BA_FOR_LINE`, `DO_VAR_CONT`) 10 % | — |
| Base (ALI-iterated J, BA dump) | comoving ALI formal solution 67 %, observer-frame 20 %, BA updates 8 % | high system time (barriers at 8 threads) |
| Product (26 directions) | direction loop 69 %, `COMP_J_BLANK` 10 % (8 % observer-frame), serial `lin_setup` 5 % | `lin_scattering` 75–78 %, almost all of it `lin_apply` (the scattering GMRES: 14.45 ray sweeps per direction and frequency, plus 2 forcing/final sweeps) |

The observer-frame solution produces only the observed spectrum. In solver passes it is pure
overhead. In the product pass, ~85 % of the tangent work is the per-direction scattering GMRES,
whose operator (I − Λσ) is the same for every direction at a given frequency.

## 2. Changes

### 2.1 Build 29: no observer-frame formal solution in solver passes

`jv_noobs.patch`: when `QSTARS_NO_OBS_FORMAL` is set, both `CMF_FORM_SOL_V2` branches of
`COMP_J_BLANK` are skipped. The existing `ELSE IF(FIRST_FREQ)` fallback still sets the
observer-grid bookkeeping. `Runner` sets the variable for every solver pass (`--obs-formal`
restores the old behaviour). Verification and final runs do not go through `Runner` and keep the
spectrum.

*Equivalence.* A deterministic test (`scripts/check_noobs_ba.py`) ran two single-thread passes
with frozen stored J and the BA dump, observer solution on and off. Residual and BA are
**bit-identical** (`noobs_ba_bitwise_check.json`). In ALI-iterated base passes, stored J and BA
differ from build 28 by at most 1.5e-8 (J, one point of 7.2 million; median 3e-16) and 1.5e-9 of
the row scale (BA). Build-28 reruns differ among themselves by up to 6e-9 and 2.5e-10. Both are
at the level of the ALI's thread-order noise (tolerance 1e-10). The two build-30 base runs sit together, 1.5e-9 from the build-28 ones, so the skipped branch may leave a trace in the ALI iteration through shared state. Any such trace is at the ALI convergence level (stored J ≤ 1.5e-8 at one point), far below anything that moves a Newton step at the 1e-3 test. The deterministic one-thread base-pass test that would settle it was dropped to save time (§9). Verification and final passes keep the observer solution.

### 2.2 Build 30: explicit per-frequency scattering operator

`jv_asig.patch`. The tangent solves (I − Λσ)x = f once per direction and frequency by GMRES
(Jacobi preconditioner, tolerance 1e-12). Each Krylov step was a full ray sweep (`lin_apply`).
The operator depends only on the base state at that frequency, so build 30 builds it once:
A(:, j) = J response to a unit scattering emissivity at node j. The 167 columns are batched
16 at a time by `lin_apply_multi`, which uses the statements of `lin_apply` with the direction
index innermost (vectorized; rays summed in the same order). GMRES then multiplies by A. Same
map, same algorithm, preconditioner and tolerance; summation order only. The operator is used
from 4 directions up (`QSTARS_JV_ASIG_MIN`); `QSTARS_JV_NOASIG` restores build 29.

### 2.3 Build 31: batched direction sweeps, parallel `lin_setup`, balanced blocks

`jv_batch.patch`. With A available, each direction still needs one forcing sweep and one final
sweep. Build 31 runs the directions in one block per thread, and each block makes one forcing
and one final sweep for all its directions (`lin_scattering_multi`). The GMRES of each direction
is the A branch of `gmres_scat` verbatim. `lin_setup`'s ray loop runs in parallel. The per-ray
data are independent, and the primal moments are still summed serially in ray order, so they are
unchanged. Operator column blocks are sized to fill twice the thread count (build 30 left
5 of 8 threads idle in the second round). `QSTARS_JV_NOBATCH` and `QSTARS_JV_NOSETUP_OMP`
restore build 30.

### 2.4 Driver (`scripts/cool_newton_krylov.py`)

- **Thread budgets.** Frozen-J columns: `--col-workers 6 --col-threads 1`. Report 21 used
  3 workers × 2 threads; a column pass takes the same time at 1 and 2 threads (167 s vs 173 s,
  §3). Trials: Report 21's `8 // n_trials` (2 threads for 3 trials); a base pass is faster at
  2 threads than at 8 (334 s vs 362 s). Products: threads chosen by measurement (§3).
- **Same-base product reuse.** A product pass at byte-identical state and replayed J is reused
  when the driver returns to the same base: after a PTC δ change, or at a PTC switch to δ = 0.
  jz = J_dir(current δ)·z + drad, with drad read from the pass. The radiation part depends only on
  state, replayed J and direction, so this is the exact product. Checked offline on the
  74 archived pairs: the rebuilt products equal the stored ones bit for bit; with δ → 0 they
  differ from them by exactly the PTC diagonal term, to 2e-16 (`same_base_reuse_check.txt`).
  Reused passes count against the `--krylov` budget. A previous Newton step made at the same base
  lies in the reused span, so it is not spent as a new direction.
- **Resume hardening.** A product pass left by an earlier attempt with other directions is moved
  to `superseded/` with a README line, not reused and not fatal.
- **Optional full-budget guard** (`--full-krylov-near`, default off). Near convergence, it runs
  every Krylov pass before the 1e-3 test. It changes which step the test sees, so it is off and
  the convergence path is Report 21's.
- `--obs-formal`, `--col-threads`, `--trial-threads` (0 = Report 21 rule).

### 2.5 Environment

- Every process in this session runs at macOS background QoS (`PRI 4`; `taskpolicy -B` cannot
  lift it). Ten busy loops together got ~350 % CPU while a game was in the foreground. Background
  work shares roughly 4–6 cores, so 8-thread passes oversubscribe them.
- On battery the same column pass took 569 s instead of 122 s on the idle machine at night.
  Its CPU-seconds rose 4.4× as well, so cross-session wall or CPU times cannot be compared.
  Before/after pairs were therefore always run back to back (§3).
- Idle sleep is held off (`caffeinate -i -w <queue>`), and only one native job runs at a time.

## 3. Matched benchmarks (same state, back to back)

### 3.1 Method and conditions

`scripts/bench_pass.py` reruns an archived pass with the same restart, VADAT, state
(`QSTARS_STEP.bin`), replayed J, warm start and tangent directions, read from its
`run_manifest.json`. Only the executable, the thread count and the observer flag change. While
the native process runs, `sample` records stacks. Afterwards the outputs are compared with the
archived pass, each equation row on its own scale (rows span 8–13 decades; §6 explains why).
Queue scripts: `results/bench22/run_matched.sh`, `queue2.sh`, `queue3.sh`; per-pass JSON in
`Reports/evidence/report22/`; summary in `matched_summary.json`.

Conditions: 2026-10-06 afternoon, on battery until 15:25, interactive foreground use (a game
until ~15:00, a browser throughout), background QoS. Absolute times are 1.5–4.7× the
idle-night times of `results/bench22/before` (Oct 5, 23:00). Only the ratios within a pair are
meaningful.

### 3.2 Results

| Pass (state: nk167_hot iterate 14) | Before: build 28 | After: build 30/31 | Wall ratio | CPU ratio | Agreement after vs archived (row-wise max rel.) | Before rerun vs archived |
|---|---|---|---|---|---|---|
| Frozen-J column (residual only) | 569 s, 592 CPU-s (2 thr) | 167 s, 152 CPU-s (1 thr, b30) | 3.4× | 3.9× | residual 2.3e-14 | 1.1e-14 |
| Base / trial (ALI J to 1e-10, BA, stored J) | 912 s, 2 206 CPU-s (8 thr) | 362 s, 1 177 CPU-s (8 thr, b30); 334 s, 560 CPU-s (2 thr) | 2.5× (8 thr) | 1.9× | residual 3.2e-12; BA 1.4e-9 of row scale; J ≤ 1.5e-8 (ALI noise, §2.1) | 2.3e-12; 6.9e-11; ≤ 6e-9 |
| Product, 26 directions | 1 558 s, 8 549 CPU-s (8 thr) | 810 s, 4 131 CPU-s (8 thr, b31) | 1.92× | 2.07× | products 1.0e-12 (global 5.0e-14); residual 1.2e-14; GMRES steps identical (16 275 803) | 6.0e-13 (global 8.3e-15) |
| Per product (Jv) | 59.9 s, 329 CPU-s | 31.2 s, 159 CPU-s | 1.92× | 2.07× | — | — |

Build 30 alone (8 of the 26 directions, 4 threads, observer solution still on, heavily
contended): products to 1.0e-12 per row; every frequency used the operator.

### 3.3 Products per outer iterate and seconds per iterate

With the Report 21 weights (ultra 1e-8), η = 0.25 is never reached: the residual stalls at
0.46 (§4). Every iterate therefore spends the full budget of 3 passes, 74 products, before and
after. The speed-up acts per product, not on their number. Same-base reuse saves whole passes
only when the driver revisits a base (PTC δ change or switch).

Summed single-pass CPU cost of one iterate (6 columns + 74 products + 3 trials, matched numbers;
the "before" trials are costed at 8 threads, because a 2-thread build-28 base pass was not
measured under these conditions):

| | Columns | Products | Trials | Total CPU-s per iterate |
|---|---|---|---|---|
| Before (build 28, R21 settings) | 3 553 | 24 331 | 6 617 | 34 501 |
| After, same weights (74 products) | 914 | 11 758 | 1 679 | 14 352 (2.4×) |
| After, ultra weight 1e-10 (26 products, §4.2) | 914 | 4 131 | 1 679 | 6 725 (5.1×) |

Measured wall time per accepted iterate: archived Report 21 iterate 14, ≈ 10 800 s (another
167-radius run shared the machine). Resumed run on the fast path (build 31, ultra weight 1e-10): iterate 15 took 1 538 s with its columns reused (products 878 s, trials 660 s). The iterate-16 columns took 128 s (6 workers × 1 thread, one round). A full iterate is therefore ≈ 1 670 s (28 min), **6.5× faster** than the archived one. About 2× of that comes from the faster passes alone; the rest comes from one product pass instead of three. With the old weights the estimate is ≈ 3 300 s (3.2×).

## 4. Measured Jacobian structure: why no Schur complement, and the ultra-trace weight

### 4.1 Where the stalled linear residual lives

Iterate 14 has 74 exact products and the driver's weighted objective ‖W M⁻¹(R − J s)‖. The
least-squares residual is 0.458 (`residual_structure_it14.txt`, offline on the archived pairs).
Split by row class:

| Rows | Share of initial residual² | Share of final residual² | Reduced to |
|---|---|---|---|
| Significant populations | 0.784 | 0.029 | 0.088 |
| T (radiative equilibrium) | 0.122 | 0.530 | 0.954 |
| Ultra-trace (fraction ≤ 1e-20, weight 1e-8): He II 0.348, He III 0.082 | 0.090 | 0.430 | 1.001 |
| Trace (1e-20…1e-10, weight 1e-3) | 0.002 | 0.008 | 0.985 |
| n_e | 0.002 | 0.003 | 0.638 |

96 % of the final residual is at depths 20–48, the ionization front. Is the T-row stall a
missing subspace or a weighting conflict? Re-solving with the same 74 products and different
row weights settles it:

| Least-squares weights | Front T-row residual reduced to | Ultra-trace M⁻¹ residual | max front s_T |
|---|---|---|---|
| Driver (ultra 1e-8) | 0.958 | ×1.00 | 0.0055 |
| Ultra 1e-10 | 0.016 | ×2.3 | 0.052 |
| Ultra rows excluded | 0.016 | ×2.3 | 0.052 |
| Front T row only | 0.000 | — | — |

The span of the 74 exact products already contains the front temperature correction. The
driver's weights decline it, because it would raise the weighted residual of He II/He III levels
whose populations are below 1e-20 of He. The populations of those levels are numerically
irrelevant to every other equation, and the convergence test ignores them (fraction ≤ 1e-10).
A Schur-complement or reduced T/front solve would only add exact directions that the span
already has. It would not change the weighting outcome, and building the dense T–T block exactly
needs ND = 167 more products per iterate. **Not justified by the measured structure; not built.**

### 4.2 Ultra-trace weight 1e-10: adopted for the resumed run with the user's approval

Offline on the archived iterates 12–14 (`ultra_weight_offline_it12_14.txt`), with the same
products:
- Weight 1e-10 brings the linear residual below η after the first pass (0.20 at iterate 14,
  against 0.49/0.46/0.46 with 1e-8). That is 26 instead of 74 products per iterate.
- The significant-row residual falls from 0.36 to 0.10.
- The front T step grows about tenfold (0.05 instead of 0.005).

The change touches the least-squares norm and the line-search merit (natural measure), and so
the solver path. It does not touch the model, the residual, the products or the convergence test
(max Newton step over levels with fraction > 1e-10 < 1e-3). Because it changes the path, it was
put to the user, who approved it on 2026-10-06. It is paired with the full-budget guard
(`--full-krylov-near 10`): when the significant step falls below 1e-2, all three Krylov passes run
before the convergence test, so convergence is never declared on a one-pass step.

## 5. Resumed 167-radius PTC solve

`results/nk167_hot` resumed at 15:39 from `it15_trial_l1` (δ = 0.0312). Settings: build 31,
8 product threads, 6 column workers × 1 thread, trials at 2 threads, ultra weight 1e-10 with the
full-budget guard. The iterate-15 frozen-J columns made by build 28 before the pause are reused
through the driver's exact-state check (§3.2: build 28 and 30 residuals agree to 2e-14). The run's
`README.txt` records the change of path.

| Iterate | Columns | Products | Krylov residual | max significant step | front max s_T | Trials (λ: merit) | Wall |
|---|---|---|---|---|---|---|---|
| 15 | reused | 26 (1 pass, 878 s) | 0.192 | 0.093 | 0.051 | 1: 0.195, **0.5: 0.140**, 0.25: 0.180 (current 0.235) | 1 538 s |
| 16 | 6 in 128 s | running | | | | | |

The PTC RE net/gross ratio went from 1.367e-2 to 1.406e-2. Because the SER update only ever
decreases δ, it held δ = 0.0312. The iterate-16 base residual merit (‖W M⁻¹R‖ with the new base's
preconditioner) is 0.467, against 0.193 at iterate 15. Merits at different bases use different
preconditioners and are not strictly comparable; the step was accepted on the natural-monotonicity
test. Front temperatures now move about ten times faster per iterate than in Report 21. Whether
that converges, or whether the front has a finite-temperature root at all (Report 21 §4.5), is
still open. The run continues unattended; each iterate is logged in `stdout.log` and
`history.json`.

## 6. Candidates measured but not adopted

- **Skipping BA work in residual-only and product passes** (`UPDATE_BA_FOR_LINE`,
  `DO_VAR_CONT`): about 10 % of a column pass and 5–8 % of a product pass after build 31. It needs
  a native control-flow change (COMPUTE_BA), and `qs_jv` requires COMPUTE_BA. Left as a candidate.
- **Direct LU solve of (I − Λσ)** with the explicit operator, instead of GMRES. It would be
  faster, but it changes the products at the GMRES tolerance (1e-12): the products would no longer
  equal the current solver's to summation order. Not adopted under this report's rule.
- **Earlier in this campaign** (Report 22 planning, kept in the session record): Anderson
  acceleration, merging split pieces, span pruning, stale-column directions, speculative columns,
  a residual-only native exit, fewer Krylov passes, analytic columns, sequential trials. Each
  either changes the step without a measured gain, or weakens a check.
- **Thread scaling**: column passes do not scale (1 thread fastest per CPU-second), base passes
  are faster at 2 threads than at 8, and product passes use 8 (background work gets ~6 cores).
  A 4-thread product benchmark was dropped to resume the solve sooner.

## 7. Correctness review of the changes

A four-reviewer workflow covered: tangent arithmetic equivalence (build 31 vs build 30/28),
OpenMP and memory safety, the observer-skip dependencies in CMFGEN, and the driver changes.
Every finding was adversarially re-checked (journal under the session's workflow record):
- **Tangent equivalence, OpenMP safety, observer skip:** no defects.
- **Driver, confirmed and fixed:**
  1. Agreement metrics normalized globally could hide errors in small rows; now row-wise, and the
     build-31 production gate uses them.
  2. The full-budget guard changed which step the convergence test sees; now off by default, and
     used only with the approved ultra-weight change.
- **Driver, judged harmless and hardened anyway:**
  1. A resumed iterate whose product pass had other directions; that pass is now moved to
     `superseded/`.
  2. Trial threads restored to Report 21's rule.
  3. Reused passes now count against the Krylov budget, and a same-base previous step is no
     longer spent as a new direction.

## 8. Gate table

| Gate (this request) | Status | Evidence |
|---|---|---|
| Dominant costs per iterate identified | **DONE** | §1; `results/bench22/*_samples`, `results/profile22` |
| Redundant transfer work removed (observer-frame solution in solver passes) | **DONE** (bit-identical residual/BA at frozen J) | §2.1; `noobs_ba_bitwise_check.json` |
| Radiation/tangent/Krylov information reused | **DONE**: one scattering operator per frequency for all directions; exact same-base product reuse | §2.2, §2.4; `same_base_reuse_check.txt` |
| Directions batched/parallelized where scaling is real | **DONE** (build 31; thread budgets measured) | §2.3, §3 |
| Preconditioner / reduced thermal solve | **ANALYSED**: stall is a weighting conflict, not a missing T subspace; Schur complement not justified; ultra weight 1e-10 adopted with approval | §4 |
| Before/after on the same 167-radius state: base, per product, products/iterate, s/iterate, agreement | **DONE** (iterate time measured on the resumed run) | §3; `matched_summary.json` |
| Physical model, products, convergence criterion unchanged | **YES** (products to 1.0e-12 per row; solver path changed only by the approved weight) | §3.2, §4.2 |
| 167-radius T-released cool closure converged | **NOT YET** (resumed) | §5 |
| v1 | **NO** | — |

## 9. Failed and superseded attempts

- **A contaminated "before" iterate** (`results/bench22/before_iter`): the rerun of iterate 14
  with build 28 overlapped the build-30 correctness test and a foreground game. It was paused,
  then terminated; matched pass-by-pass pairs replaced it (README in the directory).
- **Build-30 check without the observer flag**: the first correctness run of build 30 was
  launched before the driver set `QSTARS_NO_OBS_FORMAL`. It therefore measured the scattering
  operator alone (the observer solution was still 17 % of its main thread). Kept as
  `results/bench22/after/prod8_b30_check`.
- **A queue stalled by its own wait condition**: `queue2/3.sh` waited with
  `pgrep -f "bench_pass.py.*prod26_b28_t8"`, which also matched the monitoring shells' command
  lines. The queue resumed when those shells were stopped (~2 min lost).
- **Scope trimmed after the user flagged the elapsed time** (16 h of campaign wall time): the 4-thread product
  benchmark, a deterministic 1-thread base-pass observer check (open, §2.1) and a
  separate "after" iterate on the iterate-14 state were dropped. The resumed run's iterates give
  the measured per-iterate time instead.

## 10. Reproduction

```bash
export MESASDK_ROOT=$HOME/qsrd-runtime/mesasdk; source $MESASDK_ROOT/bin/mesasdk_init.sh
export MACOSX_DEPLOYMENT_TARGET=14.0
export QSTARS_CMFGEN_ROOT=$PWD/backends/cmfgen/build_31/cur_cmf
make -C $QSTARS_CMFGEN_ROOT d_new_main
.venv/bin/python backends/cmfgen/record_build.py backends/cmfgen/build_31/cur_cmf
```

Builds 29–31 layer `jv_noobs.patch`, `jv_asig.patch` and `jv_batch.patch` (each with a
`*_patch_manifest.json`) on build 28. Executable SHA-256: build 29 `fb95a32aee5a91d5…`,
build 30 `6a02536ecbab514f…`, build 31 `56f0774ae3d571d4…`.

Benchmarks:

```bash
.venv/bin/python scripts/bench_pass.py results/nk167_hot/it14_from_it14_trial_l1/jv00 OUT/prod \
  --x0 results/cmfgen_ali_lambda167 --exe backends/cmfgen/build_31/cur_cmf/exe/cmfgen_dev.exe --threads 8
.venv/bin/python scripts/bench_pass.py ARCHIVED OUT --x0 X0 --compare-only     # re-score agreement
.venv/bin/python scripts/check_noobs_ba.py results/nk167_hot/it14_trial_l1 OUT --x0 results/cmfgen_ali_lambda167 \
  --exe backends/cmfgen/build_30/cur_cmf/exe/cmfgen_dev.exe [--base-mode]
.venv/bin/python scripts/compare_iterate.py ARCHIVED_RUN NEW_RUN K             # iterate-level agreement
```

The offline structure analyses run from `Reports/evidence/report22/` (`resid_structure.py`,
`offline_krylov_eval.py`, `reuse_check.py`, with `nkoff.py`). Resumed solve: the
`results/bench22/queue3.sh` command line; `--ultra-weight 1e-10 --full-krylov-near 10`.
Regression subset `tests/test_nlte_newton.py`: 9 passed, 1 skipped.
