# P2: bounded comparison — deep global Krylov with the P1 corrections

Run `results/nk167_trap/v1_p2` (launch `results/nk167_trap/launch_v1_p2.sh`), base `free_r25/it16_trial_l0.5`
(frozen 2.8371e-2), stage-2 method of free_r25 **without** the dense window, plus `--struct-active --consistent-dirs
--step-solver box`; iterate 16 with 12 sequential product passes (`--krylov 12 --eta 0.01`), 34 directions per pass.

## Fit by pass count (iterate 16, 34 directions per pass, ~12.5 min per pass)

| passes | directions | raw unbounded (driver log) | raw LM-QR bounded | raw box bounded | **truncated basis (τ 1e-6) box** |
|---|---|---|---|---|---|
| 1 | 34 | 0.992 | 0.998 | 0.997 | |
| 3 | 102 | 0.945 | 0.907 | 0.977 | |
| 5 | 170 | 0.858 | 0.868 | 0.917 | |
| 7 | 238 | 0.772 | | | 0.876 (185 basis vectors) |
| 12 | 408 | 0.661 | | 0.784 (driver) | 0.781 (308 basis vectors) |

(`p2_fit_by_pass.py`, `p2_tsvd.py`, driver.log.)

## The raw least squares fits product noise

The 5-pass LM-QR step predicted 0.868 and natively gave 0.995 (λ 1) / 0.986 (λ 0.5); the gap was exactly linear in λ
and sat in the trapped-line cells (H I n=2 at d33/37/38, He I row 33 at d7/8). The shrinking-step directional
derivative (`p2_dirderiv.py`, native trials at λ = 0.02 and 0.005, one exact product pass of the step and its pieces):

| | ‖FD + J s‖ / merit |
|---|---|
| J s from a fresh exact product of the step | **0.0024** (‖FD‖ 0.298: 0.8 %) |
| J s = Σ y_i J z_i from the stored products | 0.473 |

The exact products and the native derivative agree; the model built from them does not, because the 238 directions span
only 185 independent ones and the least squares returned coefficients up to **3.8e14** whose sum cancels to
‖s‖ = 0.38: the products' own small errors, amplified 1e13–1e14, became an O(1) model. In the orthonormal basis of the
normalized directions truncated at relative singular value τ (`--span-tau`, `cool_newton_krylov.span_basis`), the
predicted fit equals the fit computed from a fresh exact product of each candidate step to four digits for every τ
from 1e-1 to 1e-8 (`p2_tsvd.py --verify`, 12 candidates in one product pass). This also explains the "fit stalls at
0.8–0.9 / predictions not realized" behaviour of Reports 24–26: part of every Krylov fit was noise.

## Native results at it16 (driver trial pipeline; local merit relative to 2.8371e-2)

| step | λ | predicted | **native** | frozen measure |
|---|---|---|---|---|
| 5 passes, raw LM-QR | 1 / 0.5 | 0.868 / 0.904 | 0.995 / 0.986 | 2.823e-2 / 2.798e-2 |
| 7 passes, τ 1e-6 box | 1 / 0.5 / 0.25 | 0.876 / – / – | **0.920** / 0.926 / 0.959 | **2.609e-2** |
| 12 passes, raw box (driver) | 1 / 0.5 / 0.25 / 0.125 | 0.784 (λ 1) | 1.197 / **0.885** / 0.927 / 0.962 | **2.510e-2** (accepted) |
| 12 passes, τ 1e-6 box | 1 / 0.5 / 0.25 | 0.781 (λ 1) | 1.301 / 0.892 / 0.927 | 2.529e-2 |

Conservation after every trial's projection ≤ 4e-6 (λ 1) and ~1e-13 otherwise (consistent directions). The driver
accepted `it17_trial_l0.5`: frozen 2.8371e-2 → **2.5102e-2 (−11.5 %)** in one iterate of ~3 h wall (12 passes 2.5 h,
column/ALI passes and trials ~0.4 h), i.e. ~3.7 %/h; the 7-pass step gave −8 % for ~1.7 h (~4.7 %/h). The dense-window
run gave −3.2 % in total over 70 h of column builds. Full steps at the |s| ≤ 0.5 bound overshoot (nonlinearity at the
fronts); λ 0.5 is the accepted length in both 12-pass variants.

## Then: three corrected iterates (plan rule: review after at most three)

`--span-tau 1e-6 --krylov 8`, same method otherwise, from `it17_trial_l0.5`; stop rules: all trials rejected, or < 1 %
cumulative frozen improvement over three accepted steps. Results are appended below by the session.

| iterate | passes | basis kept | predicted (λ 1) | trials λ 1 / 0.5 / 0.25 / 0.125 (local, relative) | accepted | frozen measure |
|---|---|---|---|---|---|---|
| 16 (raw, box) | 12 | – | 0.784 | 1.197 / **0.885** / 0.927 / 0.962 | λ 0.5 | 2.5102e-2 (−11.5 %) |
| 17 (τ 1e-6, box) | 8 | 230 / 273 | 0.851 | 1.070 / **0.925** / 0.949 / 0.972 | λ 0.5 | 2.3234e-2 (−18.1 % cum.) |
| 18 (τ 1e-6, box) | 8 | 226 / 274 | 0.918 | 1.127 / **0.964** / 0.972 / 0.984 | λ 0.5 | 2.2493e-2 (−20.7 % cum.) |

Iterates 17–18 ran while the user's own jobs loaded the machine (load ~21 on 10 cores at 00:22; the driver was then
reniced +10): passes took 12–19 min instead of 12.5.
| 19 (τ 1e-6, box) | 8 | 221 / 274 | 0.912 | 1.144 / **0.964** / 0.972 / 0.984 | λ 0.5 | 2.1853e-2 (−23.0 % cum.) |

Frozen measures from each accepted state's `_tp` residual (`frozen_scores_v1_p2.json`; it17 2.5074e-2, it18 2.3224e-2,
it19 2.2491e-2), the last from its trial pass. No stop rule fired; the run stopped at the planned review
(`--max-outer 20`). Review in [Report 26](../../../Report_26.md) §6. Resume: pointer `results/nk167_trap/v1_p2/accepted_state.json`
(it20_trial_l0.5).
