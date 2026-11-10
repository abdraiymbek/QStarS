# Report 29 — continuation across the front nonlinearity (from 2026-11-09)

**Status: QStarS is NOT v1.** The 167-radius cool closure is not converged. Live record of the continuation block
chosen by the user (option 1 of Report 28 §6, with a cheap piece of option 3 first).

## 1. Option 3, cheap piece: what the convergence scope would change (analysis only)

`evidence/report29/scope_check.py` at it35 (unbounded Newton step of that iterate, frozen it16 measure):

| threshold | levels below | their share of the residual² | radiatively coupled among them | remaining levels | remaining: median / 90 % / max step | remaining < 1e-3 |
|---|---|---|---|---|---|---|
| 1e-8 | 5 575 | 41.9 % | 805 | 4 236 | 6.4e-3 / 4.6e-2 / 98 | 23 % |
| 1e-7 | 7 914 | 55.6 % | 1 033 | 1 897 | 6.2e-3 / 7.9e-2 / 98 | 25 % |
| 1e-6 | 8 899 | 56.9 % | 1 057 | 912 | 2.0e-3 / 0.38 / 98 | 45 % |

"Radiatively coupled": the fixed-J and local-ALI diagonals differ by more than ×2, or H I n = 2 — such levels cannot be
handed to a fixed-J cleanup (Report 27 §5). **Conclusion:** moving low-population levels out of the Newton test does
not bring convergence within reach: even the 912 major variables need corrections up to 98, only 45 % meet the test,
and ~1 000 of the low-population levels are radiatively coupled anyway. No scope change is made; the decision stays the
user's and is not needed now.

## 2. What continuation has to cross

The unbounded Newton step of iterates 33–35 concentrates at **depth 33** (inner H front): H I n = 2 175–235, T 95–102,
H II and n_e 37–56, plus the near-trace He II ground at d22–25 (213 → 1 580). This is the nearly singular front mode of
Report 25 (H I n = 2 held by trapped Lyα, coupled to the ionization balance and T), now measured with correct algebra;
the λ² curvature of the trials sits in the same equations (H ionization balance at d19–47; Report 26 §6b).

## 3. Continuation: the T homotopy

`--t-homotopy --hom-mode pin`: RE row λ R_RE + (1 − λ) c (T − T_fix)/T_fix, T_fix = the start state's T. Stage λ = 0
holds T exactly (`--freeze-t`): the population problem at fixed T, without the T ↔ ionization coupling of the soft
mode. Later stages step λ → 1 with T released; new option `--hom-re-merit` makes the blended RE row count in the least
squares and the acceptance (with a homotopy it was left out, right only at fixed T). Run `results/nk167_trap/v1_p3`,
start `v1_p2/it36_trial_l0.25`, solver as Reports 27–28 (no dense columns), λ trials 1 / 0.5 / 0.25 / 0.125.

**Test first (two iterates of stage 0, ~3.5 h):** does the population problem progress at fixed T, and does the d33
mode shrink? If not, T is not the right continuation handle and the next candidate is population pseudo-time
continuation (`--pop-ptc`, built in).

## 4. Results

(appended as the block runs)

**Stage λ = 0, iterate 36 (T held):** the unbounded Newton step's largest components fall from H I 175 / H II 37 /
n_e 37 / T 98 (free T, it35) to **H I 0.75 (d16), H II and n_e 0.089 (d28), He I 0.155**; only the near-trace He II
ground at d25 remains large (51.8). The 0.25-box step equals the unbounded one (linear fit 0.918 / 0.918): the bound no
longer binds — an ordinary Newton regime. Trials λ 1 / 0.5 / 0.25 / 0.125: 0.945 / **0.945** / 0.967 / 0.982 of the stage
merit; accepted λ 0.5: 1.4937e-2 → 1.4109e-2 (−5.5 %), 2.4 h. **The d33 soft mode is the T ↔ ionization coupling:
with T held it disappears.** T is the right continuation handle.

**Stage λ = 0, iterate 37:** same picture (unbounded step: H I 0.88 at d16 row 2, H II/n_e 0.086, He I 0.16; He II d25
trace 51); 0.25-box fit 0.914 against 0.907 unbounded. Trials 0.930 / 0.944 / 0.967 / 0.982; **full step accepted**
(λ 1, ρ 0.82): stage merit 1.3971e-2 → 1.2991e-2 (−7.0 %), 1.95 h → 3.6 %/h (the free-T runs of Report 28: 0.2–0.6 %/h).

The real problem's measure (frozen it16 norm, RE row included, T unchanged) improves with it:

| state | frozen measure | vs it16 (2.837e-2) | RE net/gross max | conservation |
|---|---|---|---|---|
| v1_p2 it36 base (Report 28 end) | 1.5783e-2 | −44.4 % | 3.45e-3 | 3e-16 |
| after stage-0 iterate 36 | 1.5210e-2 | −46.4 % | 3.51e-3 | 3e-16 |
| after stage-0 iterate 37 | **1.4584e-2** | **−48.6 %** | 3.58e-3 | 3e-16 |

Groups: H I n = 1 / ionization balance 0.047 → 0.023, H I n ≥ 3 0.42 → 0.29; H I n = 2 0.128 → 0.201 and He I
0.277 → 0.316 grow (their equations are coupled to T, which is still held).

**Continuous run from iterate 38 (2026-11-09 18:03).** New driver option `--freeze-t-lam0`: T is held only while
λ = 0 and released for λ > 0, so a single run carries the fixed-T stage into the release stages. Setting: as above
without `--freeze-t`, plus `--freeze-t-lam0 --hom-re-merit --stage-max-iters 4 --max-outer 80`. Stage 0 gets two more
iterates (4 in total; the He II d25 trace keeps its step above the stage tolerance, so the cap ends it); then λ = 0.25,
0.5, 0.75 with T free, 4 corrector iterates each (earlier if the stage tolerance 1e-2 is met; the driver halves the λ step
and retreats if no trial decreases the merit); then λ = 1 (the real equations) to the convergence tolerance. Automatic
stops as before (< 0.5 %/h twice at λ = 1, all trials rejected, RE net/gross > 5e-3, the soft mode returning: significant
step > 16 with < 2 % gain).
