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

| iterate | λ | T | linear fit (8 passes) unbounded / box | accepted | trials λ 1 / 0.5 / 0.25 / 0.125 | stage merit | drop | ρ | frozen measure | wall |
|---|---|---|---|---|---|---|---|---|---|---|
| 36 | 0 | held | 0.918 / 0.918 | λ 0.5 | 0.945 / 0.945 / 0.967 / 0.982 | 1.4937e-2 → 1.4109e-2 | −5.5 % | — | 1.5210e-2 | 2.4 h |
| 37 | 0 | held | 0.907 / 0.914 | λ 1 | 0.930 / 0.944 / 0.967 / 0.982 | 1.3971e-2 → 1.2991e-2 | −7.0 % | 0.82 | 1.4584e-2 | 1.95 h |
| 38 | 0 | held | 0.954 / 0.961 | λ 1 | 0.968 / 0.973 / 0.983 / 0.990 | 1.2888e-2 → 1.2475e-2 | −3.2 % | 0.83 | 1.4524e-2 | 2.6 h |
| 39 | 0 | held | 0.970 / 0.981 | λ 0.5 | 0.996 / 0.989 / 0.993 / 0.996 | 1.2477e-2 → 1.2345e-2 | −1.1 % | 0.79 | 1.4444e-2 | 3.3 h |
| 40 | 0.25 | free | 0.965 / 0.982 | λ 0.5 | 1.000 / 0.990 / 0.993 / 0.996 | 1.2379e-2 → 1.2260e-2 | −1.0 % | 0.72 | 1.4358e-2 | 5.9 h* |
| 41 | 0.25 | free | 0.965 / 0.988 | λ 1 | 0.989 / 0.992 / 0.995 / 0.997 | 1.2288e-2 → 1.2148e-2 | −1.1 % | 0.92 | 1.4267e-2 | 2.6 h |
| 42 | 0.5 | free | 0.972 / 0.992 | λ 1 | 0.993 / 0.995 / 0.997 / 0.998 | 1.2229e-2 → 1.2140e-2 | −0.7 % | 0.93 | 1.4226e-2 | 2.2 h |
| 43 | 0.5 | free | 0.972 / 0.994 | λ 1 | 0.995 / 0.996 / 0.998 / 0.999 | 1.2193e-2 → 1.2127e-2 | −0.5 % | 0.88 | 1.4172e-2 | 3.5 h |
| 44 | 0.75 | free | 0.971 / 0.993 | λ 1 | 0.994 / 0.995 / 0.997 / 0.998 | 1.2213e-2 → 1.2145e-2 | −0.6 % | 0.76 | 1.4133e-2 | 1.85 h |
| 45 | 0.75 | free | 0.969 / 0.991 | λ 0.5 | 0.994 / 0.994 / 0.997 / 0.998 | 1.2149e-2 → 1.2081e-2 | −0.6 % | 0.87 | 1.4069e-2 | 1.9 h |
| 46 | **1** | free | 0.962 / 0.984 | λ 0.5 | 1.047 / 0.994 / 0.994 / 0.997 | 1.2310e-2 → 1.2239e-2 | −0.6 % | 0.51 | 1.3975e-2 | 1.7 h |
| 47 | 1 | free | 0.964 / 0.987 | λ 0.5 | 1.024 / 0.996 / 0.996 / 0.998 | 1.2249e-2 → 1.2199e-2 | −0.4 % | 0.47 | 1.3932e-2 | 2.25 h |

Iterate 38: the linear model stays accurate (ρ 0.83) but the 8-pass span reaches less (0.954); the unbounded step's
largest major component is H I row 2 at d19 (3.96), He II d25 trace 61.8. The real-problem measure gains only −0.4 %:
with T held, the RE part (0.121 of the frozen measure) and the T-coupled H I n = 2 (0.224) cannot improve. Stage 0 is
near what it can give; iterate 39 is its last (cap 4), then T is released at λ = 0.25.

Iterate 39 (last of stage 0): −1.1 % in 3.3 h (0.33 %/h); frozen measure 1.4444e-2 (**−49.1 %** since it16). Stage 0
has given what it can (−17 % of the stage merit, −8.5 % of the frozen measure in four iterates); the cap now releases T
at λ = 0.25 (iterate 40).

**Iterate 40 — first stage with T free (λ = 0.25).** The soft mode does not return: the unbounded step's largest
major components are H I row 2 at d19 (4.6), He I 0.34, H II / n_e 0.077, **T 3.4e-4** (it35 with T free: H I 175, T 98);
only the He II d26 trace stays large (42). Accepted λ 0.5 (ρ 0.72): −1.0 % of the λ = 0.25 merit; frozen measure
1.4358e-2 (**−49.4 %**), RE net/gross 3.64e-3. T moves very little at λ = 0.25 (the pin keeps 3/4 of the RE row's
weight), as intended: T follows the populations in small steps. \*Wall time inflated: 00:30–05:00 another workload on
the machine (load average 56–76 on 10 cores) slowed the passes ×5 (75 min instead of 13–15); progress rates for this
iterate are judged per pass, not per hour.

Iterate 41 (λ = 0.25): −1.1 % (full step, ρ 0.92), frozen measure 1.4267e-2 (**−49.7 %**); the He II d25 trace's step
fell to 5.8 (from 42–73); largest major step H I row 2 at d19 (5.7), T 1.7e-4. **Change (2026-11-10 08:50):** the
corrector iterates at intermediate λ gain ~1 % each while T moves by < 1e-3, so the stages are shortened from 4 to 2
iterates (`--stage-max-iters 2`; restart at iterate 42, completed passes reused, interrupted ones moved to
`_interrupted_it42/`). λ = 1 is now reached at iterate 46 instead of 52; the question the block has to answer is whether
the Newton step there stays in the ordinary regime.

Iterate 42 (λ = 0.5): no soft mode (largest major step H I row 2 d19 3.8, T 6.6e-4 at d33, He II d26 trace 32); full step
accepted, −0.7 %, frozen measure 1.4226e-2 (**−49.9 %**). The 0.25 box binds (box fit 0.992 against 0.972 unbounded):
the H ionization balance at d19 asks for corrections of ×4.

Iterate 43 (λ = 0.5): −0.5 %, frozen measure 1.4172e-2 (**−50.0 %**); T step 2e-4, H I row 2 d19 4.5, He II d26 trace 21.
The soft mode stays away, but the per-iterate gain shrinks (−1.1 → −0.7 → −0.5 %): the 8-pass span reaches only 0.972
even unbounded, and the box (0.994) leaves most of that. The limit is again the linear solve (span/preconditioner), now
without the soft mode in the way.

Iterate 44 (λ = 0.75): −0.6 %, frozen 1.4133e-2 (**−50.2 %**); T step 2.0e-3 at d33 (the front depth starts to ask for T
again as the pin weakens, still ×5·10⁴ below it35), H I row 2 d19 3.2.

Iterate 45 (λ = 0.75): −0.6 %, frozen 1.4069e-2 (**−50.4 %**), T step 8e-4. Iterate 46 starts the real equations
(λ = 1, T free, convergence test 1e-3).

**Iterate 46 — the real equations (λ = 1).** The front mode reappears, much weaker: the unbounded step asks T 1.10 and
H II / n_e 0.34 at **d33** (it35, same equations before the homotopy: T 98, H I n = 2 175, H II 37), H I row 2 d19 4.4,
He II d25 trace 35. The full step overshoots again (×1.047; it33–35: ×1.14–1.26), λ 0.5 accepted (ρ 0.51): −0.6 %,
frozen measure **1.3975e-2 (−50.7 % since it16; −11.5 % since the homotopy start)**, 0.34 %/h — the first λ = 1 iterate
below 0.5 %/h (stop rule: twice).

Iterate 47 (λ = 1): T 0.74 / H II, n_e 0.43 at d33 (no growth), full step ×1.024, λ 0.5 accepted: −0.4 %, frozen
**1.3932e-2**, 0.18 %/h.

## 5. Stop (2026-11-10 23:02) and what the continuation achieved

Iterates 46 and 47 ran below 0.5 %/h at λ = 1: stop rule. Best state `results/nk167_trap/v1_p3/it48_trial_l0.5`
(pointer `accepted_state.json`), frozen it16 measure **1.3932e-2**: −50.9 % since it16 (2.837e-2), **−11.7 % in this
block** (1.5783e-2 → 1.3932e-2, iterates 36–47, ~29 h of compute including ~4 h lost to machine contention).
Conservation ~3e-16, RE net/gross max 3.64e-3 (start 3.45e-3). Residual composition (frozen measure share): He I 0.31,
H I n ≥ 3 0.28, H I n = 2 0.23, RE 0.13, H⁻ 0.03, He II 0.02, H I n = 1 / ionization balance 0.008 (it36: 0.047).

**What the continuation did.** Releasing T gradually from a fixed-T solution moved the state off the branch where the
front soft mode dominated: at the same equations, the unbounded step at d33 fell from T 98 / H I n = 2 175 / H II 37
(it35) to **T 0.7–1.1 / H II 0.3–0.4** (it46–47); full steps overshoot by ×1.02–1.05 instead of ×1.14–1.26. The
fixed-T stage itself gave the largest gains (−17 % of its merit in four iterates; the H ionization balance group −83 %).

**What it did not do.** Once T is free the per-iterate gain returns to ~0.4–0.6 %: the 8-pass span reaches only
0.96–0.97 of the residual even unbounded, and the 0.25 box leaves ~0.99. With the soft mode now small, this is the
linear solve's limit (span / preconditioner quality), not the nonlinearity — the opposite of the Report 28 diagnosis
for the old state. Closing diagnostic: §6.
