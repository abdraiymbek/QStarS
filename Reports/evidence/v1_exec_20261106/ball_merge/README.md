# Ball 2011 smoke gate: the 46 Msun merge event (V1 plan P5, bullet 2)

2026-11-07. Diagnostic only. No registration, plan, amendment, run directory, or file under `src/` or `benchmarks/` was
changed. Nothing is committed. Every experiment is either a monkeypatch inside its own Python process or a scratch
DRYRUN under `results/ball_merge_20261107/` (gitignored). All jobs ran single-threaded at nice 5, at most 4 at a time;
total compute was about 65 CPU-minutes (about 50 min wall).

Code state: the stage-7 code freeze (`Reports/evidence/report23/ball/stage7_code_freeze.sha256`, 26 files including
`scripts/validate_ball_track.py`) verifies with 0 mismatches. The plan sha256 is `a10b4166…`. The ABI 7 library sha256
`f6bfdc28…` matches `run.json` of the preserved smoke run.

## 1. Reproduction (goal 1)

`scripts/repro_rejections.py` loads `benchmarks/results/2026-10-14-ball-track/20261014T051036.115106Z_smoke_C1_n200/checkpoint.npz`
(step 1255, M_BH = 46.0038, 174 zones) with its stored configuration and repeats the eight logged rejected attempts.
All eight reproduce `delta_actual` to the logged precision, and all are `profile_dlnT` rejections (> 2 δ_target = 0.01):

| dt (yr) | 172.14 | 86.07 | 43.04 | 21.52 | 10.76 | 5.38 | 2.69 | 1.34 |
|---|---|---|---|---|---|---|---|---|
| logged | 0.01525 | 0.01691 | 0.01775 | 0.01817 | 0.01839 | 0.0185 | 0.01855 | 0.01857 |
| reproduced | 0.015246 | 0.016913 | 0.017749 | 0.018172 | 0.018387 | 0.018497 | 0.018549 | 0.018569 |

Every attempt merges one pair of cells (174 → 173 zones), and Newton converges in 4–5 iterations every time. The
largest change is always in cell 167, 0.0601 Msun below the surface (`data/repro_rejections.json`, which also lists
the outer 14 cells: κ, T, Δlnκ, Γ).

The whole registered smoke run was also reproduced end to end: `E0_control` (my DRYRUN driver with the refinement
switched off) starts from the registered row 0 (`…045753.248504Z_relax_C1_n200/row0_6000rsun.npz`). Its `history.csv`
and final `checkpoint.npz` are **byte-identical** to the preserved run's. So the driver used in §3 is the registered
code path.

## 2. Error separation at the event (goal 2)

Common coordinate: mass below the surface, with ln T linear in it between cell centres. This is the convention of the
controller's mapped metric. The probe is the centre of cell 167 (0.0601 Msun below the surface). "Max outer" is the
maximum over 400 points in the outer 0.33 Msun. Script: `scripts/error_separation.py`; data:
`data/error_separation.json` (registered Newton tolerance) and `data/error_separation_tol1e-10.json` (part C with
tolerance 1e-10, see below).

### (a) Temporal error: merge suppressed (merge threshold 0.45 in memory, so no event fires)

| dt (yr) | 1 | 3 | 10 | 30 | 100 | 172.14 | 344.3 |
|---|---|---|---|---|---|---|---|
| controller δ_actual | 2.42e-5 | 7.12e-5 | 2.29e-4 | 6.55e-4 | 2.10e-3 | 3.59e-3 | 7.11e-3 |
| BE local error, front (step doubling, 2×(one step − two half steps)) | 1.5e-6 | 4.3e-6 | 1.4e-5 | 3.2e-5 | 5.3e-5 | 6.7e-5 | 9.7e-5 |

- The controller change is proportional to dt at the front: 2.1–2.4e-5 per yr (least-squares slope 2.06e-5 per yr
  for dt ≥ 10 yr). It reaches δ_target = 0.005 near dt ≈ 240 yr, which is the dt the controller was holding before
  the event (median dt 139 yr).
- Below 1 yr the change goes to a floor of about 2e-5. That floor sits at the base (dln r0, dln ρ_b), not at the front.
  At the front the change falls to 4e-7 at dt = 1e-3 yr.
- The actual backward-Euler local error at the front, from step doubling, is 1.45e-6 per yr × dt up to 30 yr. It
  grows sublinearly above that (order reduction on the stiff, quasi-static outer cells) and is 6.7e-5 at the
  rejected dt of 172 yr. That is 1.3 % of δ_target and 1/280 of the merge jump.

### (b) Remap-induced structural error: registered merge vs. the same step without it

J(dt) = ln T(merged) − ln T(no merge), same dt, common coordinates:

| dt (yr) | 1e-8 | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 | 172.14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| J at the probe | 4.4e-5 | 2.25e-3 | 5.31e-3 | 8.53e-3 | 1.286e-2 | 1.682e-2 | 1.835e-2 | 1.859e-2 | 1.863e-2 | 1.875e-2 | 1.883e-2 |
| max outer \|J\| | 1.04e-2 | 9.3e-3 | 8.5e-3 | 9.7e-3 | 1.29e-2 | 1.68e-2 | 1.83e-2 | 1.85e-2 | 1.85e-2 | 1.86e-2 | 1.87e-2 |

- J does not depend on dt across the controller's whole range (dt ≥ t_dyn = 0.83 yr): it stays at 1.86–1.88 %, which
  is 3.7 δ_target, or 1.9× the rejection threshold.
- At dt → 0 the jump does not vanish. At 1e-8 yr the remap itself already moves the merged cell by 1.0 %. The jump then
  spreads to cell 167 over the thermal time of the outer cells (≈ 1e-3 to 1e-2 yr, far below t_dyn).
- So the rejection measures a structural change, not a temporal error. Halving dt cannot reduce it.

### (c) Spatial discretization error at the probe, and the merge jump vs. resolution

No N = 400 or 800 state exists at M_BH = 46: the N = 400 and 800 relaxations are row 0 at M_BH = 5. So the outer 14
cells (160–173, the outer 0.15 Msun) were refined ×f in a one-off remesh, and the state was stepped 1 yr. That is long
compared with the thermal time of the outer cells and short compared with evolution (the temporal change is 2.4e-5),
so the result is the quasi-static structure on that mesh.

| f (outer cells refined) | 1 (174 zones) | 2 (188) | 4 (216) | 8 (272) | extrapolated |
|---|---|---|---|---|---|
| Δ ln T at probe vs. 174 zones | 0 | −1.858 % | −2.408 % | −2.548 % | −2.595 % |

- The observed orders are 1.75 and 1.98, i.e. second order. The extrapolated spatial error of the 174-zone mesh at the
  probe is **e₁₇₄ = +2.59 %**. The merged 173-zone mesh is about +4.45 %.
- **The remap jump is bounded by the spatial error of the current mesh:** J = 1.86 % = 0.72 e₁₇₄. It is the change in
  discretization error caused by coarsening one cell pair where κ changes by a factor 1.8–2.4 per cell.
- The same event at higher resolution (merging the pair of cells that straddles the old 168/169 face on the f-refined
  mesh, from the refined state, with the registered controller metric) shrinks faster than second order:

| f | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| J at probe (dt = 1 yr) | 1.86e-2 | 1.46e-3 | 1.63e-4 | 1.97e-5 |
| max outer \|J\| (dt = 1 yr) | 1.85e-2 | 1.62e-3 | 2.8e-4 | 3.7e-5 |
| registered δ_actual of the merge step, dt = 1 yr | **0.0186 (reject)** | 0.00156 | 1.1e-4 | 5.9e-5 |
| registered δ_actual of the merge step, dt = 172 yr | **0.0153 (reject)** | 0.00215 | 0.00369 | 0.00387 |

- Successive ratios at the probe are 12.7, 9.0 and 8.3, an observed order of 3.7, 3.2 and 3.05. That is third order
  asymptotically, as expected for coarsening a single cell: an O(h²) local truncation error over a cell of mass O(h).
  Second order is the global error, which sums ~1/h such cells.
- With the front refined only ×2, the same merge moves the solution by 0.3 δ_target and is accepted. At f ≥ 4 the merge
  step's δ_actual at 172 yr is the ordinary temporal change.
- **Answer to goal 2:** the merge jump is a resolution artifact, bounded by (0.72×) the spatial error of the current
  mesh. Refining the front shrinks it at third order (≥ the expected second order).

### Side finding: a Newton residual floor appears in a resolved front

On the refined meshes the registered Newton tolerance of 1e-12 cannot be met (`scripts/residual_floor.py`,
`data/residual_floor.json`):

| f | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| final max\|f\| (registered tol 1e-12) | 4.5e-13 ok | 5.5e-12 FAIL | 1.0e-11 FAIL | 2.5e-11 FAIL |
| max\|f\| under 4-ulp perturbations of x | 4.5e-13 | 0.6–2.1e-11 | 3.0–9.3e-11 | 0.9–1.6e-10 |
| last Newton step | — | 3.9e-13 | 3.5e-13 | 2.8e-13 |

- The worst rows are always the transport-closure rows of the front faces. These rows are scaled by
  max(0.25 |Δln P|, 1e-3) (`quasistatic.py:537`).
- Refining drives 0.25 |Δln P| below the 1e-3 floor: on 7 %, 20 % and 36 % of all faces at f = 2, 4 and 8. Inferred, not
  separately measured: the closure's evaluation noise (EOS and opacity roundoff) divided by that smaller scale then
  rises above 1e-12.
- The registered roundoff exit (residual ≤ 1e-10 and step ≤ 1e-13) does not fire, because the steps stall at
  3–4e-13.
- Part (c) above therefore used a Newton tolerance of 1e-10 (diagnostic only). The probe temperatures differ by at most 1.6e-11 in ln T
  between the two tolerances.
- The same floor shows up unrefined near 55 Msun. At the verifier's merge-free step-1667 state the attainable residual
  is 7.5–9.7e-13, with no margin, while Newton iterations rise 8/15/23/40 at dt = 0.25/0.5/1/2 yr.

## 3. Refining the front under the registered controller (goal 3, DRYRUNs)

`scripts/front_refine_dryrun.py` runs the registered smoke track (stored C1 configuration, amendments 3, 4, 2, 1,
registered controller and stop rules, from registered row 0) with one in-memory change, a front-refined P2 target:

- target_k × clip(C / g_k, 1/FMAX, 1), where g_k is the max |Δln κ| over the two faces of cell k of the accepted state,
  smoothed by a running max over ±W cells so that cells ahead of the moving front are refined before it arrives.
- The registered split (2.0) and merge (0.5) hysteresis acts on these targets, so every split and merge remains an
  event measured by the registered controller (no exemption, unlike Amendment 5).
- `quasistatic.targets` is monkeypatched in the DRYRUN process only.

Summaries: `scripts/summarize_dryruns.py` → `data/dryrun_summary.json`; front state at the stops:
`scripts/front_state.py` → `data/front_state.json`.

| Run | Change vs. registered | Stop M_BH | Steps | Event steps: max δ_actual (> δ_target) | Rejections | Stop cause |
|---|---|---|---|---|---|---|
| E0 | none (control) | 46.00 | 1255 | 0.0058 (1) | 8 controller | merge jump; byte-identical to the registered run |
| E1 | front C = 0.15, W = 3, FMAX = 8 | **46.74** | 1402 | 0.0046 (0) | 34 Newton, residual 1.03–2.1e-12 | t_dyn trigger: Newton floor ≈ 1e-12 |
| E4 | front C = 0.3, W = 2, FMAX = 4 | **50.45** | 1647 | 0.0096 (8, all < 2 δ_target) | 49 Newton, 1.01–2.2e-12 | t_dyn trigger: Newton floor |
| E3 | N = 400 globally (C0 resolution, C1 rules) | **48.11** | 1534 | 0.0045 (0) | 66 Newton, 1.02–2.7e-12 | t_dyn trigger: Newton floor |
| E5 | none; Newton tolerance 1e-10 (diagnostic) | 45.98 | 1254 | 0.0058 (1) | 9 controller, δ 0.0143–0.0184 | the same merge jump, one step earlier: the tolerance alone does not help |
| E2 | front C = 0.15, W = 3, FMAX = 8, plus Newton tolerance 1e-10 (diagnostic) | **60.00** | 4488 | 0.0068 (5, all < 2 δ_target) | 298 Newton (293 genuine, residual > 1e-6), 3 enforced screen | t_dyn trigger after enforced-screen rejections (Kn > 0.1 at τ ≥ 10), with dt already ≈ 9 yr from Newton growth |
| E6 | front C = 0.3, W = 2, FMAX = 4, plus Newton tolerance 1e-10 (diagnostic) | **53.03** | 1676 | 0.0088 (11) | 17 Newton, 6 controller δ 0.0125–0.0135 | a refinement split inside an already coarse front (Δln κ 0.5 per cell) jumps 1.35 % at every dt |

What the DRYRUNs show:

1. **Front refinement removes the 46 Msun merge stop under the registered controller.**
   - Every front-refined run passes 46.0 Msun with every split and merge measured.
   - E1 (C = 0.15) had no controller rejection at all; its event steps peak at δ_actual = 0.0046 < δ_target.
   - The global N = 400 control (E3) passes 46 too (event steps ≤ 0.0045). This is consistent with §2(c): the 46 Msun
     stop is an N = 200 front-resolution artifact.
   - The Newton tolerance alone (E5) changes nothing: the stop stays at 45.98 Msun on the same merge jump.
2. **Under the registered Newton tolerance, every run that passes 46 stops soon after on the Newton residual floor.**
   - Stops: 46.74 (E1), 48.11 (E3, N = 400), 50.45 (E4).
   - All their rejections are Newton stalls at max|f| = 1.01–2.7e-12, just above the registered 1e-12.
   - At E1's stop (`data/front_state.json`), repeating the last attempt at dt = 1.31, 0.66, 0.33 and 0.16 yr gives
     final residuals of 1.26, 0.96, 1.23 and 0.96e-12. These are insensitive to 4-ulp perturbations of x, with last
     Newton steps of 4e-13 to 1.7e-12. Whether a step passes is a coin toss at the tolerance and independent of dt,
     so halving dt cannot cure it, and the t_dyn trigger fires.
   - The worst rows are again front closure rows, whose scale (4.5e-3) is smaller than on the coarse front (1.1e-2).
     The floor rises as the front is resolved, by refinement or by N.
   - This is the same signature as 14 of the 20 rejections of the verifier's merge-free run near 55 Msun.
3. **The refinement must act early, or the refinement events become the jumps.**
   - With the milder setting (C = 0.3; E4, E6), splits happen only once a cell spans Δ ln κ ≈ 0.6.
   - In E4, 8 event steps exceed δ_target. E6 stops at 53.0 Msun because one split inside the front (Δ ln κ 0.5 per
     face, cell 13 from the surface) moves T by 1.35 % at every dt from 0.15 to 38 yr. The same attempt with splits
     suppressed gives 4.9e-5 (`scripts/probe_event_stop.py`, `data/event_stop_probe.json`).
   - With C = 0.15 (E1, E2) no event step exceeded 2 δ_target.
4. **With the floor removed by fiat (Newton tolerance 1e-10, which is not admissible as an amendment), the best
   variant E2 reaches M_BH = 60.0 and stops for a different reason.**
   - Past 46 Msun, the median Newton iteration count jumps from 7 to 22–23, and genuine Newton failures from the
     predictor (residual 1e-4 to 1e-3) recur every ~10 steps at dt ≈ 20–40 yr.
   - The median dt falls from 349 yr (40–46 Msun) to 134, 40, 23, 13, 10 and 9 yr in the bins up to 60 Msun.
   - It ended at step 4488 with three rejections by the enforced screen "radiation gradient scale outside the
     diffusion limit (Kn > 0.1) at τ ≥ 10". That screen is enforced even in published_assumption mode. dt halved from
     9.4 to 0.94 yr < t_dyn = 1.05 yr.
   - At the stop, Γ in the outer cells reaches 6.7 and Δ ln κ is still up to 0.47 per face, with 176 zones.
   - Rejected fraction 6.3 % (gate ≤ 20 %). Overall median dt 19.4 yr: the smoke-gate requirement of median dt ≥ 50 yr
     would fail on this trajectory even if it reached 100. Separately, the plan's dt floor (dt < 10 yr for 500
     consecutive steps) was close to firing.
   - Not classified: whether the Kn screen at 60 Msun is physical (diffusion breaking down in the super-Eddington
     front) or a mesh-dependent evaluation of the screen. Also open whether the Newton growth (the verifier's 55 Msun
     "unclassified" stop) is a predictor or conditioning problem.

## 4. What a registered amendment would have to say

**DRAFT TEXT ONLY. It needs the user's approval and is not registered.** Under the plan's rule, an amendment is filed
as `benchmarks/ball2011_track_plan_amendment_<n>.json`, registered by sha256 (`.registration.json`) before any run that
uses it, after a preserved failing registered run that shows the defect, and backed by a test that failed before and
passes after. Tolerances (including `numerics.solver.tolerance` 1e-12 and the roundoff exit), physics, microphysics and
targets are never amended.

> **Amendment 6 (draft): opacity-front refinement of the P2 mesh targets.** Kind: numerical-method bug fix
> (registration.amendment_rule).
>
> **Defect.** The registered mesh targets (STARS eq. 2.19 Q, Table 2.1 constants) contain no opacity term. In the
> super-Eddington outer hydrogen opacity front (Γ up to 3.4; κ changing 1.8–2.4× between neighbouring cells), the
> Lagrangian mesh therefore coarsens: 200 → 174 zones in the preserved smoke run, with merges only. A merge there
> changes the discrete solution by the change of its discretization error, independent of dt. At step 1255 this is
> J = 1.86 % at 0.060 Msun below the surface, flat for dt from 0.1 to 172 yr, 3.7 δ_target. Halving dt can never cure
> it, and the dt < t_dyn trigger fires at M_BH = 46.0 (preserved run
> `benchmarks/results/2026-10-14-ball-track/20261014T051036.115106Z_smoke_C1_n200`).
>
> **Change.** Each target is multiplied by clip(C / g_k, 1/F, 1), where g_k is the largest |Δ ln κ| over the faces of
> cell k of the accepted state, maximised over ±W cells. C, W and F must be fixed in the amendment before any
> registered run that uses it, and must not be tuned on registered output. Candidates measured here: C = 0.15, W = 3,
> F = 8.
>
> **Unchanged.** The split/merge thresholds (2.0 / 0.5; C7 1.5 / 0.67), the controller (every component, δ_target,
> reject ratio; every mesh event is measured and can be rejected), the Newton tolerance and roundoff exit, physics,
> microphysics, targets, and every immutable section.
>
> **Bound (why this is not a tolerance relaxation).** The accepted state is the measured state. The merge jump at the
> step-1255 front falls from 1.86 % to 0.146 %, 0.016 % and 0.002 % at ×2, ×4 and ×8 front resolution (third order).
> It is bounded by the spatial error of the mesh, which itself converges at second order (e₁₇₄ = 2.6 %).
>
> **Tests.** A failing-before / passing-after behavioural test on the preserved checkpoint: at step 1255 the registered
> targets merge cells 168+169 and the step is rejected at every dt from 1.34 to 172 yr. With the amendment, no merge
> fires in the front, and a forced merge of the straddling pair on the refined mesh has δ_actual < δ_target. Plus a
> toy-model test that the default (C = 0, off) path is bitwise unchanged.
>
> **Consequences.** Every row 0 must be re-relaxed with the new amendment list. N changes during the run, so the
> convergence campaign's resolution labels (C0/C1/C2) refer to the base N.

**What the draft would NOT achieve, on the evidence of §3:** on its own it moves the stop from 46.0 to 46.7–50.4 Msun.
There the registered Newton tolerance (1e-12) meets a residual floor of about 1e-12 in the resolved front. A global
N = 400 run meets the same floor at 48.1. The tolerance and the roundoff exit are registered tolerances and cannot be
amended. Getting further therefore needs a second, separate numerical-method fix that lowers the attainable residual of
the front closure rows: for example a better-conditioned closure formulation, or smoother microphysics evaluation. That
fix would need its own preserved failing run, which the registered smoke run with Amendment 6 would produce near
46.7 Msun. A change of the 1e-3 scale floor or of the roundoff-exit numbers would be a de facto tolerance change, like
the rejected Amendment 5. Even with both, the diagnostic E2 suggests the track then runs into growing Newton work
(dt ≈ 9–20 yr past 55 Msun, below the smoke gate's 50 yr median) and an enforced Kn screen at 60 Msun, well short of
M_BH = 100. A registered fix would also have to be developed and verified on a trajectory that is not the registered
smoke run itself, so that C, W and F are not tuned on registered output.

## Files

- `scripts/bm_common.py`: shared setup (reads the preserved run; writes only `data/` and `results/ball_merge_20261107/`).
- `scripts/repro_rejections.py` → `data/repro_rejections.json`: goal 1.
- `scripts/error_separation.py` → `data/error_separation.json`, `data/error_separation_tol1e-10.json` (`BM_TOL=1e-10`):
  goal 2.
- `scripts/residual_floor.py` → `data/residual_floor.json`: the Newton residual floor (refined front; 55 Msun state).
- `scripts/front_refine_dryrun.py NAME C W FMAX [TOL] [WALL] [START]` (env `BM_ZONES=400` for E3): goal 3 DRYRUNs, output
  in `results/ball_merge_20261107/<NAME>/` (history.csv, events.jsonl, checkpoint.npz, DRYRUN_report.json), logs
  `results/ball_merge_20261107/<NAME>.log`.
- `scripts/summarize_dryruns.py` → `data/dryrun_summary.json`; `scripts/front_state.py [RUN C W FMAX DT]` →
  `data/front_state.json`; `scripts/probe_event_stop.py RUN C W FMAX DT` → `data/event_stop_probe.json` (E6 stop).
