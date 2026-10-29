# Smoke-gate stop at M_BH = 46: diagnosis and adversarial verification

## Diagnosis

**Smoke-gate stop at step 1255: diagnosis**

The jump comes from a split/merge mesh event (two outer cells merged), not from evolution over dt, and the step controller counts it as a timestep error. This is class (a), a numerical-method defect. The draft fix gets the run past M_BH = 46, but a DRYRUN stops again at about 53.7 Msun on a different failure. So Amendment 5 alone will not get the smoke gate to M_BH = 100.

Nothing was written under `benchmarks/results`, and `src/qstars` and the driver are unchanged; the code-freeze check passes and the plan sha256 is still `a10b4166…`. Nothing was committed. All work ran single-thread at nice 15, about 15 min of CPU in total.

**1. Reproduction** (unpatched code, from the run's `checkpoint.npz` with its stored config and `build_ball_abi7`)

| dt (yr) | delta_actual (= profile_dlnT) |
|---|---|
| 1e-3 | 0.01286 |
| 1e-2 | 0.01682 |
| 0.1 | 0.01835 |
| 1 | 0.01857 |
| 10 | 0.01840 |
| 100 | 0.01664 |
| 172.144 | 0.01525 |

- The 172.144 yr value matches the run's first rejection exactly.
- Every attempt, at every dt, merges the same pair of cells (174 → 173 zones). Newton converges in 3–5 iterations every time.
- The jump does not vanish as dt → 0. Going further down, it falls to 0.0053 at 1e-5 yr and 0.0018 at 1e-8 yr. It reaches its full size over the thermal time of the outer cells (hours to days), far below anything the controller can choose (t_dyn = 0.83 yr).
- **Same step without the merge** (on the mesh before the event): the change is 1.6e-5 to 1.9e-5 for dt ≤ 0.1 yr, 2.4e-5 at 1 yr, 2.1e-3 at 100 yr and 3.59e-3 at 172 yr. It grows in proportion to dt, matching the previous accepted steps.

**2. Location**
- **The merge:** cells 168 and 169 of 174 merged. Their mass-to-target ratios were 0.662 and 0.497, against the registered merge threshold of 0.5. The merged cell sits 0.050 Msun below the surface, at T ≈ 7400 K.
- **Largest change:** cell 167, the next cell inward, which the merge itself does not touch. It is 0.060 Msun below the surface and 9949 Msun above the base, at optical depth about 2.8. At dt = 1 yr:
  - T goes from 8190 to 8343 K (+1.86%).
  - Density falls about 10% (3.66e-11 to 3.31e-11).
  - Opacity rises from 1.07 to 1.20.
  - Cells 166, 165 and 164 change by +1.50%, +0.80% and +0.44%.
  - The merged cell itself moves only +0.23% from its remapped temperature.
- **Which unknowns move:** only T jumps. Radii change by at most 1.2e-4, surface luminosity by 1.4e-5, photospheric density by +0.1%. The 33-cell fixed-mass base region does not change.
- **What the region is:** the outer hydrogen opacity peak. Density increases outward by a factor of 9, and the local Eddington factor (κL/4πcGm) exceeds 1 in cells 158–168, peaking at 3.4. Every face there is convectively unstable, but convection carries almost nothing (efficiency ≤ 2e-3), so the gradient is the radiative one. No opacity clamp applies, and all cells are above the STARS opacity-table switch (T ≥ 6500 K).

**3. Candidate causes ruled out**
- **Opacity or EOS switch, clamp or table edge:** opacity, pressure, internal energy, adiabatic gradient and heat capacity are smooth over 6400–13500 K at these densities, with no isolated jumps.
- **Convection on/off switch:** no face anywhere in the profile changes regime, with or without the merge.
- **Photosphere branch:** same branch, density change 0.1%.
- **Base region boundary:** unchanged.
- **The mapped metric itself:** the registered Lagrangian metric shows the same jump (0.0186).
- **Fold or second root:** continuing the step in dt from 1e-8 yr to 172 yr gives one smooth solution branch. The smallest singular value of the Jacobian is the same with and without the merge (4.91e-4 at 1 yr; 3.8e-3 vs 3.6e-3 at 172 yr), and at 1 yr its singular vector lies outside the merge region.
- **`fold.classify`** (run in scratch at step 1255): it finds no fold in M_BH from 44.8 to 49.2 Msun (smallest singular value flat at 1.9e-4, determinant sign constant) and returns `numerical_failure`.
- **Size check:** halving the outer 34 cells changes T at the same mass coordinate by −1.74%. That is about the size of the merge jump (+1.86%), so the jump is the N = 200 discretization error in the outer opacity peak.

**4. Classification: (a).** The physics is continuous: the same step without the merge is smooth and proportional to dt. The merge is the registered split/merge rule acting in the Lagrangian envelope. It changes the discrete solution by the change in its discretization error, the same at every dt. The controller treats that as a timestep error, so halving dt can never succeed and Ball's trigger fires. This is the same defect class Amendment 1 removed at the base; the envelope above the fixed-mass base region still has it. The 22 earlier merges in the run were small (0.0013–0.0058). Merges became dangerous once the hydrogen opacity peak moved into the outer cells as Teff fell to about 6460 K. The mesh only ever merged (200 → 174 zones, no splits), so that region keeps getting coarser.

**5. Fix, test and draft**
- **Fix:** on a step with a split or merge, the controller measures every component on the same step solved on the mesh before the event. The code already runs that counterfactual solve on every event step for its existing remesh diagnostics. It is enabled by a new setting `stepcontrol.events='counterfactual'`; the default `'measured'` is the registered rule and leaves non-event steps bitwise unchanged. The accepted state is still the solution on the new mesh, and the measured jump is still logged. The core change:
  ```python
  free=diag.get('event_free') if s.get('events','measured')=='counterfactual' else None
  if free:diag,after=free,free['observables']
  ```
- **Patch:** `amendment5_mesh_event_controller.patch` (in scratch, with a preserved copy in the evidence folder). It touches `quasistatic.py`, `stepcontrol.py`, the driver's `AMENDABLE` list and `tests/test_quasistatic.py`, and `git apply --check` passes on the frozen tree.
- **Tests:**
  - `test_registered_rule_counts_a_mesh_event_jump_at_every_dt` documents the defect on the toy model; it passes before and after.
  - `test_mesh_event_jump_is_not_a_timestep_error` is a toy-model test with a forced merge: the amended metric equals the no-merge step to 1e-12 and is ≤ 1e-6 at dt = 1e-6 yr.
  - `test_smoke_step_1255_merge_is_measured_on_the_pre_event_mesh` (MESA opt-in) uses the preserved checkpoint: the registered rule rejects at 1 and 172 yr, the amended rule accepts.
  - Both of these fail before the patch, but only because the `events` setting does not exist yet (`ValueError: Unknown stepcontrol settings ['events']`). That is the same kind of "fails before" the earlier amendments used.
- **Regression:** on a scratch copy of the tree with the patch applied, 144 passed and 6 skipped without MESA. With MESA, 156 passed and 2 skipped.
- **Patched reproduction at step 1255:** delta_actual is 2.4e-5 at 1 yr, 3.59e-3 at 172 yr and 6.2e-3 at 300 yr, all accepted.
- **Draft:** `benchmarks/ball2011_track_plan_amendment_5.draft.json` (status DRAFT). It cites the preserved smoke run and documents the defect, the change, the tests, the DRYRUNs and the consequences.
- **Registration cost:** the driver checks that a row 0 was relaxed with the same amendment list. So with Amendment 5, every row 0 must be relaxed again with amendments 3, 4, 2, 1, 5 — about 15–20 min. The results should be bitwise the same, because relaxation never reads the controller.

**6. DRYRUN** (scratch, from step 1255 with the controller restored to dt = 172 yr, smoke options, at most 400 steps; 290 s)
- 343 accepted steps in 355 attempts, median dt 127 yr.
- 7 more outer merges were accepted. Their measured jumps were up to 0.058; the counterfactual metric was 0.0036–0.0038.
- It reached M_BH = 53.75 at step 1598 (t = 0.2506 Myr), then hit the termination trigger again on a different failure:
  - Another outer merge moves T of cell 160 by +8.8% at every dt.
  - The solution on the new mesh exists and is smooth in dt: continuation from 1e-3 yr converges all the way to 77 yr.
  - But Newton started from the registered predictor stalls at residual 2.65 for every dt ≥ 1 yr, and the F1 fallback also fails ("Pseudo-time step collapsed").
  - Amendment 5 does not cover this. It would need its own preserved failing run, which the registered smoke run with Amendment 5 should produce near 53.7 Msun.
- **Exploratory extra, not proposed:** I also prototyped a Newton start that continues the step in dt. With it the run gets through step 1598 and one more merge, to M_BH = 55.22 (step 1624). There it fails on ordinary steps with no merge: continued in dt, the step becomes singular between 0.48 and 0.70 yr, below t_dyn of 0.98 yr. The singular direction is T of cell 160, again in the outer hydrogen zone. `fold.classify` finds no fold at that state. I have not classified this.

**Not established**
- The converged size of the outer discretization error: refining those cells by 3×, 4× and 8× failed in Newton from my scratch starting states.
- Whether the N = 400 and N = 800 runs see similar merge jumps; I did not run them.
- Any fix for the step-1598 failure, and the cause of the step-1624 stop.
- Whether the smoke gate can reach M_BH = 100 with any set of fixes.

Files are in `/Users/abraiymbek/Projects/QStarS/Reports/evidence/report23/ball/smoke_jump/`:
- smoke_jump_evidence.json
- amendment5_mesh_event_controller.patch
- data/
- scripts/
- DRYRUN_amendment5_from_step1255/
- DRYRUN_exploratory_ladder_from_step1598/

Draft amendment: `/Users/abraiymbek/Projects/QStarS/benchmarks/ball2011_track_plan_amendment_5.draft.json`

Scratch:
- `/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/ball_jump/amendment5_mesh_event_controller.patch`
- `/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/ball_jump/patched/`

## Verification (structured)

```json
{
 "agrees": false,
 "classification": "a",
 "amendment_admissible": false,
 "checks_run": "All checks ran single-thread at nice 15 with my own scripts in scratch (verify/), on the unpatched frozen tree. Nothing was written under benchmarks/results. src and the driver are unchanged: the stage7 freeze check gives 0 mismatches and the plan sha256 is still a10b4166... Evidence is preserved in /Users/abraiymbek/Projects/QStarS/Reports/evidence/report23/ball/smoke_jump/verifier/ (scripts/, logs/, DIAGNOSTIC_DRYRUN_registered_controller_merges_off_from_step1255/).\n(1) Reproduction from the smoke checkpoint (step 1255, 174 zones, M_BH 46.0038). With the registered merge, delta_actual = profile_dlnT is 0.0085, 0.01682, 0.01857 and 0.01525 at dt = 1e-4, 1e-2, 1 and 172.144 yr. Every attempt merges the 5th and 6th cells from the surface (174 -> 173) and the largest change is in cell 167 (0.0601 Msun below the surface, T 8190 -> 8343 K). With the merge suppressed (merge threshold 0.45, in memory only) delta_actual is 2.42e-5 at 1 yr and 3.585e-3 at 172 yr.\n(2) Outer structure at step 1255: dlnkappa between neighbouring cells is -0.87 (167/168) and -0.60 (168/169), dlnT is 8-10% per cell, and the local Gamma reaches 3.4.\n(3) Forced remesh variants at 1 yr. Merging 150+151 or 140+141 gives 3.4e-4 and 6.2e-4. Merging 169+170 gives 2.37%, and 167+168 gives 1.58%. Splitting 168 and 169 gives -0.77% at 1 yr and 1.12% at 172 yr, so even a refinement in that zone is rejected at the controller's dt.\n(4) Refinement, which the investigator did not establish. T at the mass coordinate of cell 167 changes by -1.88% (x2), -2.41% (x4) and -2.55% (x8). That is second-order convergence: the 174-zone solution is about +2.6% off and the merged mesh about +4.5%. The refined solves stop at a roundoff floor (residual 2e-12 to 2.5e-11, last Newton step 3-6e-13); they are not near-singular.\n(5) fold.classify at step 1255: numerical_failure. sigma_min stays at 1.884-1.885e-4 over M_BH 44.83-49.21 and the sign of det J does not change.\n(6) Patch, applied to a scratch copy. git apply works. In the default and 'measured' modes the new state, diag and metrics are bitwise identical to the unpatched code at the event steps (1 and 172 yr); diag gains only four extra event_free_* keys. The 'counterfactual' metrics equal exactly those of the merge-suppressed step. Tests: before the patch, 2 fail (both with ValueError for the unknown 'events' key) and 1 passes; after it, 3 pass (the MESA opt-in test included).\n(7) The investigator's states from the Amendment 5 DRYRUN, rerun on unpatched code. At step 1598 (167 zones) the opacity front lies across cells 160/161: dlnkappa = -2.03 (a factor 7.6) and dlnT = -16% per cell. The step with the registered merge fails in Newton (residual 2.65 at 1 and 30.9 yr). The same step without the merge converges: 5 iterations and delta 5.5e-5 at 1 yr, 20 iterations and delta 1.65e-3 at 30.9 yr.\n(8) DIAGNOSTIC DRYRUN (not a proposed setting): unpatched code, registered controller, merges switched off in memory, continuing from step 1255. It passed 53.75 Msun, then stopped at step 1667, M_BH = 55.01, with termination_trigger (dt 0.50 yr < t_dyn 0.97 yr): 412 accepted steps in 432 attempts. 14 of its 20 rejections are Newton stalls at residual 1.03e-12 to 2.4e-10, just above the registered 1e-12 tolerance; 6 are genuine (4.5e-6 to 3.9e-4). At step 1667 the opacity front sits at cells 158/159 (dlnkappa = -2.25, Gamma up to 4.4). Newton iterations from the predictor rise 4 / 10 / 15 / 23 / 40 at dt = 0.1 / 0.3 / 0.5 / 1 / 2 yr, and Newton fails from 5 yr up. Refining the front (x4 or x8) makes Newton harder, not easier.",
 "reasoning": "I agree with the diagnosis but not with the fix.\n\nWhat I confirm:\n- **Location:** the merge of cells 168 and 169 (dm/target 0.662 and 0.497), with the largest change in cell 167 of the outer hydrogen opacity peak.\n- **Mechanism:** the merge causes a jump that builds up over the cells' thermal time (well under 0.01 yr) and then stays at 1.5-1.9% for every dt from 0.01 to 172 yr. The registered controller counts it, so halving dt can never succeed.\n- **Classification (a):** the refinement sequence converges at second order with ordinary Newton steps. sigma_min and fold.classify show no fold, and merges outside the peak are harmless (about 5e-4). So the jump is a change of discretization error in an outer opacity peak that the registered Lagrangian split/merge mesh resolves very poorly (kappa changes by a factor 1.8-2.4 between neighbouring cells), not physics.\n\nThe defect, though, sits in the mesh rule (merges in an unresolved peak). Here the controller is right to refuse a 3.7x-tolerance change of the accepted state.\n\nWhy I reject the proposed Amendment 5:\n- It is numerical-method only, off by default and bitwise neutral in measured mode; I verified all three.\n- But it removes the per-step tolerance from mesh-event jumps entirely, with no bound and no check that the solution after the event continues the same branch. The state it accepts is not the state the controller measured.\n- In the investigator's DRYRUN it accepted jumps up to 5.84%, which is 11.7x delta_target. That is in effect a tolerance relaxation, which the amendment rule forbids, even though no number changes.\n- It also lets the opacity peak coarsen further: 7 more merges put the front across a single cell pair with a kappa ratio of 7.6.\n- The step-1598 blocker it then meets is the same root cause, not a different failure: the same step without its merge converges easily.\n- Precedent points the other way. Amendment 1 cured the identical defect class by removing events from the region, so the accepted state stayed the measured state.\n\nWhat any merge fix would and would not achieve: with merges switched off and the registered controller, the run passes 53.7 Msun but stops at 55.0 Msun anyway, on Newton growth in the same front and roundoff stalls just above 1e-12. Neither Amendment 5 nor a mesh fix gets the smoke gate to M_BH = 100. That stop near 55 Msun is unclassified (d) and must stay a numerical blocker under stop_rules.outcome_if_exhausted until it is diagnosed. I found no use of the published tables.",
 "problems": [
  "Amendment 5 takes mesh-event jumps out of error control altogether, with no upper bound. The accepted state is the solution on the new mesh, but delta_actual is measured on the counterfactual (old-mesh) solution. The DRYRUN accepted measured jumps of 0.0152, 0.0308, 0.0059, 0.0038, 0.0584, 0.0366 and 0.0027 (up to 11.7x delta_target, about 16% in total that no tolerance checked). This is in effect a relaxation of the registered tolerance, which the amendment rule forbids even though no number changes.",
  "There is no safeguard that the solution after an event continues the same branch. Continuity in dt was shown only for step 1255. Under events='counterfactual', any later merge or split would be accepted silently, whatever its jump, including a jump onto another branch.",
  "The step-1598 blocker is mislabelled as 'a different failure'. It has the same root cause: a merge inside the unresolved hydrogen opacity front. On the investigator's step-1598 state the same step without the merge converges (5 iterations, delta 5.5e-5 at 1 yr; 20 iterations, 1.65e-3 at 30.9 yr). Amendment 5 speeds up the coarsening: by step 1598 the front sits across one cell pair with dlnkappa = -2.03 (a factor 7.6) and dlnT = -16% per cell.",
  "The fix targets the symptom (the controller) rather than the defect (the registered split/merge mesh under-resolving the outer hydrogen opacity peak). A mesh-side fix, as Amendment 1 did for the base, would keep the accepted state equal to the measured one. Caveat: even splitting cells 168 and 169 at step 1255 moves T by 0.77% at 1 yr and 1.12% at 172 yr, so a refinement there is also rejected at the controller's dt.",
  "No merge-related fix reaches M_BH = 100. With merges switched off and the registered controller (diagnostic DRYRUN), the run passes 53.75 Msun but stops at 55.01 Msun (step 1667) with termination_trigger. This agrees with the investigator's exploratory stop at 55.22. Of 20 rejections, 14 are Newton stalls at residual 1.03e-12 to 2.4e-10, just above the registered 1e-12, and 6 are genuine. Newton iterations grow from 4 to 40 between dt = 0.1 and 2 yr in the front at cells 158/159, where Gamma reaches 4.4. Refining the front makes Newton harder. This stop is unclassified.",
  "The tests 'fail before' only through ValueError('Unknown stepcontrol settings [events]'), not through any behavioural assertion. The toy test, that counterfactual equals the no-merge step to 1e-12, is close to a tautology: I confirmed the counterfactual metrics are bitwise those of the merge-suppressed step. Nothing tests that the trajectory the amended controller accepts is accurate.",
  "Even in the default mode, the patch adds four event_free_* keys to every step's diag and to every accepted event in events.jsonl. States and metrics are bitwise unchanged, but the run outputs are not byte-identical to the registered code, and the code-closure digest changes for every campaign run.",
  "The investigator left the converged size of the outer discretization error open; it is now measured. At 0.060 Msun below the surface, T converges at second order to about 2.6% below the 174-zone value (x2 -1.88%, x4 -2.41%, x8 -2.55%), so the merged mesh is about 4.5% off. The refined solves fail only at a roundoff floor (residual 2e-12 to 2.5e-11, last step 3-6e-13), above the registered 1e-12 tolerance and 1e-13 roundoff step. That floor is a risk for the N = 400/800 campaign in this region and also accounts for many of the near-tolerance rejections in both DRYRUNs."
 ]
}
```
