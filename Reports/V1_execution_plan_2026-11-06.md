# QStarS v1: findings, decisions, and execution plan

Prepared 2026-11-06. **Recommendation: retain the scientific v1 definition, stop repeating the current dense-front run, and first establish a numerically reliable, conservation-consistent step on the archived problem.** Advance the independent published-evolution benchmark alongside the closure work. Treat a validated cool component as an intermediate milestone, not the release.

This is a planning review, not a new convergence result. It covers [Report 25](Report_25.md), the subsequent October 31–November 3 run, the solver implementation, registered validation plans, and the original completion ledger. No production solver, registration, or gate status was changed. Existing uncommitted work was preserved. New read-only calculations and source hashes are in [review_summary.json](evidence/v1_plan_20261106/review_summary.json) and [residual_map.txt](evidence/v1_plan_20261106/residual_map.txt).

**1. What the latest evidence actually establishes**

Report 25 is the latest written report. Its fixed-temperature residual improvement is real within its stated measure: 2.125 to 0.051, with conservation near machine precision. It does not establish statistical equilibrium, released-temperature convergence, or a passing physical audit. The fixed-temperature and released-temperature measures must not be joined into one convergence curve.

The subsequent work materially changes the next-step recommendation:

| Evidence | Finding | Consequence |
|---|---|---|
| `free_r25/README.txt` and `driver.log` | The proposed dense hydrogen-front experiment was already run, with about 1,200 exact columns at depths 19–47 per iteration. Thirteen column builds consumed about 70 hours in total, excluding other work. | Do not propose or restart the same experiment as an untried solution. |
| Read-only replay with the final base's fixed row weights/scales | The starting dense-run state scores 0.02930; the last valid base scores 0.02837, approximately 3.2% lower. | The multi-day run produced little progress and then stalled. This norm is a solver diagnostic, not a percentage physical error. |
| Final iteration 16 | The unconstrained linear residual ratio is 0.811; the bounded-step ratio is 0.998766. The largest unbounded significant-variable step is about 1,396. All eight native trials worsen the merit. | A richer linear span alone is insufficient: the useful unconstrained direction is unusably large, and its bounded replacement predicts almost no progress. |
| Final residual decomposition | Depths 62–167 contain 54.8% of the squared weighted residual; depths 19–47 contain 21.5%; outer depths 1–9 contain 12.0%. | The remaining problem is distributed across hydrogen, helium, and the interior. A hydrogen-window-only strategy misses most of it. |
| Implementation of `--dense-window` | Columns enter the least-squares span and are rebuilt at each base. | This tests fresh span enrichment. It is not the cached two-level subdomain preconditioner proposed in Report 25; that remains a distinct, unproven option. |
| Ball benchmark registration and Report 23 | Q2 already fails at row 0. The smoke track separately stops at a mesh event near 46 solar masses; a merge-free diagnostic later stalls near 55. | There are both a benchmark-definition issue and a numerical evolution issue. Neither is solved by CMFGEN convergence. |

There are two implementation findings to address before another campaign:

- `lm_step_bound` forms `A.T @ A` and the step Gram matrix before solving its damped system. This is a numerical risk in the very ill-conditioned subspace being investigated. It also searches a scalar quadratic-penalty family until a component bound holds; it is not a direct solution of the component-constrained least-squares problem. Test stable augmented QR/SVD and an explicit constrained solve on identical stored products. Neither is yet demonstrated to cure the physical soft modes.
- `history.json` writes `accepted` before checking whether the best trial decreases the merit. Consequently its final entry names `l0.125_d16` even though the driver explicitly rejects every trial and stops. The valid restart is **`results/nk167_trap/free_r25/it16_trial_l0.5`**, not the named iteration-17 trial. This is a reporting/restart hazard; it is not evidence that the rejected state actually advanced the solver.

The tiny singular values and trapped-transition balances support the soft-mode diagnosis. They do not prove that the full nonlinear problem has no equilibrium, nor that every remaining difficulty is physical. The new run makes step construction, feasible directions, and global residual coverage the immediate engineering questions.

**2. Viable routes and my decisions**

| Route | Scope and benefit | Main risk | Decision |
|---|---|---|---|
| A. Repair and extend the existing exact-action solver | Reuse the audited CMFGEN residual/tangent, conserved states, and expensive stored products. Stabilize and constrain the step, then test deeper global Krylov iteration and cached subdomain preconditioning. | Soft modes may still prevent useful global progress. | **Primary route.** Smallest evidence-based change with the best reuse of completed work. |
| B. Nonlinear population relaxation plus continuation | Solve stiff population groups more fully between radiation updates; use that as nonlinear preconditioning for the same released residual. Continue from nearby admissible states if needed. | Fixed-radiation relaxation alone can hide trapped-line errors; a continuation branch may never reach the target. | **Fallback after a bounded A experiment.** Every accepted improvement must survive a fully released radiation pass. |
| C. Move the nonlinear driver to PETSc/SUNDIALS | Retain the residual/Jv and microphysics; obtain a mature solver framework and an independent numerical comparison. | Adapter, scaling, state reuse, and native-pass cost remain; changing libraries does not remove weak physical restoring rates. | Use a small comparison adapter if A exposes driver limitations. Do not rewrite the radiation backend now. |
| D. Release a narrower research milestone first | Package the validated hot-envelope/static components, or later the validated cool component, with explicit scope. Delivers something usable sooner. | Does not satisfy the requested general-purpose scientific v1. | Reasonable optional milestone under a distinct name/version; **do not redefine it as v1**. |

PETSc documents trust-region methods and nonlinear preconditioning, making it a suitable comparison framework, not a promised cure. [PETSc SNES manual](https://petsc.org/release/manual/snes/). SUNDIALS explicitly distinguishes a small-step termination from residual convergence and warns that it may indicate a stalled invalid solution; the QStarS release verifier should preserve that distinction. [KINSOL documentation](https://sundials.readthedocs.io/en/latest/kinsol/Usage/index.html).

**3. Execution order and concrete deliverables**

The chain is: trustworthy baseline → useful constrained steps → released 167-radius equilibrium → physical audits → spatial/frequency convergence → self-selected outer solution → conservative envelope coupling. The Ball replication is independent of the CMFGEN chain. Both feed the broader v1 validation ledger.

**P0 — Freeze the baseline and make decisions reproducible. Planning allowance: one engineering day.**

- Capture the current dirty source snapshot, executable/atomic-data hashes, full launch flags, equation map, two-photon input, and last accepted state. Do not overwrite or trim historical runs.
- Correct future history records to distinguish `best_trial`, `accepted`, and terminal status. Make restart selection use an explicit accepted-state pointer. Add a regression test reproducing iteration 16's all-rejected-trials case.
- Maintain both a frozen reference norm for campaign comparison and the local norm used to construct/accept each step. Preserve per-region and per-species contributions, maxima, unfloored diagnostics, and raw physical audits. Freezing scales must not exclude newly significant equations from the physical checks.
- Reconcile the README, validation documentation, and completion ledger with the latest accepted evidence. Their current summaries range from Reports 15–24; the active gate entries include superseded statements that no selected-wind module exists. Update chronology without converting implemented modules into validated results.
- Create a gate-to-artifact inventory covering the original criteria as well as Gates 1–6. Record unresolved registration decisions explicitly.

Deliverable: one restartable baseline manifest, a trustworthy run-status record, and a consistent release checklist. Exit criterion: the accepted baseline and its recorded residual can be replayed with the same physics; no native solve is launched from a rejected trial.

**P1 — Use archived products to determine why the useful step disappears. Planning allowance: two to three engineering days, then a small native confirmation budget.**

Start with the stored iteration-4 and iteration-16 products, which span an early improvement and the final stall.

1. Compare the current damped normal-equation solve with the equivalent augmented least-squares system `[A; sqrt(mu) W S] y ≈ [b; 0]`, using rank-revealing QR/SVD and basis scaling. Record numerical rank, backward error, step bounds, predicted merit, and sensitivity to the rank cutoff. Do not change the residual weights between methods.
2. Compare that penalty family with an explicit component-constrained solve in the same span. Enforce the linearized hydrogen/helium/charge constraints during construction. Use elimination or a null-space basis where practical, retaining independent checks of the eliminated equations. Include trace-ion adjustment in the modeled update; do not score one step and apply another.
3. Compare the planned step, the complete realized update after every map/projection, and the native residual response. Verify shrinking-step directional derivatives on the actual update path, including released T and the two-photon outer update. The retained one-column finite-difference result is useful but not a complete current-state audit.
4. Separate three outcomes: unstable subspace algebra; good local derivatives but inadequate feasible span; or a poor local nonlinear model despite an adequate feasible span. This determines whether to improve algebra, add directions, or change globalization/nonlinear preconditioning.

Deliverable: an offline method-comparison table plus native confirmations of the best candidates, with no new dense builds needed for the initial comparison. The existing saved products are reusable at their own base; cross-base products are approximations and must not be called exact.

Proposed experiment promotion rule, **not a new scientific release tolerance**: advance a candidate if it predicts at least a 10% decrease with an admissible step and realizes a meaningful decrease on both archived test bases, with conservation preserved and no new native-transfer failure. If the native result contradicts the prediction, diagnose the update/Jv discrepancy before running more iterations. Failure to meet this target does not prove no equilibrium; it rejects an expensive campaign launch on that method.

**P2 — Resolve the distributed slow modes. Planning allowance: three to five engineering days for a first method comparison, conditional on P1.**

- First test 10–20 sequential global Krylov products on a fixed base, with the true scaled residual as the minimization target. Measure wall time, memory, residual by region, and the resulting bounded step. More inner iterations only help if they produce a better feasible nonlinear step.
- If the span is inadequate, build a genuine two-level preconditioner: overlapping hydrogen/helium/outer subdomains plus global interior/thermal directions. Choose coverage from residual and coupling diagnostics, not the old hydrogen window alone. Preserve global exact Jv for the outer solve.
- Reuse/factor subdomain approximations across nearby accepted states only while fresh probe products and actual/predicted reductions support reuse. Rebuild on failed probes or front motion. Account for setup cost amortized over accepted progress, rather than comparing raw iteration counts.
- Consider a reduced physical basis using ion fractions and excitation departures if constrained full-population steps remain badly conditioned. Verify the coordinate transformation and its tangent; logarithms alone do not enforce nuclei conservation or solve trapped-line coupling.
- If neither the stable constrained solve nor deeper/global preconditioning produces useful steps, switch to Route B: positivity/conservation-preserving population relaxation in overlapping regions, followed by full radiation and temperature updates. Use it as an initializer or nonlinear preconditioner. Local fixed-J convergence is never the acceptance certificate.

Deliverable: one selected solver configuration and a rejected-alternatives record. Proposed stop rules: stop after three accepted steps with less than 1% cumulative frozen-norm improvement; stop on all-trials-rejected; and review after at most three full nonlinear iterations per new expensive configuration. A longer run requires demonstrated progress per wall hour. Do not spend another 70 hours rebuilding the same front span.

**P3 — Close the component gates with the final physics. Effort estimated only after P2 establishes a convergence rate.**

- Finish the released-T 167-radius solve. Keep unshifted equations as the target, retain number/charge conservation, and distinguish approximate bounded steps from a converged unbounded Newton correction. A small bounded step is not convergence.
- Repair the closing chain before using it: `post_convergence.sh` defaults to build 31, while this solve uses build 36; its child scripts default even earlier when invoked directly. Carry the accepted executable and complete physics flags explicitly through cleanup and verification. Verify two-photon stored-J handling and self-consistency, `INC_AD/INC_ADV`, and native electron heating/balance diagnostics. Do not assume a matching executable alone reproduces the same residual.
- Run the registered trace cleanup and all-level native verification. Existing limits include last fractional correction ≤1e-3 on all levels, nuclei/charge ≤1e-6, electron balance ≤1%, and global **and cumulative** radiation-energy error ≤0.1%. Apply the remaining spectral-temperature, recombination-coverage, and material checks in the registered plans. Keep floor-limited failures visible; no automatic trace-level exemption.
- Rebuild the clean 333-radius base, transfer conservatively with verified equation mapping, solve independently, and audit it. Perform the full registered doubled-frequency recipe (expected 85,359 frequencies), auditing that state independently too. Check the two-photon field on the refined frequency grid. Comparisons need passing inputs and the registered internal-profile limits, not only similar endpoint luminosities.
- Finish atom-size and H-minus rate sensitivity and the relevant relaxed-response checks before production use. These do not disappear when the Newton solver converges.

Deliverable: a reproducible, independently audited cool closure and passing spatial/frequency evidence. This is the first major scientific milestone. Source contracts: [cool validation](../benchmarks/cmfgen_cool_validation_plan.json), [component thresholds](../benchmarks/cmfgen_validation_plan.json), and [frequency/refinement plan](../benchmarks/cmfgen_cool_frequency_plan.json).

**P4 — Select the outer solution and couple it to a real envelope. Starts only after P3.**

- Exercise the existing `cmfgen_wind.py` module using converged released states, re-solving closure after each hydro update. Check both optical-depth and radius formulations, critical regularity, base conditions, outer extent, energy and momentum, and resolution on the selected structure.
- Preserve the registered outcomes: transonic wind, bounded hydrostatic/no-wind solution, photon-tiring exit, Mach-one MLT exit, invalid reduction, or numerical failure. A legitimate no-wind state is useful evidence for that case; it does not supply missing escaping-wind validation. If the isothermal-type reduction fails its screens, the coupled characteristic/energy problem is additional work, not a tolerance adjustment.
- Construct an admissible envelope at the selected node's base. Establish the MESA-to-CMFGEN energy-zero bridge, retain convective flux at the handoff, and implement an independent Gate 4 evidence verifier; the current node interface intentionally reports `validated: false` without one.
- Build enough validated nodes to cover a declared evolution interval. Enforce no extrapolation and tested interpolation error. Check independent feed/wind streams, finite thermal storage, mass/energy accounting, retries, remapping, restart equivalence, and timestep/spatial refinement on the actual coupled model.

Deliverable: one conservatively evolved physical envelope with a justified outer boundary and a bounded validity domain. Component tests and synthetic nodes do not satisfy this milestone. See [wind plan](../benchmarks/cmfgen_selected_wind_plan.json) and [node contract](../docs/wind_node.md).

**P5 — Repair the independent Ball replication path. Can proceed while P1–P3 run.**

The published-assumption Ball model uses no CMFGEN wind, so it should not wait for cool non-LTE closure. Schedule separate work periods or compute slots; do not oversubscribe the shared machine during timing comparisons.

- Audit Q2 against its cited source and distinguish the paper's claimed onset/location of an inversion from the registration's blanket row-0 absence test. Quantify the inversion on common physical coordinates; a shrinking adjacent-cell increment is not proof that the resolved inversion disappears.
- Preserve the existing registered campaign as FAIL. If the source audit establishes a registration mismatch, prepare a **new, explicitly approved registration** that tests the supported physical claim and still reports every inversion and the old Q2 result. No silent amendment: the existing rule forbids changing targets/physics/tolerances. If Q2 is faithful, investigate model/microphysics differences instead; do not re-register merely to obtain PASS.
- Reproduce the 46-solar-mass merge event. Separate temporal error from remap-induced structural error, bound both, and refine the opacity front rather than exclude mesh events from the controller. Retain rejected Amendment 5 as rejected. Then classify and fix the separate merge-free Newton stall near 55 solar masses.
- Re-run the registered smoke target, then the complete C0–C7 convergence/restart campaign and 44 comparison cells, including endpoint classification and the existing ledger tolerances.

Deliverable: one full independently registered published evolutionary replication, or a precisely located model discrepancy. Even a passing Ball result supplies only one replication, not the complete v1 requirement.

**P6 — Complete the requirements beyond Gates 1–6 and release.**

Before declaring v1, close the original ledger explicitly:

- Validate a conservative moving physical cavity/reservoir and the claimed evolving inner-boundary families. The quasi-static Ball path and moving-face identities do not certify this for the finite-thermal evolution engine.
- Register a second independent complete evolutionary benchmark and a cross-paper prediction with no retuning. First candidate to assess: the Hassan/MESA-QUEST family already named in `benchmarks/manifest.yaml`, given the existing MESA infrastructure. Fiacconi or Dotan remain alternatives but require their own published physics/outer closures. Perform a source/input-availability audit before selecting; this review has not established that any candidate is immediately reproducible. Static profiles and the Ball thesis/paper versions do not count as two independent evolutions.
- Validate the claimed production microphysics/atomic sufficiency and model-family robustness; account for existing printed-value discrepancies without hiding them. Each discrepancy needs reproduction, attribution, and an explicit verdict under its registered scope.
- Run complete production-family resolution, timestep, remap, rejection, and restart campaigns; measure full-run Mac memory/storage/wall time. Hot controlled examples are insufficient substitutes.
- Reproduce from a clean documented installation with the native dependencies configured. Update the release evidence, examples, version, and citation only after the gate verifier is green. Confirm `qstars validate --require-v1` against the complete ledger, including any required verifier implementation work.
- Keep QSRD downstream of independently validated framework v1. Its campaign cannot be used to validate the framework that it depends on.

Deliverable: a release dossier mapping every claim to a reproducible passing artifact and every remaining limitation to the declared supported domain.

**4. Decisions that need to be settled before the affected campaigns**

These do not block the numerical diagnostics in P0–P2. This review changes none of the registrations.

| Decision | Recommended treatment |
|---|---|
| Meaning of v1 | Preserve the general-purpose scientific requirements. Name smaller deliveries explicitly as research milestones. |
| Ball Q2 | Authorize a source-faithful successor registration only if the audit demonstrates a mismatch; retain the original FAIL. A numerical restart cannot solve this decision. |
| Frequency luminosity/cumulative-heating criteria currently marked pending | Retain the proposed 1% common-grid criteria with the registered normalizations, and settle their status before the campaign. |
| Floating-point-floor exemptions | Do not exempt levels in the initial campaign. Measure reproducibility and improve the verifier first; any later exemption requires a justified, visible scope decision. |
| Wind mass-loss refinement | Retain the existing 1e-3 log-mdot criterion; do not relax it to 1% to obtain a pass. Keep the stricter dual-opacity and dual thermal-screen defaults pending formal disposition. |
| Second evolutionary target and exact supported families | Select after a bounded source/input audit, then register all primary observables, tolerances, and no-retuning rules before prediction. |

The relevant recorded constraints are in the [Ball registration rule](../benchmarks/ball2011_track_plan.json), [cool closing phase](../benchmarks/cmfgen_cool_validation_plan.json), [frequency plan](../benchmarks/cmfgen_cool_frequency_plan.json), and [wind plan](../benchmarks/cmfgen_selected_wind_plan.json). They are scientific registration decisions, not a request for permission to perform the already scoped analysis.

**5. Schedule, cost control, and confidence**

Allocate the first **one to two engineering weeks** to P0–P2, with the Ball source/mesh audit progressing independently when resources permit. This is a planning allowance for a decisive method evaluation, not a prediction that equilibrium or v1 will be reached in that time. The first checkpoint should deliver stable constrained-step comparisons; the second should either demonstrate sustained released-state progress or trigger Route B/C.

Once a candidate demonstrates progress, estimate closure cost from accepted reduction per wall hour, including setup, failed trials, cleanup, and native verification. Re-estimate after 333-radius and doubled-frequency pilots. Do not extrapolate Report 25's initial improvement into a release date: the subsequent 70-hour dense-column campaign demonstrates a different regime.

Confidence is high that repeating the current dense-front setup is a poor use of compute, and high that v1 has independent requirements beyond closure. Confidence is moderate that step algebra, conservation-consistent directions, and broader preconditioning are the best next experiments. Whether those produce an audited equilibrium, an admissible selected outer solution, and full published tracks remains unproven. A defensible v1 date comes after those feasibility gates, not before them.

Immediate execution decision: **P0 → P1 using existing products → a bounded P2 comparison; start the independent Ball audit without changing its verdict.** Keep the expensive full campaign gated on measured useful steps.
