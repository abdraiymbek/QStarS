# QStarS v1 gate-to-artifact inventory

Prepared 2026-11-06 for [the v1 execution plan](../../V1_execution_plan_2026-11-06.md) P0, bullet 5. **QStarS is not v1.** Every status below is quoted from its recorded source and none is changed here. "Implemented" never means "validated". Machine-readable copy: [gate_inventory.json](gate_inventory.json).

Latest evidence: [Report 25](../../Report_25.md) is the latest written report. The stage-2 dense-window run (2026-10-31 to 11-03, `results/nk167_trap/free_r25`; no Report 26 was written) lowered the frozen norm from 0.02930 to 0.02837. At iterate 16 every trial was rejected; the last valid base is `it16_trial_l0.5`.

Groups: **OC** are the original twelve completion criteria. **PG** are the nine prioritized gates (the user's priorities in `completion_gates.json` and the docs/validation.md matrix). **G0-G6** are the v1 chain of Report 23 §8. **X** are other open criteria.

## Original completion criteria

### OC-1 — Equations and boundary conditions sourced or derived

- **Recorded status:** PASS. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=1] (2026-10-22 assessment).
  > "status": "PASS", "scope": "Implemented closures have derivations and references; this does not supply missing production closures."
- **Artifacts:** `docs/equations.md`; `docs/sources.md`; `src/qstars/data/sources.json`; `docs/literature.md`; `docs/cmfgen_selected_wind.md`; `docs/nlte_cool_closure.md`
- **Registered plan:** none registered
- **Missing:** Derivations exist for implemented closures only. Any closure added for the selected outer boundary, moving cavity or second benchmark family must be derived and sourced when it is added.
- **Plan phase:** P4, P6

### OC-2 — Mass and energy conserved to quantified tolerance

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=2] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Tested resolved-domain equations pass; missing cavity/atmosphere prevent whole quasi-star certification."
- **Artifacts:** `benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T030639.520967Z_robustness_validation/campaign.json`; `benchmarks/results/2026-10-14-ball-track/ (Ball smoke: per-step energy <= 1.2e-14, mass <= 3.6e-16)`; `Reports/evidence/report24/ (number/charge conservation defect and --conserve)`; `Reports/evidence/report25/progress_log.txt (closure conservation 1e-15)`
- **Registered plan:** `benchmarks/ball2011_track_plan.json (conservation section)`
- **Missing:** Whole quasi-star certification needs a validated cavity/reservoir in the thermal engine and a validated atmosphere/wind boundary; the node-coupled ledger (Gate 5) has only a synthetic node.
- **Plan phase:** P4, P6

### OC-3 — Conservative retries, remeshing, restart and timesteps

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=3] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Controlled hot regressions/campaigns pass; broader physical-model coverage remains."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T030639.520967Z_robustness_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T031253.242612Z_time_method_validation/campaign.json`; `Reports/evidence/report23/ball/smoke_jump/diagnosis_and_verification.md`; `Reports/evidence/v1_plan_20261106/review_summary.json (history.json acceptance caveat)`
- **Registered plan:** `benchmarks/ball2011_track_plan.json (restart, convergence_campaign, stop_rules)`
- **Missing:** Production-family retry/remap/restart campaigns. The Ball merge event at M_BH = 46 Msun shows a remap-induced structural jump that is not yet bounded. MESA remap-threshold cases never remapped. The cool-closure driver records "accepted" before its rejection check, so restart selection needs an explicit accepted-state pointer.
- **Plan phase:** P0, P5, P6

### OC-4 — Representative timestep and mesh convergence

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=4] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Controlled hot internal profiles pass; general structures untested."
- **Artifacts:** `benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T031651.519611Z_midpoint_512_refined_time/method_comparison.json`
- **Registered plan:** `benchmarks/cmfgen_cool_frequency_plan.json`; `benchmarks/cmfgen_cool_validation_plan.json`; `benchmarks/ball2011_track_plan.json (C0-C7)`; `benchmarks/cmfgen_selected_wind_plan.json`
- **Missing:** Closure 167/333 and doubled-frequency pairs (Gates 2-3), the Ball C0-C7 campaign, wind resolution on the self-selected structure, and production-family timestep/mesh convergence.
- **Plan phase:** P3, P4, P5, P6

### OC-5 — Several independent published results reproduced

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=5] (2026-10-22 assessment).
  > "status": "FAIL", "scope": "Static subsets pass; the clarified complete-evolution requirement fails."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T031903.825105Z_validation/validation.json`; `benchmarks/manifest.yaml`; `benchmarks/results/2026-10-14-ball-track/`; `Reports/Report_23.md`
- **Registered plan:** `benchmarks/ball2011_track_plan.json (+ registration, amendments 1-4; amendment 5 rejected)`
- **Missing:** Zero complete evolutionary replicas. Gate 6 is FAIL. A second independent complete evolution and a cross-paper prediction without retuning are neither selected nor registered.
- **Plan phase:** P5, P6

### OC-6 — Prescriptions explicit and traceable

- **Recorded status:** PASS. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=6] (2026-10-22 assessment).
  > "status": "PASS", "scope": "Implemented families and disagreements are documented; missing evolving families remain failures under priority 3."
- **Artifacts:** `docs/model_families.md`; `src/qstars/data/sources.json`; `benchmarks/manifest.yaml`; `docs/sources.md`
- **Registered plan:** none registered
- **Missing:** Missing evolving inner-boundary families are carried by PG-3.
- **Plan phase:** P6

### OC-7 — Thermal quasi-star evolution through time

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=7] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Finite thermal evolution with BH growth works in restricted envelopes; full lifecycle unverified."
- **Artifacts:** `configs/hot_resolved_validation.yaml`; `benchmarks/results/2026-09-28-verified-spatial/`; `benchmarks/results/2026-10-14-ball-track/ (quasi-static track to M_BH = 46 Msun)`
- **Registered plan:** `benchmarks/ball2011_track_plan.json`
- **Missing:** No full lifecycle with a justified outer boundary. The quasi-static Ball track stops at a numerical blocker, and no node-coupled envelope has been evolved.
- **Plan phase:** P4, P5, P6

### OC-8 — Independent effective feeding and wind

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=8] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Separate stream accounting works; self-consistent wind/atmosphere coupling is missing."
- **Artifacts:** `configs/two_stream_pulse.yaml`; `benchmarks/results/2026-09-28-closure/20260928T024813.279478Z_fiacconi_wind_validation/validation.json`; `src/qstars/cmfgen_wind.py`; `tests/test_cmfgen_wind.py`; `src/qstars/wind_node.py`; `tests/test_node_wind.py`
- **Registered plan:** `benchmarks/cmfgen_selected_wind_plan.json`; `benchmarks/self_selected_plan.json (Report 12, FAIL)`
- **Missing:** A self-consistent selected wind and its coupling (Gates 4-5). cmfgen_wind.py and wind_node.py are implemented and tested on analytic and synthetic inputs only.
- **Plan phase:** P4

### OC-9 — Identify breakdown of assumptions

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=9] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Domain exits and numerical failures are distinct; current screens do not prove all production-model validity."
- **Artifacts:** `docs/validity.md`; `src/qstars/validity.py`; `benchmarks/ball2011_track_plan.json (screen registry)`
- **Registered plan:** `benchmarks/ball2011_track_plan.json (replication_mode published_assumption: screens recorded, not enforced)`; `benchmarks/cmfgen_selected_wind_plan.json (classification order)`
- **Missing:** Screens are not shown to cover production models. The selected-wind classifications have never run on real inputs.
- **Plan phase:** P4, P6

### OC-10 — Reproducibility from inputs and version

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=10] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Fresh runs freeze/verify source and hash MESA dependencies; older campaigns need chronology and fresh confirmation."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/README.md`; `Reports/evidence/v1_plan_20261106/review_summary.json`; `docs/operations.md`
- **Registered plan:** none registered
- **Missing:** P0 baseline manifest: source snapshot including the uncommitted solver edits, executable and atomic-data hashes, full launch flags, equation map, two-photon input and the last accepted state. Later, a reproduction from a clean documented installation.
- **Plan phase:** P0, P6

### OC-11 — Useful individual MacBook runs

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=11] (2026-10-22 assessment).
  > "status": "CONDITIONAL", "scope": "Controlled runs fit; full published-track costs unknown."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T031938.396825Z_performance_validation/performance.json`; `docs/performance.md`; `results/nk167_trap/free_r25/driver.log (dense builds about 5.4 h each)`
- **Registered plan:** none registered
- **Missing:** Full-track and converged-closure cost. Accepted progress per wall hour once a method makes progress. Full-run memory, storage and wall time on the Mac.
- **Plan phase:** P2, P6

### OC-12 — External QSRD configuration

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json original_completion_criteria[number=12] (2026-10-22 assessment).
  > "status": "FAIL", "scope": "Deferred by the independent-v1 prerequisite."
- **Artifacts:** none
- **Registered plan:** none registered
- **Missing:** Deliberately deferred until the framework v1 is independently validated. QSRD must stay downstream.
- **Plan phase:** P6

## Prioritized gates

### PG-1 — Spatial and nonlinear convergence

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=1]; same status in docs/validation.md "Release gate matrix".
  > "status": "CONDITIONAL", "remaining": "Production atmosphere/cavity/wind models untested. ..."
- **Artifacts:** `benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T030639.520967Z_robustness_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T031651.519611Z_midpoint_512_refined_time/method_comparison.json`; `benchmarks/results/2026-09-28-closure/20260928T031253.242612Z_time_method_validation/campaign.json`; `benchmarks/results/2026-09-28-verified-spatial/20260928T033503.024409Z_spatial_validation/previous_campaign_comparison.json`; `Reports/Report_13.md`; `Reports/Report_14.md`; `Reports/Report_15.md`; `Reports/Report_16.md`; `Reports/Report_21.md`
- **Registered plan:** `benchmarks/cmfgen_cool_validation_plan.json`; `benchmarks/cmfgen_cool_frequency_plan.json`; `benchmarks/cmfgen_validation_plan.json`
- **Missing:** Released 167-radius closure and its 333-radius and doubled-frequency pairs; production-model convergence.
- **Plan phase:** P1, P2, P3, P6
- **Notes:** Linked original criteria: OC-4

### PG-2 — Production microphysics

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=2]; same status in docs/validation.md "Release gate matrix".
  > "status": "FAIL", "remaining": "Candidate fails geometry and subsonic-convection validity; no validated conservative atmosphere/envelope/wind coupling. ..."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T021532.646459Z_mesa_microphysics_validation/validation.json`; `benchmarks/results/2026-09-28-closure/20260928T030928.629959Z_mesa_hot_current_verification/termination.json`; `benchmarks/results/2026-09-28-closure/20260928T025423.972541Z_robustness_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T030927.417120Z_mesa_data_inventory/inventory.json`; `benchmarks/results/2026-09-28-outer-boundary/20260928T042655.865567Z_atmosphere_reference/validation.json`; `benchmarks/results/2026-09-28-outer-boundary/20260928T042655.119787Z_freeeos_audit/audit.json`; `docs/atmosphere.md`; `docs/self_selected_outer.md`; `Reports/Report_13.md`; `Reports/Report_14.md`; `Reports/Report_15.md`; `Reports/Report_16.md`
- **Registered plan:** `benchmarks/cmfgen_cool_validation_plan.json (requirements: H- sensitivity, relaxed derivatives)`
- **Missing:** Validated cool non-LTE closure; H-minus rate and atom-size sensitivity; relaxed response; production atmosphere/envelope coupling; accounting for the printed-value discrepancies.
- **Plan phase:** P3, P6
- **Notes:** Linked original criteria: OC-1, OC-9

### PG-3 — Inner boundary and cavity

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=3]; same status in docs/validation.md "Release gate matrix".
  > "status": "FAIL", "remaining": "No conservative moving physical cavity reservoir or alternative evolving published inner-flow families. ..."
- **Artifacts:** `docs/equations.md`; `docs/model_families.md`; `Reports/evidence/development_20260928/final_tests/result.json`
- **Registered plan:** `benchmarks/ball2011_track_plan.json (quasi-static cavity only)`
- **Missing:** A conservative moving physical cavity/reservoir in the finite-thermal engine and the evolving inner-flow families. The quasi-static Ball path does not certify these.
- **Plan phase:** P6
- **Notes:** Linked original criteria: OC-2, OC-6

### PG-4 — Wind physics

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=4]; same status in docs/validation.md "Release gate matrix".
  > "status": "FAIL", "remaining": "No combined H/He population/energy/force closure, convective handoff or physical envelope-to-escaping-wind eigenvalue. ..."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T024813.279478Z_fiacconi_wind_validation/validation.json`; `docs/equations.md`; `Reports/evidence/development_20260928/final_tests/result.json`; `benchmarks/results/2026-09-28-spherical/20260928T172742.131636Z_spherical_wind_validation/validation.json`; `Reports/Report_10.md`; `Reports/Report_11.md`; `Reports/Report_12.md`; `Reports/Report_13.md`; `Reports/Report_14.md`; `Reports/Report_15.md`; `Reports/Report_16.md`; `Reports/Report_21.md`
- **Registered plan:** `benchmarks/cmfgen_selected_wind_plan.json`; `benchmarks/self_selected_plan.json`; `benchmarks/self_selected_refinement_plan.json`
- **Missing:** Gates 4-5 on converged inputs; independent Gate 4 evidence verifier; MESA-CMFGEN energy-zero bridge.
- **Plan phase:** P4
- **Notes:** Linked original criteria: OC-8

### PG-5 — Multiple full published evolutions

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=5]; same status in docs/validation.md "Release gate matrix".
  > "status": "FAIL", "remaining": "Zero complete published quasi-star evolutionary tracks reproduced. ..."
- **Artifacts:** `benchmarks/manifest.yaml`; `benchmarks/results/2026-09-28-closure/20260928T031903.825105Z_validation/validation.json`; `benchmarks/results/2026-09-28-closure/20260928T024813.279478Z_fiacconi_wind_validation/validation.json`
- **Registered plan:** `benchmarks/ball2011_track_plan.json`
- **Missing:** Gate 6 repair (Q2 decision, the 46 Msun merge, the 55 Msun stall, C0-C7, the 44 cells) and a second independent complete evolution.
- **Plan phase:** P5, P6
- **Notes:** Linked original criteria: OC-5

### PG-6 — Cross-validation without retuning

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=6]; same status in docs/validation.md "Release gate matrix".
  > "status": "FAIL", "remaining": "No cross-paper evolutionary prediction without retuning has passed. ..."
- **Artifacts:** `src/qstars/data/sources.json`; `benchmarks/manifest.yaml`; `docs/model_families.md`
- **Registered plan:** none registered
- **Missing:** Registered cross-paper evolutionary prediction without retuning. The target is not selected.
- **Plan phase:** P6
- **Notes:** Linked original criteria: OC-5, OC-6

### PG-7 — Numerical robustness

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=7]; same status in docs/validation.md "Release gate matrix".
  > "status": "CONDITIONAL", "remaining": "Full production-family campaigns absent. ..."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T030639.520967Z_robustness_validation/campaign.json`; `benchmarks/results/2026-09-28-closure/20260928T025423.972541Z_robustness_validation/campaign.json`; `Reports/evidence/development_20260928/final_tests/result.json`; `docs/validity.md`; `benchmarks/results/2026-09-28-outer-boundary/20260928T042828.972984Z_outer_final_verification/verification.json`; `Reports/evidence/outer_report10_20260928/hot_mesa_previous_comparison.json`; `Reports/evidence/spherical_20260928/final_tests.log`; `Reports/Report_13.md`; `Reports/Report_14.md`; `Reports/Report_15.md`; `Reports/Report_16.md`
- **Registered plan:** `benchmarks/ball2011_track_plan.json (restart, convergence campaign)`
- **Missing:** Production-family robustness campaigns; P0 restart/acceptance record fixes for the closure driver.
- **Plan phase:** P0, P6
- **Notes:** Linked original criteria: OC-3, OC-4

### PG-8 — Apple Silicon MacBook

- **Recorded status:** CONDITIONAL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=8]; same status in docs/validation.md "Release gate matrix".
  > "status": "CONDITIONAL", "remaining": "Full published-track cost unknown; Linux MESA bridge unvalidated; development high-resolution timings include concurrent workloads. ..."
- **Artifacts:** `benchmarks/results/2026-09-28-closure/20260928T031938.396825Z_performance_validation/performance.json`; `benchmarks/results/2026-09-28-spherical/20260928T172742.131636Z_spherical_wind_validation/validation.json`; `benchmarks/results/2026-09-28-outer-audit/20260928T175917.423111Z_outer_transport/campaign.json`; `benchmarks/results/2026-09-28-outer-audit/20260928T175246.266274Z_outer_discrepancy/campaign.json`; `Reports/Report_13.md`
- **Registered plan:** none registered
- **Missing:** Full-track and closure costs on the Mac.
- **Plan phase:** P2, P6
- **Notes:** Linked original criteria: OC-11

### PG-9 — QSRD readiness after independent v1

- **Recorded status:** FAIL. Source: Reports/evidence/completion_gates.json prioritized_gates[priority=9]; same status in docs/validation.md "Release gate matrix".
  > "status": "FAIL", "remaining": "Deferred until independent framework gates pass; no QSRD science campaign performed. ..."
- **Artifacts:** `docs/model_families.md`; `docs/validity.md`
- **Registered plan:** none registered
- **Missing:** Deferred by design.
- **Plan phase:** P6
- **Notes:** Linked original criteria: OC-12

## v1 chain: closure and Gates 1-6

### G0 — 167-radius T-released cool closure converged

- **Recorded status:** FAIL. Source: Reports/Report_25.md section 9 (Gates); Reports/evidence/completion_gates.json active_request_gates["167-radius T-released cool closure solved"]: "status": "FAIL".
  > | 167-radius T-released closure converged | **FAIL**: stage 1 stopped at a solver limit (§0), stage 2 started |
- **Artifacts:** `results/nk167_trap/fixedT_b35/README.txt`; `results/nk167_trap/free_r25/README.txt`; `results/nk167_trap/free_r25/driver.log`; `results/nk167_trap/free_r25/history.json`; `results/nk167_trap/free_r25/it16_trial_l0.5 (last valid base)`; `Reports/evidence/report24/`; `Reports/evidence/report25/`; `Reports/evidence/report26/`; `Reports/evidence/v1_plan_20261106/review_summary.json`; `Reports/evidence/v1_plan_20261106/residual_map.txt`
- **Registered plan:** `benchmarks/cmfgen_cool_validation_plan.json (closing_phase)`; `benchmarks/cmfgen_validation_plan.json (thresholds)`; `benchmarks/cmfgen_cool_closure_stop_rule.json + amendment_1 (scoped to results/nk167_hot)`; `benchmarks/cmfgen_cool_existence_test.json + amendment_1 (registered existence question)`
- **Missing:** A converged released state. Stage 1 stopped at a solver limit; stage 2 and its dense-window run stopped at iterate 16 with all trials rejected, with about 3.2 % frozen-norm progress in about 70 h of column builds. Needed: stable and component-constrained subspace solves, global coverage of the distributed residual (54.8 % at depths 62-167), and an explicit accepted-state pointer.
- **Plan phase:** P0, P1, P2, P3
- **Notes:** No registered stop rule names results/nk167_trap; the P2 stop rules in the v1 plan are proposals.

### G1 — Post-convergence cleanup and audits

- **Recorded status:** NOT EVALUATED. Source: Reports/Report_25.md section 9 (Gates); Reports/Report_23.md section 8 (Gate table); Reports/evidence/completion_gates.json active_request_gates["Electron balance, global and cumulative radiation energy, last correction, charge, spectral material temperature"]: "status": "NOT_EVALUATED".
  > | 1. Post-convergence cleanup and audits | NOT EVALUATED |
- **Artifacts:** `scripts/post_convergence.sh`; `scripts/trace_cleanup.py`; `scripts/native_correction_check.py`; `Reports/evidence/report23/native_correction/`
- **Registered plan:** `benchmarks/cmfgen_cool_validation_plan.json (closing_phase; floating-point-floor exemption pending_user_approval)`; `benchmarks/cmfgen_validation_plan.json (last correction <= 1e-3 all levels, nuclei/charge <= 1e-6, electron balance <= 1 %, global radiation energy <= 0.1 %)`; `benchmarks/cmfgen_cool_validation_plan.json thresholds (cumulative radiation energy <= 0.1 %)`
- **Missing:** A converged state. post_convergence.sh defaults to build 31 (its child scripts default earlier) while the solve uses build 36, so the accepted executable and full physics flags must be passed explicitly. COMP_EHB should be added to the verify pass. Stored-J two-photon handling and INC_AD/INC_ADV must be verified in the chain.
- **Plan phase:** P3

### G2 — 333 radii (167/333 spatial pair)

- **Recorded status:** NOT REACHED. Source: Reports/Report_25.md section 9 (Gates); Reports/evidence/completion_gates.json active_request_gates["Spatial (167/333) and doubled-frequency resolution"]: "status": "NOT_REACHED".
  > | 2. 333 radii | NOT REACHED |
- **Artifacts:** `scripts/prepare_mesh333.py`; `scripts/transfer_state_nested.py`; `Reports/evidence/report21/mesh333_control.json`; `Reports/evidence/report21/nk333_README.txt`
- **Registered plan:** `benchmarks/cmfgen_cool_validation_plan.json (closing_phase applies_to the 333-radius solve)`; `benchmarks/cmfgen_cool_frequency_plan.json (profile criteria shared by the spatial pair)`; `benchmarks/cmfgen_validation_plan.json (profile_refinement_relative 1 %)`
- **Missing:** A converged 167 state. The clean 333-radius base pass must be rebuilt (Report 23 found stale bootstrap traces), followed by conservative transfer, an independent solve and its audit.
- **Plan phase:** P3

### G3 — Doubled frequency and resolution comparison

- **Recorded status:** NOT REACHED. Source: Reports/Report_25.md section 9 (Gates); Reports/Report_23.md section 8 (Gate table).
  > | 3. Doubled frequency | NOT REACHED |
- **Artifacts:** `scripts/compare_resolution.py`; `Reports/evidence/report23/compare_resolution_hot_frequency_regression.json`
- **Registered plan:** `benchmarks/cmfgen_cool_frequency_plan.json (85,359 frequencies; luminosity and cumulative-heating criteria pending_user_approval)`
- **Missing:** A converged and audited 167-radius state; the doubled-frequency solve and its own audit; the two-photon field on the refined grid; settled luminosity and heating criteria.
- **Plan phase:** P3

### G4 — Self-selected CMFGEN wind (mdot and sonic point)

- **Recorded status:** NOT REACHED. Source: Reports/Report_25.md section 9 (Gates); Reports/Report_23.md section 8 (Gate table) ("NOT REACHED; module, 21 tests, registered plan"); Reports/evidence/completion_gates.json active_request_gates["Self-selected wind, coupled evolution, published evolutionary benchmark"]: "status": "NOT_ATTEMPTED" (stale label; see unresolved items).
  > | 4. Self-selected CMFGEN wind | NOT REACHED |
- **Artifacts:** `src/qstars/cmfgen_wind.py`; `tests/test_cmfgen_wind.py (21 tests)`; `docs/cmfgen_selected_wind.md`; `Reports/evidence/report23/grad_input_audit.json`
- **Registered plan:** `benchmarks/cmfgen_selected_wind_plan.json (registered before any result; four pending_user_approval items)`
- **Missing:** Implemented and tested, not validated. Needed: real inputs from a converged released closure, the outer fixed-point driver with a new closure per hydro update, non-nested state transfer, a restart writer, resolution on the selected structure, and an independent Gate 4 evidence verifier.
- **Plan phase:** P4

### G5 — Node-coupled conservative evolution

- **Recorded status:** NOT REACHED. Source: Reports/Report_25.md section 9 (Gates); Reports/Report_23.md section 8 (Gate table) ("NOT REACHED; interface, 11 tests").
  > | 5. Node-coupled conservative evolution | NOT REACHED |
- **Artifacts:** `src/qstars/wind_node.py`; `tests/test_node_wind.py (11 tests)`; `docs/wind_node.md`; `docs/model_families.md (CMFGEN wind node coupling interface)`
- **Registered plan:** none registered
- **Missing:** Implemented and tested with a synthetic node only; runs always report validated: false. Needed: a validated Gate 4 node, an admissible envelope at its base, the MESA-to-CMFGEN energy-zero bridge, an independent Gate 4 evidence verifier, nodes covering a declared interval, and refinement/restart checks on the coupled model.
- **Plan phase:** P4
- **Notes:** No benchmarks/ registration governs Gate 5; docs/wind_node.md calls its flag list "a contract proposal, pending Gate 4 and user approval".

### G6 — Published evolutionary benchmark: Ball et al. 2011 fiducial quasi-star

- **Recorded status:** FAIL. Source: Reports/Report_25.md section 9 (Gates); Reports/Report_23.md section 8 (Gate table).
  > | 6. Ball et al. 2011 benchmark | **FAIL** (Report 23 §6). It cannot pass under its registration: check Q2 ("no density inversion at row 0") fails at every resolution ... Moving it needs the owner's decision on Q2, not more computation |
- **Artifacts:** `benchmarks/results/2026-10-14-ball-track/`; `Reports/evidence/report23/ball/stage7_commands.log`; `Reports/evidence/report23/ball/smoke_jump/diagnosis_and_verification.md`; `Reports/evidence/report23/ball_final_review.md`; `scripts/validate_ball_track.py`; `tests/test_validate_ball_track.py`; `benchmarks/ball2011_track.yaml`; `benchmarks/ball2011_table1.csv`; `benchmarks/ball2012_thesis_table3_1.csv`
- **Registered plan:** `benchmarks/ball2011_track_plan.json`; `benchmarks/ball2011_track_plan.registration.json`; `benchmarks/ball2011_track_plan_amendment_{1,2,3,4}.json + .registration.json`; `benchmarks/ball2011_track_plan_amendment_5.rejected_draft.json (rejected, not registered)`; `benchmarks/ball2011_track_spec.md`
- **Missing:** Q2 disposition (source audit; any successor registration must be explicitly approved, and the original FAIL is retained). Also needed: a fix for the 46 Msun merge event, classification of the 55 Msun Newton stall, the smoke target, C0-C7, restarts and the 44 comparison cells. A PASS would be one published-assumption replica only.
- **Plan phase:** P5

## Other open criteria

### X-1 — H-minus collisional/rate data sensitivity before production use

- **Recorded status:** OPEN. Source: Reports/Report_23.md section 8 (Gate table).
  > | Other ledger criteria (multiple replicas, cross-paper prediction, H⁻ sensitivity, atom size, moving-cavity conservation) | **OPEN**, outside this chain |
- **Artifacts:** `docs/nlte_cool_closure.md`
- **Registered plan:** `benchmarks/cmfgen_cool_validation_plan.json requirements: "H- guessed collisional data must not be certified without sensitivity or independent validation" (no numeric criterion)`
- **Missing:** No sensitivity run; no registered numeric criterion.
- **Plan phase:** P3

### X-2 — Model-atom size sensitivity (atomic sufficiency)

- **Recorded status:** OPEN. Source: Reports/Report_23.md section 8 (Gate table).
  > see X-1 (same Report 23 section 8 row)
- **Artifacts:** `docs/nlte_backend.md ("Validation sequence and release boundary")`
- **Registered plan:** none registered
- **Missing:** No atom-size run and no registered criterion in benchmarks/.
- **Plan phase:** P3, P6

### X-3 — Release verifier: qstars validate --require-v1 against the complete ledger

- **Recorded status:** exit 2 (v1 not validated). Source: Reports/evidence/completion_gates.json cli_validation.
  > "cli_validation": {"checks": 31, "passed": 29, "failed": 2, "require_v1_exit_code": 2, ...}
- **Artifacts:** `src/qstars/cli.py`; `src/qstars/validation.py`; `Reports/evidence/cool_recovery_20260929/cli_release_validation/20260930T012647.271050Z_validation/validation.json`
- **Registered plan:** none registered
- **Missing:** The verifier does not read completion_gates.json: --require-v1 always exits 2, and its open-gate list is hard-coded in validation.py. A gate-reading verifier is needed before release.
- **Plan phase:** P6

## Unresolved registration decisions (plan section 4)

| Id | Decision | Registration | Recorded state | Plan recommendation | Phase |
|---|---|---|---|---|---|
| D-1 | Ball Q2 | benchmarks/ball2011_track_plan.json qualitative_checks Q2 | claim "a density inversion appears once Teff < 8000 K (paper 3.1)"; test "outward density increase > 1e-3 relative present at rows 1-4 and absent at row 0"; row-0 maxima 9.07e-3 / 5.36e-3 / 2.75e-3 at N = 200 / 400 / 800 | Authorize a source-faithful successor registration only if the audit demonstrates a mismatch; retain the original FAIL. | P5 |
| D-2 | Frequency luminosity and cumulative-heating criteria | benchmarks/cmfgen_cool_frequency_plan.json (luminosity_lab_over_Lscale, cumulative_heating_over_Lscale: pending_user_approval) | pending_user_approval | Retain the proposed 1 % common-grid criteria with the registered normalizations; settle before the campaign. | P3 |
| D-3 | Floating-point-floor exemptions | benchmarks/cmfgen_cool_validation_plan.json closing_phase (exemption: pending_user_approval) | pending_user_approval; floor-limited levels reported as FAIL-at-floor, never exempted | No exemptions in the initial campaign; measure reproducibility and improve the verifier first. | P3 |
| D-4 | Wind mass-loss refinement and the stricter defaults | benchmarks/cmfgen_selected_wind_plan.json pending_user_approval (4 items: dual-opacity Mach-one guard, relaxing ln mdot refinement 1e-3 to 1 %, thermal-relaxation screen, rule-7 enforcement) | pending_user_approval; non-weakening defaults applied meanwhile | Keep 1e-3 on log mdot; keep the dual-opacity and dual thermal-screen defaults pending formal disposition. | P4 |
| D-5 | Second evolutionary target and exact supported families | none (benchmarks/manifest.yaml lists candidates, e.g. hassan_stellar_tracks, quest_stellar_tracks) | not selected | Select after a bounded source/input audit; register every primary observable, tolerance and no-retuning rule before prediction. | P6 |

None of these was decided or changed here.

## Unresolved ledger and record items

1. completion_gates.json active_request_gates["Self-selected wind, coupled evolution, published evolutionary benchmark"] keeps "status": "NOT_ATTEMPTED". That is inaccurate since Report 23: Gate 6 was attempted (FAIL), and Gates 4-5 have implemented components (NOT REACHED). Only its scope text was corrected; the status needs an owner decision (split or relabel).
2. completion_gates.json outer_campaign_completion.outcome keeps "NUMERICAL BLOCKER (constrained re-solve of stage 1 in progress)"; the parenthetical is stale. An outcome_note was added and the value was not changed.
3. completion_gates.json assessment_date_utc (2026-10-22), automated_tests_passed/skipped (341/37) and cli_validation (29/31, last run 2026-09-30) were not re-measured. Report 25 records 357 passing tests without a skip count.
4. completion_gates.json deliverables.latest_release_report still names Reports/Report_16.md. Its meaning (last release-validation report or last report) is ambiguous, so it was left as recorded; latest_report now names Report 25.
5. No registered stop rule governs the results/nk167_trap solves. cmfgen_cool_closure_stop_rule.json is scoped to results/nk167_hot, and the P2 stop rules of the v1 plan are proposals.
6. results/nk167_trap/free_r25/README.txt has no entry for the iterate-16 stop. The stop reason is only in driver.log, and results/ was out of scope for this edit.
7. benchmarks/manifest.yaml status strings were not reviewed or changed (benchmarks/ is out of scope).
