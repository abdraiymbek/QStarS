## Stage 7 final fix round: findings verified, fixes applied, runbook

**Stage 7 can start.** I reproduced every reviewer finding that can be reproduced, and none turned out to be wrong. I fixed all the code findings except one, which I rejected (below). The rest change nothing in code: they are recorded as known outcomes or added to the amendment text.

The photosphere solve fix is now an option that a new **Amendment 4** turns on. That means four amendments must be registered in this order: **3 → 4 → 2 → 1**.

No registered run was executed and nothing was committed. Nothing was written under `benchmarks/results` (I checked before and after). The plan sha256 is still `a10b4166…cc0e`, and the target CSVs, the registration record (`2ef6f603…`) and `ball2011_track.yaml` are unchanged.

### What happened to each finding

**First reviewer (plan compliance)**

| Finding | Reproduced? | What I did |
|---|---|---|
| Blocker: the `photosphere.solve` change was on for every run | Yes. With Amendment 3 alone and the registered solve, C1, C0 and C2 all end `relaxation_blocker`. The 3464 Rsun seed fails with a false "Γ_ph ≥ 1" at every N, and at N = 400 every 6000 Rsun step is rejected at the predictor. | The registered solve (`'doubling'`) is back as the default and matches the preserved snapshot line for line. The new scan is `quasistatic.photosphere_bracket='scan'`, which only Amendment 4 turns on. It reaches the seed, the predictor and `fold.py`, and enters the row-0 signature. With the option on, the dry runs reproduce the first round's numbers to the printed digits (0.015022…, 0.0031408…, 0.00065076…, 3.8686e-8). The registered C1 relax now reproduces the preserved 2026-10-14 run exactly: same messages, same 21 rejection records, bitwise-equal state arrays. |
| Major: `relax.json` could be lost | Yes, both injected-error cases from the review. | Any unexpected error now fails only that seed and records the traceback; such a seed never counts as converged. `relax.json` is always written. After an interrupt its status is `control_stop`, and `preserved_failure` no longer accepts a `control_stop` run as amendment evidence. |
| `cells()` turned an enforced-screen row into a FAIL | Yes | Those cells are now "not evaluable" (`passed=None`), so the verdict is INCONCLUSIVE and can never be PASS. |
| Wrong row-0 label for a nonfinite start or a restart segment | Yes | A nonfinite start state is now `numerical_failure`. A restart segment gets its own label, "model_domain_exit at the restart row". |
| Inaccurate Amendment 3 text | Yes, including the file times (21:11:52 versus 21:12:15–29) | Rewritten: the real order of events, with the registered uniqueness test named as the safeguard; the MESA test moved to `tests_asserting_the_registered_defect`; a new `seed_caveats` section (only two seeds work, the two seeds both start above the relaxed radius, the narrow basin, untested configurations); the solve fix moved out to Amendment 4. |
| `row0_record` trusted any `relax.json` | Yes, a stray copy was accepted | The relax directory must now have a `run.json` of kind `relax`. For smoke, track and campaign, its `plan_registration.json` and the sha256s of its amendments (in order) must also match. For registered output it must lie under `benchmarks/results/<date>-ball-track/`. `run.json` records `row0_relax` (directory plus sha256 of `relax.json` and of the checkpoint). |
| F1 under-reported Newton iterations | Yes on MESA (reported 1; now 3+4+1 = 8 and 4+7+1 = 12) | A rescued step now counts the failed Newton, pseudo-time and final Newton iterations, so the registered 0.8 slow-step factor applies to it. |
| Source consistency was reported, not enforced | Yes | Every run directory now gets a `code_closure.json`: a digest of a fixed list of 24 Ball modules (the driver imports all of them up front), plus any qstars module loaded outside that list. `compare` requires one code (that digest, the driver and the MESA library) for the convergence pairs, Richardson, the N = 400/800 fold agreement and the restarts. C0 must also share it with its own row-0 relax run, or the verdict is INCONCLUSIVE. |
| Dry runs already reveal Q2 and restart outcomes | Yes | No change, because those plan sections are immutable. Both are listed below as known outcomes. |

**Second reviewer (MESA)**

| Finding | Reproduced? | What I did |
|---|---|---|
| Q2 fails at row 0 for certain | Yes: 9.07e-3, 5.36e-3, 2.75e-3 at N = 200, 400, 800, against 1e-3 | No change. A dry-run compare returns FAIL from Q2 alone. |
| Row-0 Richardson error for R* is at its limit | Yes: 0.00991 against 0.01 | No change. |
| N = 800 Newton tolerance sits on the roundoff floor | Yes on the current code: 5 of 65 attempts rejected with residuals 1.03–1.91e-12 | Tolerances unchanged. New per-run `rejection_causes` counting (in `track.json`, `smoke.json` and `campaign.json`), including how many Newton failures were at the tolerance floor. |
| Narrow relaxation basin, no spare seed | From the reviewer's preserved logs | Added to `seed_caveats` in Amendment 3. No radius added. |
| F1 is a formality on MESA | The registered smoke run is identical with and without F1 | Stated in the Amendment 1 text. |
| Wall times about double | Yes: 0.21, 0.41 and 1.12 s per accepted step at N = 200, 400, 800 | Used in the runbook. |
| Deep Knudsen screen has about 6% margin | Yes: maximum 0.0955 | Added as a stop condition. |
| A starting guess outside the domain raises ValueError | — | **Rejected.** No registered caller passes a guess, and the error is labelled numerical, which is the safe direction. |
| Energy budget does not check solver accuracy | — | No change; noted for reading C5. |

### Tests
- **Targeted Ball files without MESA:** 106 passed, 3 skipped.
- **Same files with MESA ABI-7:** 109 passed.
- **Full suite, run once:** 324 passed, 37 skipped.
- **`golden_outputs.py`:** 96 files, 0 differences.
- **Against the source as it was before this round:** all 15 new or changed tests fail, each at the defect it tests.

### Files
- **Code:** `src/qstars/photosphere.py`, `quasistatic.py`, `fold.py`, `config.py`, `scripts/validate_ball_track.py`.
- **Tests:** `test_photosphere.py`, `test_quasistatic.py`, `test_validate_ball_track.py`, `test_ball_integration.py`.
- **Amendments:** new `benchmarks/ball2011_track_plan_amendment_4.draft.json`. Text edited in drafts 1, 2 and 3; their settings are unchanged.
- **Evidence:** `Reports/evidence/report23/ball/stage7_final_fix.json`, plus `dryrun_logs_final/`, `dryrun_scripts_final/` and `review_evidence/` in the same folder.

---

## Stage 7 runbook

Run everything from `/Users/abraiymbek/Projects/QStarS` in zsh. Use these arrays, with the same paths in the same order on every command:
```
A3=(--amendment benchmarks/ball2011_track_plan_amendment_3.json)
A34=($A3 --amendment benchmarks/ball2011_track_plan_amendment_4.json)
A342=($A34 --amendment benchmarks/ball2011_track_plan_amendment_2.json)
A3421=($A342 --amendment benchmarks/ball2011_track_plan_amendment_1.json)
```

**How to register amendment n:**
1. Copy the draft to `benchmarks/ball2011_track_plan_amendment_<n>.json`. Leave the `.draft.json` file in place.
2. Set `status` to any text that does not start with "DRAFT", for example `"registered <UTC>"`.
3. Fill `preserved_failing_run` as a list of run-directory paths relative to the repository root, and fill `relax_result` or `smoke_gate_result`. Edit nothing else; the file is hashed after this.
4. Write the registration record next to it:
```
.venv/bin/python -c "import hashlib,json,datetime,pathlib,sys;p=pathlib.Path(sys.argv[1]);r=dict(amendment=str(p),amendment_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),plan_sha256='a10b4166222342acee53da9622885845ff7744914fb38316d8ea4f48fce8cc0e',registered_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'));p.with_name(p.stem+'.registration.json').write_text(json.dumps(r,indent=2)+'\n');print(r)" benchmarks/ball2011_track_plan_amendment_<n>.json
```
5. Check it without running anything:
```
.venv/bin/python -c "import sys;sys.path.insert(0,'scripts');import validate_ball_track as v;print(v.verify_amendment('benchmarks/ball2011_track_plan_amendment_<n>.json',v.verify_registration()['plan_sha256'])[0])"
```

The driver checks that `amendment_sha256` and `plan_sha256` match, that the status is not a draft, that the settings are amendable, and that every cited run is under `benchmarks/results/<date>-ball-track/` and records a failure. A failure means: smoke not passed (and not `control_stop`), relax `relaxation_blocker` or `not_unique`, or a track not `completed` (and not `control_stop`). In the steps below, "row0" is the file named by `row0_checkpoint` in that `relax.json` (expected `row0_6000rsun.npz`).

**Step 0 (before anything):**
- Freeze the code:
  ```
  shasum -a 256 src/qstars/{__init__,cavity,config,constants,controls,engine,eos,fold,io,ledger,mass_transport,mesa_microphysics,mesh,mlt,nonlinear,numerics,photosphere,quasistatic,state,stepcontrol,structure,transport,validity,winds}.py scripts/validate_ball_track.py backends/mesa/build_ball_abi7/libqstars_mesa.dylib > Reports/evidence/report23/ball/stage7_code_freeze.sha256
  ```
  Before every later command, `shasum -a 256 -c --quiet` on that file must print nothing.
- Optionally rerun `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 .venv/bin/python scripts/validate_ball_track.py relax --run C1` (then `C0`, then `C2`). About 4–10 s each. Expected exit 1, `relaxation_blocker`, matching the 2026-10-14 runs. Add these to Amendment 3's citations.

**Step 1.** Register Amendment 3, citing the three 2026-10-14 relax runs (plus any Step 0 reruns). Then:
```
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 .venv/bin/python scripts/validate_ball_track.py $A3 relax --run C1
```
Repeat with `--run C0` and `--run C2`. Expected: exit 1, `relaxation_blocker`, with 3464 `seed_failed` on the "(Gamma_ph >= 1); no Ball surface exists" message. About 11–17 s, 7 s and 85 s.

**Step 2.** Register Amendment 4, citing the three Step 1 runs. Then `… $A34 relax --run C1` (then C0, then C2). Expected: exit 1, `not_unique` (about 1.5e-2, 3.1e-3, 6.5e-4). About 20–25 s, 40–60 s and 200 s.

**Step 3.** Register Amendment 2, citing the three Step 2 runs. Then `… $A342 relax --run C1`. Expected: exit 0, `converged_unique` (about 3.9e-8), in about 45 s.

**Step 4.** Run both smoke gates from the Step 3 row 0:
```
… $A342 smoke --row0 <Step 3 C1 dir>/row0_6000rsun.npz
… $A342 --fallback F1 smoke --row0 <same>
```
Expected for both: exit 1, `termination_trigger` at step 0, 4 of 4 attempts rejected on `controller_metric` (`profile_dlnr`). About 3 s each.

**Step 5.** Register Amendment 1, citing both smoke directories. Then `… $A3421 relax --run X` for C1, C0, C2, C5, C6 and C7. Each must end `converged_unique` with `row0_enforced_screen_errors` empty. About 2 min, 1.5–2 min, 6 min, 1.5–2 min, 2–2.5 min and 1.5–2 min.

**Step 6.** `… $A3421 smoke --row0 <Step 5 C1 dir>/row0_6000rsun.npz`. It must pass (exit 0). About 1,700 steps, roughly 6–10 min.

**Step 7.** Run the campaign, without `--fallback`. C3 and C4 reuse C0's row 0:
```
… $A3421 campaign --row0 C0=<C0 dir>/row0_6000rsun.npz --row0 C1=… --row0 C2=… --row0 C5=… --row0 C6=… --row0 C7=…
```
Rough wall times, assuming about 3,000 steps for C0 (beyond M_BH = 8.6 this is untested):

| Run | Estimate |
|---|---|
| C0 | 20–40 min |
| C1 | 10–15 min |
| C2 | 1–1.5 h |
| C3 | 10–20 min |
| C4 | 40–80 min |
| C5, C6, C7 | 20–40 min each |
| Restarts R1–R3 | a few minutes |
| **Total** | **about 4–7 h** |

The worst case is the 24 h per-track stop rule.

**Step 8 (after the primary verdict, reported only).** Run S1, S2, S4 (only if C0 clamped any cell) and V1–V3 as `relax --run X` and then `track --run X --row0 …`. None of these has a dry run, and V3 will probably fail to build seeds.

**Branches if an outcome differs from the expectation.** Never register an amendment whose defect did not show.
- If Step 1 ends `not_unique`: skip Amendment 4 and have Amendment 2 cite the Step 1 runs.
- If Step 2 ends `converged_unique`: skip Amendment 2.
- If the Step 4 smoke passes: skip Amendment 1 and relax everything with `$A342`.

**Stop and record a numerical blocker if any of these happens:**
1. Any required configuration ends anything other than `converged_unique` with all amendments, or Step 1, 2 or 3 ends `relaxation_blocker` when the expectation was otherwise. Do not add radii or settings.
2. A row 0 has enforced-screen errors (the run ends as a model_domain_exit of row 0; INCONCLUSIVE).
3. The Step 6 smoke gate fails. F1, F2 (Amendment 1) and F3 (the automatic secant landing) are then all used up.
4. A track ends on a stop rule: 24 h, 20,000 steps, or dt < 10 yr for 500 steps.
5. The code-freeze check fails, `source_snapshot_verified` is false, or `loaded_outside_closure` is not empty. Rerun the affected runs.
6. A termination that the fold protocol does not confirm stays `numerical_failure`; never relabel it.

**Not stop rules, but report them as they fall:**
- 12 consecutive C2 retries at residuals around 1e-12.
- An exit on `radiation_knudsen_deep`; report the optical-depth width of the cell at τ = 10.

**Decide before Step 7.** The registered verdict will be **FAIL on Q2 at row 0** no matter what else happens. The restart `must_cross` "a rejected step" will probably be unmet, because the probes had 0 rejections in 100–300 steps.

**Also report these margins:**
- R* Richardson error at row 0: 0.00991 against 0.01.
- Early smoke median dt: 46.8 yr over the first 100 steps and 56.5 yr over 300, against the 50 yr gate.
- C2 tolerance-floor rejections.
- The one-sided, two-seed basin.
- F1 never rescued a step.

If C0 passes 1.1 × 1194 Msun with no fold, the endpoint is already FAIL. Stopping it at that point only saves wall time.