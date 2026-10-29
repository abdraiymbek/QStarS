**Verdict: NO-GO for Stage 7 beyond the registered `relax` runs.** The driver, thresholds, screens and conservation/restart checks hold up, but the registered relaxation cannot produce a row-0 model on MESA, so smoke, track and campaign have nothing to start from.

I changed nothing in the repository. The plan sha256 is still `a10b4166…cc0e`, `results/` is untouched, no CMFGEN was started and no registered run was executed.

## What I ran
- **Full suite, once** (single thread, nice 15): 299 passed, 35 skipped. That is one more pass than the Stage 6 fix reported; `tests/test_nlte_newton.py` was modified at 20:28, apparently by the concurrent CMFGEN workflow.
- **MESA opt-in Ball tests** (ABI-7 build): 62 passed, 1 skipped.
- **Real driver on the real plan and YAML:** `compare` with no runs gives INCONCLUSIVE; a draft amendment and a non-row-0 checkpoint are refused before any directory is created.
- **Mutation checks on the analytic toy:**
  - A 1e-3 error in the piston booking moves the energy residual by exactly the expected 8.5e-10.
  - A 1e-6 error in `mass_defect` is rejected by the registered 1e-14 budget.
  - Loosening the Newton tolerance to 1e-4 shows up as 4e-9 in the ledger.
  - A 1-ulp change in T, ledger, aux or cavity mass fails `restart_check`.
- **Dry run on MESA, not a registered run** (scratch only, under one minute in total): the driver's own relax path with the C1 mapping relabelled `DRYRUN`, one seed, plus seed-only and microphysics point probes. No step was ever accepted.

## The blocker
The first Stage 7 command will end as `relaxation_blocker` at every resolution.

| Seed radius | Outcome on MESA |
|---|---|
| 2000 R☉ | Cannot be built: seed law luminosity 6.2e8 L☉ gives Γ_ph ≥ 1.55 for every ρ_ph |
| 1000 R☉ | Cannot be built: 3.5e9 L☉, Γ_ph ≥ 9.3 |
| 12000 R☉ | Cannot be built: T_eff = 2713 K, below the 3000 K MESA floor |
| 6000 R☉ | Builds at N = 200, 400, 800 (2–5 s), but every relaxation step is rejected in the predictor at any dt |

- **Why 6000 R☉ fails:** the last cell takes T ≈ 19,700–19,900 K from the EOS at the polytrope (ρ, p). Through the grey relation that implies L_N = 8.2–8.7e8 L☉, about 21 times the law luminosity (4.0e7 L☉), so Γ_ph ≥ 2 and no photosphere exists.
- **Labelling is correct:** the driver records it as a numerical relaxation blocker, not a model domain exit.
- **Why the toys missed it:** constant electron-scattering opacity and reduced coupling keep the toy sub-Eddington.
- **Row 0 itself looks fine:** at the printed Table 1 values (3.48e8 L☉, 3030 R☉) a photosphere exists at ρ_ph = 2.4e-12 with Γ_ph = 0.92.

## Must be fixed first
1. **Seed construction** (`src/qstars/quasistatic.py` `seed()`, lines 197–247; predictor lines 530–534).
   - The surface cell must be consistent with a sub-Eddington photosphere, and at least two seeds must work for the uniqueness test.
   - The plan text fixes "T from the MESA EOS at (rho, p)" and the four radii, so this is a new Amendment 3, and `AMENDABLE` in `scripts/validate_ball_track.py:52` has no seed key.
   - By the continuum estimate, seeds between roughly 2400 and 10,900 R☉ can be built.
2. **`photosphere.solve` gives false "Γ_ph ≥ 1"** (`src/qstars/photosphere.py:54–86`).
   - The bracket scan catches only `ValueError`, so a trial that lands in Γ_ph ≥ 1 raises `DomainError`.
   - With the default guess it fails at (2.195e8 L☉, 3030 R☉) and (3.48e8 L☉, 3500 R☉), where ordinary roots exist at Γ_ph = 0.74 and 0.92.
   - It affects seeds, the first step from a seed, and the fold classification when no photosphere density is stored.
3. **Uniqueness is not enforced downstream.** `row0_for()` (driver lines 285–293) accepts any t = 0 checkpoint with a matching signature. A row-0 from a `not_unique` relax can be tracked and can reach PASS.
4. **Always-enforced screens are not enforced at row 0.** `quasistatic.run` (lines 1318–1321) only logs `enforced_screen_errors` for the start row, and `cells()` looks at report-only flags.
5. **Failing runs can lose their records.** `track()` calls `fold.classify` unguarded (driver lines 839–853), so an exception leaves no `track.json`; `smoke` crashes without `smoke.json` if `history.csv` was never written.
6. **Spec and both amendment drafts exist only in the scratchpad**, which has been purged once already.

## Checked and sound
- **Registration:** every subcommand refuses on a plan, CSV or YAML mismatch.
- **Thresholds:** the 44 cells, endpoint, Q1–Q5, convergence pairs, Richardson, conservation, restart and verdict rules match the plan values.
- **Screen policy:** `validity_warnings.json` is written before the start row, and on the MESA build the action table matches the plan's three lists exactly.
- **Amendments:** limited to base mesh, controller metric and relaxation mesh cycles; drafts are refused.
- **Default behaviour:** unchanged (suite green).

## Decisions for you
- **C6 "must converge with N"** cannot be evaluated, because the plan registers C6 at one resolution only. The driver reports it and does not enforce it.
- **Fallback F1 is skipped** on toy evidence; the plan says to apply fallbacks in order.
- **Momentum rows default to `endpoint`**, not the spec's secant form. It is recorded in every `config.yaml`.
- **Extra rejection causes:** steps are also rejected on the per-step conservation budget and on EOS or solver errors, which the plan's reject list does not name.
- **Amendment evidence:** the driver checks only that a cited run failed, not that it shows that amendment's defect. The seed-failure run must not be cited for Amendments 1 or 2.
- **No stop at 1.1 × 1194 M☉:** a track with no fold runs to the 24 h or 20,000-step stop rule.
- **Restart `must_cross`** needs a rejected step inside a segment; with the amended numerics the toy had none in 60 steps, which would make the verdict INCONCLUSIVE.
- **Source tree must stay frozen** during runs; another workflow is editing the repository now.

## Commands
Prefix every line with `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15`, from `/Users/abraiymbek/Projects/QStarS`.

Runnable now; each takes about 5 s, exits 1, and leaves the preserved failing directory under `benchmarks/results/<UTC date>-ball-track/` that an Amendment 3 must cite:
```
.venv/bin/python scripts/validate_ball_track.py relax --run C1
.venv/bin/python scripts/validate_ball_track.py relax --run C0
.venv/bin/python scripts/validate_ball_track.py relax --run C2
```

After the fixes and registered amendments (`$A` = the same `--amendment <path>` flags, same order and spelling, on every command, because they enter the row-0 signature):
```
.venv/bin/python scripts/validate_ball_track.py $A relax --run C1
.venv/bin/python scripts/validate_ball_track.py $A smoke --row0 <relax C1 dir>/row0_<R>rsun.npz
.venv/bin/python scripts/validate_ball_track.py $A relax --run C0      # also C2, C5, C6, C7
.venv/bin/python scripts/validate_ball_track.py $A track --run C0 --row0 <relax C0 dir>/row0_<R>rsun.npz
.venv/bin/python scripts/validate_ball_track.py $A track --run C3 --row0 <same C0 row0>      # C4 likewise; C1, C2, C5, C6, C7 with their own
.venv/bin/python scripts/validate_ball_track.py $A track --restart <C0 dir>/restart_bh20msun.npz
.venv/bin/python scripts/validate_ball_track.py $A track --restart <C0 dir>/restart_before_landing_500msun.npz
.venv/bin/python scripts/validate_ball_track.py $A track --restart <C0 dir>/restart_bh900msun.npz
.venv/bin/python scripts/validate_ball_track.py $A compare --run C0=<dir> ... --run C7=<dir> --restart R1=<dir> --restart R2=<dir> --restart R3=<dir>
```
`campaign --row0 C0=… --row0 C1=… --row0 C2=… --row0 C5=… --row0 C6=… --row0 C7=…` does the eight tracks, the three restarts and the compare in one sequential process.

**Time cost:** unknown beyond the relax failures. No quasi-static step has yet been accepted on MESA, so the spec's 1–1.5 s per step and 20–30 CPU-hours are untested.

Scratch outputs are in `/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/final_review/`; the dry-run directories are `dryrun_manual/` and `dryrun_driver/`.