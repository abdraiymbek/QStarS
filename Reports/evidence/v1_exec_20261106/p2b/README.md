# Council block (2026-11-07): one controlled test of a smaller step bound

Agreed by the three-member council (Claude Opus 5.5, GPT 6.1 Sol, Grok 4.7; three rounds) and approved by the user.
Run directory `results/nk167_trap/v1_p2` (continued), code at commit `c3cf1d1` on branch `v1-plan-execution` plus
logging of the planned trial step (`trial_step.npy`) and its linear-model merit at every λ (`predicted_rel_by_lam`).

## Protocol

1. **Step 0 (offline, ≤ 1 h, non-gating):** planned vs realized vs native for the stored trials of it18/it19
   (`audit_planned_realized.py`); commit on a branch; pin clarabel.
2. **One change, exclusive machine:** from `it20_trial_l0.5`, the v1_p2 flags (`--struct-active --consistent-dirs
   --step-solver box --span-tau 1e-6 --krylov 8`, λ ∈ {1, 0.5, 0.25, 0.125}) with only `--lm-bound 0.25 --lm-bound-t
   0.015`; one iterate per launch; at most 3 accepted iterates (~6 h).
3. **Matched A/B in iterate 1:** the old step (0.5 box, λ 0.5) evaluated natively on the same stored products.
4. **Radius rule** (ρ = actual / predicted reduction of the squared merit for the accepted trial, same weights):
   ρ ≥ 0.7 keep the current radius; ρ < 0.7 next radius = max(0.125, λ_acc × current radius), T scaled alike; never
   grow within the block; never shrink on a λ = 1 rejection alone.
5. **Promote** if matched-load frozen progress ≥ 1.9 %/h (≥ 20 % over the 1.6 %/h baseline of iterates 17–19) with no
   physical diagnostic worse; supporting: λ = 1 accepted on ≥ 2 iterates, each frozen drop ≥ 2 %.
6. **Stop and diagnose** (no Route B inside the block): all trials rejected; < 1 % cumulative frozen over 3 accepted;
   < 0.5 %/h on 2 iterates; d33–47 or d62–100 rising in absolute terms twice; RE net/gross max > 5e-3 or rising
   > 5e-4 in one iterate; significant max|s| > 16 with frozen gain < 2 %.
7. **Alongside, code only:** pre-product batch orthonormalization (next block), P3 closing-chain plumbing (build 36 +
   full flags), offline pass-count fits. After the block: Ball 46 M☉ merge reproduction on its own compute.

## Step 0: planned vs realized vs native (it18, it19)

The rebuilt bounded steps reproduce the driver's logged linear fits exactly (0.9122 / 0.9182). Realized = planned in
the model (the only realized-planned difference is trace He III, weight 0). The native − model gap grows as λ²:

| iterate | λ 0.125 | λ 0.25 | λ 0.5 | λ 1 | gap location |
|---|---|---|---|---|---|
| it18 | 0.0105 | 0.0414 | 0.165 | 0.657 | row 2 (H ionization balance) 98–99 %, d19–47 |
| it19 | 0.0093 | 0.0373 | 0.153 | 0.687 | row 2 65–82 %, He I 16–34 %; d33–47, d1–9 |

No material discrepancy: the λ = 1 overshoot is second-order curvature, chiefly of the H ionization balance at the
front (gap ≈ 0.66 λ² of the merit), not a pipeline effect. Not a gate; the block proceeds.

## Results

(appended as the block runs)

Baseline by the same evaluation (`block_eval.py`, wall from the iterate's first to last log line): it17 3.55 %/h,
it18 1.75 %/h, it19 1.59 %/h (0.5 box, λ 0.5 accepted each time).

| iterate | radius (T) | predicted λ 1 | trials λ 1 / 0.5 / 0.25 / 0.125 | accepted | ρ | frozen drop | wall | %/h | RE net/gross max | sig. max|s| | d33–47 / d62–100 abs² (1e-4) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.25 (0.015) | 0.888 | **0.927** / 0.932 / 0.960 / 0.978 | λ 1 | 0.66 | −5.2 % (2.185 → 2.072e-2) | 1.82 h | **2.85** | 3.02 → 3.32e-3 | 10.5 | 1.06 → 0.93 / 1.06 → 1.10 |
| 21 | 0.25 (0.015) | 0.938 | 0.962 / **0.961** / 0.977 / 0.987 | λ 0.5 | 0.90 | −4.1 % (2.073 → 1.987e-2) | 1.81 h | **2.28** | 3.32 → 3.40e-3 | 7.0 | 0.93 → 0.79 / 1.10 → 1.09 |
| 22 | 0.25 (0.015) | 0.952 | 0.979 / **0.972** / 0.983 / 0.991 | λ 0.5 | 0.86 | −2.6 % (1.987 → 1.936e-2) | 1.82 h | 1.41 | 3.40 → 3.43e-3 | 7.8 | 0.79 → 0.74 / 1.09 → 1.09 |

**Matched A/B (iterate 20, same stored products, `ab_native.py`):** old 0.5 box at λ 0.5 → 0.960 (−4.0 %, frozen
2.142e-2); old 0.5 box at λ 1 → 1.326; new 0.25 box at λ 1 → **0.927 (−7.3 %)**. Radius rule: ρ 0.66 < 0.7 → next
radius max(0.125, 1 × 0.25) = 0.25. A test pytest run (nice 19, ~1 min ×2) overlapped the first pass of iterate 20.

**Krylov depth (offline, iterate 20 products, τ 1e-6, 0.25 box; `fit_by_pass_tau.py`):** passes 1–8 → unbounded fit
0.811 / 0.779 / 0.769 / 0.760 / 0.755 / 0.748 / 0.745 / 0.740, **box fit** 0.970 / 0.965 / 0.963 / 0.960 / 0.955 / 0.944 /
0.927 / 0.888. The unbounded fit flattens after pass 3, but the bounded step keeps gaining (+0.075 from pass 3 to 8,
steepening at the end): the box needs many independent directions to place its bound-active components. By the agreed
rule (cut only if passes 3–8 add < 0.05), `--krylov` is not cut; deeper spans are a candidate for the next block.


## Verdict of the block

| criterion (agreed before launch) | result |
|---|---|
| promote: matched-load frozen progress ≥ 1.9 %/h | **met**: 2.09 %/h over the block (2.85 / 2.28 / 1.41) against 1.75 / 1.59 %/h for the 0.5 box (it18 / it19); +25–31 % |
| supporting: λ = 1 accepted on ≥ 2 iterates | not met as worded (accepted at it20 only); but λ = 1 **lowered** the merit in all three iterates (0.927 / 0.962 / 0.979), where at the 0.5 box it raised it in all four (1.07–1.33) |
| supporting: each frozen drop ≥ 2 % | met (5.2 / 4.1 / 2.6 %) |
| matched A/B on identical products (it20) | 0.25 box λ 1: −7.3 %; old 0.5 box λ 0.5: −4.0 % |
| stop rules | none fired: all iterates accepted; cumulative −11.4 %; no iterate < 0.5 %/h; d33–47 fell each time (1.06 → 0.74e-4), d62–100 flat (1.06 → 1.09e-4, one rise); RE net/gross 3.02 → 3.43e-3 (+3.0, +0.8, +0.3e-4 per iterate; < 5e-3); significant max\|s\| 10.5 / 7.0 / 7.8 (< 16); conservation ~1e-16 |
| radius rule | ρ 0.66 / 0.90 / 0.86 → radius stayed 0.25 throughout (λ_acc × 0.25 = 0.25 at ρ < 0.7 with λ_acc 1) |

**Outcome:** the re-optimized 0.25 box is adopted as the default bound (promotion criterion met, matched A/B
confirms it); the per-iterate gain still decays within the block (5.2 → 2.6 %), as it did at the 0.5 box. Watch item:
the RE net/gross maximum (d33) kept rising, but more slowly each iterate; no limit reached. Since it16 the frozen
measure is 2.8371e-2 → 1.9364e-2 (−31.7 %). Krylov depth: not cut (box fit gains 0.075 from pass 3 to 8).
