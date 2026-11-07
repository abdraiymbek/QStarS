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
