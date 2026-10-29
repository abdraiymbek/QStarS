# Gate 6 synthesized implementation spec: replicating the Ball et al. 2011 fiducial quasi-star track in QStarS

The base is the **numerics** design. Grafted onto it: from **fidelity**, the boundary state taken at the r0 face, quasi-static mechanics, the fold continuation, the STARS opacity switch and the published variant runs; from **minimal**, the ML1 closure that avoids cancellation, the exact piston work, the bookings by functional difference and the per-column tolerance derivations. Seven defects that I verified, several of them shared by all three designs, are fixed. The pre-registration JSON in §12 parses.

## 1. Scores (each criterion out of 10)

| Criterion | fidelity | minimal | numerics |
|---|---|---|---|
| Correctness against the paper and thesis | 8 | 6 | 7 |
| Conservation and ledger | 7 | 8 | 9 |
| Feasibility | 5 | 7 | 7 |
| Testability | 8 | 8 | 9 |
| Honesty of tolerances | 5 | 8 | 7 |
| **Total /50** | **33** | **37** | **39** |

**fidelity.**
- Strengths:
  - Most faithful to STARS: the base state is evaluated at r0 (thesis §2.2.2 says the innermost meshpoint lies on the boundary).
  - Uses the STARS log T = 3.75 opacity switch.
  - Diagnoses the endpoint fold by frozen-entropy continuation.
  - Lists the published variant endpoints.
- Defects:
  - Makes F0 a Newton unknown, which requires differentiating the remap.
  - The J_L proportional shrink remaps half the envelope every step.
  - W_inner omits the face-0 potential term.
  - Solver tolerance 1e-10 cannot meet a 1e-8 cumulative ledger.
  - Its 10%/30% PASS levels rest on QUEST's "within 10%", which refers only to the final BH mass (quest2025.txt:203-206) and to a model using TOV, Ledoux and X = 0.75.

**minimal.**
- Strengths:
  - ML1 constants verified against MESA.
  - The excess identity Y = Γ(Γ+1)/A² is exact (I re-derived it from the cubic).
  - Piston work is complete: p0·dV0 + G·w0·M0·(1/r̄0 − 1/r0).
  - Tolerances are derived per column.
- Defects:
  - r_B lags one step behind (d ln r_B/d ln r0 ≈ −0.6 to −0.85 gives an alternating mismatch that decays slowly).
  - The photospheric p_ext also lags.
  - Reuses the dynamic step with v/c_s unknowns at dt/t_dyn ≈ 1e3 to 1e5.
  - Boundary state taken from cell 0, which is only first order.

**numerics.**
- Strengths:
  - Correctly lists the kernel blockers.
  - Slice-split-merge mesh.
  - Tolerance 1e-12 derived from the ledger budget.
  - Step count from the path length of Table 1.
  - Finds the seed bug (I reproduced it).
  - Best test list.
- Defects:
  - Boundary state from cell 0.
  - The Schur complement of the M0 border is not a fold detector (see C3).
  - Dynamical probes at every epoch are expensive.
  - ρ0 at 30% and T0 at 10% (via 2τ_L) are loose and weakly justified.

## 2. Facts I checked, and corrections

**Confirmed against the sources and the repository**
- Thesis eq 3.15 prints 16π (ball_thesis.txt:1583-1585); the paper's eq 15 prints 4π. engine.py:18-20 `quest2025_eq11` already uses prefactor 16, so the new alias is a rename.
- Thesis eq 2.16, P_g = (2g/3κ)(1 − L/L_Edd), equals 2g/(3κ) − aT_eff⁴/6. The three designs' surface formulas are therefore equivalent.
- `grey` with external τ = 2/3 reproduces T_c⁴ = ¾T_eff⁴(τ_half + 4/3) (transport.py:47-49).
- MESA ML1 constants are ff1 = 1/8, ff2 = 1/2, ff3 = 24, ff4 = 0 (qsrd-runtime/mesa/turb/private/mlt.f90:98-102).
- The thesis eq 2.10 inconsistency holds: its own eqs 2.6-2.9 give (ℓ/H_p)².
- The bridge is ABI 6 (microphysics.f90:22). The paper-opacity routine switches handles at 7000 K (:84). Handle 1 blends at log T 3.75-3.76 (:64-65). The derivatives dkd/dkt are computed but not exported.
- `structure.initialize` crashes for the row-0 geometry at 200 zones with both meshes ("`x` must be strictly increasing", structure.py:21). I reproduced this read-only today.
- Printed variant endpoints exist: ρ ∝ r^−1/2 gives a final 1077 M☉ (paper §5.2, thesis §3.3.2); M_cav = 0 gives 800 M☉ (thesis §3.3.2). The lifetime fit is eq 51, log τ = −0.126 log M* + 7.12.
- The thesis stability check uses the volume-weighted (3γ−4)p; the paper prints 3γ−4.
- Paper §5.4 states that M* decreases by the radiated mass-energy.
- Paper vs thesis tables, as max|ln ratio| per column:

| Column | t | Ṁ | M_cav | L | ρ0 | T0 | T_eff | r_B | R* |
|---|---|---|---|---|---|---|---|---|---|
| Spread | 0.0024 | 0.0085 | 0.0159 | 0.0034 | 0.0351 | 0.0123 | 0.0141 | 0.0297 | 0.0293 |

  Rows 1-4 are nearly identical; almost all of the spread comes from row 0. The thesis Fig. 3.2 caption still says r0 = 1.66 R☉ while its table gives 1.71 R☉. Row 0 is the least reproducible row of Ball's own output.

**Corrections to the designs**

- **C1 (all three).** Every controller measures Lagrangian |Δln T| and |Δln r| over all cells.
  - The base slab swept per step (0.02-1 M☉) far exceeds the mass of the base condensation (about 1e-3 M☉ at row 0). Gas from about 6 r0 lands at r0, with Δln r ≈ 1.8 and Δln T ≈ 0.4.
  - Any 0.5-2% cap therefore forces dt down to about 0.1 yr: the cell-0 drain problem in disguise.
  - Fix: exclude gas within 2× that step's swept mass of the base, and control the base through the boundary observables (r0, ρ_b, T_b).
- **C2 (fidelity, numerics).** `MesaMicrophysics.raw` caches a single key (mesa_microphysics.py:55-62).
  - Evaluating the boundary and photosphere states in separate calls evicts the N-state cache on every colour.
  - Fix: an 8-entry LRU cache (bitwise neutral), or evaluate the concatenated arrays in one call.
- **C3 (fidelity, numerics).** After the mass-budget bookkeeping, M_cav's gravitational coupling cancels. For j ≥ 1, m_j = M0ⁿ + Δ_pred − εX_acc + Σ dm̄.
  - The only dense coupling left is εṀdt. One dense column is enough, versus three in fidelity.
  - A Schur complement on one border cannot detect the fold Ball reports, which is a fold of the whole hydrostatic problem.
  - Fix: monitor sign det J and the smallest pivot of the frozen-entropy Jacobian.
- **C4 (fidelity).** W_inner must include G·w0·m0·(1/r̄0 − 1/r0), derived from the summation-by-parts secant identity.
- **C5 (fidelity).** The run's throughput is about 2.5e56 erg, about 2000× E_int,0 ≈ 1.3e53 erg. Over about 2000 steps a 1e-8 cumulative ledger needs Newton tolerance ≤ 1e-12, as numerics derived.
- **C6 (fidelity).** A full remap every step recomputes T by EOS inversion in every cell, which adds roughly 1e-10 energy noise per step. Copy unchanged cells bitwise.
- **C7 (minimal).** Lagged r_B and p_ext are replaced by implicit rows.

**New analysis that shapes the design**
- **N1. Base inflow.** BE compresses inflowing gas (about 3× in ρ) within one step; for γ = 4/3 the rectangle pdV gives u₂ = 3u₁ instead of 1.44u₁.
  - Because dt is much longer than the convective turnover at the base, this acts as a heat source bounded by about Ṁc_s² ≈ (c_s²/εc²)·L_BH.
  - That is ~1e-4 of L_BH at row 0 and ≤ 1e-6 late in the track. It is negligible, but it does not converge over the campaign's dt range.
  - It is reported every step. Fallback F2 (ALE base) applies if it exceeds 1e-3.
- **N2. Near thermal equilibrium.** t_KH/t_evol ≈ 0.15 at row 0 and about 4e-4 at row 2, so the track is close to a sequence of thermal-equilibrium models.
  - This explains why both of Ball's tables agree on rows 1-4 despite a 3% difference at row 0.
  - It justifies steps of 1e2-1e4 yr.
- **N3. Near-Eddington row-0 photosphere.** Γ_ph ≈ 0.90-0.93 at row 0, so P_g ∝ (1 − Γ_ph) turns a 3% difference in κ_ph into about 30% in P_g.
  - Row-0 T_eff and R* are sensitive; this is consistent with Ball's own 3% spread.
  - Γ_ph ≥ 1 is a hard domain exit, so row 0 might not exist at all with QStarS microphysics (risk R2).
- **N4. Implied atmosphere column.** 4πR²P_ph/g grows to about 1e2 M☉ late in the track (κ_ph ≈ 3e-4) while the published model has no mass above τ = 2/3.
  - QStarS replicates this and reports it as a diagnostic. The P_ph·ΔV work is about 4e-6 of L·dt.

## 3. The synthesized model

### 3.1 Equations mapped onto QStarS

| Ball (B = paper eq, T = thesis eq) | QStarS realization |
|---|---|
| r0 = r_B/b, r_B = 2GM_BH/c_b², c_b² = Γ1p/ρ at r0 (T3.2, T3.4, b = 1; B1-2) | Piston row at face 0: ln r0 − ln(2GM_BH^{n+1}/c_b²) = 0. The b = r_B/r0 option exists only as a config switch. The state (ρ_b, T_b) comes from `cavity.boundary_state(mode='face')`. |
| M_cav = (8π/3)ρ_b r0³ (T3.7, B5) | `State.cavity_mass`. M_cav^{n+1} is evaluated locally inside the residual. |
| M0 = M_BH + M_cav (T3.10, B8) | `State.enclosed = bh_mass + cavity_mass + cumsum(dm)` |
| L_BH = εṀc² = ε′Ṁ_BH c², Ṁ_BH = (1−ε)Ṁ (T3.11-3.12) | `face_luminosity[0] = εc²X_acc/dt` (X_acc = rest mass swallowed this step); ledger `central_generated`, `mass_defect` |
| Ṁ_BH = 16πη(GM_BH)²ρ_b/(ε′Γ1 c_b c²), η = ε = 0.1 (T3.15) | Engine `ball2012_eq3_15`. Global row: M_BH^{n+1} − M_BH^n − dt·A_b·(M_BH^{n+1})² = 0. |
| Hydrostatics (T2.2, T2.22); dP/dr at r0 (B16) | Quasi-static face rows: area_j·(p_{j−1} − p_j)/w_j = G·m_j/(r_j·r̄_j). B16 enters through the half-cell extrapolation. |
| Energy, ε_grav (T2.3, T2.23) | Backward-Euler first law per cell with secant volumes |
| ∇ from Böhm-Vitense MLT, α = 2 (T2.4-2.10, §2.1.1) | `newton_flux` inverse closure: contrast − Y(∇rad(L))·Δln p = 0, ML1 |
| Surface (T2.15-2.16) | `stars_photosphere`: grey L_N; implicit ρ_ph row; p_ext = P_g + aT_eff⁴/3 |
| Mesh (T2.19, Table 2.1; 399 meshpoints) | Q-function targets + split/merge; N = 400 |
| ε mass loss (paper §5.4) | m_j −= εX_acc at fixed r̄, booked as central_potential_work |
| Stability check (thesis §3.2.1) | ⟨(3Γ1 − 4)p⟩_V reported |

### 3.2 One quasi-static step: `quasistatic.bondi_step(old, dt, eos, opacity, config, mode)`

`mode` is one of `evolve`, `relax`, or `land(M_T)`.

**P0. Predict the swept mass.** Δ_pred = dt·Ṁ_rest + ΔM_cav, extrapolated to second order from the last two accepted rates stored in `state.aux` (and checkpointed).

**P1. Slice the base.** Remove the bottom Δ_pred of envelope gas with decimal overlap integrals.
- The cut cell is split at fixed (ρ, u), with r_cut³ = r_k³ + 3fV_k/(4π).
- The slab joins a transit central mass, so the enclosed masses above the cut do not change.
- Untouched cells are copied bitwise.
- E_slab = E(n) − E(post), computed with `energies()`.

**P2. Re-split the base and adapt the surface.**
- Targets: dm*(μ) from Q of state n mapped in mass-above-base μ. The base structure is quasi-steady in μ.
- Split cells above 2× target and merge pairs below 0.5× target.
- Each changed cell is remapped conserving E and V. The EOS-inversion defect is booked to `remap_energy_defect`.

**P3. Newton solve.** Unknowns, 3N+3 in total:
- y_j = ln(r_j/r̄_j) for faces 0..N;
- ln T_k for cells;
- L_j/L_scale for interior faces;
- ln ρ_ph;
- δ (cell-0 sink correction);
- X = ln(M_BH^{n+1}/M_BH^n). In landing steps X is replaced by ln dt, so the landing is exact.

Rows:
- (a) Momentum, j = 1..N: [area_j(p_{j−1} − p_j)/w_j − G·m_j/(r_j·r̄_j)]/g_j, with p_N = p_ext.
  - The masses are end-of-step values: m_0 = M_BH^{n+1} + M_cav^{n+1}, and m_j = m_0 + Σ_{k<j} dm_k with dm_0 = dm̄_0 − δ.
  - Optional flag `inertia=backward_euler` adds w_j(v − v̄)/dt with v = (r − r̄)/dt. It is used for identity tests and the probe.
- (b) Energy, per cell: [dm_k·u_k − Ū′_k + p_k(V_k − V̄_k) + dt(L_{k+1} − L_k)]/(dm̄_k·ū_k).
  - Ū′_0 = dm̄_0·ū_0 − δ·h̄_0 (enthalpy sink, as in `exchange`).
  - L_0 = εc²X_acc/dt when evolving, the law value when relaxing.
  - L_N from the grey relation with τ_half + 2/3.
- (c) ML1 closure at interior faces.
- (d) Piston row.
- (e) δ − (X_acc + (8π/3)ρ_b r0³ − M_cavⁿ − Δ_pred) = 0, scaled by dm̄_0.
- (f) [P_g,EOS(ρ_ph, T_eff) − (2g_N/(3κ_ph) − aT_eff⁴/6)] / (2g_N/(3κ_ph)), with T_eff⁴ = L_N/(4πr_N²σ) and g_N = G·m_N/r_N².
- (g) The global accretion row.

**Boundary state (face mode).**
- p_b = p_0 + G·m_0·dm_0/(8π·r0⁴).
- ln T_b = ln T_0 + ∇_01·ln(p_b/p_0).
- ln ρ_b = ln ρ_0 + (Δln ρ/Δln p)_01·ln(p_b/p_0).
- The EOS at (ρ_b, T_b) then gives Γ1 and c_b.
- Cell mode (cell-0 values) is variant C6.

**Jacobian.** Complex-step `pattern_jacobian` with greedy distance-2 colouring of the declared pattern: about 10 colours, one of them for the dense X (or dt) column. Then SuperLU. Damped Newton (`nonlinear.newton`) to 1e-12 plus one polishing iteration.

**Predictor.** State n mapped in μ (r and T log-linear in μ; r0 scaled to the predicted r_B).

**P4. Predictor-corrector passes.** If |δ| > 0.05·dm_0(post), set Δ_pred ← Δ_pred + δ and repeat P1-P3 from state n, at most 3 passes. Mass is exact at every pass because δ leaves cell 0 implicitly. The step is therefore fully implicit without differentiating the remap. This replaces fidelity's F0 unknown.

**P5.** Bookings (§4), controller metrics (§5), screens (§8) and the ledger check.

## 4. Conservation and ledger for the moving cavity

Ball has no energy equation for the swept gas. QStarS books every transfer exactly: each sub-operator is applied either at fixed geometry or through the secant work identity.

**Bookings per accepted step:**
1. **Slice + δ sink.**
   - E_sw = E_slab + E_δ, where E_δ = −[E(after sink at r̄) − E(post)].
   - The split by mass fraction is a declared labelling choice: X_acc/(X_acc + ΔM_cav) goes to `advected_energy`, the rest to `cavity_swept_energy`.
   - Entropy is split the same way, into `advected_entropy` and `cavity_swept_entropy`.
   - Masses: `cavity_swept_mass += Δ_pred + δ`; `accreted_rest_mass += X_acc`; `mass_defect += εX_acc`.
2. **ε defect at fixed r̄.** `central_potential_work += G·εX_acc·Σ_j w′_j/r̄_j` (exact, because Ω is linear in M0).
3. **Geometry/thermal update at the new masses (secant identity).**
   - `inner_boundary_work += p_0·(4π/3)(r0³ − r̄0³) + G·w_0·m_0·(1/r̄0 − 1/r0)`;
   - `boundary_work −= p_ext·(4π/3)(r_N³ − r̄_N³)`;
   - `radiated_energy += dt·L_N`;
   - `central_generated += dt·L_0` when evolving; `relaxation_injected_energy += dt·L_0` when relaxing.

**Identities.**
- Mass: bh + cav + Σdm + mass_defect − initial = 0.
- Energy: E − E0 + radiated + advected + cavity_swept − central_generated − relaxation_injected − boundary_work − inner_boundary_work − central_potential_work − remap_energy_defect = 0. This holds to the Newton residual.
- Engine: mass_defect·c² = central_generated holds exactly by construction.

**Not claimed.** No energy balance inside the cavity or the BH.

**Diagnostics reported:**
- cavity_swept_energy / ∫L_BH dt (expected 1e-5 to 1e-4);
- (swept − inner work) / (L_BH·dt);
- the inflow heating ratio (N1);
- the implied atmosphere column (N4).

The ε share leaves from the centre. If STARS removed it at the surface instead, the structural effect is negligible because t_KH ≪ t_evol (N2). This is declared.

## 5. Time stepping over 4.23 Myr

**Controller (`stepcontrol.py`), δ_target = 0.005.** The metrics are:
- Lagrangian |Δln T| and |Δln r| per cell and face, excluding μ < 2·F0;
- |Δln r0|, |Δln ρ_b|, |Δln T_b|;
- |Δln M_BH|;
- |ΔM_cav|/M0;
- |Δln L*|, |Δln R*|, |Δln T_eff|.

**Rules:**
- dt_next = dt·clip(0.9·δ_target/δ_actual, 0.5, 1.5);
- no growth for 3 steps after a rejection;
- an extra factor 0.8 if Newton needs more than 8 iterations;
- reject if δ_actual > 2·δ_target, Newton fails, or more than 3 passes are needed; then halve dt (up to 12 retries).

**dt range.** Starts at 1 yr after relaxation; maximum 2e4 yr. Expect about 100 yr early (Ṁ/M_BH = 4.3e-5 per yr), about 470-7300 yr later, and about 1900 accepted steps (r_B limits about 87% of them).

**Landings.** M_BH = 100, 200, 500 and 1000 are hit exactly, to ≤ 1e-6, by solving for ln dt. Full profiles and checkpoints are saved at each landing.

**Termination trigger.** dt < t_dyn = √(R*³/(G·M*)) (Ball's operational criterion), or 12 consecutive retries. Then the protocol in §12 (endpoint) runs:
- frozen-entropy pseudo-arclength continuation in M_BH (`fold.py`, cavity and photosphere BCs active);
- the fold must be within 2% of the terminal M_BH and agree between N = 400 and N = 800 to 2%;
- then the run is labelled `model_domain_exit` ("no quasi-static continuation"); otherwise `numerical_failure`, which is never relabelled.
- An optional BE-inertia probe at dt = 0.05·t_dyn for 100·t_dyn is reported only.

**Fallbacks, in order:**
1. F1: pseudo-transient continuation inside a step.
2. F2: ALE base region (cells at fixed μ with donor-cell advection), introduced by amendment.
3. F3: secant landing.

**Stop rules.** Projected time above 24 h per track, more than 2e4 steps, or dt < 10 yr for 500 consecutive steps before M_BH = 1000 records a numerical blocker.

## 6. Initial model (row 0)

1. **Seed.** An n = 3 loaded polytrope on the small-ξ0 branch for M_BH/M* = 5e-4 (paper eqs 24-33, `reference.ball_loaded`), scaled to R_seed ∈ {2000, 6000} R☉, with fallbacks {1000, 12000}. Ball's own 3030 R☉ is deliberately not used. T comes from the MESA EOS at (ρ, p), and the mesh is built from Q targets.
2. **Relaxation.**
   - M_BH is frozen.
   - L_BH from the law is on, without rest-mass consumption.
   - The cavity and photosphere BCs are on. F0 = ΔM_cav may be negative; returned gas carries the cell-0 state (declared).
   - dt starts at 0.01 yr and doubles up to 1e8 yr.
3. **Converged** when dt ≥ 1e6 yr and two consecutive steps have changes ≤ 1e-10 and |L*/L_BH − 1| ≤ 1e-6.
4. **Uniqueness.** The two seeds must agree to ≤ 1e-6 on every row-0 observable. This is repeated per N (200, 400, 800).
5. The ledger is re-initialized; the relaxed model defines t = 0.

## 7. Microphysics and expected code-to-code differences

**Primary:**
- MESA 24.08.1 eosDT (FreeEOS-dominated in this regime; frac_ideal = 0 on all probed states; first-law defect ≤ 2.4e-4);
- new ABI 7 `qstars_opacity_stars`: OPAL Type 1 above log T 3.76, AF94 below 3.75 (kh_paper(1) at all T), with exported derivatives;
- ML1 with α = 2.

ABI 7 needs a rebuild; the only current ABI-6 build is build_outer_2.

**Expected differences, none of them tunable:**

| Source | STARS | QStarS | Main columns affected |
|---|---|---|---|
| EOS | Eggleton/Pols 1995 | MESA blend | Interior: ≤ 1e-4 in Γ1 at every printed base state. H/He zones: 1-3%. Affects T_eff, R*. |
| Opacity | Eldridge & Tout 2004, nearest value off-table | Same OPAL/AF94 families, MESA interpolation | ≲ 3% on L, Ṁ, t; κ^−6/5 on ρ0 |
| MLT | α = 2, constants unpublished | ML1 | Largest effect on T_eff and R*; S2 measures it |
| Constants | Older L☉ | L☉ = 3.828e33 | ≤ 0.5% |

## 8. Validity screens (declared replication mode)

**Report-only:**
- convective Mach > 1;
- Kn > 0.1 and free-streaming at τ_surface < 10;
- outer-face Kn, replaced by an Eddington-relation check;
- opacity clamp (STARS nearest-value rule; clamped mass must be ≤ 0.01 M☉ at the compared rows);
- density inversion;
- local Eddington factor.

**Not applicable:** Pe < 100 (ML1 includes radiative losses).

**Always enforced:**
- nonfinite values;
- T within 3000 K to 1e8 K;
- MESA first-law, backstop, dilute-gas and conduction screens;
- compactness;
- flow Mach > 0.2;
- P_g,ph ≤ 0;
- Kn and free-streaming at τ ≥ 10;
- geometry and conservation.

**Recording:** `validity_warnings.json` before step 1, per-row flags, `screen_violations.json`, and a `qstars_production_valid` flag on every comparison cell.

## 9. Staged implementation

Stages 0-5 do not touch config.py, io.py, evolution.py, mass_transport.py or winds.py, which the wind-node workflow is editing. Stage 6 waits for that commit. Stage 7 runs only on integrated, committed code.

### Stage 0: pre-registration and records (≈0.5 day; before any Ball run)

**Files:**
- New: /Users/abraiymbek/Projects/QStarS/benchmarks/ball2012_thesis_table3_1.csv (thesis lines 1590-1606)
- New: /Users/abraiymbek/Projects/QStarS/benchmarks/ball2011_track_plan.json (§12)
- Edit: /Users/abraiymbek/Projects/QStarS/src/qstars/data/sources.json — ball2011 entry: thesis 3.15 (16π), 2.15-2.16, 2.19, the eq 2.10 note, the opacity switch
- Edit: /Users/abraiymbek/Projects/QStarS/benchmarks/manifest.yaml — `ball_stellar_tracks` plan and command; `fiacconi_evolution` → "Section 3.3, Figure 6"
- New: /Users/abraiymbek/Projects/QStarS/tests/test_ball_tables.py

**Test U1:**
- eq 3.15 vs Table 1 using ideal-gas+radiation Γ1 at the printed (ρ0, T0): within 2%;
- 4π law off by ≥ 3.5×;
- r_B within 3.5%;
- M_cav from the printed (ρ0, r_B) within 2.5%;
- T_eff from (L, R) within 0.5%.

**Accept:** tests pass; the plan commit precedes every Ball results directory.

### Stage 1: MESA ABI 7 and cache

**Files:**
- Edit: /Users/abraiymbek/Projects/QStarS/backends/mesa/microphysics.f90 — `qstars_api_version = 7`; `qstars_opacity_stars(n,rho,T,X,clamp,kappa,dlnk_dlnrho,dlnk_dlnT,flags,ierr)`. The clamp maps log R to the table edge at fixed T and sets dlnk_dlnrho = 0. Off-table T sets ierr.
- Build: `scripts/build_mesa_microphysics.py --output backends/mesa/build_ball_abi7`
- Edit: /Users/abraiymbek/Projects/QStarS/src/qstars/mesa_microphysics.py — accept ABI 7; `opacity_stars(rho, T, clamp=False)` with complex-step propagation; 8-entry LRU cache in `raw()`.
- New: /Users/abraiymbek/Projects/QStarS/scripts/ball_domain_sweep.py
- New: /Users/abraiymbek/Projects/QStarS/tests/test_mesa_abi7.py

**Tests:**
- bitwise equal to `qstars_opacity_paper` at T ≥ 7000 K and log T ≤ 3.74;
- complex-step vs finite differences ≤ 1e-6;
- error at log R = −7.5 without clamp; clamp flags set;
- an ABI-6 library still loads;
- the LRU cache is bitwise neutral.

**Accept:** the domain sweep (after Stage 0), covering table base states, scaled loaded-polytrope profiles and eq-2.16 photospheres, shows all states inside the EOS domain with frac_ideal = 0 and defect ≤ 1%. The minimum log R is reported (estimate about −6.2).

### Stage 2: ML1 closure

**Files:**
- New: /Users/abraiymbek/Projects/QStarS/src/qstars/mlt.py
  - `ml1()`: Λ = αH_p, ω = Λρκ, a0 = 9/4, A = 4c_p√(ff1·p·Q·ρ)·α·ω/(ff3·a·c·T³), B³ = (A²/a0)(∇rad − ∇ad). Solve the cubic with MESA's closed form plus 2 Newton polishes, branching on real parts. Then Y = Γ(Γ+1)/A² and v = α√(Qp/(8ρ))·Γ/A. On radiative faces, Y = ∇rad − ∇ad, which is continuous across the boundary.
  - `thesis_lossless()` for S2.
  - `forward_luminosity()` for predictors.
- Edit: /Users/abraiymbek/Projects/QStarS/src/qstars/transport.py
  - `mlt_ml1` and `mlt_thesis_lossless` branches in `luminosities` (diagnostics) and in `inverse_transport_residual`: contrast − Y·log1p(Δp/p), with ∇rad = 3κ̄p̄L/(16πacG·m·T̄⁴).
  - `stars_photosphere` thermal branch.
- New: /Users/abraiymbek/Projects/QStarS/tests/test_mlt_ml1.py

**Tests:**
- against `MesaMicrophysics.mlt('ML1', 2)`: ∇ and v ≤ 1e-10 on 1e3 states (∇rad/∇ad from 1+1e-6 to 1e4, ω from 1e-4 to 1e4);
- the efficient and radiative limits;
- the Y identity;
- complex-step vs finite differences;
- a manufactured inverse.

**Accept:** the `mlt_lossless` path is bitwise unchanged.

### Stage 3: core extensions (bitwise-neutral defaults)

**Files:**
- /Users/abraiymbek/Projects/QStarS/src/qstars/state.py — `cavity_mass = 0.`, `photosphere_density = None`, `kinetic_tracked = True`, `aux = {}`; `enclosed` adds cavity_mass (+0.0 is exact); `energies()` drops kinetic energy when not tracked.
- /Users/abraiymbek/Projects/QStarS/src/qstars/ledger.py — `CAVITY_FIELDS` (`cavity_swept_mass`, `cavity_swept_energy`, `cavity_swept_entropy`, `inner_boundary_work`, `remap_energy_defect`, `relaxation_injected_energy`), added only when `cavity=True` and read with `.get(k, 0.)`, so schema-6 checkpoints still load; `initial_mass` includes cavity_mass.
- /Users/abraiymbek/Projects/QStarS/src/qstars/engine.py — alias; complex-safe `rates_from(rho, p, cs, gamma1, bh_mass, config, time)`.
- /Users/abraiymbek/Projects/QStarS/src/qstars/mesh.py — `remap(..., kinetic=True, offset=0., copy_unchanged=False)` returning slab integrals; `stars_q()`: Q = 0.5 ln p + 0.45 ln(T/(T+2e4)) − ln((m/M*)^{2/3}/0.02 + 1) − 0.05 ln(r²/1e18 + 1) (c2 = 0; ct1 → ∞ per the thesis text, so the c5 term is constant); `q_targets()`; `split_merge()`.
- /Users/abraiymbek/Projects/QStarS/src/qstars/nonlinear.py — `colouring(pattern)` (cached) and `pattern_jacobian(fun, x, pattern)`.
- /Users/abraiymbek/Projects/QStarS/src/qstars/structure.py — enforce strictly increasing qm before the Pchip interpolation, only when it is violated.

**Tests:**
- U2: alias bitwise equal to `quest2025_eq11`; `rates_from` derivatives correct.
- U7: slice/remap/split/merge — exact mass and volume; functional-difference booking to 1e-13; offset 0 gives a bitwise copy; merge entropy ≥ 0.
- U9: `pattern_jacobian` equals the dense complex-step Jacobian at N = 16 to 1e-13.
- U12: the full existing suite plus golden `qstars validate --quick` and hot-run checkpoints are byte-identical.
- U15: the Ball geometry initializes at 200 and 400 zones; 96-zone output is unchanged.

**Amendment A1 to Stage 3 / U7 (Stage 3 review, 2026-10-07; spec test only, the registered plan is unchanged).** Pending explicit user confirmation; implemented as the recommended default because the plan's conservation block is immutable.
- "Merge entropy ≥ 0" cannot hold for the real event once total energy is conserved. Splitting or merging re-faces the face-average (trapezoid) potential, and conserving u + potential moves that change, dPhi, into u, so dS = mixing − Σ dPhi/T + O(dPhi²). On the Ball seed every split cools (m/r is convex in mass inside a uniform-density cell) and 89 of 100 pairwise merges cool. Booking dPhi to `remap_energy_defect` instead would cost 7.6e-6 of E_int for the 31 seed splits and 1.3e-4 for merging all 100 seed pairs, far beyond the plan's 1e-12 per event and 1e-8 cumulative.
- U7 now requires: each split or merged cell's T equals eos.temperature(rho, u_mix − dPhi/M) to roundoff; entropy_change − merge_mixing_entropy agrees with −Σ dPhi/T through second order; mixing entropy ≥ 0 (thermodynamics of the overlap shares only). `split_merge` returns `info["potential_reallocation"]` (Σ dPhi over the events, erg).
- Stage 4 must report, per accepted step, the split/merge entropy change and potential_reallocation/(L_BH dt).
- P1/P2: only the cut remainder (new cell 0 of a slice or of returned gas) keeps a fixed (rho, u); every other changed cell conserves E and V. Returned gas cannot exceed `cavity_mass`.
- P2 merge rule (reading of the plan's "merge a pair below 0.5 × target"): a cell with dm < 0.5 × target joins the neighbour with the smaller dm/target; a merged group still below 0.5 absorbs the next cell; a merged group above 2.0 targets (x = Σ dm/target) is cut into 2^k equal-mass parts with overlap shares of volume and internal energy. A slice remainder or returned gas therefore never survives as a sliver.

### Stage 4: cavity, photosphere, step, screens

**Files:**
- New: /Users/abraiymbek/Projects/QStarS/src/qstars/cavity.py — `boundary_state`, `bondi_radius`, `cavity_mass(profile ∈ {r-3/2, r-1/2, none})`, `slice_base`, `bookings`
- New: /Users/abraiymbek/Projects/QStarS/src/qstars/photosphere.py — residual, p_ext, Γ_ph, implied column
- New: /Users/abraiymbek/Projects/QStarS/src/qstars/quasistatic.py — `bondi_step`, as in §3.2
- New: /Users/abraiymbek/Projects/QStarS/src/qstars/stepcontrol.py
- Edit: /Users/abraiymbek/Projects/QStarS/src/qstars/validity.py — `check(..., policy=None, face_luminosity=None)` with a screen registry; `policy=None` is bitwise-current behaviour.
- New tests: /Users/abraiymbek/Projects/QStarS/tests/test_cavity.py, /Users/abraiymbek/Projects/QStarS/tests/test_photosphere.py, /Users/abraiymbek/Projects/QStarS/tests/test_quasistatic.py

**Tests:**
- U5: constant-κ grey atmosphere — exact P_g = (2g/3κ)(1−Γ); Γ ≥ 1 raises DomainError; order ≥ 1.
- U6: second-order face extrapolation on an analytic loaded polytrope.
- U8: each booking equals the `energies()` difference to 1e-13; swept volume vs `outer_interface.moving_flux` ≤ 2e-7.
- U14: screens are report-only in replication mode and enforced otherwise.
- I1: well-balanced — from a relaxed equilibrium with growth frozen, a 1e6 yr step changes the state ≤ 1e-10.
- I2: the quasi-static vs BE-inertia difference scales as (t_dyn/dt)² at dt = 10, 30 and 100·t_dyn.
- I3: with the cavity frozen and inertia on, the step equals `thermal_mechanical_step(backward_euler)` to ≤ 1e-9.
- I5: a 64-zone analytic cavity run over 1e4 yr — ledger ≤ 1e-10 per step, geometry ≤ 1e-8, exact mass.

### Stage 5: driver, relaxation, termination

**Files:**
- `quasistatic.run` (transactional, modelled on `evolution.run`) and `relax`
- New: /Users/abraiymbek/Projects/QStarS/src/qstars/fold.py — frozen-entropy residual, det sign and minimum pivot, continuation, classification
- New: /Users/abraiymbek/Projects/QStarS/scripts/validate_ball_track.py — verifies the plan hash; relax, smoke, track, campaign, comparison; writes `comparison.json` and `campaign.json`
- New tests: /Users/abraiymbek/Projects/QStarS/tests/test_fold.py, /Users/abraiymbek/Projects/QStarS/tests/test_ball_restart.py

**Tests:**
- I4: with an isentropic Γ = 4/3 EOS, the continuation folds at M_BH/M = 0.1015 ± 0.5% (`reference.ball_limit`).
- I6: bitwise restart in an analytic cavity run (stop at step 20, resume to 40) across a slice, a split/merge, a rejection and a landing.

### Stage 6: integration (after the wind-node commit)

**Files:**
- /Users/abraiymbek/Projects/QStarS/src/qstars/config.py
  - New choices: `bondi_cavity`, `stars_photosphere`, `mlt_ml1`, `mlt_thesis_lossless`, `quasi_static`, `ball2012_eq3_15`, `opacity.tables`.
  - New sections: `cavity`, `quasistatic`, `replication`, `relaxation`.
  - Cross-checks: cavity requires quasi_static, newton_flux, wind off and zero feeding.
- /Users/abraiymbek/Projects/QStarS/src/qstars/io.py — checkpoint schema 7 (cavity_mass, photosphere_density, kinetic_tracked, aux JSON; schema ≤ 6 gets defaults); Ball history columns; `photospheric_radius_cm = r[-1]`.
- /Users/abraiymbek/Projects/QStarS/src/qstars/cli.py — dispatch on `mechanics == 'quasi_static'`.
- /Users/abraiymbek/Projects/QStarS/benchmarks/ball2011_track.yaml — the driver generates the variants from the plan.

**Tests:** U13 (config cross-validation), U16 (schema-7 round trip). **Accept:** the existing suite and golden outputs are unchanged.

### Stage 7: registered execution

Order:
1. Domain sweep.
2. Relaxation ×2 seeds per N.
3. Smoke gate.
4. C0.
5. Comparison.
6. C1-C7, S1, S2 (S4 only if any clamping occurred).
7. Restarts R1-R3.
8. Fold classification at N = 400 and 800.
9. V1-V3, only after the primary verdict.

**Accept:** the verdict rules in §12.

### Stage 8: records

- completion_gates.json:
  - `moving_inner_boundaries` true only if U7, U8, I5, I6 and the C0 ledger pass;
  - `published_full_evolution_replicas` incremented only on PASS, labelled "published-assumption replica" if screens fired;
  - priority 5 stays FAIL until there are multiple replicas, and priority 6 (cross-paper prediction) is still needed.
- `validation.py` pending list, the docs/validity.md:58 cavity line, the Report 23 section.

## 10. Convergence/restart campaign and cost

**Runs and thresholds:**
- C0: N = 400, δ = 0.005, tolerance 1e-12. Spatial runs C1/C2 at N = 200/800; temporal C3/C4 at δ = 0.01/0.0025.
- Pass: C0 vs C2 and C0 vs C4 within ≤ 1% on table observables and on mass-weighted L1 profile norms at 20 exterior-mass fractions.
- Richardson error of C0 ≤ min(1%, τ/3) per cell.
- C5 (tolerance 1e-10) ≤ 0.1%.
- C6 (cell-mode boundary state) reported; must converge with N.
- C7 (split/merge thresholds) ≤ 1%.
- Endpoint converged to 2%.
- N = 1600 only if its projected wall time is ≤ 8 h.
- Restarts: 3 points × 100 steps, bitwise.

**Cost:**
- Development: about 2,000-2,500 lines including tests, about 35-50 agent-hours.
- Per step at N = 400: about 1-1.5 s (MESA 13-36 ms per call, about 10 colours, 5-8 iterations).
- Wall time with 3 niced runs beside nk167_hot:

| Item | Time |
|---|---|
| C0 | 1-2 h |
| C2 | 3-4.5 h |
| C4 | 2-4 h |
| Each other run | 1-2 h |
| Campaign total | ≈ 20-30 CPU-h, ≈ 8-12 h wall after C0 |
| Pessimistic (dt ≈ 10 yr) | > 24 h per track, which hits the stop rule |

## 11. Risks and decision points

- **R1. Screens will fire** in the super-Eddington cool surface. A PASS then means a published-assumption replica, not physics within the QStarS domain.
- **R2. Row-0 photosphere may not exist.** Γ_ph ≈ 0.9-0.93 at row 0. If Γ_ph ≥ 1 with MESA opacity, row 0 is a model_domain_exit, which is a FAIL that cannot be tuned away.
- **R3. Newton stiffness** at the H-ionization opacity peak, the density inversion, and the M_BH < 120 transition.
- **R4. ML1 vs the true STARS MLT.** S2 measures the effect. Changing MLT after a run is forbidden; any change is a new plan.
- **R5. Microphysics.** A FAIL with converged numerics is recorded as FAIL. Upgrades (the STARS source, an Eggleton EOS port) require new pre-registration.
- **R6. Merge timing** with the wind-node workflow.
- **R7. CPU contention** with nk167_hot.
- **R8. Gate 6 alone does not satisfy v1** (completion_gates.json:112, :124, :384).

## 12. Pre-registered `benchmarks/ball2011_track_plan.json` (content; file not created)

```json
{
  "plan_id": "ball2011_track_plan",
  "plan_version": 1,
  "gate": "Gate 6 of the QStarS v1 chain: reproduce the Ball, Tout, Zytkow & Eldridge 2011 (MNRAS 414, 2751) fiducial quasi-star evolution, with numerical convergence and restart certification",
  "registration": {
    "status": "registered_before_any_ball_configuration_run",
    "commit_rule": "Commit this file before any execution that uses the Ball configuration (domain sweep, relaxation, smoke, track, campaign). scripts/validate_ball_track.py stores git rev-parse HEAD:benchmarks/ball2011_track_plan.json and the file sha256 in every run directory, and refuses to run if the working-tree file differs from HEAD.",
    "immutable_sections": ["model", "replication_mode", "comparison", "endpoint", "qualitative_checks", "conservation", "restart", "convergence_campaign.thresholds", "verdict_rules"],
    "amendment_rule": "Only numerical-method bug fixes may be amended, each backed by a test that failed before and passes after, recorded as benchmarks/ball2011_track_plan_amendment_<n>.json citing the preserved failing run directory. Tolerances, primary physics, microphysics and targets are never amended. Sensitivity and secondary runs never select or replace a primary result.",
    "preservation": "Every run directory, including failed, rejected and exited runs, is kept under benchmarks/results/<date>-ball-track/."
  },
  "sources": {
    "paper": {"cache": "sources/cache/ball2011.txt", "table": "Table 1, lines 483-492", "termination": "sections 3.2-3.3, lines 550-566: timestep falls below the dynamical time at M_BH = 1194 Msun"},
    "thesis": {"cache": "sources/cache/ball_thesis.txt", "table": "Table 3.1, lines 1590-1606", "equations": "2.15-2.16 surface; 2.19 and Table 2.1 mesh; 3.2-3.4 inner radius; 3.7 cavity mass; 3.10 inner mass; 3.11-3.12 luminosity; 3.15 accretion law with 16 pi"},
    "targets": {"primary": "benchmarks/ball2011_table1.csv", "secondary_reported": "benchmarks/ball2012_thesis_table3_1.csv"}
  },
  "model": {
    "total_mass_msun": 10000.0,
    "initial_bh_mass_msun": 5.0,
    "composition": {"X": 0.7, "Y": 0.3, "Z": 0.0},
    "accretion_law": {"engine_model": "ball2012_eq3_15", "equation": "thesis 3.15: Mdot_BH = 16 pi eta (G M_BH)^2 rho_b / (eps_prime Gamma1_b c_b c^2), eps_prime = eps/(1-eps)", "eta": 0.1, "epsilon": 0.1, "printed_paper_eq15_4pi": "not run; contradicted by the printed table by a factor of 4"},
    "luminosity": "thesis 3.11-3.12: L_BH = eps Mdot c^2 and Mdot_BH = (1-eps) Mdot, where Mdot is the rest-mass flow across r0; eps Mdot leaves the gravitating mass and is booked as mass_defect",
    "inner_radius": {"equation": "thesis 3.2 and 3.4: r0 = r_B / b with r_B = 2 G M_BH / c_b^2 and c_b^2 = Gamma1 p / rho at r0", "b": 1.0, "boundary_state": "face: half-cell hydrostatic extrapolation from cells 0 and 1 to r0, because STARS evaluates the state on the innermost meshpoint (thesis 2.2.2)"},
    "cavity": {"profile": "rho proportional to r^(-3/2) inside r0", "equation": "thesis 3.7: M_cav = (8 pi / 3) rho_b r0^3", "inner_mass": "thesis 3.10: M0 = M_BH + M_cav"},
    "surface": "Outer face at tau = 2/3 with no mass above it (STARS: surface at m = M*). Eddington grey T^4 = (3/4) Teff^4 (tau + 2/3) between the last cell centre and the face. L* = 4 pi R*^2 sigma Teff^4 (thesis 2.15). P_g = (2 g / 3 kappa_ph)(1 - L*/L_Edd) with L_Edd = 4 pi G c M*/kappa_ph (thesis 2.16). P_ph = P_g + a Teff^4 / 3. kappa_ph = kappa(rho_ph, Teff), with rho_ph solved implicitly. P_g <= 0 is a model_domain_exit and is never clipped.",
    "convection": {"closure": "mlt_ml1", "alpha": 2.0, "constants": "MESA 24.08.1 turb/private/mlt.f90 ML1: ff1 = 1/8, ff2 = 1/2, ff3 = 24, ff4 = 0, a0 = 9/4", "criterion": "Schwarzschild", "velocity_cap": "none, as in STARS", "basis": "Thesis 2.1.1 cites Boehm-Vitense 1958 with ALPHA = 2, and QUEST used ML1 with alpha = 2 to match Ball. Printed thesis eq 2.10 lacks the square on l/Hp that its own eqs 2.6-2.9 require; the lossless form is sensitivity run S2 only."},
    "eos": "MESA 24.08.1 eosDT through the QStarS bridge, X = 0.7, Z = 0",
    "opacity": {"tables": "stars2004_switch", "description": "MESA OPAL Type 1 (gs98, Zbase = 0) for log T >= 3.76, Alexander & Ferguson 1994 (lowT_af94_gn93) for log T <= 3.75, MESA blend in between: the Eldridge & Tout 2004 STARS switch (thesis 2.1.2)", "domain": "log R in [-7, 1] and T in [3000, 1e7] K", "off_table": "Replication mode follows the STARS nearest-table-value rule (log R clamped to the table edge at fixed T) and reports it per step as screen opacity_clamp. Production mode exits."},
    "excluded": ["nuclear burning", "neutrino losses", "rotation", "winds", "surface accretion", "general relativity"],
    "initial_model": "The envelope is relaxed to thermal equilibrium with M_BH frozen at 5 Msun, L_BH taken from the accretion law, and the cavity and photosphere boundaries active (paper section 3). t = 0 and the evolution ledger are defined at the relaxed model."
  },
  "replication_mode": {
    "mode": "published_assumption",
    "label": "Published-assumption replica of Ball et al. 2011; QStarS applicability screens violated in the listed intervals",
    "report_only_screens": ["convective Mach > 1 (ML1 velocity, Cox & Giuli 14.86b)", "radiation Knudsen > 0.1 at tau_from_surface < 10", "radiative flux above free streaming at tau_from_surface < 10", "outer-face Knudsen, replaced by an Eddington-relation consistency check", "opacity_clamp (clamped mass must be <= 0.01 Msun at every compared row)", "density inversion (diagnostic)", "maximum local radiative Eddington factor (diagnostic)"],
    "not_applicable": ["Pe < 100 screen of the lossless MLT; ML1 includes radiative losses"],
    "always_enforced": ["nonfinite state", "EOS temperature domain 3000 K to 1e8 K", "MESA first-law defect > 1%", "MESA neutral ideal backstop", "MESA dilute-gas pressure defect > 1%", "omitted conduction > 1%", "Newtonian compactness > 1e-3", "quasi-static flow Mach (|r - rbar|/dt)/c_s > 0.2", "photospheric P_g <= 0, i.e. Gamma_ph >= 1", "radiation Knudsen > 0.1 and free streaming at tau_from_surface >= 10", "geometry error > 1e-8", "conservation tolerance"],
    "recording": "validity_warnings.json lists this policy before step 1. Every history row carries screen values and flags. screen_violations.json gives, per screen, the first (t, M_BH), the number of affected rows and the extreme value. Every comparison cell carries qstars_production_valid."
  },
  "numerics": {
    "mechanics": "quasi_static: no radial inertia, as in the STARS hydrostatic equations (thesis 2.2, 2.22)",
    "time_integration": "fully implicit backward Euler (thesis 2.3)",
    "unknowns": "ln(r/rbar) at faces 0..N, ln T per cell, face luminosity at interior faces, ln rho_ph, delta (cell-0 sink correction), ln(M_BH_new/M_BH_old); landing steps replace the last unknown by ln dt",
    "swept_mass": "Before the solve, slice the predicted swept mass off the base. Outer predictor-corrector passes, at most 3. The residual delta leaves cell 0 implicitly as an enthalpy sink, so mass is exact at every pass.",
    "mesh": "Lagrangian interior. The base is sliced and re-split every step. Cells are split or merged against targets equally spaced in the STARS eq 2.19 function with Table 2.1 constants, c2 = 0 and ct1 -> infinity (thesis text). Hysteresis: split at 2.0 x target, merge a pair below 0.5 x target. Cells that do not change are copied bitwise.",
    "solver": {"method": "damped Newton; complex-step Jacobian coloured from a declared sparsity pattern, with one dense column; SuperLU", "tolerance": 1e-12, "polish_iterations": 1, "roundoff_exit": "accepted if max residual <= 1e-10 and ||dx||_inf <= 1e-13; counted and reported per step"},
    "controller": {
      "delta_target": 0.005,
      "metrics": ["Lagrangian |dln T| per cell and |dln r| per face, excluding gas within 2 x (that step's swept mass) of the base", "|dln r0|, |dln rho_b|, |dln T_b|", "|dln M_BH|", "|dM_cav| / M0", "|dln L*|, |dln R*|, |dln Teff|"],
      "dt_next": "dt * clip(0.9 * delta_target / delta_actual, 0.5, 1.5); no growth for 3 steps after a rejection; extra factor 0.8 if Newton needed more than 8 iterations",
      "reject": "delta_actual > 2 * delta_target, Newton failure, or more than 3 predictor passes; then halve dt, up to 12 retries",
      "dt_initial_years": 1.0,
      "dt_max_years": 20000.0,
      "termination_trigger": "dt < t_dyn = sqrt(R*^3 / (G M*)) (Ball's operational criterion), or 12 consecutive retries at one step"
    },
    "landing": {"bh_masses_msun": [100, 200, 500, 1000], "method": "the landing step fixes M_BH_new at the target and solves for ln dt (dense dt column); secant on dt is the fallback", "tolerance_relative": 1e-6, "outputs": "full profile and checkpoint at every landing"},
    "determinism": "OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = OPENBLAS_NUM_THREADS = 1, recorded in metadata; colouring cached per sparsity pattern"
  },
  "initial_model": {
    "seed": "n = 3 loaded polytrope on the small-xi0 branch for M_BH/M* = 5e-4 (paper eqs 24-33, reference.ball_loaded), scaled to R_seed; T from the MESA EOS at (rho, p)",
    "seed_radii_rsun": [2000, 6000],
    "fallback_seed_radii_rsun": [1000, 12000],
    "relaxation": "M_BH frozen; dt starts at 0.01 yr and doubles up to 1e8 yr; injected L_BH dt is booked as relaxation_injected_energy",
    "converged_when": "dt >= 1e6 yr, and two consecutive accepted steps have max |dln T| and |dln r| <= 1e-10 and |L*/L_BH - 1| <= 1e-6",
    "uniqueness": "two seeds agree on every row-0 observable to <= 1e-6 relative",
    "per_resolution": "repeated for N = 200, 400 and 800",
    "failure": "If no two seeds converge, record a relaxation blocker. It is numerical unless a converging path reaches Gamma_ph >= 1, which is reported as a model_domain_exit of row 0."
  },
  "comparison": {
    "metric": "d = |ln(Q_qstars / Q_ball)| at the landed M_BH, with no interpolation. Row 0 is the relaxed model. t is t(M_BH), measured from t = 0 at the relaxed model.",
    "units": "QStarS constants: MSUN = 1.98847e33 g, RSUN = 6.957e10 cm, LSUN = 3.828e33 erg/s. Unit-convention differences <= 0.5% are inside the allowances.",
    "observables": {"t": "state.time", "Mdot_BH": "eq 3.15 rate at the landed state", "M_cav": "state.cavity_mass", "L*": "luminosity at the photosphere face", "rho0": "rho_b at r0 (face state)", "T0": "T_b at r0", "Teff": "(L*/(4 pi R*^2 sigma))^(1/4)", "r_B": "r0 (= r_B for b = 1)", "R*": "radius of the photosphere face"},
    "cells": {"row_0": ["Mdot_BH", "L*", "rho0", "T0", "Teff", "r_B", "R*", "M_cav < 0.005 Msun"], "rows_1_to_4": ["t", "Mdot_BH", "M_cav", "L*", "rho0", "T0", "Teff", "r_B", "R*"], "count": 44},
    "effective_tolerance": "tau_model + realization_spread + rounding, where rounding = half a unit in the last printed digit of that cell divided by |value|, computed from the CSV string",
    "realization_spread_definition": "max over rows 0-5 of |ln(paper/thesis)| for the column, measured from Table 1 and thesis Table 3.1",
    "allowance_components": {"numerical": 0.01, "opacity_realization": 0.03, "surface_closure_mlt_eos": 0.03, "interior_adiabat": 0.01, "unit_conventions": 0.005},
    "columns": {
      "t": {"tau_model": 0.05, "realization_spread": 0.0024, "basis": "t = integral of dM_BH/Mdot_BH, so it inherits the Mdot allowance"},
      "Mdot_BH": {"tau_model": 0.05, "realization_spread": 0.0085, "basis": "Mdot_BH = L_BH/(eps_prime c^2) with L_BH ~ L* (Ball: 0.999-1.011); L ~ L_Edd at the innermost convective boundary scales as 1/kappa: 0.03 opacity + 0.01 numerical + 0.005 units, rounded up"},
      "L*": {"tau_model": 0.05, "realization_spread": 0.0034, "basis": "same as Mdot_BH"},
      "Teff": {"tau_model": 0.05, "realization_spread": 0.0141, "basis": "0.03 surface closure (ML1 vs unpublished STARS constants; EOS in H/He ionization zones) + 0.01 numerical; Teff is pinned near H recombination and only weakly opacity-dependent"},
      "R*": {"tau_model": 0.10, "realization_spread": 0.0293, "basis": "identity R* = (L/(4 pi sigma Teff^4))^(1/2): 0.5 x 0.04 + 2 x 0.04, added linearly"},
      "T0": {"tau_model": 0.05, "realization_spread": 0.0123, "basis": "radiation-dominated interior with identical ideal-gas+radiation limit (MESA vs ideal Gamma1 agree to 1e-4 at all printed base states); loaded-polytrope scaling T0 ~ kappa^(-2/5) gives 0.012, plus 0.01 adiabat and 0.01 numerical, rounded up"},
      "r_B": {"tau_model": 0.08, "realization_spread": 0.0297, "basis": "r_B ~ M_BH beta / T0 at fixed M_BH: tau_T0 + 3 x 0.01 for beta"},
      "rho0": {"tau_model": 0.15, "realization_spread": 0.0351, "basis": "rho ~ beta T0^3 / (1 - beta): 3 x tau_T0"},
      "M_cav": {"tau_model": 0.15, "realization_spread": 0.0159, "basis": "thesis 3.9: M_cav ~ (M_BH/K)^3 at matched M_BH, with the adiabat allowance set equal to tau_T0 = 0.05"}
    },
    "strict_tier_reported": "0.01 + realization_spread + rounding, with no physics allowance; reported per cell and never required",
    "row_1_note": "M_BH = 100 lies in the interval where radiative zones disappear (M_BH < 120, paper section 3.2). It is flagged as sensitive but its tolerance is not relaxed."
  },
  "endpoint": {
    "ball": {"M_BH_msun": 1194, "t_myr": 4.23, "M_BH_over_initial_Mstar": 0.1194},
    "classification_protocol": [
      "the termination trigger fires",
      "from the last accepted state, continue the frozen-entropy hydrostatic problem in M_BH with the cavity and photosphere boundaries active (pseudo-arclength); a fold is a sign change of det J with the smallest pivot vanishing as a square root",
      "the fold must lie within 2% of the terminal M_BH and agree between N = 400 and N = 800 to 2%",
      "optional dynamical probe, reported only: same residual with backward-Euler inertia, dt = 0.05 t_dyn for 100 t_dyn"
    ],
    "labels": {"fold_confirmed": "model_domain_exit: no quasi-static continuation (Ball's operational criterion)", "otherwise": "numerical_failure, never relabelled"},
    "criteria": {"M_BH_end": "|ln(M_end / 1194)| <= 0.10", "t_end": "|ln(t_end / 4.23 Myr)| <= 0.10"},
    "outcomes": {"PASS": "fold confirmed and both criteria met", "FAIL": "fold confirmed outside tolerance, or the track continues past 1.1 x 1194 Msun with no fold", "INCONCLUSIVE": "termination classified as numerical_failure"},
    "row_5_columns": "reported only"
  },
  "qualitative_checks": [
    {"id": "Q1", "claim": "Mdot_BH between about 1.8e-4 and 3.7e-4 Msun/yr (paper 3.1)", "test": "every accepted step after t = 0 lies in [1.6e-4, 4.0e-4] Msun/yr", "required": true},
    {"id": "Q2", "claim": "a density inversion appears once Teff < 8000 K (paper 3.1)", "test": "outward density increase > 1e-3 relative present at rows 1-4 and absent at row 0", "required": true},
    {"id": "Q3", "claim": "dynamical stability: volume-weighted (3 Gamma1 - 4) p > 0 (thesis 3.2.1; the paper prints 3 gamma - 4)", "test": "integral of (3 Gamma1 - 4) p dV > 0 at rows 0-4 and at the last accepted model", "required": true},
    {"id": "Q4", "claim": "all but the outermost few 10 Msun are convective (paper 3.1)", "test": "mass above the outermost convective boundary < 100 Msun at rows 2-4", "required": true},
    {"id": "Q5", "claim": "envelope adiabatic index about 1.34 (paper 3.1)", "test": "mass-weighted Gamma1 over convective cells in [1.33, 1.36] at rows 1-4", "required": true},
    {"id": "Q6", "claim": "|dR*/dt| at most about 0.1 Rsun/yr (paper 3.2)", "test": "maximum reported", "required": false},
    {"id": "Q7", "claim": "L* close to L_BH (table: 0.999-1.011)", "test": "|L*/L_BH - 1| reported at rows 0-4", "required": false},
    {"id": "Q8", "claim": "Mdot kinks below M_BH ~ 120; a late Mdot maximum coinciding with a Teff minimum and an R* maximum (paper 3.2)", "test": "reported", "required": false}
  ],
  "conservation": {
    "mass_per_step_relative": 1e-14,
    "mass_cumulative_relative": 1e-12,
    "energy_cumulative_over_initial_internal_energy": 1e-8,
    "energy_also_reported_against": ["initial binding energy", "throughput"],
    "engine_identity_relative": 1e-12,
    "geometry": 1e-8,
    "swept_volume_unit_test": 2e-7,
    "reconstruction_defect_per_event_relative": 1e-12,
    "inflow_numerical_heating_over_L_BH": {"reported_every_step": true, "flag_if_above": 1e-3},
    "booking": "Swept gas (base slice plus cell-0 sink) leaves with exactly the energy difference of the envelope functional at fixed geometry: cavity_swept_energy plus advected_energy, split by the mass fractions dM_cav and Mdot dt. The eps share of the rest mass leaves the central mass: mass_defect, and +G eps Mdot dt sum(w/rbar) in central_potential_work. Piston work p0 dV0 + G w0 M0 (1/rbar0 - 1/r0) goes to inner_boundary_work. Photospheric -P_ph dV_N goes to boundary_work. Reconstruction defects go to remap_energy_defect. No energy balance is claimed inside the cavity or the black hole."
  },
  "restart": {
    "points": ["first accepted step with M_BH >= 20", "last accepted step before the M_BH = 500 landing", "first accepted step with M_BH >= 900"],
    "segment": "restart from the checkpoint and run 100 accepted steps, or to termination",
    "requirement": "Bitwise equality with the uninterrupted C0 run for every history row and for the end-of-segment checkpoint: r, T, dm, v, cell_volume, radius_compensation, face_luminosity, bh_mass, cavity_mass, photosphere_density, aux (controller and predictor state), ledger",
    "must_cross": ["a base slice", "a split/merge event", "a rejected step", "a landing (restart point 2)"]
  },
  "convergence_campaign": {
    "runs": {
      "C0": "production: N = 400, delta_target = 0.005, tolerance 1e-12, face boundary state, split/merge 2.0/0.5",
      "C1": "N = 200",
      "C2": "N = 800",
      "C3": "delta_target = 0.01",
      "C4": "delta_target = 0.0025",
      "C5": "tolerance 1e-10",
      "C6": "boundary_state = cell",
      "C7": "split/merge 1.5/0.67"
    },
    "thresholds": {
      "C0_vs_C2_and_C0_vs_C4": 0.01,
      "metrics": "every compared Table 1 observable at the landed rows, plus mass-weighted L1 norms of ln rho, ln T, L/L_BH and ln r at 20 fixed exterior-mass fractions",
      "richardson_error_of_C0": "<= min(0.01, tau_model/3) per compared cell, from C1/C0/C2 and C3/C0/C4",
      "observed_order": "reported; expected >= 0.8 and monotone",
      "C0_vs_C5": 0.001,
      "C0_vs_C7": 0.01,
      "C6": "reported; must converge with N",
      "endpoint_M_BH_between_N400_and_N800": 0.02
    },
    "escalation": "If C0 vs C2 exceeds 1%, run N = 1600 only if its projected wall time is <= 8 h; otherwise the verdict is INCONCLUSIVE for numerical reasons."
  },
  "sensitivity_runs_reported_only": [
    {"id": "S1", "change": "MESA default opacity handle (OPAL Type 2 + Ferguson 2005 low T)"},
    {"id": "S2", "change": "mlt_thesis_lossless: thesis eqs 2.6-2.9, (l/Hp)^2, delta = 1, alpha = 2, no cap"},
    {"id": "S4", "change": "strict opacity domain, run only if C0 clamped any cell"}
  ],
  "secondary_published_checks_after_primary_verdict": [
    {"id": "V1", "change": "rho ~ r^(-1/2) cavity, M_cav = (8 pi/5) rho_b r0^3 (thesis 3.33, paper eq 46)", "target": "final M_BH 1077 Msun", "tolerance_ln": 0.10},
    {"id": "V2", "change": "M_cav = 0 (thesis 3.3.2)", "target": "final M_BH 800 Msun", "tolerance_ln": 0.10},
    {"id": "V3", "change": "M* = 1e5 Msun", "target": "lifetime from paper eq 51: log10 tau = -0.126 log10 M* + 7.12, i.e. 3.09 Myr", "tolerance_ln": 0.10},
    {"rule": "Reported separately. They never change the primary verdict, and each needs a fold-classified endpoint."}
  ],
  "verdict_rules": {
    "PASS": "All 44 cells within effective tolerance in C0; endpoint PASS; required qualitative checks pass; convergence thresholds met; conservation thresholds met at every accepted step of C0-C7; restarts bitwise; plan hash verified",
    "FAIL": "Any cell outside tolerance while convergence is met, an endpoint FAIL, or a failed required qualitative check",
    "INCONCLUSIVE": "Convergence, conservation or restart not met, or an endpoint numerical_failure while all cells pass",
    "reporting": "Cell by cell with no averaging. Each failed cell is listed with its deviation, effective tolerance and strict-tier value, and every screen-violation interval is listed. A replica is labelled published-assumption if any report-only screen fired."
  },
  "stop_rules": {
    "per_track_wall_hours": 24,
    "max_accepted_steps": 20000,
    "dt_floor": "dt < 10 yr for 500 consecutive accepted steps before M_BH = 1000 is a numerical blocker",
    "smoke_gate": "Before C0, an N = 200 run must reach M_BH = 100 with median dt >= 50 yr and <= 20% rejected attempts; otherwise apply the fallbacks in order",
    "fallbacks_in_order": ["F1 pseudo-transient continuation inside a step", "F2 ALE base region (cells at fixed mass above the base with donor-cell advection), introduced by amendment", "F3 secant landing instead of the dt-unknown landing"],
    "outcome_if_exhausted": "a numerical blocker is recorded; it is never relabelled as physics"
  },
  "cost_estimate": {
    "per_step_seconds_N400": "1.0-1.5 (MESA call 13-36 ms per 400 states measured; about 10 colour evaluations per Jacobian; 5-8 Newton iterations)",
    "accepted_steps_C0": "about 1900 (path length of Table 1 at delta = 0.005: r_B dominates), plus about 20% retries",
    "C0_wall_hours_contended": "1-2",
    "campaign_cpu_hours": "20-30",
    "campaign_wall_hours_3_niced_cores": "8-12 after C0",
    "pessimistic": "dt collapse to about 10 yr gives > 24 h per track and triggers the stop rule"
  }
}
```

## 13. Decisions that need the user (recommended default for each)

1. **Replication-mode screen policy.** Report-only screens are convective Mach > 1, Kn and free-streaming at τ < 10, and the STARS opacity clamp (≤ 0.01 M☉). This determines what a Gate 6 PASS means in the v1 ledger.
   - **Default:** approve. Label results "published-assumption replica" and count them in completion_gates.json only under that label.
2. **Tolerance table in §12, before it is committed.**
   - **Default:** commit as written. It lies within the range of the three independent designs and is tighter than fidelity's for every column except r_B and R*, where all values sit close together.
3. **Compute scheduling beside nk167_hot.**
   - **Default:** run the smoke gate and C0 now, niced (about 2-3 h wall). Start the full campaign (about 8-12 h wall on 3 niced cores) automatically only if C0 passes the comparison; otherwise stop and write the Report 23 blocker.
4. **Exact STARS MLT and EOS.** This would mean downloading the public STARS source and pre-registering a new plan.
   - **Default:** no. ML1 stays primary, S2 is reported, and the upgrade is proposed only if failures concentrate in T_eff and R*.
