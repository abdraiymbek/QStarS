# Synthesis: RE-row heating vs GENCOOL cooling at the H I front (Report 23)

**Verdict.** There is no sign or unit error. The RE row is positive because radiation is pumping an H I n=2 level whose rate equation has never been balanced, and three quarters of that pumping is a two-photon approximation artifact. The electron thermal balance at the same state is cooling, and CMFGEN's own electron-energy-balance equation agrees with GENCOOL to four digits. The closure failure is numerical/formulation; a finite-temperature equilibrium is expected near 7–9 kK but is not yet demonstrated.

Nothing was started or modified in the repo. Paths: `Q=/Users/abraiymbek/Projects/QStarS`, `C=Q/backends/cmfgen/build_31/cur_cmf`, `S=/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/synth` (my check scripts `v1.py`–`v9.py`).

## 1. Established facts (each re-verified by me)

**Sign, units, content of the RE row**
- Trace row 134 is raw `STEQ_T` (`C/new_main/subs/generate_full_matrix_v3.f:186`; `C/new_main/subs/cmf_blkband_v3.f:635,730`).
- Continuum term: `STEQ_T += FQW*(CHI_NOSCAT*RJ-ETA_NOSCAT)` (`C/new_main/cmfgen_sub.f:2427-2429`). Positive means net radiative absorption.
- `CHI_NOSCAT` is frozen after bf+ff+two-photon and before scattering (`C/new_main/mod_subs/comp_opac.f:121-127,166-167`). Lines go only into CHI/ETA (`cmfgen_sub.f:2311-2318`) and enter RE as `-ETAL*ZNET` (`:2555,2588-2591`).
- No collisional term enters `STEQ_T`; collisions go only to `STEQ_T_EHB` (`C/new_main/subs/steq_multi_v10.f:194-196`).
- The Newton convention matches the driver's: "STEQ is +ve when the populations are to decrease" (`C/new_main/subs/solveba_v13.f:195,344`).
- Conversion: erg cm⁻³ s⁻¹ = 4π·1e-10·R_T. A direct integral of `QSTARS_SPECTRAL.bin` reproduces the trace to 0.1% (5.9644e-9 vs 5.9588e-9 at d35).

**GENCOOL and the native electron balance**
- GENCOOL is the electron-pool balance only: bf, ff, collisional and adiabatic terms, no line or two-photon rows (`cmfgen_sub.f:3208-3217`; file `Q/results/gencool_167/pin100_it05/GENCOOL`). The "percentage" is 200·net/Σ|terms| (`C/subs/wrcoolgen_v3.f:178-186`).
- The `COMP_EHB=T` pass `Q/results/gencool_167/pin100_it05_ehb` finished at 20:01 during this analysis. Its `STEQ_VALS` "Electron Energy Balance Equation" (`cmfgen_sub.f:4127-4132`) at d31–40 equals minus GENCOOL net to four digits (d35: −1.1734e-9 vs 1.1734e-9). Its RE row is unchanged (4.7419 at d35).

**Numbers at the fixed-T state** (`Q/results/nk167_pin/theta1.00/it05_trial_l1`)

| d | T [K] | R_T | R_T/l1 | 4π·1e-10·R_T | native EHB | R_3 (n=2) [cm⁻³ s⁻¹] | hν_Lyα·R_3 / (R_T·c) | 2γ share of R_T | 2γ up-rate, OLD / exact | defect D |
|---|---|---|---|---|---|---|---|---|---|---|
| 33 | 10279 | 1.994 | 7.2e-3 | 2.506e-9 | −4.341e-10 | 147.7 | 0.963 | 76.0% | 81.3 / 0.80 | 5.76e-10 |
| 35 | 10387 | 4.742 | 8.75e-3 | 5.959e-9 | −1.173e-9 | 356.9 | 0.979 | 74.9% | 186.6 / 3.09 | 1.41e-9 |
| 37 | 10250 | 10.925 | 6.0e-3 | 1.373e-8 | −2.659e-9 | 806.8 | 0.960 | 73.1% | 401.7 / 12.2 | 3.47e-9 |

- **Energy identity.** 4π·1e-10·R_T = EHB + h·Σ(E1−En)·R_n + D, the form of Hillier's own check (`C/new_main/subs/chk_cool_rates_consist.f:39-42`) plus the two-photon defect D. It closes to 3.5e-4: 5.9567e-9 vs 5.9588e-9 at d35. Without D it is 24% short.
- **Two-photon method.** Every pass sets `COMP_F=T` (`Q/backends/cmfgen/formal_ali_trial.py:11-13`), so `comp_opac.f:123-124` forces `OLD_DEFAULT` (`it05_trial_l1/OUTGEN:2599`). That mode uses exp(−hν_B/kT_e) as the partner-photon occupation (`C/new_main/subs/two/two_phot_opac_v3.f:118-130`). At 4500 Å that is 0.046 against an actual 1.5e-4 at d35. My reconstructed up-rate matches `TWO_PHOT_SUM` (1.5599e-6 vs 1.5601e-6 s⁻¹ per atom).
- **Temperature derivative.** The driver's frozen-J difference dR_T/dlnT is +3.408, +7.383, +14.997 at d33/35/37. The derivative of the OLD_DEFAULT term alone is 3.407, 7.380, 14.985. The native analytic entry is −4e-3 because `C/new_main/subs/two/two_phot_var_opac.f:71-77` has no dT term.
- **n=2 is frozen.** n2(d35) stays between 4.3649e-5 and 4.3754e-5 cm⁻³ over 60 passes spanning 7.7–12.9 kK. The accepted it05 Newton step for n2 at d35 is −9.2e-5.
- **The exact Jacobian shows why.** In the archived products `it05_from_it05_trial_l1/jv00` (direction 6, d35), the n=2 row's fixed-J action is +0.41615 and its radiative response is −0.41549. The net is +6.6e-4, so 99.8% of the frozen-J correction is cancelled by Lyα coupling.
- **Sinks and feed for n=2.** Photoionization from n=2 is 0.2442 cm⁻³ s⁻¹ at d35 (5.58e3 s⁻¹ per atom, `HIPRRR`). The feed into excited levels is 351 cm⁻³ s⁻¹, or 164.5 without two-photon. Balance therefore needs n2 raised by about 675× at d35 (180–720× over d31–37, under 10× at d≥41).
- **Lyman pumping.** The Lyα zone absorbs a net 51.6 cm⁻³ s⁻¹ at d35, 90% of it at one frequency point at +4 Doppler widths. Other Lyman lines add 8.6e-10 erg. My Sobolev photon-supply estimate 4πJ_c|dv/dl|/hc is about 49 cm⁻³ s⁻¹.
- **Cool restart.** At x0 (`Q/results/nk167_hot/it00_base`, 3806 K at d35) GENCOOL net is −4.5e-13 (slight heating) while 4π·1e-10·R_T is +3.19e-9, of which hν·R_3 is 2.88e-9. x0 is a Λ-iterated state (`Q/results/cmfgen_ali_lambda167/closure_control.json`, `IN_ITS`).

## 2. Ranked explanations

1. **n=2 rate equation unbalanced and invisible to the solver (numerical: convergence definition and metric).**
   - For: the identity; R_T·c ≈ 1.02·hν_Lyα·R_3; frozen n2; the exact-product cancellation; native EHB = −GENCOOL.
   - Against: nothing found.
   - This alone produces the sign conflict. RE equals the thermal balance only when every level row is zero.
2. **`OLD_DEFAULT` two-photon absorption (formulation artifact of single-iteration passes).**
   - It supplies 73–76% of R_T and 187 of the 357 cm⁻³ s⁻¹ feeding n=2 at d35, and overestimates the rate by 33–101× at d33–37.
   - It is the only temperature-dependent part of the row, so it explains the growth from 7.9 to 12.3 kK and the absence of a root at frozen populations.
   - It also explains the frozen-J thermal-curve root. That solve (`Q/scripts/thermal_curve.py:55-60`) zeroes the level rows at frozen J, leaving RE = EHB + D. At d35 this is +3.5e-4 of l1 at 10.39 kK and crosses zero near 10.5 kK. The "10.2 kK root" is real Lyα collisional cooling balancing the spurious D.
   - It would persist even with n=2 balanced.
3. **Continuum photons redshifting into trapped Lyman lines (real within the model; excitation, not heat).**
   - It supplies 24–29% of R_T and is temperature-independent.
   - In statistical equilibrium it is disposed of by Balmer photoionization of n=2, so it drives ionization and about 1 eV of photoelectric heating per event, not 10.2 eV.
   - Caveat (speculation): pure Doppler profiles, Rayleigh wings and no metals may misstate the UV supply.
4. **Excluded.** Sign or unit error; state mismatch (the pass's own GENCOOL is identical); `SCL_LN`; adiabatic/advection (5e-13 against 6e-9); non-coherent scattering; quadrature.

**Corrections to the investigators**
- "GENCOOL is insensitive to n2" is only partly true. At balanced n2, Balmer photoheating offsets about 15–25% of net cooling at d33–37 and roughly all of it at d≤30 (my first-order estimate).
- The proposed frozen-J replay with `NOSTIM` is uninformative: with J frozen it can only reproduce arithmetic already done.

## 3. Consequence for the existence question

- The observed failure is numerical. "No root up to 12.3 kK" and "RE heating at 10 kK" are properties of an unconverged n=2 plus the OLD_DEFAULT term, not of the H/He physics.
- The fixed-T "statistical equilibrium" state is not one. n2 is short by a factor of several hundred at the front. Because n=2 photoionization sets the ionization there (0.244 of about 0.5 cm⁻³ s⁻¹ recombination at d35), H II and n_e are wrong too. The existence-test readout (sign of RE/c) is therefore invalid.
- The 10.3 kK front is not an equilibrium. The electron balance cools at 1.2e-9 erg cm⁻³ s⁻¹ at d35. Raising n2 cannot reverse that, since each event adds at most about 1 eV of heating against 10.2 eV of cooling.
- Inference, not demonstrated: the electron balance is heating at 3.8 kK and cooling at ≥7.7 kK at frozen populations. With n2 balanced, pumped Balmer photoheating is about 1e-10 erg cm⁻³ s⁻¹. A root must lie in between.
- Speculation: my rate estimate puts it at about 7–9 kK with H II fraction 0.1–0.3 at d35.

## 4. Most decisive next test

One native pass at the fixed-T state with a consistent two-photon rate, the native electron balance, and n=2 multiplied by its first-order balance factor at d26–46. I prepared the step file but did not run it: `S/n2boost_state/QSTARS_STEP.bin` (factors in `S/n2boost_state/F_by_depth.txt`; only row 3 at 21 depths differs from `it05_trial_l1/QSTARS_STEP.bin`).

```
cd /Users/abraiymbek/Projects/QStarS && .venv/bin/python scripts/verify_cool_state.py \
  /private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/synth/n2boost_state \
  results/gencool_167/pin100_it05_nostim_n2bal \
  --x0 results/cmfgen_ali_lambda167 --exe backends/cmfgen/build_31/cur_cmf/exe/cmfgen_dev.exe \
  --replace results/nk167_adv/transfer/replace_167.txt --free-t --threads 4 \
  --vadat INC_AD T --vadat INC_ADV T --vadat TWO_METH NOSTIM --vadat COMP_EHB T
```

- **Cost.** About 11–20 min wall on 4 threads; the four earlier passes of this kind took 672–1205 s. The wrapper labels every such pass `numerical_failure`; that is routine and the outputs are still written.
- **Why `NOSTIM` works here.** `comp_opac.f:125-126` honours it, and the rate and variation routines follow the same module switch (`two_phot_opac_v3.f:56-58`). Check that OUTGEN says "Using NOSTIM method" and that the `TWO_PHOT_SUM` up-rate is 0.
- **Read.** Trace rows 3–31, 2 and 134; the last "Electron Energy Balance Equation" block of `STEQ_VALS`.

| Quantity | d33 | d35 | d37 | If the frozen-J view were right |
|---|---|---|---|---|
| Σ R_n, n≥2 [cm⁻³ s⁻¹] (now 145 / 351 / 793) | within about ±10 | within about ±25 | within about ±60 | about −1e7 |
| Ion row magnitude (now 1e-3 to 8e-3) | about 64 | about 165 | about 391 | unchanged |
| EHB [erg cm⁻³ s⁻¹] | −3.7e-10 to −3.3e-10 | −1.0e-9 to −0.9e-9 | −2.1e-9 to −1.9e-9 | — |

- I also predict EHB near zero or positive at d≤30 (T ≤ 8.1 kK). That would place the thermal balance of the n2-balanced front near 8 kK before ionization relaxes.
- The RE row is not a clean prediction in this state, because the ionization and n≥3 rows become unbalanced.
- Control without the boost, fixed by existing data: R_T = 0.478 / 1.192 / 2.94 and R_3 = 66 / 170 / 405 at d33/35/37, EHB unchanged. Run it in parallel if a measured control is wanted: same command with state `results/nk167_pin/theta1.00/it05_trial_l1` and destination `results/gencool_167/pin100_it05_nostim`.
- The definitive follow-up is to iterate n2 and ionization to balance at fixed T (about 5–8 such passes), then read EHB at θ=1.0 and at a cooler θ.

## 5. Driver and report interpretations that are wrong

- **"Net heating raises T"** (`Q/scripts/cool_newton_krylov.py:306-315,359-365`). The sign algebra is right, but R_T is net radiative absorption by matter, which equals electron heating only when all level rows vanish. At the front PTC raised T from 7.7 to 10.5 kK while the electron balance was cooling throughout.
- **Positive T diagonal.** The completed frozen-J T diagonal at the front is positive (+7.38 at d35, entirely the OLD_DEFAULT factor). A plain Newton step and the PTC step therefore point in opposite directions there.
- **Trace-level assumptions.** "Numerically irrelevant to the other equations" (`:383-384`) and "trace levels do not feed back on J" (`:50-55`) are false for H I n=2. Its fraction is 2.6e-13 to 8.6e-13, yet it carries 96–98% of the RE residual, sets the Lyα source function and controls front ionization.
- **Convergence test.** n=2 has weight 1e-3 (`:389-400`) and is outside the test (`:442-445,820-823`). "Significant step 0.09 and falling" says nothing about it.
- **Frozen-J metric.** The `precond` merit W·M⁻¹R (`:372-380`) understates the n=2 error by at least 600×. This is Λ-iteration-type false convergence, and raising the weight alone will not fix it. A Lyα-aware preconditioner or a dedicated direction is needed.
- **RE/c is not net/gross heating.** The scale c is the Lyα gross scattering rate (273.4 + 267.2 of 542.1 at d35). The whole electron budget is only 1.7e-3 of c, so "|net/gross| ≤ 1e-3" cannot certify thermal balance near the front.
- **"Lagged two-photon J"** (`:16`; comment at `C/new_main/qs_jv.f90:612-613`) is not what runs. The branch taken is OLD_DEFAULT (`qs_jv.f90:623-628`).
- **Report 23.** In §2.3 the 10.2 kK local root is artifact-set. The §2.5 existence readout is invalid. The §2.4 statement that the two diagnostics "must agree in statistical equilibrium" is correct, but the state is not in statistical equilibrium.