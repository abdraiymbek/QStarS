# RE-row puzzle at the H I front: numbers lens

Read-only analysis; nothing in the repo was changed and no CMFGEN was started. Path prefixes: `SRC` = `/Users/abraiymbek/Projects/QStarS/backends/cmfgen/build_31/cur_cmf`, `Q` = `/Users/abraiymbek/Projects/QStarS`, `S` = `/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad`.

## Verdict

- **No sign or unit error.** R_T > 0 means the gas absorbs more radiative energy than it emits; the driver's convention matches CMFGEN's.
- **R_T and GENCOOL measure different balances.** They have opposite signs because the H I n=2 rate equation is not in statistical equilibrium at the front; n=2 is outside the driver's convergence test and its population has not moved during the campaign.
- **About 75% of R_T at the front is an artefact.** It is CMFGEN's `OLD_DEFAULT` two-photon 1s→2s absorption, which every QStarS pass selects. It grows as exp(-hν_B/kT_e), which is why R_T grows with front temperature.

## 1. Units and sign of the stored residual

- **Trace content.** Row 134 is the raw `STEQ_T`: `STEQ_STORE(:,K)=STEQ(:,K)` is taken right after `GENERATE_FULL_MATRIX_V3`, before the Krylov hook and row scaling (`SRC/new_main/subs/cmf_blkband_v3.f:290,358,490,635`; writer `:776-799`). `STEQ_VEC(NT)=STEQ_T(K)` (`SRC/new_main/subs/generate_full_matrix_v3.f:186`); the electron-heat-balance branch (`:270-273`) is off (`USE_EHB` default F, `SRC/new_main/mod_subs/rd_control_variables.f:980`).
- **Definition.** `STEQ_T` = +∫(χ_noscat J − η_noscat)dν (`SRC/new_main/cmfgen_sub.f:2427-2429`) − line net emission (`:2578-2591`) − adiabatic/advective work (`SRC/new_main/subs/eval_adiabatic_v3.f:234-249`).
- **Units.** erg cm⁻³ s⁻¹ × 10¹⁰/4π (`eval_adiabatic_v3.f:132-138,145-158`; `SRC/new_main/subs/chk_cool_rates_consist.f:39-42`). So erg cm⁻³ s⁻¹ = 1.2566e-9 × R_T.
- **Newton convention.** BA columns are scaled by POPS (`generate_full_matrix_v3.f:395`) and `POPS=POPS*(1-STEQ)`, "+ means decrease" (`SRC/new_main/subs/solveba_v13.f:344,359,377,385`). The driver reads the row unmodified (`Q/src/qstars/nlte_newton.py:12-13,25-39`; `Q/scripts/cool_newton_krylov.py:306-317,359-365`).
- **Unit test.** `Q/tests/test_nlte_newton.py:218-243` checks only the algebra of the added PTC term; it assumes rather than tests that R_T > 0 is heating.
- **GENCOOL.** "Net Cooling" is the electron thermal balance only: bound-free + free-free + collisional (+ adiabatic), positive = cooling (`SRC/subs/wrcoolgen_v3.f:94-101,181-186`; `SRC/subs/wr_ad_cool.f:28`).

The two are linked by an energy identity, where R_j are the rate-row residuals and e_j the excitation energies:

4π·10⁻¹⁰·R_T = −GENCOOL_net − adiabatic + h·Σ_j e_j R_j + M_2γ

Row 2 is the ionization equation and row 32 number conservation (`generate_full_matrix_v3.f:263,350-352`; `QSTARS_REPLACE_TRACE`). M_2γ is the two-photon term of section 4.

With M_2γ omitted, the identity closes to ≤1e-6 at depths 5–15, which confirms sign, units and row meaning. With M_2γ included it closes to ≤1e-3 at depths 26–46 in the four states analysed (`S/ident_pin.txt`, `S/budget.txt`, `S/budget_T.txt`).

## 2. Magnitudes at depths 30–41

State `Q/results/nk167_pin/theta1.00/it05_trial_l1`. The l1 scale is from the raw native BA blocks; the driver's `ptc_c` uses completed blocks and differs by about 1%. The pass's own GENCOOL agrees with the new `Q/results/gencool_167/pin100_it05/GENCOOL` (finished 19:33) to the fourth digit.

| d | T [kK] | HI/H | R_T raw | l1 | R_T/l1 | 4πe-10·R_T | GENCOOL net | −GENCOOL as raw R_T | R_3 (n=2) | R_3/l1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 8.11 | 0.985 | 0.364 | 845 | 4.3e-4 | +4.57e-10 | +1.06e-11 | −0.008 | 24.5 | 3.8e-4 |
| 31 | 8.92 | 0.990 | 0.668 | 411 | 1.6e-3 | +8.39e-10 | +3.93e-11 | −0.031 | 42.4 | 1.3e-3 |
| 32 | 9.82 | 0.990 | 1.200 | 336 | 3.6e-3 | +1.51e-9 | +1.76e-10 | −0.140 | 83.6 | 3.2e-3 |
| 33 | 10.28 | 0.990 | 1.994 | 275 | 7.2e-3 | +2.51e-9 | +4.34e-10 | −0.345 | 147.7 | 6.9e-3 |
| 34 | 10.33 | 0.989 | 3.059 | 386 | 7.9e-3 | +3.84e-9 | +7.07e-10 | −0.563 | 228.3 | 7.6e-3 |
| 35 | 10.39 | 0.988 | 4.742 | 542 | 8.7e-3 | +5.96e-9 | +1.17e-9 | −0.934 | 356.9 | 8.4e-3 |
| 36 | 10.35 | 0.986 | 7.181 | 991 | 7.2e-3 | +9.02e-9 | +1.82e-9 | −1.44 | 539.6 | 7.0e-3 |
| 37 | 10.25 | 0.984 | 10.93 | 1811 | 6.0e-3 | +1.37e-8 | +2.66e-9 | −2.12 | 806.8 | 5.7e-3 |
| 38 | 9.73 | 0.980 | 15.57 | 4774 | 3.3e-3 | +1.96e-8 | +2.59e-9 | −2.06 | 1074 | 2.9e-3 |
| 39 | 9.23 | 0.972 | 22.69 | 12621 | 1.8e-3 | +2.85e-8 | +2.66e-9 | −2.12 | 1472 | 1.5e-3 |
| 40 | 9.06 | 0.953 | 33.49 | 42255 | 7.9e-4 | +4.21e-8 | +4.66e-9 | −3.71 | 2314 | 7.1e-4 |
| 41 | 8.42 | 0.924 | 48.11 | 141610 | 3.4e-4 | +6.05e-8 | +3.72e-9 | −2.96 | 2902 | 2.7e-4 |

- **Size and sign.** At depths 33–37 the converted R_T is +5.0 to +5.8 times GENCOOL's net cooling, with opposite sign.
- **Hot it28.** `Q/results/nk167_hot/it28_trial_l1` is the same within 3%: R_T/l1 = 7.1e-3, 7.7e-3, 8.6e-3, 7.1e-3, 5.9e-3 at depths 33–37, and GENCOOL net 3.86e-10 to 2.42e-9. Its parsed GENCOOL values match `Q/results/gencool_167/hot_it28/GENCOOL` to the printed digits (`S/hot28.txt`).
- **What l1 is.** The row scale is the Lyα gross rate: at depth 35, dR_T/dln n1 = +273.4 and dR_T/dln n2 = −267.2 out of l1 = 542. R_T/l1 and R_3/l1 are therefore the same Lyα-referenced net/gross ratio.

## 3. Hydrogen rate rows: statistical equilibrium does not hold for n=2

- **n=2 (row 3).** R_3 = +148 to +807 cm⁻³ s⁻¹ at depths 33–37, 0.6–0.8% of its row scale. In energy, hν_Lyα·R_3 is 96–98% of 4πe-10·R_T there (5.83e-9 of 5.96e-9 at depth 35).
- **Other rows at depths 33–37.**
  - Ionization equation (row 2): |R/l1| ≤ 6e-4 — balanced.
  - n=3 (row 4): −2.0e-2 to −2.4e-2.
  - Rows 5–31: up to 3e-2.
  - H⁻ (row 1): +0.10 to +0.12.
  - These are energetically negligible next to n=2 (−1.1e-10 erg at depth 35).
- **n=2 is frozen.** n2(depth 35) is the same to 0.2% across the campaign while T went from 8.87 to 12.32 kK:

| State | T(d35) [kK] | n2(d35) [cm⁻³] |
|---|---|---|
| hot it04 | 8.87 | 4.365e-5 |
| hot it28 | 10.34 | 4.371e-5 |
| pin it05 | 10.39 | 4.375e-5 |
| theta1.15 it02 | 12.32 | 4.372e-5 |

- **Fixed-T iterates.** Over it01–it05, R_3(d35) went 354.0 → 356.9, i.e. it is not falling (`S/hist_pin.txt`). The it04 Newton step for n2 at depths 33–37 is −2.3e-4 to −3.4e-4, against 1.2–1.7e-2 needed even at frozen J.
- **Why the driver misses it.** The n=2 fraction is 2.6e-13 to 8.6e-13 at depths 33–37, below `sig_frac` = 1e-10 until depth 43. It is excluded from the "significant step" (`cool_newton_krylov.py:820-823,442-445`) and weighted 1e-3 in the least squares (`:389-400`). So "significant step 0.09" says nothing about the level that carries the energy; the all-level step is 0.75–0.90 (`history.json`).
- **Consequence.** GENCOOL's collisional cooling (about 72 excitations cm⁻³ s⁻¹ × 10.2 eV at depth 35) is real at this state. It is invisible to the RE row until n=2 is allowed to radiate that energy.

## 4. Two-photon term: about 75% of R_T, and why it grows with T

- **Which method runs.** Every pass sets `COMP_F=T` (`Q/backends/cmfgen/formal_ali_trial.py:11-13`). With `TWO_METH` at its default `USE_RAD` (`rd_control_variables.f:728-729`), `SRC/new_main/mod_subs/comp_opac.f:123-124` then forces `OLD_DEFAULT`. OUTGEN confirms: "unable to open EDDFACTOR_INFO -- will compute new f … Using OLD_DEFAULT method for two-photon decay" (`it05_trial_l1/OUTGEN:2590-2599`).
- **What it does.** `OLD_DEFAULT` gives the companion photon the occupation exp(−hν_B/kT_e) instead of the actual J (`SRC/new_main/subs/two/two_phot_opac_v3.f:118-131`, `two_phot_rate.f:68-76`). The exact form uses J for both photons (`two_phot_rate.f:77-84`).
- **Reproduction.** From the pass's stored `QSTARS_RJ.bin` I reproduce CMFGEN's `TWO_PHOT_SUM` up-rate to 4 digits (1.5599e-6 vs 1.560e-6 s⁻¹ per atom at depth 35).

Budget at the pin state, erg cm⁻³ s⁻¹ unless noted. "RE at SE" is what the row would read with all R_j = 0 at the same T, n1, n_e and J.

| d | 4πe-10·R_T | −GENCOOL | h·Σe_jR_j | M_2γ | 2γ up, OLD [cm⁻³ s⁻¹] | 2γ up, exact [cm⁻³ s⁻¹] | collisional ≤ [cm⁻³ s⁻¹] | 2γ share of R_T | RE at SE, OLD | RE at SE, exact |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 4.57e-10 | −1.06e-11 | 3.48e-10 | 1.20e-10 | 11.8 | 0.11 | 0.65 | 0.68 | +1.1e-10 | −1.1e-11 |
| 33 | 2.51e-9 | −4.34e-10 | 2.36e-9 | 5.76e-10 | 81.3 | 0.80 | 26.5 | 0.76 | +1.4e-10 | −4.3e-10 |
| 35 | 5.96e-9 | −1.17e-9 | 5.72e-9 | 1.41e-9 | 186.6 | 3.09 | 71.7 | 0.75 | +2.4e-10 | −1.17e-9 |
| 37 | 1.37e-8 | −2.66e-9 | 1.29e-8 | 3.47e-9 | 401.7 | 12.2 | 162 | 0.73 | +8.1e-10 | −2.66e-9 |
| 39 | 2.85e-8 | −2.66e-9 | 2.33e-8 | 7.90e-9 | 713 | 48 | 163 | 0.69 | +5.2e-9 | −2.66e-9 |
| 41 | 6.05e-8 | −3.72e-9 | 4.72e-8 | 1.69e-8 | 1353 | 184 | 237 | 0.65 | +1.3e-8 | −3.72e-9 |

- **Overestimate.** `OLD_DEFAULT` gives 60 times the exact rate at depth 35 (33–140 times over depths 28–37).
- **Energy inconsistency.** The RE row charges the absorbed photon's mean energy (0.73 hν₀) while the rate equation charges hν₀/2. That difference is M_2γ; it is larger than GENCOOL's net cooling at every depth 26–44.
- **T-derivative.** The driver's own frozen-J central differences (`Q/results/nk167_pin/theta1.00/cols_it05_trial_l1/col_T_p`, `col_T_m`) give dR_T/dlnT = +3.408, +7.383, +15.00 at depths 33, 35, 37. The T-derivative of the `OLD_DEFAULT` term alone predicts 3.407, 7.380, 14.99. The native analytic BA entry is −4.2e-3. At fixed J and populations the RE row rises as roughly T^1.6, entirely from this term, so it has no root in T.

Depth 35 across states:

| T [kK] | 4πe-10·R_T | −GENCOOL | h·Σe_jR_j | M_2γ | R_3 | 2γ up, OLD | collisional ≤ | RE at SE, OLD [raw/l1] | RE at SE, exact [raw/l1] |
|---|---|---|---|---|---|---|---|---|---|
| 8.87 | 5.35e-9 | −1.4e-10 | 4.07e-9 | 1.42e-9 | 252 | 136 | 8.3 | +1.9e-3 | −2.0e-4 |
| 10.39 | 5.96e-9 | −1.17e-9 | 5.72e-9 | 1.41e-9 | 357 | 187 | 72 | +3.5e-4 | −1.7e-3 |
| 12.32 | 7.91e-9 | −5.82e-9 | 1.29e-8 | 8.7e-10 | 774 | 339 | 356 | −7.2e-3 | −8.5e-3 |

The "RE at SE, OLD" column sits on the same curve as the frozen-J thermal curve in `Q/Reports/evidence/report23/thermal_curve_d35_it25.txt` (+1.47e-3 at 9.47 kK, −1.21e-3 at 10.89 kK). The thermal curve's root near 10.2 kK is therefore where real Lyα collisional cooling balances the spurious M_2γ.

## 5. Best-supported explanation

1. **Established.** No sign or unit error; RE and GENCOOL are different balances linked by the identity in section 1.
2. **Established.** n=2 is frozen and out of equilibrium at the front; R_T is its energy image, so the "close to statistical equilibrium" premise is false for the level that matters.
3. **Established.** 64–76% of R_T at depths 30–41 is the `OLD_DEFAULT` two-photon absorption, forced by `COMP_F=T` in single-iteration passes. Its exp(−hν_B/kT_e) factor explains the growth with T from 7.9 to 12.3 kK and the absence of a coupled root.
4. **Inference.** With a consistent two-photon rate and n=2 actually equilibrated, R_T = −GENCOOL net, i.e. −1.2e-3 to −1.7e-3 of l1 at depths 33–37. The front would then cool, as GENCOOL says, and the 10 kK front would not be a solution. This assumes GENCOOL's net cooling is insensitive to n2; bound-free heating from n=2 is currently 1e-15 to 7e-13, so n2 would have to rise by more than 10³ to matter.
5. **Not measured directly.** The remaining 24–32% of R_T (about 98 of 357 cm⁻³ s⁻¹ at depth 35) is net H-line absorption, attributed to Lyα because n2/n1 has a minimum at depths 33–35. I inferred it from the budget and did not compute it from the line radiation field.

## 6. Single most decisive cheap test

Run one frozen-J replay pass at `Q/results/nk167_pin/theta1.00/it05_trial_l1` with the single extra setting `--vadat TWO_METH NOSTIM`. Use the same mechanism as the column passes (`Q/scripts/cool_newton_krylov.py:657`, `replay=base/'QSTARS_RJ.bin'`). Cost is about 2 CPU-minutes on one thread. `NOSTIM` removes two-photon absorption and keeps spontaneous emission, and it does not trigger the `comp_opac.f:123` fallback. Confirm that OUTGEN prints "Using NOSTIM method" before reading the result.

Predictions at unchanged populations, depths 33–37:

| d | R_T raw now | R_T raw predicted | R_T/l1 predicted | R_3 now | R_3 predicted |
|---|---|---|---|---|---|
| 33 | 1.994 | 0.48 | 1.7e-3 | 147.7 | 66 |
| 34 | 3.059 | 0.75 | 1.9e-3 | 228.3 | 106 |
| 35 | 4.742 | 1.19 | 2.2e-3 | 356.9 | 170 |
| 36 | 7.181 | 1.87 | 1.9e-3 | 539.6 | 266 |
| 37 | 10.93 | 2.94 | 1.6e-3 | 806.8 | 405 |

An unchanged row 134 would falsify explanation 3. A drop as predicted confirms it. The follow-up, which is not cheap, is then to converge n=2 at depths 25–45 with full weight; the prediction is that R_T turns negative there.

Scratch scripts and outputs are in `S` (`ana1.py`–`ana6.py`; `pin.txt`, `hot28.txt`, `hist_pin.txt`, `ident_pin.txt`, `twophot_pin.txt`, `twophot_recalc_pin.txt`, `budget.txt`, `budget_T.txt`).