# Ball 2011 Q2 source audit (V1 execution plan P5, bullet 1)

Date: 2026-11-06. Read-only audit: no files changed under `benchmarks/`, `src/`, `scripts/`, or in any existing
results or report file. The only files written are in this directory:

- `quantify_row0_inversion.py`: the analysis. It reads the preserved profiles and writes `row0_inversion.json`.
  Run with `.venv/bin/python -I Reports/evidence/v1_exec_20261106/ball_q2/quantify_row0_inversion.py`; it takes about 2 s.
- `row0_inversion.json`: all numbers quoted below.

**The registered Ball campaign verdict stays FAIL.** Nothing here changes it. Section 5 is a draft for the user to
review. It is not registered and has not been approved.

## 1. What the sources say

Ball, King & Begelman 2011, MNRAS 414, 2751 (`sources/cache/ball2011.txt`). The PDF could not be rendered here
because poppler is not installed. Text and figure-axis labels come from the PDF text layer.

- Sec. 3.1 "Structure", PDF p. 5 (MNRAS p. 2755). This is the only statement about the inversion's onset:
  > "The second feature, apparent in all but the first density profile in Fig. 1, is the density inversion in the
  > outer layers. It appears once the photospheric temperature T_surf drops below about 8000 K. From then, the
  > surface opacity increases owing to hydrogen recombination. The Eddington luminosity falls and the star's
  > luminosity apparently exceeds the Eddington limit. It is well known that hydrostatic models can sustain this
  > super-Eddington luminosity through a density inversion (Langer 1997)."
- Fig. 1 caption, PDF p. 4 (MNRAS p. 2754):
  > "Plot of density against radius for models in the fiducial run with M_BH/M_sun = 5, 100, 200, 500, 1000 and 1194
  > (see Table 1). ... In the outer layers the density is inverted, as discussed by BRA08. Note that the initial model
  > has inner radius of 1.66 R_sun which is too small to be seen."
  
  The axis labels in the text layer are: rho/g cm^-3 on a log axis with ticks 10^-12, 10^-9, 10^-6, 10^-3 (nine
  decades), and r/10^4 R_sun on a linear axis from 0 to 4.
- Related statement, PDF p. 9 (MNRAS p. 2759): "Once the surface temperature of the quasi-star decreases below
  about 8000 K, the Eddington limit becomes much larger owing to the decreasing opacity." This is about the
  Eddington limit, not the inversion.
- The paper gives no amplitude, no threshold, no location (beyond "in the outer layers"), and no resolution study
  for the inversion. It says nothing about an inversion anywhere other than the outer layers at row 0.

Ball PhD thesis (`sources/cache/ball_thesis.txt`, secondary source), Sec. 3.2.1 "Structure", printed pp. 33-34
(PDF pp. 51-52):

> "Secondly, the density in the outer layers is **usually** inverted as in all but the first density profile in
> Fig. 3.2. The density inversion appears once the photospheric temperature Teff drops below about 8000 K.
> Thereafter, the surface opacity increases owing to hydrogen recombination ... Density inversions are possible
> when the temperature gradient becomes strongly superadiabatic (Langer 1997), as is the case in regions where
> convection occurs but is inefficient."

The thesis Fig. 3.2 log-rho axis runs from 10^-12 to 10^-4. Its Table 3.1 gives Teff = 14.1, 5.22, 4.77, 4.55,
4.49, 4.51 kK. The thesis adds "usually" and the superadiabatic mechanism. Otherwise it repeats the paper's
claim, so it is not independent evidence.

## 2. How the registration maps the claim

The registered check (`benchmarks/ball2011_track_plan.json`, Q2):

- claim: "a density inversion appears once Teff < 8000 K (paper 3.1)"
- test: "outward density increase > 1e-3 relative present at rows 1-4 and absent at row 0"

The evaluator is `scripts/validate_ball_track.py:714`. It reads `max_density_inversion` from `src/qstars/quasistatic.py:1193`:

```python
max(np.diff(state.rho)/state.rho[:-1])
```

This is the largest relative increase between adjacent cells, taken anywhere in the envelope. The evaluator
requires it to be <= 1e-3 at M_BH = 0 (taken as 5) and > 1e-3 at M_BH = 100, 200, 500 and 1000.

| Q2 row | Table 1 M_BH (M_sun) | Table 1 Teff (kK) | Paper Fig. 1 | Q2 requires |
|---|---|---|---|---|
| 0 | 5 | 14.3 | first profile, no inversion "apparent" | adjacent-cell max <= 1e-3 |
| 1 | 100 | 5.22 | inverted | > 1e-3 |
| 2 | 200 | 4.77 | inverted | > 1e-3 |
| 3 | 500 | 4.55 | inverted | > 1e-3 |
| 4 | 1000 | 4.49 | inverted | > 1e-3 |
| 5 | 1194 | 4.51 | inverted | reported only |

The row assignment matches the source: row 0 is the only row above 8000 K, and it is the "first density profile".
The absence requirement at row 0 is grounded in the source sentence "all but the first", although the registered
claim field quotes only the onset sentence.

## 3. Row-0 inversion in the existing runs

Inputs are the relaxed row-0 profiles `profile_row0_{6000,3464}rsun.npz` of the last relax campaign, step 5 of
Report 23 sec. 6, in `benchmarks/results/2026-10-14-ball-track/`. The runs are 045753 C1 N=200, 045958 C0 N=400,
050203 C2 N=800, and the N=400 variants 050616 C5, 050739 C6 and 050913 C7. Both seeds give the same numbers to
about 1e-9 (relax status `converged_unique`).

Columns in the table below:

- **Adjacent max**: the registered Q2 metric.
- **Resolved rise A**: max over i<j of rho_j/rho_i - 1. This is the total outward rise across the inversion,
  measured from the density trough to the peak.

  The same value comes out when the search is limited to a fixed physical window, either T in [20, 40] kK or
  exterior mass M*-m in [1e-3, 1e-1] M_sun. The rise lies entirely inside both windows and is the only inversion
  in the model.
- **tau**: approximate. It is integrated inward from the outer face, plus 2/3.
- **Gamma_rad**: kappa L_r/(4 pi c G m) at the trough.

| run | N | adjacent max (Q2) | cells with a rise | resolved rise A | trough to peak r/R* | exterior mass (M_sun) | T (kK) | tau | Gamma_rad | R* (R_sun) |
|---|---|---|---|---|---|---|---|---|---|---|
| C1 | 200 | 9.07e-3 | 2 | 0.01553 | 0.9656 to 0.9731 | 0.0175 to 0.0129 | 28.6 to 26.5 | 22 to 16 | 0.995 | 3197 |
| C0 | 400 | 5.36e-3 | 5 | 0.01885 | 0.9643 to 0.9738 | 0.0169 to 0.0116 | 29.1 to 26.5 | 22 to 15 | 0.996 | 3114 |
| C2 | 800 | 2.75e-3 | 11 | 0.01944 | 0.9649 to 0.9750 | 0.0162 to 0.0107 | 29.1 to 26.2 | 22 to 14 | 0.996 | 3092 |
| C5 (tol 1e-10) | 400 | 5.36e-3 | 5 | 0.01885 | same as C0 | same as C0 | same as C0 | same as C0 | 0.996 | 3114 |
| C6 (cell boundary state) | 400 | 5.42e-3 | 6 | 0.01974 | 0.9640 to 0.9752 | 0.0165 to 0.0105 | 29.3 to 26.2 | 22 to 14 | 0.996 | 3078 |
| C7 (split/merge 1.5/0.67) | 400 | 5.36e-3 | 5 | 0.01885 | same as C0 | same as C0 | same as C0 | same as C0 | 0.996 | 3114 |

`ln rho` on a common temperature grid from 26 to 30 kK converges in the same way. The N=400 and N=800 profiles
differ by 2.3e-3 to 3.1e-3 in `ln rho`. That difference is mostly a uniform offset that follows the 0.7 %
change in R*, and the shape of the bump is the same (`row0_inversion.json`, `ln_rho_on_T_grid`).

**Convergence.** On common physical coordinates, the row-0 inversion is resolved and converging. It is not a
grid artifact:

- **Amplitude**: A = 0.01553, 0.01885 and 0.01944 at N = 200, 400 and 800. The successive differences are
  3.3e-3 and 5.9e-4, a ratio of 5.6 (observed order about 2.5). The Richardson estimate is A of about 0.0196,
  that is a 2.0 % rise, or 0.0085 dex.
- **Location**: the inversion stays fixed at T of about 29 to 26 kK, 0.016 to 0.011 M_sun below the surface,
  tau of about 22 to 14, and r/R* of about 0.965 to 0.975.
- **Q2 metric**: the adjacent-cell value roughly halves with each doubling of N (9.07e-3, 5.36e-3, 2.75e-3) only
  because the same rise is spread over 2, 5 and then 11 cells. Extrapolating that 1/N trend, it would drop below
  1e-3 at roughly N of 2300 while the physical inversion stays at 2 %. In that case the row-0 outcome of Q2 would
  change with mesh resolution, not with the physics.

Physically, the rise sits on a small opacity bump: kappa goes from 0.362 to a peak of 0.371 cm^2/g near 29 kK
(N=800), against an electron-scattering floor of about 0.34. Gamma_rad is about 0.996 there, so the radiative
layer is at the Eddington limit, and a 2.7 % opacity excess produces a 2 % density rise. The row-0 Teff is
13.96 kK at N=200 and 14.20 kK at N=800, against 14.3 kK in Table 1. This is not the hydrogen-recombination
inversion the paper describes. It is a different, hot, near-Eddington feature of about 0.01 dex, sensitive to
small opacity differences such as MESA tables versus Ball's STARS tables.

**Along the track.** The smoke run (C1, N=200, `20261014T051036.115106Z_smoke_C1_n200/history.csv`, 1255 steps
to M_BH = 46) shows the following:

- The hot inversion fades as Teff falls. The adjacent metric drops below 1e-3 at M_BH = 6.59 (Teff = 12.34 kK).
  From M_BH = 6.85 (Teff = 12.12 kK) the model has no outward density increase at all.
- An outward increase reappears at M_BH = 23.47 (Teff = 7.60 kK). The adjacent metric exceeds 1e-3 at
  Teff = 7.60 kK, 1e-2 at Teff = 7.57 kK and 0.1 at Teff = 7.38 kK.
- At the last accepted model (M_BH = 46.0, Teff = 6.46 kK), the adjacent metric is 0.60 and the resolved rise is
  A = 8.27, a 9.3-fold increase. It runs from T = 13.1 to 7.0 kK, tau of about 22 to 1.2, and exterior mass 0.25
  to 0.04 M_sun. This is the hydrogen-recombination, super-Eddington outer-layer inversion. Its onset at about
  7.6 kK agrees with the paper's "below about 8000 K".

This is one resolution only, and the run did not reach row 1 (M_BH = 100).

## 4. Verdict

**Mismatch: the registered test makes a stronger claim than the source (category b), around a core the source
does support.**

- **The source does support absence at row 0.** The paper says the outer-layer inversion is "apparent in all but
  the first density profile" and appears once T_surf is below about 8000 K. Row 0 (Teff 14.3 kK) is that first
  profile.
- **The source does not support the registered threshold.** Its evidence for absence is visual: Fig. 1 plots
  log rho over nine decades. The row-0 star occupies r of 0 to 0.3 on an axis that runs from 0 to 4 (×10^4 R_sun).
  The resolved row-0 rise in our model is 0.0085 dex over about 1 % of R* (about 30 R_sun). That is about 0.1 % of
  the figure's height over less than 0.1 % of its width, far below what the figure can show.
- **The paper is silent on amplitude.** It gives no amplitude or threshold. It ties "the" inversion to hydrogen
  recombination and super-Eddington surface layers at T_surf below about 8000 K. It neither says nor implies that
  a 14.3 kK model has no outward density increase above 1e-3 anywhere.
- **The registered metric is not a physical amplitude.** "Any adjacent-cell increase > 1e-3 anywhere" scales
  like 1/N for a resolved feature. Section 3 shows the Q2 row-0 outcome would flip at about N = 2300 with no
  change in the converged 2 % inversion. A test whose verdict depends on N while the physics is converged is not
  a faithful encoding of a qualitative claim.
- **The thesis weakens the claim further.** It says the outer layers are "usually inverted", which softens the
  statement rather than sharpening it.

**Why not (c) ambiguous.** The source is clear about what it observes, a visibly inverted outer layer below about
8000 K that is absent in the first profile. The problem is the registration's added quantitative strength, not
an unclear sentence.

**Why not (a) faithful.** The row-0 failure comes from a converged 2 % feature that the source could not have
shown. It is not a contradiction of anything the source asserts. Whether Ball's own row-0 model contains a
similar sub-visual bump cannot be decided from the paper.

**The FAIL stands.**

- The original Q2 fails at row 0 at every resolution (9.07e-3, 5.36e-3 and 2.75e-3 against 1e-3), and that
  registered result is preserved.
- Gate 6 also fails independently: the smoke numerical blocker at 46 M_sun, and the campaign was never run.
- The rows 1-4 part of Q2 has not been evaluated, because no run reached M_BH = 100.
- Under the plan's decision table, this audit only establishes the precondition for asking the user to authorize
  a successor registration. It authorizes nothing by itself.

## 5. DRAFT successor check, Q2-S (not registered, needs explicit user approval)

> Status: draft text for review. It is not a file under `benchmarks/` and is not registered. Adopting it requires
> the user's explicit approval and a new registration (sha256 plus UTC time) made before any run that it scores.
> The original Q2 and its FAIL are kept and reported alongside it. No existing target, tolerance, physics or
> microphysics changes.

```json
{
  "id": "Q2-S",
  "supersedes": null,
  "reports_alongside": "Q2 (original, unchanged; its value and verdict are always printed)",
  "claim": "The outer-layer density inversion is apparent in all but the first density profile of Fig. 1 and appears once T_surf drops below about 8000 K (paper 3.1, p. 2755; thesis 3.2.1 'usually inverted')",
  "metric": "resolved outward rise A = max_{i<j, r_j > r_i} rho_j/rho_i - 1 over the envelope, with trough/peak located in T, exterior mass, tau and r/R*; evaluated at N = 200, 400, 800 and Richardson-extrapolated; the adjacent-cell metric of Q2 is also printed",
  "visibility_threshold": "A_vis = 10^0.1 - 1 = 0.259 (0.1 dex): proposed as the smallest rise legible on Fig. 1's nine-decade log-density axis. USER DECISION REQUIRED (see hazard note)",
  "test": [
    "rows 1-4: converged A >= A_vis with the peak at tau_from_surface < 100 (outer layers)",
    "row 0: converged A < A_vis",
    "onset (reported, required only if the track reaches row 1): first accepted model with A >= A_vis has Teff within [7000, 9000] K. USER DECISION REQUIRED on the bracket"
  ],
  "convergence_requirement": "the Q2-S verdict at each row must be the same at N = 400 and 800 and for the Richardson value; otherwise INCONCLUSIVE, not PASS",
  "always_reported": "every outward density rise with A > 1e-4 at every compared row and at the last model: amplitude, location (T, exterior mass, tau, r/R*), local Gamma_rad and kappa, and its N-dependence; the original Q2 values and verdict"
}
```

**Hazard note for the decision.**

- The row-0 amplitude (A of about 0.020) and the N=200 smoke onset (Teff of about 7.6 kK) are already known, so
  any threshold or bracket chosen now is chosen after seeing the data.
- A_vis = 0.1 dex and the onset bracket of 7000 to 9000 K are proposed from the source, from Fig. 1's axis span
  and the paper's "about 8000 K". The user should judge whether those are genuinely source-derived.
- An alternative that does not depend on the figure's resolution is to require only the onset and rows 1-4
  presence, and to report the row-0 amplitude without a verdict. That still drops a source-supported claim (row-0
  absence), so it is the weaker option.
- If neither is acceptable, the plan says to keep Q2 as registered and investigate the model and microphysics
  differences instead. Here that means the 2.7 % opacity bump near 29 kK at Gamma_rad of about 0.996 in MESA
  versus the STARS opacities. Do not re-register merely to obtain a PASS.

## User decision (2026-11-07, after the three-member council review of Report 26)

- The registered Q2 result stays **FAIL**; Gate 6 stays FAIL (completion_gates.json: the composite active-gate entry was
  split by component, Gate 6 = FAIL).
- A source-faithful successor check (Q2-S) is **authorized in principle only**. No numeric threshold is approved: the
  draft's 0.1 dex amplitude and 7000–9000 K onset bracket were chosen after seeing the data. A threshold must first be
  derived from the paper's claim ("apparent in all but the first density profile … once T_surf drops below about
  8000 K") and written down, with its derivation, before any scored run; Q2-S is then registered as a new plan, never
  as an amendment of Q2.
