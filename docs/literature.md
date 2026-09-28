# Literature audit and interpretation corrections

Thirteen primary PDFs were downloaded and inspected in relevant sections.
Hashes and page counts are in `sources/acquisition.json`; the full texts are
cached locally and not redistributed. The source registry states which
prescriptions are implemented and the exact reference locations. The complete
145-page Ball thesis was not audited line by line.

## Accretion convention

BRA08 Eq 3 is an explicit luminosity law. Ball 2011 Eq 15 prints `4π`, while
substitution of its Bondi radius into its earlier convection limit suggests
`16π`. The later Santarelli paper uses `16π` and a gravitational-growth
efficiency correction. The early QUEST note prints a different convention.
QStarS exposes named laws instead of adjusting their coefficients to force
agreement. See [Ball](https://arxiv.org/abs/1102.5098),
[early QUEST](https://arxiv.org/abs/2507.12651), and
[Santarelli](https://arxiv.org/abs/2510.11772).

## Boundary-dependent mass limits

The BH-only Bondi radius plus cavity gives Ball's n=3 polytropic limit near
0.102. Including the cavity in the Bondi radius gives a limit near 0.017 for
N=1. These are distinct static models, not universal destruction thresholds.

Hassan's Appendix A cubic has a small positive root up to the quoted critical
mass. At the critical mass the inspected displayed Eq 6 returns twice the
positive double root when substituted into the cubic. QStarS solves the cubic
directly and checks it by substitution. This is an equation-level discrepancy,
not evidence that the authors' MESA implementation used a wrong root.
[Hassan et al.](https://arxiv.org/abs/2510.18301)

## Saturated-convection cutoff

Integrating the printed ODE gives an entropy minimum at
`x_tilde=1.01483543`, consistent with the stated `≈1.015`, with
`Ki=0.46835699` and `MBH/Mi=0.62510156`. The paper's quoted
`Ki≈0.46554`, `MBH/Mi≈0.6321` instead occur at **x_tilde=1**.
Matching the outer envelope at that printed cutoff gives
`MBH/Mstar=0.617828` and adiabatic mass fraction `0.02256`, reproducing
the quoted ~62% maximum and ~2% outer mass. Both interpretations are retained;
the entropy root is not shifted to make the scalars match. Integrations in raw
radius and log radius agree, and tighter tolerances preserve the discrepancy.
[Coughlin & Begelman](https://arxiv.org/abs/2405.00084)

## Winds and recent work

Dotan's launch/tiring prescription is an explicit optional boundary closure.
Fiacconi & Rossi's radiation-enthalpy advection cannot be reproduced by a
tiring cap alone. Their transonic solution remains an unimplemented benchmark.
Cantiello's modes require nonadiabatic GYRE analysis: QStarS domain exits are
not evidence of these instabilities. Roman-Garza uses cold net mass gain,
nuclear evolution, GR, and late analytic extrapolation; this does not validate
a thermodynamic two-stream boundary. Begelman–Dexter and the 2026 CMFGEN
atmosphere work further motivate separating interior and spectral predictions.
These works were audited; their full tracks have not been replicated.

## Backend decision

The conservative kernel is useful for independent bookkeeping and equation
checks, and restricted hot-envelope evolution. It does not replace validated
ionization EOS, opacity blends, or stellar atmospheres. MESA remains the
preferred route for complete stellar models; see [MESA](mesa.md).
