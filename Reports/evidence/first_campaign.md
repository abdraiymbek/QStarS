# QStarS validation evidence

**QStarS v1 is not validated.** This campaign separates numerical identities from publication replication.

| Check | Result | Measured | Reference | Absolute error |
|---|---|---:|---:|---:|
| Lane-Emden n=1 radius | PASS | 3.1415926535873955 | 3.141592653589793 | 2.397637643980488e-12 |
| Ball 2011 n=2 maximum BH fraction | PASS | 0.12675259050436447 | 0.127 | 0.0002474094956355366 |
| Ball 2011 n=3 maximum BH fraction | PASS | 0.10151516249158155 | 0.102 | 0.0004848375084184481 |
| Ball 2011 n=4 maximum BH fraction | PASS | 0.08852070726685707 | 0.089 | 0.00047929273314292753 |
| CB24 corrected Bondi N=1 | PASS | 0.017040027221003805 | 0.017 | 4.002722100380357e-05 |
| CB24 corrected Bondi N=2 | PASS | 0.005755921466625292 | 0.0058 | 4.4078533374707594e-05 |
| CB24 corrected Bondi N=0.333333 | PASS | 0.1263147141037867 | 0.13 | 0.0036852858962133084 |
| CB24 saturated entropy minimum radius | PASS | 1.0148354305817882 | 1.015 | 0.0001645694182117463 |
| CB24 printed Ki_max | FAIL / investigate | 0.4683569905792735 | 0.46554 | 0.0028169905792734817 |
| CB24 printed BH/inner mass | FAIL / investigate | 0.6251015602357203 | 0.6321 | 0.006998439764279718 |
| Two-stream evolution completes | FAIL / investigate | model_domain_exit | — | — |
| mass_residual_relative | PASS | 1.1756803922309254e-16 | 0.0 | 1.1756803922309254e-16 |
| energy_residual_relative | PASS | 2.4005670206736495e-12 | 0.0 | 2.4005670206736495e-12 |

## Convergence measurements

| Axis | Zones | dt [yr] | Solver tolerance | Status | Radius difference | Luminosity difference | Energy residual / U₀ |
|---|---:|---:|---:|---|---:|---:|---:|
| mesh | 16 | 0.0001 | 1e-10 | numerical_failure | 0.00633 | 0.305 | -8.81e-14 |
| mesh | 32 | 0.0001 | 1e-10 | numerical_failure | 0.00456 | 0.199 | -5.96e-15 |
| mesh | 64 | 0.0001 | 1e-10 | model_domain_exit | 0.00237 | 0.0667 | -4.41e-13 |
| mesh | 128 | 0.0001 | 1e-10 | numerical_failure | 0 | 0 | 4.57e-14 |
| timestep | 32 | 0.0004 | 1e-10 | numerical_failure | 0.000118 | 0.0039 | -5.19e-15 |
| timestep | 32 | 0.0002 | 1e-10 | numerical_failure | 6.03e-05 | 0.00221 | -1.03e-14 |
| timestep | 32 | 0.0001 | 1e-10 | numerical_failure | 2.14e-05 | 0.000857 | -5.96e-15 |
| timestep | 32 | 5e-05 | 1e-10 | model_domain_exit | 0 | 0 | -5.91e-15 |
| solver | 32 | 0.0001 | 1e-07 | model_domain_exit | 1.45e-09 | 2.84e-09 | 5.35e-11 |
| solver | 32 | 0.0001 | 1e-09 | model_domain_exit | 1.45e-09 | 2.45e-09 | 2.4e-12 |
| solver | 32 | 0.0001 | 1e-11 | numerical_failure | 0 | 0 | -1.61e-14 |

The finest run is a comparison reference, not an exact solution. Small conservation residuals do not establish physical convergence.

## Open gates

- Full Ball 2011 Table 1 evolutionary replication
- MESA-QUEST thermal tracks
- Fiacconi & Rossi transonic advective wind solutions
- Cantiello nonadiabatic stability
- Validated ionization EOS and Rosseland opacity
- Moving cavity/Bondi/saturated inner evolution
- Long-duration quasi-static sequences and QSRD demonstration
