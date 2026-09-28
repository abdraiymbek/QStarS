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
| CB24 printed cutoff Ki (x_tilde=1) | PASS | 0.46554529552993895 | 0.46554 | 5.295529938942334e-06 |
| CB24 printed cutoff BH/inner mass | PASS | 0.6320874468538957 | 0.6321 | 1.2553146104266055e-05 |
| CB24 composite maximum at printed cutoff | PASS | 0.6178279952051948 | 0.618 | 0.0001720047948051695 |
| Two-stream evolution completes | PASS | completed_requested_interval | — | — |
| mass_residual_relative | PASS | 1.1756803922309254e-16 | 0.0 | 1.1756803922309254e-16 |
| energy_residual_relative | PASS | 1.2029733266271373e-12 | 0.0 | 1.2029733266271373e-12 |

## Convergence measurements

| Axis | Zones | dt [yr] | Solver tolerance | Status | Radius difference | Luminosity difference | Energy residual / U₀ |
|---|---:|---:|---:|---|---:|---:|---:|

The finest run is a comparison reference, not an exact solution. Small conservation residuals do not establish physical convergence.

## Open gates

- Full Ball 2011 Table 1 evolutionary replication
- MESA-QUEST thermal tracks
- Coupling stationary Fiacconi & Rossi winds to a conserved evolving envelope
- Cantiello nonadiabatic stability
- Production atmosphere and stellar evolution with audited MESA microphysics
- Moving cavity/Bondi/saturated inner evolution
- Long-duration quasi-static sequences and QSRD demonstration
