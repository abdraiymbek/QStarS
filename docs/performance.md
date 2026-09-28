# Performance on the available MacBook

Measured on Apple M5, 16 GiB RAM, macOS 27 arm64, Python 3.12.14. Each example used a fresh CPU process. No GPU, remote compute, or MESA run was used. Solver timing includes initialization and output I/O but excludes the separate plotting process. Peak RSS is measured per fresh process.

| Case | Simulated yr | Zones | Accepted / retries | Runtime s | Peak MiB | Final energy residual / U0 |
|---|---:|---:|---:|---:|---:|---:|
| closed_relaxation | 1 | 32 | 107 / 0 | 0.696 | 73.7 | -4.41e-15 |
| hot_accreting | 10 | 32 | 214 / 0 | 5.064 | 73.4 | 2.86e-15 |
| two_stream_pulse | 10 | 32 | 219 / 0 | 5.374 | 73.4 | 4.3e-15 |
| bra_luminosity | 0.1 | 32 | 52 / 0 | 0.576 | 73.2 | -1.2e-13 |
| dotan_launch | 0.1 | 32 | 25 / 0 | 0.253 | 73.1 | -5.29e-17 |
| pressure_confined | 0.01 | 64 | 279 / 81 | 20.052 | 73.9 | 1.99e-11 |

All six completed the configured interval. The ten-year two-stream case transported 0.16 solar masses in, 0.08 out, and 0.002 into the BH engine. BH gravitational mass increased by 0.0018 solar masses. Its energy residual was 4.30e-15 of initial internal energy, 1.16e-13 of initial binding energy, and 1.93e-12 of boundary/engine energy throughput.

The 64-zone concentrated, pressure-confined case required 81 retries and had an energy residual about 8.5e-6 of the transported energy, despite only 2.0e-11 relative to the much larger initial internal energy. This is why both normalizations are reported.

## High-resolution exception

A 128-zone core-concentrated thermal case failed with the default nonlinear method. An experimental trust-region fallback was rejected after 275 seconds: it advanced only 2.68e-5 years, with 504 accepted steps and 171 retries. Its checkpoint, source snapshot, and interruption context are preserved. The fallback was not retained in the release. This is an unresolved conditioning/timestep problem, not a physical instability.

Useful restricted individual runs are practical locally. These measurements do not establish the cost or correctness of million-year stellar evolution or a MESA backend. Independent sweeps can use `--workers`; retain one process per modest CPU allocation, then measure throughput rather than assuming linear scaling. Cloud/HPC is optional for independent ensembles, not a cure for unconverged physics or solver failures.

[Raw per-process measurements](../Reports/evidence/performance_20260927T201211.json) · [Rejected fallback](../Reports/evidence/rejected_trust_region.json) · [Saved runs](../benchmarks/results/2026-09-28/)
