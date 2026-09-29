# Apple Silicon performance report

Report 13 adds native CMFGEN on the same Apple M5 / 16 GiB MacBook, using
GNU Fortran 13.3.0, Accelerate and one OpenMP thread per process. The full
metal/clumped reference takes **883.66 s** for two initialization iterations
and reaches **9.78 GB peak footprint**. Reduced H/He single-transfer cases at
84/167/333/665 radial points take about **26.3/53.6/148.9/727.9 s**;
the initial 25-iteration thermal run takes **919.57 s / 0.57 GB**. These are
measured component costs with concurrent work, not isolated scaling or the
cost of a physical wind. The full native-enabled regression suite has
**108 passes in 135.30 s**. Exact per-case wall, resident and peak-footprint
measurements are retained in
[the campaign summary](../Reports/evidence/backend_20260928/final_audit/campaign_summary.json)
and [Report 13](../Reports/Report_13.md).

The reduced component is useful locally. Repeated full-metal models and global
wind-response campaigns should be evaluated on CPU machines with at least
32 GiB, preferably 64 GiB RAM; independent cases can use separate workers.
This is a resource recommendation based on measured memory, not a demonstrated
cloud speedup. Linux numerical agreement and within-model OpenMP scaling remain
unmeasured. No coupled non-LTE quasi-star evolution cost is established.

Report 12's main component campaign took **105.08 s / 762.42 MiB** and its
refinement **156.50 s / 936.75 MiB**. The 1601-point two-band thermal solve took
**51.39 s**; distant native diffusion matching with rejected-trial recovery took
**17.65 s**. The complete native-enabled suite passed **97 tests in 50.52 s**.
These are stationary reference components on the Apple M5 / 16 GiB MacBook,
not physical quasi-star tracks. Ray storage scales approximately as bands times
radial-points squared; independent larger frequency/grid campaigns can use cloud
CPU/memory resources. MESA calls require process isolation, not concurrent
threads. No cloud scaling measurement or full non-LTE model runtime is claimed.
See [Report 12](../Reports/Report_12.md).


Report 11's sequential grey radiation/stationary-control campaign took
**55.42 s / 745.03 MiB**; its outer-extent follow-up took **68.63 s / 1213.56 MiB**.
Individual 401/801-node stationary controls took **3.95/8.13 s**; the extended
2001-node, 100,000R control took **25.95 s**. These are useful individual local
CPU workloads on the same Apple M5 / 16 GiB MacBook. The full native-enabled
85-test suite took **34.29 s**. No coupled stellar evolution or cloud scaling
was measured. Angular-ray storage grows roughly quadratically with radial
resolution; independent larger grids can benefit from cloud CPU/memory capacity.
See [Report 11](../Reports/Report_11.md) for the narrow physical domain and failures.


Report 10 measured **205.23 s / 676.58 MiB** for 42 stationary attribution
runs and **31.74 s / 555.00 MiB** for the moving-face, MLT-exit and coverage
audit on the Apple M5 MacBook. The full 76-test suite took 32.69 s. These
are component measurements; no coupled evolution was run because Stage A
failed. Independent cases can use separate processes locally or on cloud CPU;
MESA calls are not thread-safe. No parallel scaling measurement is claimed.
See [Report 10](../Reports/Report_10.md) and its frozen drivers.


Report 9 adds a spherical wind component campaign: **170.41 s**, **569.16 MiB**
peak RSS on the same Apple M5 / 16 GiB machine. The two baseline low-gravity
controls take 0.77 and 0.64 s including diagnostics; base root searches dominate
total cost. These timings do not imply a validated stellar atmosphere or an
escaping wind. Cloud CPUs can parallelize future independent benchmark and
coupled-refinement cases, but cannot resolve the physical closure uncertainty.
See [Report 9](../Reports/Report_9.md) for the scope and evidence.


The Report 8 sequential verification with the new dilute-gas validity screen
completes the same 128-zone MESA ten-year envelope in **30.830 s**, at
**1045.0 MiB** peak RSS, with 500 accepted steps and no retries. The 14-case
uncoupled atmosphere campaign takes **7.159 s**, with **527.4 MiB** peak RSS.
Those atmosphere cases fail physical validity and are not production-track
performance measurements. Full quasi-star evolutionary-track costs remain
unknown; extra cloud CPUs cannot repair the demonstrated physical closure
failures. See [Report 8](../Reports/Report_8.md).

Measured on **Apple M5, 16 GiB RAM, macOS 27 arm64**, Python 3.12.14. Final cases ran sequentially in fresh CPU processes after the other QStarS campaigns finished. No GPU, cloud resources or compiled MESA redistribution is required for the analytic-EOS solver. MESA is a separate optional installation.

Model runtime includes initialization, evolution and model-output I/O. Fresh-process wall time additionally includes interpreter startup and provenance setup. Plotting is excluded. Peak RSS is per fresh process.

| Case | Zones | Simulated years | Accepted / retries | Model s | Process wall s | Peak MiB | Energy residual / throughput |
|---|---:|---:|---:|---:|---:|---:|---:|
| closed_relaxation | 32 | 1 | 107 / 0 | 0.585 | 0.823 | 73.8 | 0 |
| hot_accreting | 32 | 10 | 214 / 0 | 1.642 | 1.879 | 73.8 | -8.74e-11 |
| two_stream_pulse | 32 | 10 | 219 / 0 | 1.709 | 1.949 | 74.2 | -1.11e-10 |
| bra_luminosity_transient | 32 | 0.1 | 52 / 0 | 0.643 | 0.900 | 74.0 | 6.74e-10 |
| dotan_launch_closure | 32 | 0.1 | 25 / 0 | 0.162 | 0.410 | 73.7 | 8.09e-08 |
| pressure_confined_hot_envelope | 64 | 0.01 | 100 / 0 | 0.786 | 1.031 | 74.2 | -2.26e-10 |
| hot_resolved_validation | 512 | 10 | 500 / 0 | 24.215 | 24.463 | 88.9 | 3.15e-13 |
| mesa_hot_ten_year | 128 | 10 | 500 / 0 | 31.293 | 31.556 | 1045.4 | -3.45e-11 |

All eight complete their requested interval and pass unchanged-source provenance verification. Coarse exploratory examples are not evidence of research spatial convergence. The 512-zone `hot_resolved_validation` row has the bounded convergence evidence described in the validation report. The MESA row uses 128 zones and direct EOS/opacity calls; it is not a complete MESA stellar track.

## High-resolution and native stellar workloads

The ten-year backward-Euler surface-refined campaign took 12.5, 22.3, 51.0 and 513.1 seconds at 256, 512, 1024 and 2048 zones. Those were development measurements with other campaigns active, so they are not isolated scaling benchmarks. The 2048-zone run needed 2,372 accepted steps because boundary replacement limited the timestep; it is not a fixed-step factor-of-two scaling experiment. All reached the common ten-year endpoint.

A subsequent sequential repeat, with no other QStarS campaigns active, took
12.627, 22.752, 52.916 and 408.516 seconds at those four resolutions. All states
and ledgers reproduce the earlier results bitwise, and every run passes the
unchanged-source guard. These cases share one Python process, so its peak RSS
is cumulative rather than an independent per-case memory measurement. The
adaptive-step difference still prevents interpreting the timings as pure mesh
scaling. The roughly 6.8-minute 2048-zone case is a useful local workload for
this controlled ten-year problem.

The symmetric cosine mesh overconcentrates both boundaries and can impose very small source-cell timesteps. The repaired 2048-zone version was deliberately stopped after verifying the earlier geometry failure was gone. The surface-refined alternative and conservative implicit mechanics are substantially more useful for the tested secular hot-envelope problem. Midpoint with a 0.005-year maximum step required about 121 seconds at 512 zones; its interior profiles agree within one percent with backward Euler at the matched timestep.

The ordinary native MESA smoke test used about 85 seconds and 1.86 GiB peak child RSS to reach model 100, including pre-main-sequence construction. QUEST attempts stopped during initialization/inner relaxation; their wall budgets do not establish the cost of a successful quasi-star track. The installed microphysics dataset inventory occupies 15.54 GB on disk; the measured direct-library envelope run peaks near 1.05 GiB RSS.

## Local and cloud use

Useful individual controlled 1D runs fit this MacBook. Prefer local exploratory models and individual resolved hot-envelope runs; inspect their validity exits and convergence evidence before interpreting the science. Independent mesh/timestep/physics grids, many pinned stellar families and large parameter ensembles are reasonable CPU-cloud workloads. Use isolated processes and allocate memory per worker; the MESA interface is not thread-safe across multiple installations. Measure aggregate throughput instead of assuming linear scaling.

The Python solver has no GPU dependency. The MESA shared-library build has only been validated on Apple Silicon; a Linux build must be verified before calling a cloud MESA workflow supported. Cloud CPUs cannot fix an invalid closure, missing cavity accounting, unconverged stellar initialization or unavailable published reference data. Million-year production-track performance remains unknown.

[Validation scope](validation.md) · [Machine-readable measurements](../benchmarks/results/2026-09-28-closure/20260928T031938.396825Z_performance_validation/performance.json) · [Historical baseline](../Reports/evidence/performance_20260927T201211.json)
