# Report 3 — Verification, performance, and remaining scientific gates

Date: 2026-09-28 UTC.

The release has 29 passing automated tests and six completed example simulations.
CLI coverage includes run, restart, validation, saved-output plotting, independent
CPU sweeps, and configurable static published-reference families. Primary source
metadata is included in the installed package and every new simulation output.

The final validation campaign retains two source-cutoff discrepancies and one
128-zone concentrated-mesh nonlinear failure as failed checks. Default-mesh
time and solver-tolerance comparisons complete at a common final time; spatial
internal-energy errors remain too large for research certification. A trial
trust-region fallback was rejected after slow, retry-heavy progress; its full
checkpoint, source, event log, and actual interruption context are preserved.

Configuration testing caught YAML 1.1 scalar issues (`off` and positive scientific
exponents). The safe scientific loader now handles both predictably. Additional
checks reject unsupported incoming thermodynamic states and test exact pulse
mass delivery. An explicit wall-clock budget and interruption handler preserve
the accepted checkpoint. Saved output plots were rendered and visually inspected.

On the Apple M5/16 GiB machine, ten-year hot-envelope and pulsed two-stream
examples took approximately 5–6 seconds with about 74 MiB peak process memory.
The concentrated 64-zone case was slower and less accurate relative to energy
throughput. These results do not predict the cost of full stellar lifetimes or
a MESA backend.

**QStarS v1 remains incomplete.** The implementation is a useful, working
development framework, not a completed generally validated scientific simulator.
Full published evolutionary/wind replication, production microphysics, moving
inner boundaries/cavity accounting, and stronger convergence remain required.
The requested QSRD demonstration is deferred until those gates are satisfied.
No science campaign, remote deployment, Git push, or upstream communication was
performed. Evidence and original code are committed locally.
