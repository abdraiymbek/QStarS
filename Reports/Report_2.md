# Report 2 — Conservative kernel and failed campaign

Date: 2026-09-28 UTC.

The first implementation included EOS, opacity, hydrostatic initialization,
coupled evolution, independent mass streams, engine, wind, ledger, remapping,
checkpointing, controls, CLI, and static reference solvers. Tests covered the
first law, mass-energy conversion, equal feed/wind, immutable rejected trials,
remapping, and exact restart continuation.

The first short run conserved energy to about 4e-15 of initial internal energy.
Longer runs reached the Mach guard; the first convergence campaign also suffered
nonlinear failures at tiny retry timesteps. It is preserved in
`Reports/evidence/first_campaign.json` and `.md`. Its models ended at unequal
times, so those comparison numbers are NOT valid convergence evidence.

Investigation found artificial convection: an independently differenced pressure
scale height and temperature gradient did not preserve the initial polytropic
entropy contrast. Computing `Delta ln T - nabla_ad Delta ln P` from the same
states fixed the cause. A regression now requires zero initial convective flux
in the n=3 polytrope. No physical coefficient was adjusted.

The corrected two-stream model completed 0.1 years without retries, then ten
years in about five seconds. Original failed run directories remain untouched.
Other repairs addressed overlap-remap endpoint roundoff, scalar opacity-table
shape, and the plotting scale setter. None changed physical prescriptions.
