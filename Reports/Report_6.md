# Report 6 — Fine-mesh arithmetic, stiff time integration and wind audit

QStarS remains a development build. This report records further repairs and
counterexamples; Report 7 gives the final gate matrix for this campaign.

The 2048-zone investigation separated four effects. Dense hydrostatic
initialization was replaced, on large meshes, by a sparse solve of the same
discrete equations. A fixed 40-iteration Newton limit stopped a still-converging
convective front; extending the iteration allowance reached the original
residual tolerance. An attempted linear initial flux predictor failed four
regressions and was removed. A previous accepted flux can now be reused only
when the inverse closure admits it; otherwise the original predictor is used.

Repeated remaps exposed cancellation in huge cumulative CGS mass coordinates.
The overlap sweep now calculates endpoints with decimal arithmetic and applies
conservative cell-integral weights. Forty remaps of 2048 cells preserve volume,
energy and radial geometry in an automated regression. Compensated radial
updates are checkpointed. Failed pre-repair runs remain intact.

A surface-refined mesh avoids excessive refinement of the first accreting
cell. Midpoint integration still developed grid-scale oscillations: its
2048-zone trial exited the lossless-convection domain near 0.214 years.
That is not a physical destruction result. A selectable backward-Euler
mechanical method now deposits its analytically derived kinetic dissipation
as internal heat. It preserves the whole-energy identity and records that
conversion explicitly. No residual compensation, coefficient fit or tolerance
relaxation is used. Its ten-year 256/512/1024/2048-zone campaign completed.
The kinetic-heating diagnostic alone cannot prove profile accuracy, so mesh,
timestep, remesh and cross-method comparisons are retained separately.

The stationary Fiacconi–Rossi radiation-dominated wind equations are now an
independent executable reference family. Inward integration includes radiation
enthalpy advection, kinetic and gravitational energy, the published thin
outer boundary, a photospheric diagnostic and the explicitly imposed base Mach
number. This gas-pressure-free model has no gas sonic eigenvalue selection.
It is not yet coupled to the evolving envelope. Figure 4 was digitized from
its vector paths with saved axis calibration, source hash and line-width
uncertainty. Not all plotted quantities agree with the exact caption setup.
In particular, optically thin sound-speed diagnostics depend on whether one
continues p=U/3 or uses the paper's g(tau)U. Both diagnostics are saved. No
parameter was fitted to the figure. Percent-level luminosity differences are
retained rather than attributed conclusively to a paper error.

The source-provenance audit found a reproducibility hazard: a long-running
Python campaign can keep imported modules while later snapshots read edited
files. New runs freeze package bytes at import and verify that the package
stayed unchanged. A changed tree is explicitly unverified. Earlier exploratory
campaigns lack this guard; their initial source snapshot and chronology must
be retained, and fresh unchanged-source verification is required for release
claims. A fresh MESA ten-year run matches the earlier state within 4e-11 in
relative density. A separate inventory hashes 8,135 installed microphysics
files totaling 15.54 GB without editing them.

These repairs do not supply a moving physical cavity, a production atmosphere,
a coupled stellar wind, full published evolution replication or cross-paper
predictive validation. Those remain scientifically meaningful requirements in
1D and cannot be waived as implementation inconveniences. QSRD is deferred.
