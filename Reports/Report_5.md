# Report 5 — Convergence repair and direct MESA microphysics

This is an ongoing v1 closure campaign. QStarS is not yet v1.

The previous statement that MESA was absent was incorrect. A compiled MESA
24.08.1 and ARM SDK 24.7.1 already existed under the earlier QSRD runtime.
Read-only inspection found local build changes (disabled graphics and optional
pulsation tools, version-number bookkeeping); QStarS records those changes and
reuses the libraries without changing the other project's runtime or results.

The fixed-excision solver's concentrated-mesh failure had multiple causes:

- subtracting nearby radial faces to reconstruct tiny cell volumes introduced
  density noise; stored cell volumes now advance by exact secant face motion;
- finite-difference nonlinear Jacobians and forward evaluation of almost
  adiabatic MLT amplified that noise; a sparse complex-step Newton solve now
  treats face luminosities as unknowns and inverts the same transport law;
- a depleted first mass cell was never replenished by mesh redistribution;
  conservative overlap remapping now restores the reference mass fractions
  when boundary-cell masses change by the configured fraction;
- the old invalid-geometry penalty had a spurious zero far outside the
  admissible domain. Invalid trial geometry now rejects the trial.

The first three fixes retain the documented equations and physical
coefficients. Stored volumes and solved luminosities are checkpointed; old
checkpoint schemas remain readable. The old nonlinear choices remain
selectable for reproduction. Failed forward-flux solves, wall-budget stops,
and the invalid-geometry exception remain in distinct result directories.

A separate continuum n=3 boundary-value solution provides an independent
hydrostatic spatial reference. Radially designed Lagrangian meshes address
poor resolution near fixed mass-grid boundaries. The new campaign measures
density, temperature, pressure, specific energy, radius, and entropy profiles,
including whole-mass overlap norms and explicitly reported interpolation
coverage. Small global energy residuals do not certify these profile errors.

The MESA interface links its public EOS and opacity APIs through an explicit
versioned C ABI. It currently supports fixed H1/He4 composition only. It does
not supply a production atmosphere or MESA stellar evolution. A 128-zone hot
two-stream run completed using direct MESA calls with independently accounted
energy residual about 2.4e-14 of initial U, and 1.0e-8 of energy throughput.

A wider 153-point microphysics audit found large first-law defects in some
cool, dense EOS blend states (up to about 0.31 with the stated normalization).
This is an observed limitation to investigate and restrict, not a passing
production-microphysics gate. Low-density ionization examples have much
smaller defects. MESA documentation explicitly discusses inconsistencies
between blended sources; its neutral ideal backstop must not be silently
accepted as ionization physics.

Sources: [MESA EOS overview](https://docs.mesastar.org/en/24.08.1/eos/overview.html),
[public MESA source](https://github.com/MESAHub/mesa/tree/24.08.1), the existing
QStarS equation derivations and source registry. Subsequent evidence reports
will give final campaign thresholds and results. No QSRD experiment is run.
