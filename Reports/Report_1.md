# Report 1 — Scope, source audit, and architecture

Date: 2026-09-28 UTC; task began September 27 in America/Chicago.

The repository started with a placeholder README and initial Git commit.
Work uses `codex/qstars-framework`. No previous private QSRD material was
supplied. Thirteen primary PDFs, MESA documentation, and pinned QUEST code
were inspected; source hashes, locations, and implementation status are recorded.
Only relevant portions of the 145-page Ball thesis were inspected.

The development kernel uses consistent thermodynamics and compatible mechanical
work to conserve internal+kinetic+gravitational energy. Independent boundary
streams transport mass, enthalpy, potential energy, and entropy. Radial inertia
is retained; excessive Mach number is a domain exit, not physical destruction.

MESA is preferred for full production microphysics and complete stellar tracks.
No MESA installation or full stellar replication is claimed. The dedicated
kernel and its static reference ODEs do not close the user's v1 objective.

The audit identified Ball/QUEST normalization differences, a saturated-convection
cutoff inconsistency, and a displayed Hassan radius differing from its defining
cubic. These are preserved explicitly rather than absorbed into tuned coefficients.
