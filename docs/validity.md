# Validity domain and release limitations

**This is not validated general-purpose QStarS v1.** Current time-dependent
results concern smooth, Newtonian, fixed-excision envelopes with fixed
composition and explicitly chosen outer boundary conditions. Some hot
pressure-confined models have quantitative convergence evidence. They are
not complete quasi-star evolutionary tracks.

## Active checks

| Check | Threshold and purpose |
|---|---|
| Mass and energy ledger | Default residual limit 1e-7 of initial mass/internal energy; validation additionally reports binding-energy and throughput normalizations |
| Radial motion | Configured maximum Mach 0.2; the run stops within 0.01% of that boundary to avoid endless bisection masquerading as a solver failure |
| Newtonian gravity | Maximum compactness 2GM/(rc²) = 0.001 by default |
| Ideal EOS | T=1e5..1e8 K, H and He ionization screening >=99%, electron occupation <=0.1 |
| MESA interface | H1/He4 only; reject neutral ideal backstop, nonfinite values, unsupported query coordinates and first-law defect >1% |
| MESA dilute gas pressure | Reject >1% discrepancy from rho kT/(mu m_u) when T>=40000 K, eta<-10 and conservative Coulomb Gamma<0.01; radiation must not hide gas-pressure error |
| Omitted conduction | Reject conductive contribution >1% of diffusive heat transport |
| Radiative diffusion | Radiation-gradient Knudsen number <=0.1 and flux below free streaming; outer flux is checked separately |
| Efficient convection | Pe>=100 on faces carrying >1% of the transported flux in convection; no supersonic blobs |
| Geometry | Positive ordered faces and volumes; independent geometry/volume agreement within 1e-8 |
| Boundary replacement | <=10% of the affected cell per substep; retries preserve the accepted ledger |
| Nonlinear solve | Recomputed momentum, energy and transport residuals; no energy residual reset |
| Restart | Versioned, finite physical arrays, complete finite ledger, mesh/volume consistency and stored solver state |

Kn=0.1 and Pe=100 are conservative asymptotic screening choices: diffusion
and lossless convection need scale separation. They do not establish LTE
thermalization, a realistic photosphere or one-percent evolutionary accuracy.
Cell optical depth is saved but no longer used as a rejection rule: refining
a smooth diffusive field must not make its physical validity disappear.

The default conservation limit is an execution safeguard. Passing it is not
a sufficient research accuracy criterion, especially when transported energy
is tiny relative to stored radiation energy. Reports retain all three energy
normalizations and require convergence of interior profiles.

## Limits that invalidate a general v1 claim

- A moving cavity is absent. Mesh redistribution conserves resolved-domain
  mass and energy but does not create the missing physical cavity reservoir.
- A validated coupled production atmosphere and full radiatively lossy convection are absent.
  Effective temperature from luminosity and outer radius is a diagnostic;
  the envelope solver has no validated photospheric radius.
  Native full MLT and explicit-pressure atmosphere references now execute,
  but the 6000 K candidate requires Mach 2.34 convection and the 8000 K
  candidate is geometrically extended. Numerical convergence does not make
  those assumptions physical; [Report 8](../Reports/Report_8.md) records the
  specific release blocker and its sensitivity tests.
- The wind ODE component is not coupled to the envelope. A prescribed wind or
  launch formula does not prove an escaping transonic solution exists.
- MESA EOS blending has retained first-law counterexamples. Rejecting them
  is necessary, but it leaves holes in the usable thermodynamic domain.
- The nuclear/composition evolution, stellar backend and several full
  published quasi-star evolutionary tracks remain unvalidated.
- No cross-paper evolutionary predictive test without retuning has passed.

These requirements are scientifically meaningful in 1D and cannot be waived
merely because implementation is difficult. In contrast, rotation, discs,
magnetic geometry and simultaneous angular inflow/outflow are unresolved by
this spherical model. Its effective streams must never be represented as a
resolved multidimensional solution. Full GR/ADM-mass accounting is also
outside the documented Newtonian convention.

`model_domain_exit`, `numerical_failure`, `control_stop` and
`completed_requested_interval` have distinct meanings. None establishes
physical destruction, instability or a quasi-star lifetime. QSRD remains
blocked by the independent framework-validation prerequisite.
