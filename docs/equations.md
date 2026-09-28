# Equations, closures, and discrete conservation

Internal units are CGS. Composition is fixed. There is no nuclear evolution,
rotation, magnetic field, relativity, or resolved accretion disc. The fixed
excision does not represent an evolving cavity reservoir.

## Continuum and thermodynamics

For envelope mass coordinate m:

\[
\partial r/\partial m=(4\pi r^2\rho)^{-1},\quad
Dv/Dt=-4\pi r^2\partial P/\partial m-G(M_{BH}+m)/r^2,
\]
\[
Du/Dt+P D(1/\rho)/Dt=-\partial L/\partial m=T Ds/Dt.
\]

The last equality is for fixed composition. The initial model is hydrostatic;
subsequent evolution retains inertia. Engine luminosity is an inner boundary
flux. Mass transport is a separate conservative source-cell substep.

With `Rgas=k_B/(mu m_u)` and `1/mu=2X+3Y/4+Z/2`:

\[
P=R_{gas}\rho T+aT^4/3,\quad u=3R_{gas}T/2+aT^4/\rho,
\]
\[
s=R_{gas}[\tfrac32\ln(T/1K)-\ln(\rho/1\mathrm{g\,cm^{-3}})]
+4aT^3/(3\rho).
\]

The entropy zero is fixed consistently for all streams. `cv=3Rgas/2+4aT³/rho`.
With `beta=Pgas/P`, `chi_rho=beta`, `chi_T=4-3beta`, derive
`Gamma3-1=P chi_T/(rho T cv)`, `Gamma1=chi_rho+chi_T(Gamma3-1)`,
`nabla_ad=(Gamma3-1)/Gamma1`, and `cp=cv+P chi_T²/(rho T chi_rho)`.
Pressure, sound speed, energy, and entropy use this same EOS. Monotonic
log-Newton inversion determines T from density and pressure or internal energy.
Saha ionization screening and nondegeneracy checks supplement the temperature
window; they are validity tests, not an ionization EOS. Metals use a trace-mixture
mean-molecular-weight approximation.

## Radiation and convection

Default `kappa=0.2(1+X)` is electron scattering, not an opacity table spanning
recombination or a thermalization model. External log-rho/log-T opacity tables
require provenance and reject extrapolation. No imported table is yet validated.

\[
L_{rad}=-4\pi r^2 ac/(3\kappa\rho)\,\partial T^4/\partial r.
\]

One face flux is shared by adjacent zones. Optional lossless MLT uses
`l=alpha Hp`, `Hp=P/(rho g)`, `D=max(nabla-nabla_ad,0)`,
`vconv²=g delta l² D/(8Hp)` and `Fconv=rho cp T vconv lD/(2Hp)`.
The numerical excess is `max[-Hp(ΔlnT-nabla_ad ΔlnP)/Δr,0]`, preserving
the entropy contrast of the same thermodynamic cell states. Flux is capped at
`P cs`; the uncapped ratio is saved. Blob radiative losses, overshoot,
composition gradients, and time-dependent convection are absent. This is not
full MESA ML1, and requires the efficient-convection limit.

The grey boundary is
`Lout=4πR² ac(Tlast⁴-Tbath⁴)/[3(tau_half+tau_external+2/3)]`.
It extrapolates constant opacity and plane-parallel transfer. The last resolved
zone must satisfy the EOS domain; an inferred cool atmosphere is not validated.
No photospheric radius is solved.

## Discrete work identity

Cells have mass `dm_i` and volume `V_i=4π(r_{i+1}³-r_i³)/3`. Nodal weights
`w_j` are half the neighboring cell masses. Define

\[
\Omega=-\sum_j G w_j(M_{BH}+m_j)/r_j,\quad
K=\tfrac12\sum_jw_jv_j^2,\quad U=\sum_i dm_i u_i.
\]

This is trapezoidal mass quadrature for spherical binding energy. During the
fixed-mass substep use midpoint velocity/pressure, gravitational acceleration
`G(MBH+m)/(rnew rold)`, and the secant area

\[
\bar A_j=\frac{4\pi}{3}(r_{new}^2+r_{new}r_{old}+r_{old}^2).
\]

Then `Δ(4πr³/3)=Abar Δr`. The cell equation is
`dm Δu=-Pbar ΔV+dt(Lleft-Lright)`. Multiplying the nodal momentum equation
by midpoint displacement and summing cancels interior pressure work, giving

`Δ(U+K+Omega)=dt(Linner-Louter)-Pext ΔVouter`.

Heat flux uses backward Euler; mechanics uses implicit midpoint. A coupled
MINPACK hybrid solve advances T and velocity. Residuals are recomputed after
the solve; the optimizer success flag alone is insufficient. Accepted time
steps have no energy-residual compensation. Rejected candidates never mutate
the accepted state or ledger.

## Independent streams and engine

Delivered matter joins the outer source cell at its resolved velocity with
specified density and temperature, or explicitly selected surface thermodynamics.
It adds `dm_feed h_feed`, with `h=u+P/rho`. Wind removes
`dm_wind h_surface`; BH consumption removes `dm_rest h_inner` at the inner cell.
Each mass and entropy flux remains separate when feed and wind happen to balance.
This is comoving boundary mixing, not resolved multidimensional inflow/outflow.

The potential per unit transported cell mass is the derivative of the discrete
binding energy:

\[
\phi_i=-G\left[\tfrac12(M_i/r_i+M_{i+1}/r_{i+1})+
\sum_{j>i}w_j/r_j\right].
\]

Evaluate at midpoint masses. Binding energy is quadratic in mass, so this
derivative reproduces its exact change. Kinetic transport per mass is
`(v_i²+v_{i+1}²)/4`. These source-cell averages require spatial convergence.
Changing the central mass gives the separate term
`Wcentral=-G ΔMBH sum(wmid/r)`; it is not fitted to the energy residual.

`ΔMBH=(1-epsilon)dm_rest`; `Egenerated=epsilon dm_rest c²`. Swallowed
non-rest energy is recorded in an unresolved engine reservoir, not added to
`MBH` in this Newtonian convention. Central-potential work is likewise an
explicit exchange with that prescribed mass source. A relativistically
consistent BH+envelope ADM-mass budget is not implemented.

If launch enthalpy plus potential and kinetic energy is insufficient for the
specified terminal wind speed, the deficit is debited from outward radiation
and credited to escaping wind exactly once. Insufficient luminosity rejects
the trial. A prescribed wind rate does not prove that a transonic solution
exists. The Dotan option implements launch/tiring formulas without a resolved
porous atmosphere or radiation-advecting transonic wind.

## Ledger and remap

Mass residual:
`Menv+MBH+Mwind+epsilon Mrest-Msupplied-Minitial`.

Energy residual:
`E-Einitial+Eradiated+Ewind+Eadvected-Esupplied-Egenerated-Wboundary-Wcentral`.

All terms are signed and saved. Energy residuals are normalized by initial
internal energy, initial binding energy, and transported energy throughput.
Entropy inventory and supplied/removed/mixing/remap entropy are exposed;
entropy is not assumed conserved during mixing or heat transport.

Remap conserves mass-coordinate overlap integrals of volume and total energy.
New internal energy equals remapped total energy minus new gravity and kinetic
energy. Negative internal energy is rejected. Entropy and hydrostatic balance
can change and must be assessed. Refinement is optional, disabled by default,
and doubles zones on a cadence; it is not feature-sensitive research AMR.
