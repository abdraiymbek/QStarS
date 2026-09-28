# Spherical atmosphere and wind components

Status: research component; **not a production outer boundary**. Reports 9–10 and
their machine-readable evidence distinguish equation verification from physical
validation and published-model replication. The Report 8 plane-parallel MLT
candidate remains rejected. Neither solver failure nor a local wind branch
proves the absence of other stationary solutions, time-dependent outflows or
porous states.

## Method audit and choice

| Primary source / implementation | Relevant method | Decision for this campaign |
|---|---|---|
| [Nakauchi et al. 2017](https://arxiv.org/abs/1611.09601), equations 1–9, Figure 2 | Spherical LTE diffusion with gas pressure, ionization energy, radiation enthalpy and gas sonic regularity | Implement this as a separate optically thick component. Its explicit gas sonic condition is useful where gas pressure becomes significant during recombination. |
| [Fiacconi & Rossi 2016](https://arxiv.org/abs/1509.07511), sections 2–3 | Gas-pressure-free radiation wind; advected radiation energy; prescribed base Mach; porous intermediate atmosphere | Preserve the existing alternative. Section 3's temperature-dependent opacity and porosity do not supply the gas critical-point eigenvalue or a general non-LTE boundary. |
| [Dotan, Rossi & Shaviv 2011](https://arxiv.org/abs/1107.3562), sections 3.3–3.4 | Porous atmosphere, wind function and tiring prescription | Keep its coefficients explicit. The existing launch formula is not a solved porous atmosphere. Do not merge it with the smooth wind to manufacture escape. |
| [Owocki, Gayley & Shaviv 2004](https://arxiv.org/abs/astro-ph/0409573) | Effective opacity from porosity length / distribution | Scientifically motivated alternative with unresolved structure parameters. Those parameters require independent calibration; no opacity reduction is fitted here. |
| [Quataert et al. 2016](https://arxiv.org/abs/1509.06370) | Heating-driven, radiation-enthalpy wind with photon trapping | Supports testing outflow rather than forcing MLT. Its adiabatic limit must not be extended automatically to a recombining photosphere. |
| [Owocki, Townsend & Quataert 2017](https://arxiv.org/abs/1708.07790) | Unified diffusive and advected radiation energy | Use its distinction between comoving diffusion and observer-frame radiative energy transport. Emergent light alone is not the wind energy supply. |
| [Hassan et al. 2025](https://arxiv.org/abs/2510.18301), section 2.4 | Cheng et al. eruptive mass removal in inefficient convective layers | A selectable future eruptive alternative, not a stationary transonic solution. The paper explicitly allows a turbulent bound outcome instead of escape. |
| [MESA-QUEST](https://github.com/andysantarelli/MESA-QUEST), cached commit `671eb9b898b28c0d02a6a2455dacc0f5b269e280` | `qs_satcon_winds_ApJ2026/src/run_star_extras.f90:333` applies `1.4e-4 Mstar^0.96 MBH^0.17` Msun/year; the LRD inlist selects TDC | The wind routine removes stellar mass; it does not solve an outer atmosphere, a sonic eigenvalue or a wind energy ledger. Merely selecting TDC is not validation in this regime. The Bondi project is a separate source directory. |
| [Jiang et al. 2015](https://arxiv.org/abs/1509.05417) | 3D radiation transport and convection at an opacity peak | Demonstrates why radiative leakage and density/flux correlations matter. Its iron-opacity composition is different from the H/He controls here. No universal turbulence coefficient is imported. |
| [Martins, Schaerer & Marques-Chaves 2026](https://arxiv.org/abs/2605.12141) | Spherical, expanding CMFGEN non-LTE atmosphere | Suitable direction for an independent escaping-radiation calculation. Prescribed velocity/density structure is not a self-consistent mass-loss eigenvalue. |

The adopted component is a smooth, stationary **outflow**, with no MLT flux.
This removes supersonic convective blobs from its equations; it does **not**
establish that neglecting turbulent convective transport is valid at every
possible envelope matching point. A decreasing entropy in a subsonic part
must not be silently declared stable. Stationary existence and stability are
separate questions.

## Region boundaries and matching

1. **Deep evolving envelope.** Hydrostatic or quasi-hydrostatic structure;
   convection requires its existing subsonic and radiative-loss validity tests.
   This campaign does not replace its transport with an unvalidated cap.
2. **Subsonic spherical component.** Radial gravity, flux and geometry vary.
   At a fixed matching sphere, supply radius, mass, composition, density,
   temperature, and total outward energy flux in a consistent energy gauge.
   Pressure and specific entropy follow from the same EOS on both sides.
3. **Gas sonic point.** Solve regularity and choose a real accelerating root
   crossing from subsonic to supersonic. A radiation sound-speed crossing in
   the Fiacconi alternative is not this mathematical critical point.
4. **Optically thick supersonic component.** Continue continuity, momentum,
   total energy and diffusion with realistic recombination opacity. Stop
   before an unresolved second sonic crossing or invalid diffusion.
5. **Photosphere and escaping radiation.** A production implementation still
   needs transport beyond the diffusion domain with an independently checked
   gas/radiation energy exchange and flux-mean opacity. This region is not
   supplied by a Rosseland table or an imposed blackbody temperature alone.

The stationary matching helper solves for sonic radius, density and
temperature from the three base constraints; `mdot` is an output. Recovery
from slices of a known solution is a numerical matching test, not independent
stellar validation. The Report 8 bases are also tried, with failures retained.
Their plane-parallel profiles are failed candidates, so even a successful
match would not validate the complete stellar boundary.

`boundary_flux()` returns mass, momentum stress, separate radiative and
advected gas/radiation internal energy, pressure work, kinetic and potential
energy, and advected entropy. Since Report 10 it accepts `face_velocity=u`.
For area A, its relative mass flux is A*rho*(v-u). Relative enthalpy advection
is supplemented by `moving_pressure_work=A*(Pg+Pr)*u`. Thus a comoving face
retains diffusive radiation and pressure work. Equivalently, energy flux is
the fixed-sphere flux minus u times area times lab-frame conserved energy
density. Momentum likewise subtracts swept momentum; entropy uses relative
mass transport. These formulas assume LTE isotropic radiation and a static
central potential. Independent swept-volume and energy-gauge checks are in
[Report 10](../Reports/Report_10.md). A moving stellar reconstruction is not
implemented or claimed. Feeding remains an independent channel in the core.

## Equations and independent derivation

Use CGS, (P=P_g+P_r), (P_r=aT^4/3), (h=u_g+P_g/\rho+4P_r/\rho).
The model has a fixed central mass (M), no distributed heating and no
convective luminosity. Define

\[
\dot M=4\pi r^2\rho v,\qquad
\mathcal E=L_d+\dot M(v^2/2+h-GM/r),
\]

\[
v\frac{dv}{dr}=-\frac{1}{\rho}\frac{d(P_g+P_r)}{dr}-\frac{GM}{r^2},\qquad
\frac{dT}{dr}=-\frac{3\kappa_R\rho L_d}{16\pi ac r^2T^3}.
\]

Here (L_d) is the **comoving diffusive** luminosity. The radiation enthalpy
is separately advected; using it again as diffusion would double count it.
Gas ionization is retained in (u_g). Pressure work is counted once inside
enthalpy and displayed separately in the boundary ledger.

For (x=\ln(r/r_s)), define (g_T=d\ln T/dx),
(c_T^2=(\partial P_g/\partial\rho)_T),
(b=(\partial P_g/\partial\ln T)_\rho/\rho). Continuity gives
(d\ln\rho/dx=-2-d\ln v/dx), so

\[
\frac{d\ln v}{dx}=\frac{N}{D},\quad
N=2c_T^2-b g_T+\frac{GM}{r}(\Gamma_d-1),\quad
D=v^2-c_T^2,\quad \Gamma_d=\frac{\kappa_R L_d}{4\pi cGM}.
\]

At the critical point (D=N=0). Thus

\[
L_{d,s}=\frac{4\pi cGM}{\kappa_s}
\frac{1-2c_{T,s}^2r_s/(GM)}{1+\rho_s b_s/(4P_{r,s})}.
\]

Given (r_s,\rho_s,T_s), this fixes (v_s=c_{T,s}), mass flux and total
energy. For the sonic slope (q=d\ln v/dx), l'Hopital's rule gives

\[
D_zq^2+(D_x+D_y g_T-N_z)q-(N_x+N_y g_T)=0,
\]

where (z=\ln(v/v_s),y=\ln(T/T_s)) and the partial derivatives keep the
mass/energy constants fixed. Only a real positive root with positive (D')
is used. Ambiguous or complex roots are reported as inadmissible **local
critical states**; they are not evidence that the global problem has no root.
Centered derivative step and sonic start offset are refinement variables.

The algebraic energy invariant is supplemented by independent quadratures:

\[
\Delta(v^2/2)+\int\left(\rho^{-1}dP+GM r^{-2}dr\right)=0,
\]
\[
\Delta P_r+\int \frac{\kappa_R\rho L_d}{4\pi c r^2}dr=0,
\qquad \Delta L_d+\dot M\int T\,ds=0.
\]

The last integral uses integration by parts with EOS entropy. It can detect
an inconsistent EOS/energy partition that exact algebra alone would miss.
The integrator's adaptive points are retained and subdivided for these
quadratures; a uniformly sampled recombination layer gave falsely large
integration defects in the first retained attempt.

## Microphysics, geometry and validity

`FreeEOSOpacity` uses direct FreeEOS 3.0.0 EOS1 derivatives, gas internal
energy and ionization, with the **same** composition and electron state passed
to MESA 24.08.1 radiative opacity. No MESA neutral fallback EOS is used in this
component. Build manifests and library hashes are saved with each campaign.
The previous FreeEOS upstream byte-test differences still apply. The
rectangle (10^{-22}\le\rho\le10^{-4}), (3000\le T\le10^7) is a query
bound, not a physical LTE certificate.

FreeEOS uses the diatomic hydrogen ground state for its energy zero
(`free_eos_detailed.f90` near line 3488). A constant shift cancels within the
stationary ODE, but **does not** cancel from a supplied total base energy
flux. Matching and an eventual evolving ledger must use the same reference
on both sides. No energy offset is fitted to reproduce a paper.

Gravity is (GM/r^2); diffusion flux is (L_d/(4\pi r^2)), and (L_d)
itself varies as mechanical and enthalpy transport change. Integrated shell
mass is reported; the point-mass approximation is screened at 1% of (M).
The separate spherical hydrostatic reference can include enclosed mass
self-consistently and has an independent mass integral test. A steady wind
with non-negligible accumulated mass cannot be promoted by keeping (M)
constant. An infinitely old constant outflow would also have unbounded wind
mass; no such global self-gravitating stationary configuration is claimed.

The geometric diagnostic is (1-(r_{\rm inner}/r_{\rm outer})^2).
For a 1% gravity/area-flux goal, a thin layer needs approximately
(|\Delta r|/r<0.005); actual luminosity variation is additional.

The strict endpoint is the gradient Knudsen number

\[
R_E=\frac{|d(aT^4)/dr|}{\kappa_R\rho aT^4}
=\frac{3L_d}{4\pi r^2caT^4}\le0.1.
\]

Optical depth is integrated only **to the finite computed endpoint**. It is
not labeled optical depth to infinity. `paper_photosphere` extends the
reference to (T^4=L_d/(4\pi r^2\sigma)), where (R_E=3/4), for the
published comparison. It necessarily fails the stricter diffusion screen.
It supplies neither an independently calculated optical photosphere nor a
terminal velocity. No automatic continuation through that domain boundary
is made. Cool reference integration also stops at 3100 K before the EOS
query limit, with no imposed temperature floor.

A decreasing velocity approaching a second sonic point is a resolved local
branch property. The code stops at Mach-squared 1.01 on the downward crossing,
not at a numerically singular division. A global fallback history requires
time-dependent dynamics and remains uncomputed.

## Reproduce

Build the independently installed libraries using `docs/atmosphere.md`.
Then set `QSTARS_FREEEOS_LIBRARY`, `QSTARS_MESA_LIBRARY` and `MESA_DIR` to
those recorded installations and run:

```sh
.venv/bin/qstars reference benchmarks/spherical_recombining_wind.yaml
.venv/bin/python scripts/digitize_nakauchi_reference.py
.venv/bin/python scripts/validate_spherical_wind.py
.venv/bin/python scripts/plot_spherical_validation.py results/CAMPAIGN
.venv/bin/pytest -q
```

The campaign includes four unaltered published sonic states, two low-gravity
component controls, tolerances, radial step/output refinement, sonic offset,
critical derivative step, independent integrator, outer placement, and
matching-radius recovery tests. Full stellar evolution is deliberately gated
on independent outer-boundary validation, as requested; no coupled parameter
sweep is performed.

## Report 10 attribution and validity enforcement

The historical-opacity adapter selects `paper_gs98` or `paper_gn93` through
`FreeEOSOpacity(eos, mesa, opacity_family=...)`. Both require ABI 6 and use
separate MESA handles, pure OPAL Type1 above 7000 K and AF94 below. No table
or sonic coefficient is fitted. The common table rectangle is checked.
The paper's stalled cutoff can be reproduced with `stalled_mach=1.5`;
the default remains sqrt(1.01). The independent second-sonic cutoff remains active even when a higher comparison Mach is never reached. This changes only stopping conditions.

42 runs show that the tested opacity families, abundance convention and
Euler discretization do not account for the published red branch. See
[Report 10](../Reports/Report_10.md) for profile comparisons, conditional
input sensitivity, improved inner coverage and retained failure evidence.

The separate grey/MLT atmosphere now raises `AtmosphereDomainExit` at the
first convective Mach-one event. This implements the explicit validity-exit
option; it is not a turbulence model or evidence that the wind's zero
convective flux assumption is valid at a physical base. Escaping radiation
and physical matching remain unresolved; Stage B was not run.
