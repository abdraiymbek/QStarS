# Self-selected outer boundary: equations, evidence and remaining closure

Report 12 outcome: **FAIL**. The production quasi-star still cannot select a validated thermodynamic escaping wind. This campaign adds and tests necessary components; it does not enable a new stellar mass-loss prescription. The Report 10 Mach-one MLT guard remains active. No physical envelope attachment, evolving-star coupling, complete published track or QSRD calculation is authorized by a component pass.

## What the installed microphysics can determine

MESA 24.08.1 supplies a [Rosseland radiative opacity](https://docs.mesastar.org/en/24.08.1/kap/overview.html). FreeEOS supplies equilibrium gas pressure, ionization, energy, entropy and derivatives. Their combination is useful in the already screened LTE diffusion domain. Neither interface supplies frequency-dependent absorption/emission, externally irradiated level populations, or the flux-weighted force in a recombining escaping spectrum.

The means required by a static continuum transfer problem are distinct:

\[
\kappa_P(T_g)=\frac{\int\kappa_\nu^a B_\nu(T_g)d\nu}{\int B_\nu(T_g)d\nu},\qquad
\kappa_E=\frac{\int\kappa_\nu^a J_\nu d\nu}{\int J_\nu d\nu},
\]
\[
\kappa_F=\frac{\int(\kappa_\nu^a+\kappa_\nu^s)F_\nu d\nu}{\int F_\nu d\nu},\qquad
\frac1{\kappa_R}=\frac{\int(\kappa_\nu^a+\kappa_\nu^s)^{-1}\partial_T B_\nu d\nu}{\int\partial_T B_\nu d\nu}.
\]

These scattering expressions assume isotropic coherent scattering. The first two means enter thermal exchange; the flux mean enters force. A zero net flux has no defined flux mean. Signed frequency fluxes are retained in the force integral. The radiation energy temperature \(T_E=(E/a)^{1/4}\) is a diagnostic, not an assertion of a blackbody spectrum or material LTE.

[Malygin et al. (2014)](https://arxiv.org/abs/1408.3377) demonstrate the need for two-temperature weighting and the possibility of several thermal equilibria. Their published metallicities are -0.3, 0 and +0.3 dex; these tables cannot be relabeled primordial. Their level populations also assume LTE. [Mayer & Duschl (2005)](https://arxiv.org/abs/astro-ph/0411613) supply primordial mean-opacity calculations including D/Li and distinguish Planck continuum and line sets. The inspected paper offers the full tables on request; a complete dataset was not acquired here. Single-temperature means would still not close arbitrary irradiated spectra or population kinetics. No external requests were sent.

[SIROCCO](https://sirocco-rt.readthedocs.io/en/latest/operation.html) solves ionization and thermal balance on a specified/imported flow and remains a serious existing backend candidate. Its [line treatment uses the Sobolev approximation](https://sirocco-rt.readthedocs.io/en/latest/physics/radiative_transfer.html), requiring a separate assessment near a slow or nonmonotonic launch. The pinned star example from Report 11 still prescribes mdot and a beta velocity law. It was not installed, run, or coupled in this campaign; an executable thermal solver alone would not validate an envelope eigenvalue.

[Sander et al. (2017)](https://arxiv.org/abs/1704.08698) provide a relevant established route: non-LTE statistical populations, comoving-frame transfer and hydrodynamic iteration with an inner atmospheric constraint. Their sections 5.2–5.3 explain why the critical point alone does not set a consistent stellar solution and why force changes with the mass flux must be recomputed. Their O/WR-star results use different compositions and model assumptions. None is represented as a reproduced quasi-star benchmark here.

The minimum defensible production extension is thus a spectrum-dependent absorption/emission/force and population calculation, with sufficient frequency groups or monochromatic resolution demonstrated by convergence for the selected H/He family. Non-LTE can be omitted only where collision/radiation equilibration and flow-time checks establish the relevant equilibrium domain. A universal number of bands cannot be inferred from the current synthetic tests. Advection of ion populations must be retained if their relaxation is slow. FreeEOS chemical energy must not be combined with non-LTE populations by simply adding recombination energy a second time.

## New thermal transfer component

`thermal_transfer.py` extends the existing spherical rays to specified frequency bands with true absorption, Kirchhoff emission and coherent isotropic scattering:

\[
\boldsymbol n\cdot\nabla I_g=\chi_g^a(B_g(T_g)-I_g)+\chi_g^s(J_g-I_g),\qquad
Q=4\pi\sum_g\chi_g^a(J_g-B_g).
\]

Each band coefficient is constant in frequency within its declared edges. The coefficients are frozen on the radial grid. Band Planck intensities and their temperature derivatives use a 64-point quadrature, separately checked against adaptive integration. Outgoing core light is a specified blackbody; inward core light is absorbed; the exterior has zero incoming intensity. The **net** core luminosity is an output and is not confused with the emitted outward luminosity.

For fixed material temperature, each scattering source equation is solved by GMRES. A local monotonic Planck root and repeated formal transfer solve \(Q=0\) for a static material. This assumes no material advection or mechanical/chemical heating. It does not impose \(T_g=T_E\). Both the local thermal residual and the independent integral

\[
L_{\rm out}-L_{\rm in}+\int4\pi r^2Q\,dr=0
\]

are reported, with no flux renormalization. Direct formal transfer can have a small algebraic residual and an unacceptable spatial flux error; coarse failures remain in the campaign.

Independent controls are exact emitting/absorbing spherical chords and an optically thin opacity \(\kappa_\nu\propto\nu\). In the latter, \(\int\nu B_\nu d\nu\propto T^5\), so a uniformly luminous core gives \(T_g=T_*W^{1/5}\), while \(T_E=T_*W^{1/4}\). The group sequence tests this known spectrum rather than only comparing two implementations of a grey mean. The frequency-independent limit and Planck derivatives have independent regression checks.

**This component does not provide realistic H/He opacities, ionization feedback, moving frequency shifts, non-LTE populations, an accelerating material energy solution, or mdot.** Extinction optical-depth radii are reported separately for each band; none is certified as a thermal photosphere or color radius. Radiative equilibrium with frozen coefficients has a monotonic local root; that does not establish uniqueness for a real opacity that depends on temperature and populations.

## Material energy and the critical condition

`material_wind.py` exposes the gas-only first-law operator using actual FreeEOS derivatives when supplied. Define

\[
c_T^2=\rho^{-1}(\partial P_g/\partial\ln\rho)_T,\quad
b=\rho^{-1}(\partial P_g/\partial\ln T)_\rho,\quad
U=(\partial u_g/\partial\ln\rho)_T-P_g/\rho,\quad
C_T=(\partial u_g/\partial\ln T)_\rho.
\]

For stationary continuity, \(y=d\ln v/d\ln r\) and \(t=d\ln T/d\ln r\),

\[
U(-2-y)+C_Tt=\frac{rQ}{\rho v},\qquad
(v^2-c_T^2)y=2c_T^2-bt-\frac{GM}{r}+rg_{\rm rad}.
\]

Eliminating the temperature derivative gives

\[
a_g^2=c_T^2-\frac{bU}{C_T},\quad
(v^2-a_g^2)y=2a_g^2-\frac{GM}{r}+rg_{\rm rad}-\frac{brQ}{\rho v C_T}.
\]

At a regular critical point both denominator and numerator vanish. The new operator assumes finite local source functions with no explicit velocity-gradient dependence. In particular, \(Q\) is thermal exchange, excluding the mechanical power \(\rho v g_{\rm rad}\). The radiation equations must lose exactly the corresponding thermal and work energy.

For ideal gas, this critical speed is adiabatic; the native FreeEOS gas gamma1 independently checks the derivative combination in recombining states. This is **not** a replacement for the isothermal critical condition when temperature is prescribed, for the LTE diffusion equations that constrain its gradient differently, or for a line force explicitly depending on dv/dr. A fully coupled transfer problem has additional characteristic and nonlocal constraints. Reusing any one of these critical formulae outside its declared equations is not justified.

For two material ODEs with specified local radiation sources, base radius, mass, composition, density and temperature set the thermodynamic boundary; continuity relates base velocity to the one unknown mdot. A regular sonic crossing supplies an eigencondition. Sonic radius is located by the solution rather than independently prescribed. When transfer and populations are solved, their incoming-characteristic/rate conditions must also be supplied. The envelope must provide a consistent net energy flux and inner radiation condition; assigning arbitrary intensity, net luminosity, radiation energy, base velocity and mdot simultaneously would overconstrain the problem. A moving matching radius can only be an extra unknown if an additional physical condition determines it. The current diffusion matcher retains its own three sonic unknowns and three base constraints; it is not the missing non-LTE boundary problem.

## Independent base-selected wind limit

The calorically ideal adiabatic Parker reference has zero radiation, zero external heating, fixed point gravity, \(u=P/[(\gamma-1)\rho]\), and \(1<\gamma<3/2\). This is an explicit mathematical fluid family, not ionized H/He or trapped-radiation-plus-escaping-transfer physics. [Waters & Proga (2012)](https://arxiv.org/abs/1207.7348) discuss the broader Parker topology; this implementation restricts itself to the accelerating subsonic-base family.

Let \(a_b^2=\gamma P_b/\rho_b\), \(\lambda=GM/(r_b a_b^2)>2\), \(h=1/(\gamma-1)\), \(p=h-3/2\), and \(x=a_s^2/a_b^2\). The sonic condition, constant entropy and Bernoulli integral reduce the eigenproblem to

\[
\frac{\lambda^4}{32}x^{2p}+h-\lambda-px=0,
\quad 0<x<\left(\frac4{\lambda^2}\right)^{1/p},\quad
\frac{r_s}{r_b}=\frac{\lambda}{2x},\quad
\frac{\rho_s}{\rho_b}=x^h.
\]

For \(2<\lambda<h\), the residual is convex, starts positive and ends negative in this entire admissible interval: one accelerating root is bracketed. It selects mdot=4 pi rs² rhos as and the base velocity. Two-sided integration uses the independently derived regular derivative; an implicit Bernoulli root at multiple radii checks its density profile. Rematching from another subsonic slice recovers mdot without passing it as input.

For \(\lambda\ge h\), no accelerating root lies inside this reference interval. The compatible hydrostatic solution has

\[
a^2(r)=a_b^2\left[1-(\gamma-1)\lambda(1-r_b/r)\right].
\]

It reaches vacuum at \(r/r_b=(1-h/\lambda)^{-1}\) for strict inequality, or only at infinity for equality. This is a scoped hydrostatic/no-wind classification with an analytic argument. It is not proof of no quasi-star wind, uniqueness in other polytropic families, or stability of a real atmosphere. Inputs outside the implemented range return a domain exit or an explicit invalid-domain exception. A numerical integration failure is kept distinct from this root-existence argument.

## Matching solver robustness

The FreeEOS/MESA diffusion matcher now uses damped Newton steps with a finite trust radius, centered admissible Jacobian stencils and backtracking. Invalid sonic branches, integration failures and parameter-box exits are recorded as rejected trials; they no longer necessarily abort a valid search. Its residual still matches base density, temperature and total outward non-gravitational energy flux. Bounds and tolerances are numerical controls, not physical fit coefficients.

The initial point must be admissible. Loss of Jacobian rank, exhausted evaluations and failed line searches return `WindMatchError` with trial history. None means physical nonexistence. A successful local root does not enumerate all branches or establish escape. Report 12 broad-guess, matching-radius and budget-failure cases use slices of a known wind as numerical witnesses; they are expressly not physical stellar-envelope tests.

## Conservative interface implemented

For interface speed \(u\), area A, common potential Phi, material speed v and realizable lab radiation E,F,Prr, `outer_interface.py` returns

\[
\dot M=A\rho(v-u),\qquad
\dot\Pi=\dot Mv+A(P_g+P_{rr}-uF/c^2),
\]
\[
\dot E=\dot M(u_g+v^2/2+\Phi)+AP_gv+A(F-uE).
\]

The ledger separates gas internal energy, relative gas pressure work, moving-face gas pressure work, kinetic energy, potential energy, lab radiation flux and swept radiation energy. Radiation enthalpy advection is already in F; it is not added again. The lab radiation momentum density F/c² contributes to the swept momentum. For a comoving face, gas mass advection vanishes but pressure work remains. A gas-energy gauge shift changes the flux by exactly mdot times that shift.

Thin, source-free matching equates mass, momentum and energy at identical radius, enclosed mass and potential. Distinct thermodynamic states can satisfy these jumps: the normal-shock regression has a temperature discontinuity. Smooth density, temperature, entropy and enthalpy continuity are additional assumptions, not extra conservation laws. Entropy production/stability and radiation spectral compatibility would need their own admissibility checks; a flux pass is not a physical solution certificate. This interface enforces an explicit common gas/potential energy gauge.

A finite layer additionally has momentum geometry/gravity sources and mass/energy storage. Its two faces cannot be compared directly with the zero-thickness jump helper. A time-dependent central potential requires its source term. No moving cavity or evolving envelope update is implemented by this interface. Outward F at a finite face is not certified escape.

## Convective handoff and execution gates

The existing [literature audit](spherical_wind.md) remains binding. Standard MLT stops at Mach one. Capping its speed or flux would leave an unaccounted energy excess, not solve the handoff. Saturated-convection scalings supply a transport bound but not, alone, a momentum stress, dissipation law or mass-loss eigenvalue. Porosity-length families require explicit structural parameters and independent calibration. Turbulent pressure requires a stress and energy closure; it cannot be added only as extra support. The cited 3D iron-opacity results have a different composition from a primordial atmosphere. Heating-driven trapped-radiation winds require verified trapping and an independently specified energy-deposition mechanism. No coefficient from those families was fitted to the present controls.

A future physical run must jointly pass atomic/thermal/force closure, subsonic convective handoff, sonic branch/matching, exterior transfer and all conservative fluxes. Multiple roots require continuation and branch retention; a failed solver is an unresolved numerical result. Only after a physical static envelope selects a validated branch may the evolving star use its rate. Restart/retry/remesh and branch-hopping campaigns then apply to that coupled model. The current component outcomes do not permit those later stages.
