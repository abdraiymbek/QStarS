# Supported model families and implementation boundaries

QStarS remains a development framework. A named prescription identifies its
source; it does not claim reproduction of the source's entire stellar model.

| Selectable capability | What executes | What is still missing |
|---|---|---|
| `ideal_gas_radiation` | Consistent hot, fully ionized ideal gas plus radiation thermodynamics | Ionization, molecules, degeneracy, evolving composition |
| `mesa` EOS and opacity | Direct MESA 24.08.1 H1/He4 EOS, ionization and radiative opacity; independent QStarS ledger | Coupled atmosphere/full MLT, nuclear/composition evolution and production stellar backend; excluded inconsistent blend and dilute-gas states |
| Native MESA atmosphere/MLT reference | Fixed/varying Eddington atmosphere, pressure inversion, Henyey/ML1 and automatic derivatives | Native starting-pressure prescription is not a validated low-gravity vacuum boundary |
| Direct FreeEOS 3.0.0 reference | Gas EOS1, ionization, derivatives and explicit blackbody radiation; MESA opacity accepts the same electron state | Global EOS certification, arbitrary composition, safe stitching and evolution coupling |
| Finite-pressure grey/MLT atmosphere reference | Hydrostatic layers, full local MLT, optical-depth profiles and mass/energy inventories | Supersonic cool candidate and extended warm candidate fail validity; spherical transfer and conservative envelope/wind coupling absent |
| `fixed_excision` | Fixed inner radius with separate swallowed enthalpy, BH rest-mass conversion and central-potential work | Moving physical cavity, its mass/energy reservoir and pressure work |
| `bra2008` | BRA08 equation 3 luminosity law | Complete BRA structure and evolution |
| `ball2011_printed` | Ball's printed equation 15 with its stated coefficient | Moving Bondi boundary and complete Ball tracks |
| `quest2025_eq11` | Later QUEST/Santarelli equation 11, with distinct rest accretion and BH growth | Upstream MESA cavity and saturated inner-flow evolution |
| `prescribed` engine | Explicit rest-accretion control and efficiency | A self-consistent inner-flow solution |
| `off`, `prescribed`, `dotan2011` wind | Separate mass/enthalpy removal and an independent energy-to-escape debit; optional published launch/tiring formula | Porous atmosphere, resolved sonic transition and radiation-enthalpy wind coupling |
| `fiacconi2016_wind` reference | Stationary constant-opacity wind ODEs with radiation advection, photosphere diagnostic and explicit base Mach closure | Coupling to the evolving envelope; gas-pressure sonic eigenvalue selection; full evolution replication |
| `nakauchi2017_lte_diffusion` reference | Spherical gas-sonic regularity, FreeEOS/MESA recombination, radiation/gas enthalpy, component base matching and flux interface | Escaping radiation, validated envelope match, turbulent transport, published escaping branch and evolution coupling |
| Spherical grey angular transfer (Python reference API) | Static VEF and moving coherent-scattering transfer; E, F, Prr, angular closure and radiation work | Non-grey absorption/flux means, material energy/populations and coupled stellar boundary |
| Constructed continuum wind (Python reference API) | Pressure-free fixed-grey stationary momentum plus escaping radiation, supplied mdot/base velocity | Recombination, gas sonic eigenvalue, physical envelope matching; illustrative dimensional scaling is not a star |
| Ball loaded-polytrope reference | Published n=2,3,4 dimensionless structure | Thermal evolution |
| Corrected Bondi reference | Cavity-inclusive dimensionless equilibrium limits | Thermal evolution and moving cavity accounting |
| Saturated-convection reference | CB24 inner/outer ODE matching, with both documented cutoff interpretations | Full evolutionary family |

The three accretion luminosity laws remain separate. The factor-of-four
coefficient disagreement is not averaged or silently corrected. The static
Bondi and saturated limits are boundary-dependent results, not universal BH
mass limits or destruction criteria.

The standalone wind alternative follows Fiacconi & Rossi's gas-pressure-free
model. Its base Mach number is visible in the YAML (default 0.8); the paper's
range 0.75–0.85 is propagated in figure comparisons. Its sonic crossing is not
a gas sonic critical point that independently determines the wind rate.
Radiation, advected enthalpy, kinetic energy and gravitational energy are all
retained. No fitted coefficient is used. The two Figure 4 configurations are
in `benchmarks/fiacconi2016_*_wind.yaml`; the benchmark retains discrepancies
between the ODE solution, plotted diagnostics and caption parameters.

Separate effective feeding/wind streams are supported by the evolution
kernel. Combining those streams into one net rate would lose energy and
entropy information, so that reduction is not performed. Their simultaneous
multidimensional geometry cannot be resolved in spherical 1D. This is a
physical validity boundary, not a software defect that a finer grid fixes.

The spherical gas-sonic component is selected with
`qstars reference benchmarks/spherical_recombining_wind.yaml`; it is kept
separate from Fiacconi and Dotan. See [its domain and equations](spherical_wind.md).
A completed reference integration is never converted automatically into
stellar mass removal. A local critical point, a resolved decelerating branch,
a photosphere comparison endpoint and a validated terminal wind are distinct.

## Report 10 additions and restrictions

The spherical wind's optional `paper_gs98` and `paper_gn93` opacity adapters
use MESA's Type1 OPAL and AF94 with the paper's 7000 K switch. They are
attribution references, not new production wind families. The mass-number
abundance experiment lives only in the campaign driver. None reproduces the
published escaping branch. All original alternatives retain their own assumptions.

Moving-face fluxes now include relative advection and explicit surface work;
no evolving cavity or physical matching result follows from that interface.
The grey/MLT atmosphere raises `AtmosphereDomainExit` at Mach one. No invalid
supersonic profile is returned. See [Report 10](../Reports/Report_10.md).

## Report 11 radiation scope

[The new angular transfer](escaping_radiation.md) is a separate explicitly grey
family. It does not replace FreeEOS/MESA recombination in the gas-sonic wind.
Pure scattering returns no gas temperature, entropy, enthalpy or ionization;
its effective temperature and tau=2/3 radius are radiation diagnostics. The
constructed stationary model prescribes mass flow and base velocity. No new
production evolution boundary is enabled by these component results.

## Report 12 component scope

- Frequency-band static absorption/emission/scattering with a solved material radiative-equilibrium temperature and four distinct mean opacities. Coefficients are synthetic and frozen; H/He populations and an accelerating energy equation are not coupled.
- A base-selected adiabatic Parker verification family, with no mdot/base-speed input and a proven hydrostatic alternative in its restricted interval. Its calorically ideal fluid is not a primordial quasi-star wind.
- Gas-only FreeEOS derivatives in a local continuum momentum/energy operator; a general moving-interface ledger using lab E/F/Prr without an LTE assumption.
- More robust diffusion matching that retains rejected trials and still does not certify a physical envelope or escaping branch.

All are separate Python reference APIs. No new production stellar wind family is enabled. See [the equations and validity boundary](self_selected_outer.md).

## CMFGEN diagnostic component (Report 13)

The optional external CMFGEN interface supplies imported H/He populations,
spectral coefficients, moments, heating and force for a **prescribed trial
flow**. The supplied metal/clumped reference, reduced H/He thermal iterations,
and frozen-state relativistic transfer tests are separate configurations.
They are not selectable production mass-loss families. No opacity/force
coefficient was fitted to a quasi-star target. Standard and advective transfer
options remain explicitly distinct; their different energy residuals are
retained. [Selection, units and validity](nlte_backend.md).
# Report 14 cool component scope

The external CMFGEN controls now explicitly include H−, HI, HeI, HeII and
stripped ions with two-photon processes and corrected neutral-H Rayleigh
scattering. T9 and T6 are prescribed primordial stress cases at literature
input scales, not new selectable production wind families. Their cool
nonlinear/thermal/energy validation fails. The resolved-atom material API
provides population energy and frozen-state partial derivatives; it does not
provide a relaxed EOS or global sonic closure. See
[the cool contract](nlte_cool_closure.md). The existing physical model families
and the Mach-one MLT guard retain their earlier validity restrictions.
