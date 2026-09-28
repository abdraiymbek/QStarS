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
