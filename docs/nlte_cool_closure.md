# Cool non-LTE material and wind closure

Report 14 extends the CMFGEN route into deliberately prescribed cool H/He
controls. This document defines the implemented material differential and the
remaining physical coupling requirements. **It does not authorize a production
wind.** Numerical results and release decisions are in
[Report 14](../Reports/Report_14.md).

## Component definition and reference boundary

The [predeclared plan](../benchmarks/cmfgen_cool_validation_plan.json) uses the
T9 and T6 luminosity, radius and mass-flux scales in
[Martins, Schaerer & Marques-Chaves (2026)](https://arxiv.org/html/2605.12141v1).
Their models prescribe the flow. Their metal-bearing atmospheres and adjusted
deep scale height are not reproduced by these primordial controls. Their inner
flux temperature is distinct from the emergent effective temperature. Agreement
with those input scales is not validation against their output profiles.
An independent cool reference reproduction remains required.

The controls include H−, 30 HI levels, 69 HeI levels, 30 HeII levels and the
stripped ions, with native full-level statistical equilibrium. Two-photon
processes and Rayleigh scattering are enabled. Electron scattering is coherent
in this energy ABI. X=.7/Y=.3 uses native integer atomic masses. Metals,
molecules, clumping and X-rays are absent. Level dissolution remains enabled;
its thermodynamic effect is not silently represented by an ideal gas.
No temperature floor, opacity multiplier or fitted mass flux is introduced.
The explicit beta-one flow and its 0.01 km/s base speed are component inputs,
not a hydrostatic inner boundary or a wind solution.

The official H− data include a guessed collisional-ionization input
(`HYD/MI/13apr15/col_guess`). These rates require a sensitivity study or
independent validation. A small H− number fraction alone does not establish
negligible opacity. The 50 km/s line-width parameter follows the reference
calculation convention; it supplies no turbulent stress, pressure, energy or
convective transport law.

## Two native defects found and repaired

The [18-Jun-2025 public CMFGEN source](https://sites.pitt.edu/~hillier/web/CMFGEN.htm)
places H0 (H−) before HI when enabled. Its RVTJ output uses fixed ion indices,
so hydrogen and helium population labels become wrong. The read-only
[population patch](../backends/cmfgen/population_export.patch) selects by name
and exports H− separately. The QStarS reader rejects mismatched ion labels.

The Rayleigh opacity and its radiation-field linearization still check the
obsolete HMI name. They therefore use H− instead of neutral HI with H0 enabled.
Even the old HMI branch passes the first ion's level count instead of the
selected HI count. The [physical repair](../backends/cmfgen/hminus_rayleigh.patch)
corrects the name, array stride and population equation index. It changes no
cross-section coefficient. A native test compiles the actual call sites and
checks nonuniform depth populations against the independent low-frequency HI
formula and its derivative. Native build capabilities and the runner reject
the affected configuration without the verified repair. This is a bounded
call-site test, not independent certification of the Rayleigh fit at every
frequency. The standalone upstream display program is not used by QStarS and
is not repaired by this patch.

## Material state and derivatives

[`nlte_material.py`](../src/qstars/nlte_material.py) treats the population vector
as an independent coordinate. Let Y_a=n_a/rho be the number of particles in
atomic level/ion a per gram, Z_a its charge, and E_a its atomic energy. Charge
neutrality gives Y_e=sum(Z_a Y_a). The neutral H and He ground levels define
the atomic energy zero; H− has energy minus its electron affinity.

For the explicitly ideal, nondegenerate resolved-atom component,

```
A = sum[(1 + Z_a) Y_a]
Pgas = rho k T A
ugas = (3/2) k T A + sum[E_a Y_a]
```

These are gas quantities. There is no blackbody radiation energy in ugas,
no radiation pressure in Pgas, and no turbulent energy. Radiation E, F and Prr
remain independent transferred moments. At fixed populations per gram,

```
P_lnrho = P_lnT = P
u_lnrho = 0
u_lnT = (3/2) k T A
dP/dY_a = rho k T (1 + Z_a)
du/dY_a = (3/2) k T (1 + Z_a) + E_a
```

`differential` checks element-nucleus conservation for every population
direction. Electrons follow charge neutrality, including the negative H−
charge. It returns dP, du and the first-law heat `du - (P/rho) dlnrho`.
Tests independently perturb density, temperature, ionization, excitation and
H− detachment, and verify cancellation of constant element energy-gauge shifts.
A ground-state Saha control checks the population chain rule; it does not
replace the non-LTE calculation or certify its rates.

These partials are **not** relaxed derivatives of a two-variable equilibrium
EOS. In particular, the heat capacity through recombination includes
`sum(u_Y dY/dlnT)` and can be much larger than the frozen-population value.
The code provides no automatic conversion of these partials into the local
equilibrium wind operator. Entropy, occupation-probability free-energy terms,
molecular contributions and Coulomb corrections are not supplied by this
ideal component. Native level dissolution affects populations; its consistency
with material thermodynamics remains to be bounded or implemented.

FreeEOS and MESA use their own mass and energy conventions. A constant energy
offset per conserved element cancels in the local first law, but mass-convention
differences, nonideal terms and population relaxation cannot be repaired by
such an offset. No LTE pressure/heat-capacity graft is made onto the CMFGEN
non-LTE state. A production envelope bridge still requires a validated shared
material model in the overlap region.

## Independent energy accounting

For a steady spherical flow with mean density rho and a smooth filling factor
of unity in these controls,

```
rho v [du/dr - (P/rho) dlnrho/dr] = Q0
du/dr = u_lnT dlnT/dr + sum(u_Y dY/dr)
```

Q0 is total comoving absorption minus emission, including atomic transitions.
Do not add an extra ionization heat term when the latent energy already appears
in u. Radiation mechanical work v f0 appears in the lab-frame four-force,
not as a second comoving gas heat source. The independent audit integrates both
the radiation luminosity/four-force identity and the gas first-law identity,
and retains cumulative interior errors as well as endpoint errors.

The initial native controls omit material and population advection. Their
radiative equilibrium is therefore not a steady wind first-law solution.
The upstream `STEQ_ADVEC_V4` and `EVAL_ADIABATIC_V3` contain those terms; they
must be enabled together, audited at relaxation parameter one, and checked
against QStarS's resolved-atom energy before coupling. A prescribed velocity
profile also requires a momentum residual; a constant Bernoulli flux cannot
be assumed without satisfying momentum balance.

## Global response and convection handoff gate

Write the transfer/rate residual as C(x,z)=0, where x contains the radiation
and population state at all radii, and z contains QStarS's trial density,
temperature, velocity and boundary variables. A response to a trial wind
perturbation must satisfy

```
C_x dx = -C_z dz
d(force,heat) = A_z dz + A_x dx
```

This is nonlocal. A diagonal opacity derivative or a pointwise tabulated force
cannot stand in for it. Converged positive/negative profile perturbations,
step-size agreement and finite-difference residual checks are prerequisites
for using the response in a coupled Newton/eigenvalue solve. Nonlinear residuals
larger than the perturbation invalidate that experiment. The full problem must
include mass flux as an unknown, momentum, the material first law, rate/transfer
equations, radiation boundaries and critical regularity. The present local
synthetic-source sonic formula is not a substitute for that global operator.

The matching surface must satisfy the existing independent moving-surface
mass, energy and momentum ledger. Convective enthalpy flux and stress cannot
be discarded at the interface. A radiation-only atmosphere can join an
envelope only through a demonstrated overlap where those omitted terms are
negligible to a declared error budget, or through an independently validated
transport/stress/dissipation closure. No such overlap has yet been established.
The Mach-one MLT exclusion remains active; a flux cap or line-width parameter
does not remove it. Failure to converge the atmosphere is a numerical failure,
not evidence that a physical quasi-star has no wind.
