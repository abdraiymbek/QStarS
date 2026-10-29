# Cool non-LTE material and wind closure

Report 14 extends the CMFGEN route into deliberately prescribed cool H/He
controls. This document defines the implemented material differential and the
remaining physical coupling requirements. **It does not authorize a production
wind.** Numerical results and release decisions are in
[Report 14](../Reports/Report_14.md) and
[Report 15](../Reports/Report_15.md) and
[Report 16](../Reports/Report_16.md). Report 15 repairs native H-minus FF
initialization/Jacobian defects. Report 16 identifies the signed previous-frequency
H/K/N forcing of a retained negative-J seed and tests a shared-face reconstruction
with a matching native nonlocal eta/chi response. That response passes only a
fixed-formal-factor fixture; the cool reconstruction remains negative and is
rejected. No relaxed population/temperature/formal-factor response or joint cool
closure has passed.

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

Report 15 also audits the native H-minus free-free contribution. Its tabulated
coefficient gives `chi = T exp(ell(log T, log nu)) ne nHI_ground`; `ell` is the
unchanged native bilinear log-table interpolation. Within a temperature table
interval, `dln chi/dln T = 1 + d ell/dln T`. Independent population/electron
partials follow the two density factors. `eta = chi Bnu(T)`, so its temperature
partial includes both this coefficient slope and
`dln Bnu/dln T = x/(1-exp(-x))`, with `x=h nu/(k T)`.

The repaired continuum and electron-heating Jacobians include these terms.
Cooling uses the same coefficient and ground-H population, with native
`4pi * 1e-10 * chi * (Bnu-Jnu)` per frequency: positive cooling, opposite to
QStarS's positive material heating. The table and its stimulated-emission
convention are unchanged. No additional stimulated factor is multiplied into
chi. The derivative is piecewise; at interpolation knots a single two-sided
derivative is not promised. Four-frequency native checks at 3300, 6100 and
11300 K pass, including independent first-entry cooling and the last-depth
path. These are FF component derivatives, **not** a relaxed non-LTE
population-energy, radiation-force or wind-eigenvalue response.

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

## Exact Jacobian action and Newton–Krylov closure (Report 20)

The cool closure residual R(x) is the native rate/radiative-equilibrium system
with the converged formal-ALI J(x). Its Jacobian action along a fractional
direction s (dPOPS = POPS·s) is assembled as

```
J·(x∘s) = J_dir·s + Σ_ν (∂R/∂J_ν) δJ_ν(s)
```

* δJ_ν is the analytic, depth- and frequency-nonlocal tangent of the formal
  transfer (Report 19), driven by the exact directional changes of extinction,
  emissivity, coherent scattering and the diffusion boundary, with blue-ray
  history. ∂R/∂J_ν are the coefficients of the native accumulation (bound-free
  QFV, low-T recombination, line ZNET, RE row, two-photon).
* J_dir is the fixed-J Jacobian. It is block-diagonal in depth. The native
  analytic blocks omit the dependence of level dissolution / occupation
  probabilities on n_e, total ion density and T; those columns are completed
  by frozen-J central differences on replay passes (no transfer noise).

Validated against exact-residual finite differences on the 84-radius state
x0 (Report 18 oracle and new T/n_e/population-direction oracle): agreement is at
or below the oracle's own step-to-step consistency in every tested direction
(≈1e-6 for T/n_e/population directions; 1.3e-3 along the native fixed-J Newton
step). Caution: a single-iteration CMFGEN restart inherits FIXED_T from the
restart's BAMATPNT, so `FIX_T=F` alone does NOT release T (the radiative-
equilibrium row stays the identity). Those T-direction runs validated the T
column only; T-released passes must remove BAMATPNT/BAMAT
(`cool_newton_krylov.release_t`).

The nonlinear solve (`scripts/cool_newton_krylov.py`) uses left-preconditioned
GMRES (preconditioner = complete J_dir, which weights every rate equation by its
own fixed-J rates), positivity-preserving log updates x←x·exp(−λs), a per-component
log-step limit (|λ s_i| ≤ ln 10⁴; affects only ultra-trace levels) and a parallel
line search on ‖M_k⁻¹R‖. Convergence is declared on the unlimited exact Newton
step over levels holding more than 1e-10 of their element. The equation
selection is forced to a fixed table so the residual definition never changes.

With T released, the RE row of the action is accurate in the cool layer (~0.4 % of
the net response). In the thick interior, where the RE row cancels 30–300×, the net
response is off by up to 16–25 %. The cause is that the native opacity derivatives
omit the level-dissolution dependence; with dissolution off (diagnostic only) the RE
row agrees to 2e-5–7e-5 (Report 20 §4.3).

Solver options added for the T-released problem (Report 20 §4.4), all recorded per
run and none changing the convergence criterion:
- `--augment` (GCR-type subspace with previous Newton steps);
- `--lm-trust` (Levenberg–Marquardt trust region on |s_T| inside the Krylov
  subspace; trials along the LM path);
- `--ultra-frac/--ultra-weight` (third significance tier for levels ≤ 1e-20 of
  their element);
- `--natural` (Deuflhard monotonicity test on the simplified Newton correction,
  estimated in the Krylov subspace).

`--merit rowscale` (with `--eq-floor`) and `--trace-sweep` are retained for
provenance. The trace sweep is wrong for CMFGEN, because trace rows can carry the
species-conservation equation.

Status (Report 20): the fixed-T problem converges linearly, with the largest
significant frozen-J correction at 1.4e-3. The T-released problem does not converge
at 84 radii. The coupled Newton direction prefers a cool neutral layer whose
ionization front (d11/d12) falls between grid points, and that cell has no
consistent RE temperature (`scripts/thermal_curve.py`: the frozen-J local balance at
d12 is net heating from 1.1 to 7.8 kK). The next steps are the 167-radius released
solve, T continuation, and the omitted wind energy terms (`INC_AD`, `INC_ADV`).

**Report 21 additions.**
- Builds 25–27 complete the level-dissolution derivatives in the tangent opacities and
  emissivities. They cover W in the bound-free LTE populations, the dissolved fraction T2 below
  each edge, and the line ratio `L_STAR_RATIO = W_up/W_low`, all chained through
  `B/A/X_LEV_DIS(n_e, POPION, T)`. With them, every row of the complete action matches
  T-released exact-residual oracles at the oracle floor, the RE row included.
- `INC_AD=T` and `INC_ADV=T` are on. The direct Jacobian is block-tridiagonal (native
  off-diagonal blocks kept; frozen-J columns corrected for them), validated by an odd-depth
  perturbation test (7.8e-6).
- T continuation (`--t-homotopy`) has three modes:
  - `--hom-mode pin`: R_T(λ) = λ·R_RE + (1−λ)·c·(T−T_fix)/T_fix. Badly scaled, because c is the RE
    row's gross scale and T stays pinned until λ ≈ 1. Kept as a failed attempt.
  - `--hom-mode shift`: R_T(λ) = R_RE − (1−λ)·R_RE(start). The Jacobian is unchanged, and the
    predictor–corrector advances λ after `--stage-max-iters` corrector iterates.
  - `--hom-mode ptc`: pseudo-transient continuation, G_T = R_RE − δ·c·(T−T_k)/T_k with T_k the
    current base T. This is the physical heating direction, and the term vanishes at fixed points.
    δ0 is automatic (`--ptc-target-dT`) or reset on restart (`--ptc-delta-reset`). SER on the
    RE net/gross is decrease-only, and δ → 0 when the regularized step falls below `--ptc-switch`.

  Frozen-J columns are shared per state across stages. `--pc-t-floor` adds a fraction of the RE
  row scale to its T diagonal, in the preconditioner only. `--split` gives T its own Krylov
  directions per depth region, which is needed for accurate T steps.
- Report 21 result: the shift continuation reaches λ = 1 at 84 and 167 radii. Free Newton then
  slides the H I front down the cool branch on both grids. At the front the coupled net heating
  is ∝ T at frozen populations (exact probe of the T columns,
  `scripts/probe_columns.py`), so Newton targets T → 0. PTC heats the front from 5.5 kK and from 7.9 and
  10.5 kK hot starts (`scripts/equilibrate_at_T.py --vadat`), but slowly. The iterates remain far
  from statistical equilibrium: at the front, CMFGEN's GENCOOL thermal balance and the RE row
  disagree by orders of magnitude (`scripts/gencool_breakdown.py`). The released closure has not
  converged.
- Build 28 (`jv_nonest.patch`) runs the tangent's direction loop outside any OpenMP construct
  when it cannot occupy every thread. This removes a nested-team libgomp deadlock and an
  intermittent SIGSEGV; products are unchanged to 1.2e-16.
- A 333-radius nested restart (`scripts/prepare_mesh333.py`) and its base pass are ready for the
  167/333 spatial check. `scripts/transfer_state_nested.py` carries a state onto it, and
  `scripts/post_convergence.sh` chains trace cleanup, verification and audits.
- Completed passes are compacted by `scripts/slim_pass.py`. It removes only source-identical
  copies and outputs that nothing reads, and logs the removals in `QSTARS_SLIM.json`.

**Report 22 additions (solver speed; same model, products and convergence test).**
- Builds 29–31:
  - The observer-frame formal solution is skipped in solver passes (`QSTARS_NO_OBS_FORMAL`, set by
    the driver; verification and final runs keep it).
  - One explicit per-frequency scattering operator serves every tangent direction's GMRES.
  - Forcing/final sweeps are batched per thread block, and `lin_setup` runs in parallel.
  - Products agree with build 28 to ≤ 1e-12 per equation row, with identical GMRES step counts.
  - Matched on the same 167-radius state: base pass 2.5×, Jv product 1.9×, column pass 3.4× faster.
- Driver:
  - Thread budgets: 6 column workers × 1 thread; trials `8 // n_trials`.
  - Exact reuse of same-base products after a PTC δ change or switch, counted against
    `--krylov`.
  - Resume hardening; an optional full-budget guard (`--full-krylov-near`).
  - `scripts/bench_pass.py` reruns any archived pass and compares it row by row.
  - `scripts/compare_iterate.py` compares whole iterates.
- Measured structure of the stalled Krylov residual (0.46 after 74 products at 167 radii): 53 %
  is in the front RE row and 43 % in ultra-trace He II/He III rows. The product span already
  contains the front T correction: weighted on the T row alone, it reduces that residual to 0.016.
  The 1e-8 ultra-trace weight declines it. A Schur-complement T solve is therefore not justified.
  Ultra-trace weight 1e-10 (user-approved, with `--full-krylov-near 10`) reaches η in one pass,
  26 instead of 74 products.
- Background QoS: processes launched from the agent session share roughly 4–6 cores. Battery
  and foreground use slow them 1.5–5×, so benchmark only in back-to-back pairs.

**Report 23 additions (what the solver could not see, and the physics it was solving).**
- **Two-photon (build 32).** Every QStarS single-iteration pass sets `COMP_F=T`, which made CMFGEN replace
  `TWO_METH=USE_RAD` by `OLD_DEFAULT` (`comp_opac.f`). That mode gives the partner photon the occupation
  exp(−hν_B/kT_e): 33–101× too much 1s→2s absorption at the H I front, about 75 % of the RE row there, energy
  inconsistent. With `QSTARS_TWO_PHOT_J` (driver `--two-phot-j`) `USE_RAD` is kept and the conjugate-frequency J
  comes from a stored-J file. The driver replays each base with its own stored J (`<base>_tp`), so a converged state
  satisfies the equations with its own J; columns, products and trials of that base use the same file.
- **Trapped transitions.** The fixed-J preconditioner M holds untrapped line and continuum rates. In optically
  thick transitions the radiation response cancels them (H I n = 2 at the front: exact diagonal 3e-4–6.5e-4 of the
  fixed-J one), so ‖W M⁻¹ R‖ understates the error of every such level (H I n = 1–9, He I ground and ¹P series,
  He II n = 1–3; `Reports/evidence/report23/trapped_levels/`). States the old measure called converged were
  not: the PTC states had b₁ up to 7e10 at the front. Convergence in the fixed-J measure is not convergence.
- **Driver options** (`scripts/cnk_trap.py`, `tests/test_cnk_trap.py`; all off by default):
  - `--ali-precond`: M + (BA[`QS_LIN_FAC 1`] − BA[`QS_LIN_FAC 0`]), the native local-ALI response of every column,
    from one frozen-J replay per base (`<base>_ali`). Build 33 makes that pass cost ~1.25 frozen-J passes. Best
    measured preconditioner (median q = 0.40 against 0.88), nearly singular at the two ionization fronts.
  - `--probe-rows/--probe-stride`: exact block-tridiagonal columns of chosen levels from one batched product pass;
    preconditioner only (in the least-squares span the colored probe directions produced a depth sawtooth).
  - `--sig-rows`: levels significant at every depth regardless of population fraction (H I n = 2).
  - `--update-map linear`: x(1 − s) for increases (exp(−s) overshoots large increases).
  - The coupling-significance rule (`--sig-coupling`) is kept for the record and must not be used: it promotes
    ultra-trace levels through ionization-equation rows.
- **Native stops.** `--vadat STOP_IF_MW F` keeps CMFGEN's flux-accuracy and inner-optical-depth checks as
  warnings instead of ending a trial pass; a final state must still satisfy them.
- **Fixed-T versus T-released.** At the transferred 84-radius fixed-T state the trapped-level solver lowers the
  statistical-equilibrium measure threefold in two iterates, after which the fixed-T populations can no longer carry
  the luminosity (L(r)/L(ND) 0.04–0.34 in every trial): the populations must be relaxed with T released.

**Report 24 additions (constraints, far-field leakage, trial masers).**
- **Exact linear constraints.** The H II row is hydrogen number conservation, the He III row helium number
  conservation and the n_e row charge balance; their residuals are ΣH − H₀, ΣHe − He₀ and the charge sum − n_e
  with the totals of x0. The per-component limiter (`--component-log-cap`) clips populations separately and breaks
  them. Every Report 23 solve drifted off them (He up to 8×, charge up to 34×). `--conserve` projects every trial state
  back (H and He totals per depth scaled to x0's, n_e = charge sum; ratios within a species, T and H⁻ kept);
  `--lams ... 0` adds the projected base as a trial. Always check conservation of a reported state
  (`Reports/evidence/report24/collect_page_data.py` prints it).
- **Preconditioner floor.** `--pc-svd-floor 1e-2`: singular values of each equilibrated diagonal block below 1 % of the
  largest are raised to that floor (`cnk_trap.floor_blocks`). Without it, the nearly singular local-ALI blocks at the
  ionization fronts amplify the far-field action of distant directions, and the least squares uses them as blunt
  instruments (iterate 1: a deep-layer direction with coefficient 172, 97 % of its action at the H front).
- **Fixed T.** `--freeze-t` removes T from every direction in the `--pin-all` stage (otherwise the weakly weighted pin
  row lets the least squares move T: coefficient 358 at iterate 4, trial T at 0.50–1.22 of T_fix).
- **Native trial failures** are supercritical infrared masers in the trial states (inverted H I lines under the
  `SRCE_CHK` floor, CHI < CHI_SCAT, Λ spectral radius > 1), not a transfer-solver defect. Build 36
  (`ali_newton.patch`) solves the scattering fixed point by Newton where the local ALI stalls; it converges the
  frequencies that have a solution and leaves a `QSTARS_ALI_NEWTON` record (and, with `QSTARS_ALI_NEWTON_DUMP`, the
  Jacobian) for those that have none.
- **Measures.** The driver's natural measure is dominated by the largest step components and accepted states that the
  fixed measure (`scripts/fixed_metric.py ... --extra`) calls worse; judge progress in the fixed measure, by depth region,
  together with conservation.
