# CMFGEN-selected mass-loss rate and sonic point (Gate 4)

Status: **component implemented and tested; gate not closed.** `src/qstars/cmfgen_wind.py` solves a
self-selected hydro step for mdot and the sonic point from a frozen CMFGEN radiation field. It runs no
native code. Its inputs must come from a converged and audited T-released closure, and none exists yet
(Reports 21–22). Nothing here selects a stellar mass-loss rate, opens `cmfgen_backend.production_outer`,
or adds a `winds.rate` mode. Every classification carries `wind_certified=False` and
`physical_impossibility_claimed=False`.

The plan `benchmarks/cmfgen_selected_wind_plan.json` was registered before any result it governs. It
fixes the equations, inputs, classification order, thresholds and outer fixed-point policy.

## Equations of one hydro step

Continuity and momentum for a steady spherical flow with gas pressure \(P=\rho a^2\):

\[
\dot M=4\pi r^2\rho v,\qquad
(v^2-a^2)\frac{d\ln v}{d\ln r}=N\equiv 2a^2-\frac{da^2}{d\ln r}-\frac{GM(1-\Gamma)}{r},\qquad
\frac{d\tau_R}{dr}=-\kappa_R\rho .
\]

- \(a^2=(n_{\rm nuclei}+n_e)kT/\rho_{\rm clump}\) comes from the CMFGEN populations. With \(T\) frozen as a
  function of the coordinate this is the isothermal-type critical speed. It is **not** the adiabatic speed
  of `material_wind.continuum_coefficients`, and the material first law is not solved here
  (`docs/self_selected_outer.md`, "Material energy and the critical condition").
- \(\Gamma=r^2g_{\rm rad}/(GM)\), where \(g_{\rm rad}\) is the local comoving force of the spectral export,
  `SpectralSnapshot.moments(coherent_scattering=True)['radiative_acceleration']`.
- \(\kappa_R=\chi_R/\rho_{\rm clump}\) from RVTJ. \(\tau_R\) of a state is integrated independently: a trapezoid
  in \(\ln r\) of \(\chi_R f r\), plus an outer power-law tail.

At the critical point, \(v_c=a_c\) and \(N=0\). The branch slope comes from `spherical_wind.regular_sonic_slope`
in \((\ln r,\ln v,\ln\tau)\), with \(d\ln\tau/d\ln r\) supplied. A unique positive slope is required, with
\(D=v^2-a^2\) increasing outward.

### Frozen representations

| Mode | Frozen on | Role |
|---|---|---|
| `tau` | \(\ln\tau_R\) | **Primary.** \(da^2/d\ln r=(da^2/d\ln\tau)(d\ln\tau/d\ln r)\) depends on \(v\) and \(\dot M\). At \(v=a_c\), \(r_c=[GM(1-\Gamma_c)-(da^2/d\ln\tau)_c\kappa_c\dot M/(4\pi a_c\tau_c)]/(2a_c^2)\) in closed form. Unknowns \((\ln\dot M,\ln(\tau_c-\tau_o))\) are solved by admissible damped Newton from an explicit seed (`tau_seed`). Every trial is kept. |
| `radius` | \(\ln r\) | **Cross-check and existence proof.** \(N(r)\) is independent of \(v\) and \(\dot M\). Every sign change of \(N\) on the whole base-to-outer domain is enumerated (32 samples per table interval, then `brentq`). \(\dot M\) follows linearly from the base condition. |

The two modes coincide at a self-consistent fixed point, where the closure's \(\tau(r)\) equals the hydro
\(\tau(r)\). Their agreement there is a declared acceptance check. Tables are monotone PCHIP of
\(\ln a^2\), \(\Gamma\) and \(\ln\kappa_R\). Evaluation outside a table raises an error, and no value is
extrapolated or clipped. A trial that leaves the table is a rejected local trial, never a result.

Subsonic layers with \(\lambda=GM/(ra^2)\sim10^4\) give \(\dot M\sim e^{-10^4}\), and the supersonic optical depth
\(\tau_c-\tau_o\) can be many orders below \(\tau_o\). Both modes therefore carry \(\dot M\) only in logarithms:
the radius mode's inner \(\tau\) integral and shell mass, and in the τ mode \(\ln[\kappa\dot M/(4\pi rv)]\).
The τ mode's outer branch carries \(w=\ln(\tau/\tau_o)/\ln(\tau_c/\tau_o)\), with \(\ln\ln(\tau_c/\tau_o)\) in closed
form, and its inner branch carries \(\ln\tau\).

- **Where the table coordinate is a dependent variable** (τ mode), a table exit at a trial Runge–Kutta stage
  makes DOP853 reject and shrink that step. No accepted step uses an out-of-table value.
- **Newton iteration.** Newton steps on the smooth event-radius residuals \(\ln r(\tau_o)-\ln r_o\) and
  \(\ln r(\tau_b)-\ln R_b\). It stops on the reported forms: \(\ln\tau(R_b)-\ln\tau_b\) to first order (larger by
  \(|d\ln\tau/d\ln r|\sim\lambda\)), and \((1+z_c)w(r_o)\), which bounds \(\ln\tau(r_o)-\ln\tau_o\).
- **Tests.** Tests at \(R_b=r_c/40\) (\(\tau_c-\tau_o\approx10^{-30}\), \(\tau_o=10^{-3}\)) and at \(R_b=r_c/7000\)
  (\(\lambda(R_b)=1.4\times10^4\), \(\ln\dot M\approx-1.39\times10^4\), \(\dot M\) underflows to 0) require the two modes to
  agree with each other and with an independent quadrature of the implicit Parker root to \(10^{-7}\) in
  \(\ln\dot M\).

**Non-regular sonic approach.** A branch can approach \(v=a\) with \(N\ne0\) without crossing it. This is the fold
\(D\to0\), the only singularity of the step's ODE, and it makes DOP853 stall. A stall at
\(|v^2/a^2-1|<10^{-3}\) rejects that crossing, as a sonic re-encounter already did. Any other integration failure
leaves the candidate `candidate_integration_failure`. Enumeration then finishes, every candidate is recorded,
and the step is `numerical_failure` (rule 7). The regression test checks such a fold independently with the
Bernoulli criterion \(\int_{r_X}^{r_o}N\,d\ln r<0\).

### Constraint count

There are three unknowns: \(\dot M\) and the integration constants of \(\ln v\) and \(\tau\). They are fixed by
three constraints:

1. a regular accelerating critical crossing;
2. \(\tau_R(r_o)=\tau_o\) at the outer condition depth;
3. \(\tau_R(R_b)=\tau_b\) (primary), or \(\rho(R_b)=\rho_b\) (rematch cross-check).

The sonic radius and \(\dot M\) are outputs, and no velocity or \(\dot M\) is input. This is the counting of
`docs/self_selected_outer.md` lines 53–82, with the transfer/population and first-law constraints moved
into the CMFGEN re-convergence.

### Inputs and provenance

- **Mass.** \(M\) is a declared component input on every call. The current CMFGEN component uses
  VADAT `[MASS]` = \(10^6M_\odot\) (`backends/cmfgen/prepare_cool.py`), which is not physically sourced.
  \(\dot M\) is exponentially sensitive to it.
- **Base and outer depths.** \(R_b,\tau_b\) are taken at 1-based depth ND−2 of the first converged closure,
  and \(r_o,\tau_o\) at depth 2. Both are held fixed through the outer iteration. Depth ND, the diffusion
  inner boundary, is excluded from the tables. Depths 1 and ND−1 are interpolation margins only.
- **Composition.** \(X=0.7015\), from the CMFGEN H/He number ratio 9.333.
- **Not inputs.** The β-law, VTURB, any force multiplier and any floor are excluded.

## Classification (first matching rule; plan `classification_order`)

1. `reduction_invalid`: an input check fails (\(|v|/c\ge0.01\), spectral vs RVTJ state above \(10^{-3}\),
   or guard composition). A guard "pass" also counts here unless it guarded at least one depth with
   both opacities (pending amendment; stricter default enforced).
2. `mach_one_mlt_exit`: the existing Mach-one MLT guard trips.
3. `reduction_invalid`: MLT would carry more than 1 % of the flux, which CMFGEN does not carry.
4. `reduction_invalid`: radius-mode \(N(R_b)\ge0\), no admissible X-type crossing (including a branch that
   re-encounters or folds back to the sonic surface), or a tangency.
5. `hydrostatic_no_wind`: \(N r/(GM)<-10^{-6}\) on the whole domain. Every regular stationary solution is
   then subsonic and accelerating, so \(d\ln P/dr\le -g_{\rm eff}/a^2\). With the same base normalization,
   whether by \(\tau_b\) or by \(\rho_b\), \(P(r_o)\le P_{\rm hs}(r_o)\) and \(\dot M<4\pi r_o^2P_{\rm hs}(r_o)/a(r_o)\).
   This is scoped to that frozen state and must be confirmed on a re-converged structure.
6. `multiple_admissible_branches`: requires continuation.
7. `numerical_failure`: an integration, Newton or budget failure, or an unresolved critical-point
   candidate. A transonic root also lands here if it is above the declared tolerances:
   - \(|\ln\tau(R_b)-\ln\tau_b|>10^{-8}\);
   - regularity above \(10^{-6}\);
   - pending amendment, stricter default enforced: the τ-mode \(|\ln\tau(r_o)-\ln\tau_o|>10^{-8}\).

   This class is never physical.
8. `photon_tiring_exit`: \(\dot MGM/(R_bL)\ge1\).
9. `reduction_invalid`: a critical-point screen fails. The screens are thermal relaxation over flow time,
   photon trapping \(v_c\tau_c/c\), the nonrelativistic profile, shell-mass fraction and GENCOOL closure,
   each at most \(10^{-2}\). State quantities are read at the step's own critical point, never clipped:
   - in τ mode at \(\ln\tau_R=\ln\tau_c\), where \(a^2,\Gamma,\kappa\) were frozen (before the fixed point the
     closure's \(r(\tau_c)\) is not the hydro \(r_c\));
   - in radius mode at \(\ln r_c\).
10. `stationary_transonic_wind`: a component outcome only.

**Thermal relaxation screen.** The screen uses \((\tfrac32P+E_{\rm atomic})/(\tfrac12\sum|\text{GENCOOL rows}|)\).

- **A proxy, not a bound.** GENCOOL rows are per-process nets, so their absolute sum is a lower bound on
  gross exchange. The screen is still only a **gross-exchange proxy**, not a conservative bound, because
  the rate that actually relaxes \(T\), \(T\,d\Lambda_{\rm net}/dT\), can be well below the gross exchange.
- **Adiabatic rows (pending amendment).** The two `Adiabatic cooling rate` rows are flow work under
  `INC_AD`, not radiative exchange. Including them biases \(t_{\rm th}/t_{\rm flow}\) toward passing. The plan
  amendment drops them from the gross-exchange estimate (`radiative_gross_exchange`) and keeps them in the
  channel-closure sum. This amendment is pending user approval. Meanwhile the screen must pass with **both**
  the registered and the radiative-only estimate.
- Without the resolved atomic energy (`load_wind_inputs(levels=...)`), the screen cannot pass.

**Mach-one guard wrapper.** `mach_one_guard` applies `atmosphere_reference.AtmosphereDomainExit` semantics
at subsonic CMFGEN depths with \(\tau_R\ge2/3\). It uses LTE FreeEOS, native MESA MLT (Henyey or ML1,
\(\alpha=2\)), \(g=GM/r^2\) and the local comoving flux.

- The opacity choice on a non-LTE structure is **pending user approval**. Meanwhile the default is `both`
  (CMFGEN Rosseland and MESA LTE), and the guard trips if either trips, so it cannot weaken the guard.
- The guard is tested with stub EOS/MLT objects. It has not yet been run with native libraries on a CMFGEN
  state.

**Ledger.** `wind_ledger` evaluates the steady total energy flux through static base and outer faces with
`outer_interface.moving_flux`, using export lab moments. The energy gauge is neutral H/He ground. The gas
Bernoulli constant and \(v/v_{\rm esc}\) are reported, and escaping radiation is not certified.

## Outer fixed point (policy fixed in the plan; driver not yet written)

Each iteration takes \(x=\ln\dot M\) of the closure and computes \(G_k=\hat x_k-x_k\) from the τ-mode step.

- **Step rule.** The first step is damped by 0.5. Later steps use the secant, unless it leaves the latest
  sign-change bracket, in which case they fall back to the damped step.
- **Trust region.** Each step is limited to \(\ln 10\), and any truncation is recorded.
- **Convergence.** \(|G_k|\le3\times10^{-4}\) and \(|\Delta\ln r_s|\le3\times10^{-4}\), with the closure's
  momentum residual at most \(10^{-3}\).
- **Budget.** 12 iterations. Running out is a `numerical_failure`.

`outer_update` implements exactly this rule.

- **Inputs.** It takes the four-column history \((x_k,\hat x_k,\ln r_{s,k},\ln\hat r_{s,k})\) and the latest
  closure's `momentum_passed`. Rows with fewer columns raise.
- **Status.** It returns `converged` only when all three registered criteria pass. If \(\dot M\) and \(r_s\)
  are at the fixed point but the momentum criterion fails, the status is `mdot_converged_pending_checks`,
  which is not converged.
- **Tests.** It is tested on oscillating and expanding maps and on a synthetic density-dependent-force
  feedback, where it matches a direct `brentq` root to \(10^{-8}\).

Repeats at the fixed point:

- restarts from 3× and 1/3× \(\dot M\);
- the base moved to two other subsonic depths (\(|\Delta\ln\dot M|\le10^{-5}\));
- 333 radii and doubled frequency on the **self-selected** structure. The bound is
  \(\max|\Delta\log|\le10^{-3}\), the existing spherical-wind refinement gate applied unchanged. Relaxing it
  to 1 % is listed as pending approval and is not applied.

## Input-chain audit (Report 23 evidence, read-only)

`scripts/audit_cmfgen_wind_inputs.py` writes `Reports/evidence/report23/grad_input_audit.json`.

**Native definition.** Native HYDRO `g_RAD` is \(\kappa_FL(\mathrm{ND})/(4\pi r^2c)\). It uses the **inner**
comoving luminosity (`hydro_terms_v5.f:214`), not the local force.

**Raw comparison.** On 8 same-pass archived exports (T9 84/333-radius cool controls, T6, and hot 84/167/333
meshes), the raw ratio of spectral to native force differs by 14–124 %. This is a definitional difference.

**Corrected comparison.** After scaling by \(L(r)/L(\mathrm{ND})\), they agree to at most
\(5.3\times10^{-5}\), against a pre-registered threshold of \(10^{-3}\). This is ES14.4 print precision. The
agreement also confirms that the export density equals native POP_ATOM·μ̄·amu. Other checks on the same
exports:

| Quantity | Agreement |
|---|---|
| Flux-mean opacity (export vs RVTJ) | \(\le5\times10^{-8}\) |
| GENCOOL channel closure | \(\le7.2\times10^{-5}\) |
| Independent \(\tau_R\) vs native Tau (which uses max(Ross, e.s.)) | within 0.9 % |

The diffusion boundary depth agrees to \(2.3\times10^{-5}\) on these states. It stays excluded as declared and is
reported separately.

**Negative control.** A `USE_FIXED_J` run with stale HYDRO is rejected as `not_same_pass`.

**Illustrative only.** The archived states are prescribed-flow, non-converged controls, so these figures
are not evidence:

- On the cool controls the momentum residual has a median of 0.58–0.62 and a maximum of up to 1.0 of the term-magnitude sum (native % error median ≈198), as expected for a β-law.
- The frozen radius-mode \(N\) is negative at every T9 depth, with no critical point in the domain.
- The T6 state has two sign changes.

## What remains

1. A converged, audited 167-radius T-released closure (Gate 1). The first force must come from its
   verify pass (`RUN/verify`: QSTARS_SPECTRAL.bin, HYDRO, GENCOOL).
2. The resolved atomic energy (`levels`) for the thermal screen and ledger. The guard must also be run with
   native FreeEOS/MESA (bridge ABI ≥ 5) at the closure's X.
3. Driver tooling:
   - the outer fixed-point driver;
   - non-nested state transfer (departure coefficients and T in \(\ln\tau_R\));
   - a restart writer for RVSIG_COL/VADAT that strips `QSTARS_*TRACE*` files.

   Every outer iteration needs a new converged T-released closure.
4. Limitations:
   - The τ mode is a local Newton from an explicit seed. It enumerates nothing, so existence and branch
     counting come from the radius mode.
   - The τ-mode stage rejection relies on DOP853 rejecting a step whose error norm is NaN. This is the
     behaviour of the installed SciPy (1.18.1), and it is exercised by the \(r_c/40\) and \(r_c/7000\) tests.
   - If the isothermal-type reduction fails its screens, a coupled characteristic analysis is required
     (`docs/nlte_backend.md`), and it is not designed here.
5. Physical base data from the QStarS envelope and the convective handoff (Gate 5).

## Tests

`tests/test_cmfgen_wind.py`, about 10 s on one thread, covers:

- **Isothermal Parker, constant Γ.** \(r_c=GM(1-\Gamma)/(2a^2)\), slope 1, the implicit solution to
  \(5\times10^{-8}\), and \(\dot M\) against an independent quadrature of the implicit root to \(10^{-6}\).
- **τ vs radius mode.** Agreement to \(10^{-7}\), with base rematch by \(\tau\) and by \(\rho\) to \(10^{-5}\). The
  \(R_b=r_c/40\) and \(r_c/7000\) cases are also checked against an independent quadrature to \(10^{-7}\).
  Radii outside \([R_b,r_o]\) are rejected in both modes.
- **Variable \(a^2(r)\) and \(a^2(\tau)\).** An independent momentum integral to \(10^{-6}\) and a τ quadrature.
- **Γ crossing 1.** The analytic root and slope, a Bernoulli integral, and Owocki et al. eq. 24 plus explicit
  isothermal work.
- **Hydrostatic bound.** Analytic for both normalizations, and above the true transonic rate. Inverted
  optical depths are rejected as input.
- **Non-regular sonic approach.** A Γ bump gives \(N\) signs \(-,+,-\). The result is rule 4, with both
  critical points recorded and the X-point's fold checked by Bernoulli. Without the stall diagnosis the
  result is rule 7, with the unresolved candidate recorded.
- **Parsers and residual.** A fixed-width HYDRO parser with overflow and 3-digit exponents, GENCOOL channels,
  and the momentum residual on an exact hydrostatic layer.
- **Guard, classification and policy.** Mach and flux exits of the guard, every classification rule, the
  declared residual and regularity tolerances (which must equal the plan's), vacuous and single-opacity guard
  passes, the outer policy including the momentum criterion, the ledger and the g_rad audit definition.
- **Screens.** The τ-mode lookup at \(\tau_c\) on a closure whose \(\tau(r)\) differs from the hydro one, raising
  outside the table, and the radiative-only gross exchange.
- **Archived export.** Loads an archived export, skipped if absent.
