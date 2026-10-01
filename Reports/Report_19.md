# Report 19 — an analytic tangent of the native formal transfer map

**Outcome: TRANSFER COMPONENT PASS; cool closure remains a NUMERICAL BLOCKER.**

Starting from Report 18, I implemented a directional derivative of the actual
`QS_FRM_ALI` characteristic map and its converged coherent-scattering source.
It carries the full ray response across depth and frequency, including the
inserted-point interpolation and thick outer boundary. It passes the 1e-6
fixture target on 20- and 39-radius grids. A matrix-free scattering solve agrees
with the dense tangent reference and needs at most seven analytic ray products
per frequency in these fixtures.

**This component is not yet connected to the population Newton solve.** The
released 84-radius cool state still has not converged. No energy/electron audit,
physical resolution pair, wind solve or production mass-loss rate is delivered.
QStarS remains **0.1 development, not v1**.

## 1. Delivered code and scope

- `backends/cmfgen/formal_tangent.f90`: reusable forward-mode tangent component.
  `tangent_map` differentiates one formal pass; `tangent_scattering` closes the
  implicit scattering response; `tangent_rsqh` differentiates the native
  midpoint interpolation of r²H.
- `backends/cmfgen/formal_tangent_probe.f90`: finite-difference oracle that calls
  the existing native `CMF_FORMAL_REL_V4` and scattering ALI to convergence.
- `scripts/audit_formal_tangent.py`: builds, runs, hashes and gates the fixtures.
- `tests/test_formal_tangent.py`: a native regression covering all seven fixture
  families, frequency-history ablation and the narrow face-flux linear regime.

The component links **build 17** libraries. There is no new CMFGEN executable
build and no new executable capability flag. The existing local-ALI and moment
population Jacobians are unchanged. This avoids assigning component validation
to an executable that does not use it.

Supported: fixed R, V, angular/ray geometry and frequency grid; coherent
scattering; DIFFUSION and ZERO_FLUX inner boundaries; thin or extrapolated thick
outer boundaries; inserted ray points; native nonnegative linear short
characteristics. Geometry, velocity, Doppler-grid and frequency-grid derivatives,
HOLLOW boundaries and incoherent-scattering redistribution are not implemented.
The wind response will need additional velocity/geometry terms.

The differentiated CMFGEN characteristic/interpolation adaptation carries the
upstream GPLv3 provenance described in the backend README.

## 2. What is differentiated

At one ray point the native code uses

    a = b / dlog(nu), or zero at the first frequency
    chi_eff = chi + (1/dlog(nu) + 3) * b
    S_eff = (eta_total + a * I_blue) / chi_eff

The tangent carries dI_blue explicitly. It also differentiates the trapezoidal
optical-depth increments, the attenuation exp(-dtau), both linear source
weights, and the diffusion boundary B + mu*DBB/chi_inner. The small-dtau
polynomial is differentiated as written, including its truncation. Nothing is
clipped or removed from the primal transfer.

Inserted ray points use the differentiated native Steffen interpolation of
log(chi) and log(eta_total), including its active slope-limiter branch. Thick
outer-boundary extrapolation includes derivatives of its fitted powers, the
scattering/absorption split and active floors. The primal replay agrees with the
native J/H/K/N and face flux to below 1e-10 in every fixture.

For scattering, let A be the derivative of the formal map with respect to total
emissivity at the converged state. The component solves

    (I - A diag(sigma)) dJ = forcing

where forcing includes thermal emissivity, d(sigma)*J, total extinction, blue
ray history and boundary derivatives. A is a derivative at the current state:
with inserted log-interpolated points, the finite-amplitude source map need not
be affine. The dense reference constructs A analytically, not by differences.
The matrix-free option applies the same A with one differentiated ray sweep per
Arnoldi product. The local L* is only its diagonal preconditioner; it does not
replace the nonlocal response.

The returned ray derivatives become the blue-side state **only after** this
scattering solve. They remain fixed during current-frequency products. The
component leaves the native primal ray buffers untouched.

## 3. Fixture validation

Evidence: `Reports/evidence/formal_tangent_20261001/nd20/audit.json` and
`nd39/audit.json`. Each grid uses 24 descending frequencies, six core rays and
six perturbation directions:

1. emissivity at one depth at frequency 9 only;
2. absorption at that same depth/frequency only;
3. electron scattering, including its contribution to total extinction;
4. mixed, signed opacity/emissivity changes over a finite frequency profile;
5. Planck and diffusion-gradient boundary changes;
6. simultaneous thermal-emission and absorption changes.

The ZERO_FLUX boundary perturbation is exactly null and is tested as such.
All nonzero responses are compared over the complete frequency/depth field.
The error for each moment is max|tangent-FD| / max|FD|, including signed H/N.
It is a response-norm error, not a promise of relative accuracy in every
negligibly small cell. The table gives the worst of J, H, K, N and midpoint r²H
at the declared gate perturbations, over all six directions.

| Fixture | 20 radii: worst error | 39 radii: worst error |
|---|---:|---:|
| Smooth, thin outer boundary | 2.11e-7 | 3.89e-7 |
| Sharp opacity front | 1.53e-8 | 1.99e-8 |
| Opacity front + inserted ray points | 1.18e-8 | 2.69e-8 |
| Thick outer boundary + inserted points | 1.17e-8 | 3.68e-8 |
| Strong absorption, thick + inserted | 2.93e-7 | 2.82e-7 |
| Scattering-dominated, thick + inserted | 1.98e-7 | 4.04e-7 |
| ZERO_FLUX + inserted points | 1.30e-8 | 1.56e-8 |

**All fourteen component gates pass.** The largest matrix-free/dense disagreement
is 6.28e-14. The scattering tangent fixed-point defect stays below 2.2e-12.
The native primal scattering residual stays below 1e-13. Maximum ray lengths
are 55 points at 20 radii and 67 at 39 radii, so inserted and extended points are
actually exercised.

Matrix-free scattering needs at most 7 products per frequency, rather than
constructing 20 or 39 emissivity columns. A product is a ray sweep at one
frequency. **These are not the full-model Newton–Krylov products of Report 18**;
no 2–4-vector population GMRES result or full-model speedup is claimed.

### Frequency history is measured, not assumed

A second tangent sweep deliberately resets the blue derivative at every
frequency. In the 20-radius inserted-front blue-emissivity test, this raises the
J response error to 2.25e-2 and the face response error to 9.67e-2. Other gated
fixtures reach 0.164. The full-history response remains below 1e-6. The bluer
profile therefore matters even in these bounded models.

The existing local-ALI `TX(:,:,:)=0`/`TVX(:,:,:)=0` defect from Report 18 is
confirmed but **not repaired in the native Newton branch by this component**.
Simply preserving the old moment matrices would not retain the directional
angular intensities required by the characteristic equation. The native
adapter must retain the ray response as well as update TX/TVX correctly.

## 4. The face-flux linear regime can be much narrower

The initial standard steps were epsilon = 1e-3, 1e-4 and 1e-5. Larger-step
failures remain in the evidence. In the strong-absorption mixed-direction
fixture, J/H/K/N pass at 1e-5, but the face-flux error is **0.153 at 20 radii**
and **0.104 at 39 radii**. Nearly cancelled H is about 1e-8–1e-9 near the inner
boundary. Perturbations change its shape enough to switch the active
`MON_INTERP` slope limiter. A small fractional coefficient perturbation is
therefore not necessarily a linear perturbation of the face flux.

At epsilon = **1e-8**, the same analytic derivative agrees with the face
finite difference to **2.93e-7 / 2.82e-7**. At 1e-9, roundoff starts to dominate;
for the 20-point mixed fixture its face error increases to 4.68e-6. The
strong-absorption boundary direction also needs the smaller step.

Conversely, the blue-only thermal responses in the scattering-dominated
fixture, and the scattering response in the strong-absorption fixture, are too
weak for the smallest standard differences. Their error grows as epsilon
shrinks. Their gate uses 1e-3, which resolves the response and passes 1e-6.

The driver has an explicit per-case/direction gate-step table, fixed from the
20-radius study and then used unchanged at 39 radii. It does **not** choose the
lowest measured error at runtime. It retains the complete step-size series,
including the limiter-crossing and roundoff failures. This is a local tangent
validation. It is not evidence that a large nonlinear Newton step will work.

## 5. Gate decision and next work

| Gate | Status |
|---|---|
| Formal-transfer tangent, nonlocal in space and frequency, fixed geometry | PASS in these fixtures, with the response-resolving gate steps above |
| Analytic tangent scattering closure | PASS; dense and matrix-free references agree |
| Native population/rate-equation Jacobian | NOT INTEGRATED. Report 18's e >= 0.94 and uphill directions remain the latest native measurements |
| Released cool populations + T + radiation | NOT CONVERGED; no new physical release run |
| Electron and radiation-energy audits | NOT EVALUATED on a converged cool state |
| 84/167 radii and doubled-frequency physical resolution pair | NOT ATTEMPTED. The 20/39 fixtures are component tests |
| Global mass-loss / sonic / self-selected wind | NOT ATTEMPTED; no production mdot |

The narrowest blocker is now **the native population-to-radiation-to-rate
adapter**, rather than an absent characteristic tangent implementation. Next:

1. Supply directional thermal emissivity, total extinction, scattering and
   diffusion-boundary changes from the existing material/atomic derivatives,
   in CMFGEN's population and temperature conventions. Retain a complete
   directional blue-ray state through the native frequency loop.
2. Accumulate the resulting rate-equation product with the direct material
   terms. Either use one global direction at a time for matrix-free Newton,
   or batch the active population/depth directions to fill TX/TVX. The latter
   must account for current coefficient slots, continuum/line history,
   clumping and the native face-flux/Gamma conventions. No native memory or
   full-spectrum runtime claim has been measured here.
3. Test that complete product along the existing frozen-J step using
   `jacobian_direction_check.py`. It must reproduce the Report 18 cancellation
   and achieve e much less than one at the measured native noise floor.
4. Only then use it in the outer Newton–Krylov solve with damping in log
   variables and a cool-layer merit function. The narrow face-limiter regime
   must be handled in that line search; this report does not change interpolation
   or declare a smooth global Jacobian at limiter ties.
5. After an actual released 84-radius convergence, run the electron/energy
   audit, the physical resolution pair and the global wind eigenvalue solve
   with the existing Mach-one MLT guard.

Nothing here indicates that the cool solution fails to exist physically.
The new evidence is about the discrete transfer derivative and its local
linear regime, not a converged stellar or wind state.

## 6. Reproduction and housekeeping

From the repository root, with the recorded build-17 libraries present:

```sh
export MESASDK_ROOT="$HOME/qsrd-runtime/mesasdk"
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MACOSX_DEPLOYMENT_TARGET=14.0
.venv/bin/python scripts/audit_formal_tangent.py --nd 20 --output results/new_tangent_nd20
.venv/bin/python scripts/audit_formal_tangent.py --nd 39 --output results/new_tangent_nd39
.venv/bin/python -m pytest -q
```

The driver links actual native libraries with bounds checking in the new
component/probe. It records compiler, library, source and native executable
hashes. Archived evidence contains every input, stdout and compressed numerical
response, plus the full audit records. Build products remain under ignored
`results/formal_tangent`. No earlier report's run directory was changed.

**Tests: 101 passed, 23 skipped**, with the compiler initialized. The skips are
the existing opt-in MESA/FreeEOS integrations; the new native CMFGEN regression
ran and passed. Nothing was committed to git. No native model run is left
running. No energy/electron, wind or mass-loss claim is inferred from these tests.
