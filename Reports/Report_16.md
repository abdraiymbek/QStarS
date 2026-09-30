# Report 16 — signed moment propagation and rejected closure repairs

**Outcome: NUMERICAL BLOCKER. The Report 15 cool radiation instability is not
resolved. No validated closure, self-selected wind or production mdot is delivered.**

The first retained He I-edge negative state is now decomposed exactly. Its
forcing is dominated by the divergence of **previous-frequency H**, with signed
third-moment contributions. The spatial J-advection diagonal is positive and
near one at that event. Independent linear algebra reproduces the negative
solution. Replacing the lagged J derivative alone cannot repair it.

A shared-face, separate-node K/N reconstruction and its native nonlocal VAR
response were implemented and tested. The response passes a bounded native
finite-difference fixture. The fully requested cool trial still produces a
negative radiation state, before reaching the He I edge. The repaired frozen
edge operator also stays negative. A second, unclipped N/H representation fails
its one-way projected-coefficient requirement. Neither is accepted as a closure.

The remaining numerical requirement is **realizability-preserving coupled
angular/spatial/frequency transport, with its relaxed formal-factor response**.
The present projected moment propagation does not preserve nonnegative angular
radiation. This is not a physical nonexistence result, a proof that all moment
methods fail, or a claim that the He I atomic data are wrong. QStarS remains
**0.1 development, not v1**.

## Scope, inputs and controls

The [declared plan](../benchmarks/cmfgen_moment_closure_plan.json) retains Report
15's 0.1% nonlinear, 1% electron balance, 0.1% global/cumulative radiation-energy,
1% internal-profile and 1e−6 derivative targets. All original H-minus/HI/HeI/HeII
levels, native rates/dissolution, two-photon, Rayleigh, coherent scattering,
radiation advection and relativistic terms remain. No physical coefficient,
frequency endpoint, composition, velocity law, radius, luminosity or component
mass flux was tuned. Prescribed component mass flux remains a diagnostic input,
not production mdot. Material advection is still outside the initial equilibrium
component's scope; a wind would require it and independent accounting.

Read-only reproductions retain the failed Report 15 implementation, including
its native J/H replacement options, solely to diagnose that failure. New closure
candidates use `CHK_J=CHK_H=CHK_XM=NONE`. The shared-face candidate additionally
records and stops on a negative raw J before any native replacement. No clipped
state, frozen radiation solution or relaxed tolerance is accepted. Failed restart
checkpoints are explicitly diagnostic states, never equilibria.

## Exact He I-edge seed

Build 9 captures all coefficients and previous moments before native J
replacement, at the first negative raw solution in the declared
5.9–6.0e15 Hz window. Build 10 captures first moment iterations across
5.94–6.10e15 Hz. Both are read-only source instrumentation; their patch and
executable hashes are retained.

Write `x = gamma r² J` and `y = gamma r_mid² H`. At
**nu = 5.941825295672274e15 Hz**, row 3 of the captured system `A x = b` has:

| RHS contribution | Value in native scaled-equation units |
|---|---:|
| Local material source | +0.2110393421 |
| Previous-frequency H divergence | −165.6515943230 |
| Previous-frequency J term | +36.5389173649 |
| Left previous-frequency K/N term | +171.1471510519 |
| Right previous-frequency K/N term | −119.0654945162 |
| Total | **−76.8199810803** |

The raw solution is **x₃ = −14.0257423695**; this is the scaled variable
`gamma r² J`, not an unscaled intensity value. The independent band solve
agrees to **1.84e−16**, with componentwise backward residual **1.81e−16**.
The J-advection factor there is **P_J = 1.0013885816**, so a negative lagged
advection diagonal does not seed this event. The preceding half-face H is
−139.1651, while its neighbouring preceding x values are +244.7606 and +189.6213.

The native elimination explicitly supplies
`−HS[i−1] y_prev[i−1] + HS[i] y_prev[i]` and the signed EPS differences to
this RHS. The linear solver is accurately solving an already negative forcing
problem. This is the exact local mechanism for this retained seed, rather than
an inference from fixed global outputs. See the
[complete independent attribution](evidence/moment_closure_20260930/seed_analysis_v2/analysis.json).

A diagnostic substitution using positive formal angular H removes the negative
solution. It is **only a sensitivity experiment**, not a permitted overwrite or
an accepted repair. Making just J advection implicit via a shared upwind beta*x
flux still gives minimum x = **−13.78749**. That diagnostic retains native
boundaries and does not certify a whole-domain conservative replacement.

The edge-window sequence also shows inconsistent projected angular states
before this negative J event. Under a monotone reconstruction of x at a shared
face, the implied H/K/N at **6.0996488714e15 Hz**, face 4, require
x_face ≥ **4216.216**, while the neighbouring x values are **2185.382 / 4.0983**.
This follows from the positive-semidefinite conditions for the moment matrices
of `(1 ± mu) I(mu)`. No monotone x interpolation between those nodes can satisfy
that projected state. This diagnoses the staggered reconstructed moments; it
is not a claim that the native formal-ray angular intensities themselves obey
that particular face interpolation. The first frame in the selected window
already fails this test, so the absolute earliest loss over the entire spectrum
has **not** been located. Association with the He I edge is not proof that the
edge is the only source of loss.

Once native replacement injects the formal-ray value into a negative moment
state, subsequent frequency propagation is no longer the unmodified discrete
solution. Report 15's retained growth operator shows huge inherited-frequency
forcing. The current evidence identifies the signed seed and rejects candidate
repairs; it does not yet reconstruct every earlier angular/scattering iteration
or prove a unique full-spectrum runaway mode.

## Implemented shared-face candidate and matching response

The [experimental source patch](evidence/moment_closure_20260930/node_moment.patch)
reconstructs `gamma r² K` and `gamma r² N` at a midpoint using common linear
weights and **distinct left/right node contributions**. It replaces the
single coefficient multiplying `(x_i + x_{i+1})` by separate coefficients for
`x_i` and `x_{i+1}`. Each spatial face supplies one H flux to both neighbouring
zeroth-moment balances, so its interior divergence telescopes.

The native `VAR_JREL_V5` operator is changed with the same separate node
coefficients. Existing `EDD_JREL_VAR_V3` and `UP_TX_TVX` propagate their opacity,
emissivity and preceding-frequency variations without collapsing them back into
one coefficient. This is a numerical reconstruction experiment, not an
upstream bug allegation, changed atomic physics or a supported model family.
Boundary work, integrated energy and full relaxed formal-factor derivatives
remain separate required gates.

An actual native five-depth/two-frequency fixture retains coherent scattering,
relativistic terms and radiation advection. Its formal angular factors come
from a nonnegative two-angle distribution. Holding those factors is **only a
linearization test**, not a frozen-radiation production closure. Central finite
differences of eta/chi at one depth test the nonlocal response at all five depths
and across both frequencies:

| Actual native response | Maximum relative FD error | Target | Result |
|---|---:|---:|---|
| dJ/deta | 3.22e−8 | 1e−6 | PASS, fixed-factor fixture only |
| dJ/dchi | 8.23e−8 | 1e−6 | PASS, fixed-factor fixture only |

[Fixture inputs, compiler command and results](evidence/moment_closure_20260930/native_response_replayed/validation.json)
are preserved. [Source reapplication](evidence/moment_closure_20260930/node_source_replay.json)
reproduces both compiled modified routines exactly. This validates a bounded
Newton response; it supplies neither temperature/population relaxation nor
velocity/boundary/formal-factor sensitivities for the wind operator.

## Candidate failures and remaining requirement

Replaying the retained edge state with the new node reconstruction still gives
minimum x = **−5.30101**. It improves that particular value but does not remove
the signed moment-propagation failure.

The requested full-cycle cool trial starts from Report 15's settled fixed-J/T
material initialization, asks for fresh formal factors and releases J/T in the
inputs. Native startup still holds temperature for the first restart cycle.
That cycle stops after **6.42 s**, at
**nu = 86.1308836907e15 Hz**, row 18, with raw
**x = −2.7413667287e−244** and RHS **−5.0048712925e−240**. An independent solve
reproduces the negative sign/value with componentwise residual **1.36e−16**.
The very small absolute value is retained; no floor, clipping, spectral-tail
truncation or waiver is applied. The record demonstrates a remaining
non-positive discrete operator in that tail. It does **not** establish a
macroscopic energy defect or a converged cool closure. There is no fresh final
RVTJ/spectral pair or completed released temperature/population cycle.

The separate N/H diagnostic retains the actual formal ratio without the native
G bounds. It fails after **0.72 s** at nu = **999.6664359048e15 Hz**, face 1:
N/H = **1.0399201853**, projected H frequency coefficient W = **−0.8643643071**.
The one-way projected representation cannot simply be assumed damped there.
This is a coordinate/closure failure, not proof that physical photons reverse
frequency direction. No bound on G is used to hide it. That route is rejected,
not offered as a model family.

Together, the unchanged seed, failed implicit-J control, failed node
reconstruction, correctly differentiated bounded fixture and failed unclipped
N/H representation isolate the unresolved requirement: the **coupled moment
flux must preserve angular realizability across spatial and frequency faces**.
A fully implicit angular frequency transport and its nonlocal response are a
scientifically motivated direction; the
[Hauschildt & Baron primary reference](https://arxiv.org/abs/astro-ph/0401164)
discusses such characteristic discretizations and instability tradeoffs. It is
not an independent validation of this CMFGEN mechanism, and no mixing
coefficient from that paper is fitted here. No replacement meeting the full
requirement has passed in this campaign.

## Gate decision

| Gate | Status | Evidence boundary |
|---|---|---|
| Exact retained negative-J forcing attribution | PASS | Complete scaled recurrence and independent solve; absolute earliest full-spectrum ancestor remains open |
| Conservative shared-face K/N implementation | CONDITIONAL | Interior face cancellation; full energy/boundary/realizability validation fails |
| Matching native nonlocal eta/chi response | PASS | Fixed-factor two-frequency fixture only |
| Complete relaxed nonlocal Newton response | FAIL | Population/T/velocity/formal-factor response not validated |
| Stable finite nonnegative released cool closure | FAIL | New native candidate rejects negative J before edge; edge replay remains negative |
| Joint population/T/radiation/electron/energy convergence | FAIL | No complete synchronized equilibrium or energy certificate |
| Two converged spatial/frequency resolutions | FAIL | No converged base solution; no convergence order inferred |
| Global mdot/sonic, handoff and physical selected wind | FAIL | Invalid prerequisite closure; not attempted; Mach-one MLT guard retained |
| General v1 / QSRD | FAIL | Independent gates remain open |

The [campaign matrix](evidence/moment_closure_20260930/gate_matrix.json) and
[general matrix](evidence/completion_gates.json) retain these scope boundaries.
No physical nonexistence, one-assumption CONDITIONAL completion or
SELF-SELECTED WIND PASS is claimed.

## Reproducibility and bounded performance

Six native trials are preserved, including an IN_ITS `COMP_F` placement that
was a no-op; only its corrected VADAT trial is labelled a fresh-factor control.
Original manifests are unchanged. A response-driver metadata path failure is
also preserved, with a corrected replay in a fresh directory. Executed
runner/adapter versions remain available; the newer explicit-failure guard is
not retroactively assigned to those runs.

The [native archive](../benchmarks/results/2026-09-30-moment-closure/verification.json)
verifies **615** selected files; **942** intermediate files remain hashed local
assets. Exact read-only patches, experimental patches, native build hashes,
source snapshots and independent analyses are retained. The final native-enabled
suite passes **123 tests**, zero skips. The prior CLI validation's 29/31 result
and two CB24 discrepancies remain unchanged; no new broad CLI/science campaign
is performed for this narrow operator task. The source registry has **59** entries.

The six native processes total **96.94 s** of measured wall time on the Apple
M5 / 16 GiB MacBook. This is a process-time sum, not total campaign elapsed time,
a scaling test or a physical wind runtime. No memory/runtime limit ends these
trials. Cloud/Linux compute would not change the reproduced signed forcing;
no cloud execution or Linux agreement is claimed. The
[release verification](evidence/moment_closure_20260930/release_verification.json)
links the retained failures, final source, build, archive and regression evidence.
