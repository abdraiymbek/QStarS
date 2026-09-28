# Outer-boundary validation and the v1 release blocker

**The current atmosphere candidate fails the general-use gate.** It is an
executable diagnostic, not an enabled evolution boundary. The failure is a
specific physical validity conflict: a converged cool model requires supersonic
mixing-length convection, while a warmer model is geometrically too extended
for the plane-parallel transfer assumption. No failure is called a stellar
disruption, and no impossibility theorem for spherical 1D is claimed.

## Physics and implementation

The reference uses direct **FreeEOS 3.0.0 EOS1** for H1/He4, with gas-only option
`(ifoption, ifmodified, ifion)=(3,11,-2)`. Ionization and thermodynamic
derivatives come from the upstream free-energy calculation. Blackbody
radiation is added with the same radiation constant used by MESA. There is
no fitted EOS, Saha replacement, or interpolation between these EOS libraries.
Composition is X=0.7, Y=0.3 in the campaign. The evolution interface retains
its independently screened MESA EOS; FreeEOS is presently a reference only.

MESA 24.08.1 supplies radiative opacity using electrons and degeneracy from
that same FreeEOS state. It also supplies full Henyey or ML1 convection, with
native automatic derivatives. Alpha=2 is explicit; the Henyey constants are
nu=8 and y=1/3. Conductive opacity is not photon opacity. These calls are
process-local and must not be run concurrently in threads.

For constant gravity and outward total flux F=a c Teff^4/4, the grey radiative
layer obeys

\[
 T^4={3\over4}T_{\rm eff}^4(\tau+2/3),\qquad
 {dP\over d\tau}={g\over\kappa},\qquad
 {dP_g\over d\tau}={g\over\kappa}-{F\over c}.
\]

Below the match at tau=2/3, the native MLT gradient replaces the imposed grey
temperature gradient:

\[
 {dT\over d\tau}={T\over P}\nabla_{\rm MLT}{g\over\kappa},\qquad
 {dP_g\over d\tau}={g\over\kappa}-{4aT^3\over3}{dT\over d\tau}.
\]

The explicit outer gas pressure is an input at the stated optical depth.
It is never obtained by clipping density. Negative gas pressure rejects the
branch. Depth, column mass and LTE energy are integrated as

\[
 {dz\over d\tau}={1\over\kappa\rho},\quad
 {dy\over d\tau}={1\over\kappa},\quad
 {dE_A\over d\tau}={u\over\kappa}.
\]

These are **inventories**, not a moving-boundary conservation ledger. Evolving
this layer requires its mass, energy, gravity and boundary work to enter the
envelope/wind accounting. The small layer mass is not permission to omit its
energy. Total-flux and actual radiative Eddington factors are saved separately;
the latter uses F_rad/F = nabla_MLT/nabla_rad in the diffusion region.

## Why merely calling MESA's native atmosphere is insufficient here

The native boundary audit uses X=0.75, Y=0.25. The unmodified varying-opacity
T(tau) routine starts at tau_base*1e-5 with gas
pressure calculated from a fixed density `RHO_OUTER=1e-10`. This is a finite
starting pressure. With constant opacity 0.34, Teff=8000 K and g=3 cm/s²,
its photospheric gas pressure is **56.6504 dyn/cm²**, whereas the native
fixed-opacity analytic vacuum solution gives **0.717479 dyn/cm²**: a factor
**78.9576**. At 6000 K and g=1 the factor is **129.449**; at solar gravity it is
only **1.000742**. This is a boundary-condition difference, not a convergence
failure or a claim that all MESA atmospheres are defective. Finite external
pressure is legitimate when physically specified.

The separate 77-state low-density audit has 60 neutral ideal-backstop states.
At rho=1e-16 g/cm³ and T=6000 K, the native EOS returns free_e=1.66e-28 and
opacity 6.65e-29 cm²/g. Direct FreeEOS returns free_e=0.68785, hydrogen ionization
0.99033 and, through the same MESA opacity library, opacity **0.27555 cm²/g**.
A small total first-law error did not detect this lost ionization.

There is also a retained discrepancy inside the nominal MESA FreeEOS domain.
At rho=1e-9 g/cm³ and T=1e6 K, MESA gas pressure is 145307.689 dyn/cm²; direct
FreeEOS gives 134193.864. The dilute classical-gas screen differs by **8.282%**
while MESA's total first-law defect is **1.04e-8**. The pressure offset scales
approximately as 4.407e-6 times radiation pressure across the tested density
sequence. The historical table builder subtracts a radiation pressure from
the library's total pressure; a historical constant mismatch is a plausible
explanation, **not proven by available generation metadata**. No fitted
subtraction is applied. Evolution now rejects dilute states exceeding a 1%
gas-pressure discrepancy when eta<-10, T>=40000 K and the conservative
charge-2 Coulomb parameter is below 0.01. In this regime omitted gas
corrections are much smaller than the rejection threshold.

MESA's table loader also explicitly floors free_e at 1e-4. The largest
electron discrepancy in the 153-state overlap campaign is at a cool state
where direct FreeEOS gives 6.71e-8. This floor and the neutral fallback are
distinct mechanisms.

## Candidate atmosphere results

The 14-case campaign uses a 10000-solar-mass enclosed mass. Radius follows
R=sqrt(GM/g); it is not adjusted to force the thin-layer check to pass.
Tau ranges from 1e-6 to 100; baseline outer Pg=1e-6 dyn/cm².

| Quantity | Teff=8000 K, g=3 | Teff=6000 K, g=1 |
|---|---:|---:|
| Radius for geometry screen, cm | 6.651e14 | 1.1520e15 |
| Photospheric Pg, dyn/cm² | 1.26632 | 142.148 |
| Photospheric density, g/cm³ | 1.3644e-12 | 3.6855e-10 |
| Full integrated depth / radius | 0.08568 | 0.01079 |
| Photosphere-to-tau=100 depth / radius | 0.06676 | 0.001574 |
| Maximum convective Mach | 0.17563 | **2.33672** |
| Optical depth at peak convective Mach | 100 | 19.7878 |
| Blob optical depth at peak Mach | 206.0 | 1596.2 |
| Maximum actual radiative Eddington factor | 1.11745 | 7.75768 |
| Layer mass / enclosed mass | 7.565e-5 | 2.148e-4 |

The 8000 K layer's 6.68% depth below the photosphere alone implies a variation
in g and flux of order 13%; constant-g/constant-flux plane-parallel transfer
cannot certify one-percent structure there. In the 6000 K case, both Henyey
and ML1 yield Mach 2.33672. Standard MLT's pressure-equilibrated, subsonic blobs
are not a defensible physical solution at that speed. The failure is in an
optically thick region, not just the optically thin boundary.

Changing rtol from 1e-7 to 1e-9 leaves maximum gas-pressure/density log-profile
differences below 2.12e-5. It leaves these physical failures intact. Reducing
outer optical depth or pressure tenfold changes photospheric density by less
than 2e-5 relative, but can materially change the outermost density and depth;
those whole-profile differences are retained. A physical external Pg=1
dyn/cm² makes the cool layer thinner, but it still requires Mach 2.33224.
Changing to constant opacity removes that cool convective failure while
strongly changing the photospheric density; doing so would discard the
recombination opacity being tested, not validate the production model.

The constant-opacity analytic gas-pressure error is below 9.3e-9 of the
photospheric value. Tightened-tolerance hydrostatic integral residuals are
9.45e-13 and 9.03e-10, respectively. These numerical checks do **not** override
the geometric and convective validity failures.

## Release decision and what would resolve it

**FAIL — production outer stellar boundary.** Promoting this candidate would
require either ignoring a demonstrated validity violation or adding an
unvalidated closure. Neither is allowed. The strongest supported claim is a
screened component/reference calculation plus restricted hot-envelope evolution.
This cannot be narrowed to a general-purpose quasi-star v1.

Resolution requires a validated spherical outer-layer treatment and a
conservative dynamical/wind or explicitly sourced alternative closure where
local MLT ceases to apply, followed by the actual coupled evolution and
cross-paper campaigns. Published porous-wind and radiation-advecting wind
families remain alternatives, not interchangeable parameter fixes. LTE grey
transport also does not establish non-LTE ionization or flux-mean opacities
near a scattering photosphere. The available 1D non-LTE atmosphere literature
demonstrates that spherical 1D itself is not disproved by these failures.

## Reproduction and provenance

Install FreeEOS 3.0.0 separately from its
[official release](https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/).
The downloaded tarball SHA-256 is
`4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09`.
It is a GPL dependency and is not redistributed as MIT QStarS source or binary.

```sh
python -m pip install -e '.[test,freeeos-build]'
python scripts/build_freeeos.py --source /path/to/free_eos-3.0.0 \
  --sdk-dir /path/to/mesasdk --output backends/freeeos/build_local --arm-kind-fix
python scripts/build_freeeos_bridge.py --freeeos-build backends/freeeos/build_local \
  --output backends/freeeos/build_bridge_local
python scripts/build_mesa_microphysics.py --mesa-dir /path/to/mesa \
  --sdk-dir /path/to/mesasdk --output backends/mesa/build_outer_local
python scripts/audit_freeeos.py \
  --library backends/freeeos/build_bridge_local/libqstars_freeeos.dylib \
  --mesa-library backends/mesa/build_outer_local/libqstars_mesa.dylib --mesa-dir /path/to/mesa
python scripts/validate_atmosphere_reference.py \
  --freeeos-library backends/freeeos/build_bridge_local/libqstars_freeeos.dylib \
  --mesa-library backends/mesa/build_outer_local/libqstars_mesa.dylib --mesa-dir /path/to/mesa
```

On this ARM compiler, double/extended/quadruple kinds are 8/16/16. The explicit
portability patch removes only duplicate generic declarations for extended
precision; default double precision and all physics formulae remain unchanged.
The compiler-kind probe is a prerequisite to that patch. Earlier failed builds
are preserved. All five upstream byte comparisons still fail; the first
upstream EOS physical-state fields agree to 1.44e-14 scaled error. The full
upstream derivative suite is not silently marked passing.

The direct reference's 153-state gas-only first-law defect is at most 1.01e-7,
gas derivative error 1.99e-6 and pressure-density inverse error 2.24e-14.
Checking the gas alone prevents radiation from masking the test. These are
component checks, not an independent accuracy measurement of the entire EOS.

Evidence, compiler logs, exact profiles, failed attempts and source snapshots
are indexed in [Report 8](../Reports/Report_8.md) and the
[outer-boundary archive](../benchmarks/results/2026-09-28-outer-boundary/archive_manifest.json).
