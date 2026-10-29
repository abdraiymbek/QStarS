"""Self-selected mdot and sonic point from a frozen CMFGEN radiation field.

Gate-4 component reduction, not a certified wind and not an evolution
boundary. Plan, registered before any result: benchmarks/
cmfgen_selected_wind_plan.json. One hydro step freezes a^2=P_gas/rho,
Gamma=r^2 g_rad/(GM) and kappa_R of a converged CMFGEN state on ln tau_R
(primary) or ln r (cross-check and existence proof) and solves

    (v^2-a^2) dlnv/dlnr = N = 2a^2 - da^2/dlnr - GM(1-Gamma)/r,
    dtau/dr = -kappa_R rho,   mdot = 4 pi r^2 rho v,

with mdot unknown, a regular critical crossing, tau(r_o)=tau_o and one base
condition. The material first law and the transfer/population response are
NOT solved here; they enter only by re-converging CMFGEN on the new flow
(outer fixed point). Tables are never extrapolated or clipped. No solver
failure is labelled physical, and no classification certifies a wind.
"""
from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from .constants import G, C, KB, MU
from .cmfgen_backend import read_rvtj, read_spectral, vector_section, _numbers, sha256
from .spherical_wind import regular_sonic_slope, WindBranchError
from .atmosphere_reference import AtmosphereDomainExit
from .outer_interface import InterfaceState, moving_flux

PLAN = 'benchmarks/cmfgen_selected_wind_plan.json'
CLASSIFICATIONS = ('reduction_invalid', 'mach_one_mlt_exit', 'hydrostatic_no_wind',
                   'multiple_admissible_branches', 'numerical_failure', 'photon_tiring_exit',
                   'stationary_transonic_wind')
M_H, M_HE = 1.00794, 4.002602  # amu; composition bookkeeping only
# Plan acceptance_thresholds enforced by classify; mirrored from PLAN (tested), not overridable.
THRESHOLDS = dict(hydro_step_base_condition_log_residual=1e-8, critical_regularity_relative=1e-6)
SONIC_STALL = 1e-3  # |v^2/a^2-1| at a DOP853 stall below which the stall is the sonic fold (D->0, N!=0)


# ---------------------------------------------------------------- native files

def _fortran(field):
    """One fixed-width Fortran field; '*' overflow is NaN, never a neighbour."""
    t = field.strip()
    if not t: raise ValueError('Empty fixed-width HYDRO field')
    if '*' in t: return np.nan
    x = _numbers(t)
    if x is None or len(x) != 1: raise ValueError(f'Malformed HYDRO field {field!r}')
    return x[0]


def read_hydro_terms(path):
    """Native HYDRO_TERMS_V5 table in cgs, rows in ascending radius.

    Fields are cut at the native FORMAT widths (hydro_terms_v5.f), so an
    overflowed field is NaN and counted. g_RAD is kappa_F L(ND)/(4 pi r^2 c)
    with the INNER comoving luminosity (line 214), not the local force.
    Gamma is F9.2; recover g from g_RAD-g_TOT (ES14.4) instead.
    """
    lines = Path(path).read_text().splitlines()
    head = [i for i, s in enumerate(lines) if '% Error' in s and 'g_RAD' in s]
    if len(head) != 1: raise ValueError('Expected exactly one HYDRO header')
    turbulent = 'dTPdR/ROH' in lines[head[0]]
    spec = [(None, 1), ('radius', 15 if turbulent else 12), ('velocity', 17), ('percent_error', 9),
            ('vdvdr', 14), ('dpdr_on_rho', 14)] + ([('dpturbdr_on_rho', 14)] if turbulent else []) + [
            ('g_tot', 14), ('g_rad', 14), ('g_elec', 14), ('gamma', 9), ('depth', 7), (None, 2),
            ('tau', 10), ('sound_speed', 10), ('error_over_grad_percent', 12)]
    width = sum(w for _, w in spec); rows = []
    for s in lines[head[0]+1:]:
        if not s.strip(): break
        if len(s.rstrip()) > width: raise ValueError('HYDRO row wider than the native format')
        s = s.ljust(width); i = 0; row = []
        for name, w in spec:
            if name is not None: row.append(_fortran(s[i:i+w]))
            i += w
        rows.append(row)
    names = [n for n, _ in spec if n is not None]
    a = np.array(rows)
    if a.ndim != 2 or len(a) < 2: raise ValueError('Truncated HYDRO table')
    if not np.array_equal(a[:, names.index('depth')], np.arange(1, len(a)+1)):
        raise ValueError('HYDRO depth column is not 1..ND')
    q = {n: a[::-1, k] for k, n in enumerate(names) if n != 'depth'}
    q['radius'] = q['radius']*1e10; q['velocity'] = q['velocity']*1e5; q['sound_speed'] = q['sound_speed']*1e5
    if not np.all(np.diff(q['radius']) > 0): raise ValueError('HYDRO radius must decrease in native order')
    q['gravity'] = q['g_rad']-q['g_tot']
    text = '\n'.join(lines[head[0]:])
    def footer(label, scale=1.):
        m = re.search(re.escape(label)+r'\s+(\S+)', text)
        if not m or '*' in m[1]: return np.nan
        x = _numbers(m[1].split('(')[0])
        return x[0]*scale if x else np.nan
    q.update(turbulent_pressure=turbulent, overflowed_fields=int(np.isnan(a).sum()),
             photospheric_radius=footer('Photospheric radius is:', 1e10),
             stars_mass_Msun=footer('Stars mass is:'),
             gamma_es_photosphere=footer('Edd. fac. (e.s.) at photosphere:'),
             gamma_es_inner=footer('Edd. fac. (e.s.) at inner bound:'),
             native_grad_luminosity='inner_comoving_L(ND)')
    return q


def read_gencool_channels(path, n):
    """All GENCOOL blocks as (rows, ND) arrays in ascending radius.

    Channel rows are per-process (per-level) cooling minus heating, erg/cm^3/s.
    Their sum is checked against the native net; per-channel nets
    underestimate gross exchange, so 0.5 sum|rows| is a LOWER bound on it.
    The 'Adiabatic cooling rate' rows (V, dTdR) are flow work under INC_AD,
    not radiative exchange: kept in the closure sum, dropped from
    radiative_gross_exchange (plan amendment pending user approval).
    """
    chunks = []; cur = None; label = None
    for s in Path(path).read_text().splitlines():
        t = s.strip()
        if not t: label = None; continue
        row = _numbers(t)
        if row is None:
            label = t
            if label == 'Depth': cur = {'Depth': []}; chunks.append(cur)
            elif cur is None: raise ValueError('GENCOOL block before a Depth header')
            elif label in cur: raise ValueError(f'Duplicate GENCOOL block {label!r}')
            else: cur[label] = []
            continue
        if cur is None or label is None: raise ValueError('Unlabelled GENCOOL numbers')
        cur[label].append(row)
    if not chunks or any(set(c) != set(chunks[0]) for c in chunks):
        raise ValueError('Inconsistent GENCOOL chunk structure')
    blocks = {}
    for key in chunks[0]:
        parts = []
        for c in chunks:
            width = len(c['Depth'][0])
            if len(c[key]) < 1 or any(len(r) != width for r in c[key]): raise ValueError(f'Malformed GENCOOL block {key!r}')
            if key != 'Depth' and len(c[key]) != len(chunks[0][key]): raise ValueError(f'GENCOOL row count varies for {key!r}')
            parts.append(np.array(c[key]))
        blocks[key] = np.hstack(parts)
    depth = blocks.pop('Depth')[0]
    if not np.array_equal(depth, np.arange(1, n+1)): raise ValueError('GENCOOL depths are not 1..ND')
    blocks = {k: v[:, ::-1] for k, v in blocks.items()}
    if not all(np.all(np.isfinite(v)) for v in blocks.values()): raise ValueError('Nonfinite GENCOOL value')
    adiabatic = lambda k: k.startswith('Adiabatic cooling rate')
    channel = lambda k: ('Cooling' in k or adiabatic(k)) and not k.startswith(('Net Cooling', 'Ratio of'))
    rows = np.vstack([v for k, v in blocks.items() if channel(k)])
    net = blocks['Net Cooling Rate [ergs/cm**3/sec]'][0]
    total = abs(rows).sum(axis=0)
    radiative = sum((abs(v).sum(axis=0) for k, v in blocks.items() if channel(k) and not adiabatic(k)), np.zeros(n))
    return dict(blocks=blocks, channels=[k for k in blocks if channel(k)], net_cooling=net,
                adiabatic_channels=[k for k in blocks if channel(k) and adiabatic(k)],
                channel_sum=rows.sum(axis=0), channel_abs_sum=total, gross_exchange_lower_bound=.5*total,
                radiative_gross_exchange=.5*radiative,
                channel_closure_relative=np.divide(abs(rows.sum(axis=0)-net), total, out=np.full(n, np.inf), where=total > 0))


def read_vadat(path):
    """'value [KEY]' entries; Fortran exponents become floats, other text is kept."""
    out = {}
    for s in Path(path).read_text().splitlines():
        m = re.match(r'\s*(\S+)\s+\[([A-Za-z0-9_]+)\]', s)
        if m:
            x = _numbers(m[1]); out[m[2]] = x[0] if x and len(x) == 1 else m[1]
    return out


def rosseland_depth(radius, chi):
    """Independent tau_R on ascending radius: trapezoid in ln r of chi r.

    The outer tail assumes chi ~ r^-p beyond the grid (tau=chi r/(p-1)),
    p from the two outermost points; p<=1 has no finite tail and raises.
    """
    r = np.asarray(radius, float); c = np.asarray(chi, float); x = np.log(r)
    if np.any(c <= 0) or not np.all(np.diff(r) > 0): raise ValueError('Positive opacity on ascending radius required')
    p = -np.log(c[-1]/c[-2])/(x[-1]-x[-2])
    if not p > 1: raise ValueError('Outer opacity tail has no finite optical depth (p<=1)')
    y = c*r; seg = .5*(y[1:]+y[:-1])*np.diff(x)
    return c[-1]*r[-1]/(p-1)+np.r_[np.cumsum(seg[::-1])[::-1], 0.]


def load_wind_inputs(run_dir, *, spectral=True, hydro=True, cooling=True, levels=None, hashes=True):
    """One CMFGEN pass for the hydro step; arrays in ascending radius.

    The spectral loop precedes the final population correction, so spectral
    and RVTJ material states are compared and returned, never assumed equal.
    a^2=(n_nuclei+n_e)kT/rho_clump. levels (read_atomic_levels per ion)
    adds the resolved atomic energy of cmfgen_backend.hhe_material.
    """
    d = Path(run_dir); path = d/'RVTJ'; lines = path.read_text().splitlines()
    s = read_rvtj(path); n = len(s['radius'])
    atoms = vector_section(lines, 'Atom Density', n)[::-1]
    pressure = (atoms+s['electron_density'])*KB*s['temperature']
    nh, nhe = s['hydrogen_density'], s['helium_density']
    q = dict(run_dir=str(d), state=s, atom_density=atoms, sigma=vector_section(lines, 'dlnV/dlnr-1', n)[::-1],
             gas_pressure=pressure, sound_speed_squared=pressure/s['density_clump'],
             tau_rosseland=rosseland_depth(s['radius'], s['chi_rosseland']*s['filling_factor']),
             hydrogen_mass_fraction=float(np.median(nh*M_H/(nh*M_H+nhe*M_HE))),
             translational_energy_density=1.5*pressure, files={})
    def record(name): q['files'][name] = dict(sha256=sha256(d/name) if hashes else None, mtime=(d/name).stat().st_mtime, bytes=(d/name).stat().st_size)
    record('RVTJ')
    if (d/'VADAT').exists(): q['vadat'] = read_vadat(d/'VADAT'); record('VADAT')
    if spectral:
        m = read_spectral(d/'QSTARS_SPECTRAL.bin').moments(coherent_scattering=True); record('QSTARS_SPECTRAL.bin')
        if len(m['radius']) != n: raise ValueError('Spectral and RVTJ depth counts differ')
        q['moments'] = m
        q['spectral_sync'] = {k: float(np.max(abs(m[a]/s[b]-1))) for k, a, b in
                              [('radius', 'radius', 'radius'), ('temperature', 'temperature', 'temperature'), ('density', 'density', 'density')]}
        q['spectral_sync']['velocity'] = float(np.max(abs(m['velocity']-s['velocity'])/s['velocity']))
    if hydro:
        h = read_hydro_terms(d/'HYDRO'); record('HYDRO')
        if len(h['radius']) != n: raise ValueError('HYDRO and RVTJ depth counts differ')
        q['hydro'] = h; q['hydro_sync_radius'] = float(np.max(abs(h['radius']/s['radius']-1)))
    if cooling: q['cooling'] = read_gencool_channels(d/'GENCOOL', n); record('GENCOOL')
    if levels is not None:
        from .cmfgen_backend import read_hhe_populations, hhe_material
        mat = hhe_material(s, read_hhe_populations(path), levels)
        q['material'] = mat
        q['thermal_energy_density'] = mat['translational_energy_density_clump']+mat['atomic_energy_density_clump']
        q['specific_internal_energy'] = mat['gas_internal_energy']
    return q


# ----------------------------------------------------------- frozen sources

@dataclass(frozen=True)
class FrozenSources:
    """a^2, da^2/dx, Gamma=r^2 g_rad/(GM), kappa_R on x=ln r or x=ln tau_R.

    function(x)->(a2, da2dx, gamma, kappa), vectorized. Outside the domain
    a ValueError is raised: no extrapolation and no clipping.
    """
    coordinate: str
    function: object
    domain: tuple
    nodes: np.ndarray = None
    provenance: dict = None

    def __post_init__(self):
        if self.coordinate not in ('radius', 'tau'): raise ValueError('Frozen coordinate must be radius or tau')
        if not np.all(np.isfinite(self.domain)) or not self.domain[0] < self.domain[1]: raise ValueError('Invalid frozen domain')

    @classmethod
    def tabulated(cls, coordinate, x, sound_speed_squared, gamma, kappa, provenance=None):
        """Monotone PCHIP of ln a^2, Gamma and ln kappa on the given nodes."""
        x = np.asarray(x, float); o = np.argsort(x); x = x[o]
        y = np.column_stack([np.log(np.asarray(sound_speed_squared, float)[o]), np.asarray(gamma, float)[o],
                             np.log(np.asarray(kappa, float)[o])])
        if len(x) < 3 or not np.all(np.diff(x) > 0) or not np.all(np.isfinite(y)):
            raise ValueError('Frozen table needs >=3 strictly ordered finite nodes with a^2,kappa>0')
        p = PchipInterpolator(x, y, axis=0, extrapolate=False); dp = p.derivative()
        def f(z):
            v, d = p(z), dp(z); a2 = np.exp(v[..., 0])
            return a2, a2*d[..., 0], v[..., 1], np.exp(v[..., 2])
        return cls(coordinate, f, (float(x[0]), float(x[-1])), x, provenance)

    def __call__(self, x):
        z = np.asarray(x, float)
        if not np.all(np.isfinite(z)) or np.any(z < self.domain[0]) or np.any(z > self.domain[1]):
            raise ValueError('Outside frozen source table; no extrapolation or clipping')
        a2, da2, gam, k = self.function(z)
        if not np.all(np.isfinite([a2, da2, gam, k])) or np.any(a2 <= 0) or np.any(k <= 0):
            raise ValueError('Nonphysical frozen source value')
        return a2, da2, gam, k


def frozen_sources(inputs, mass, coordinate, *, force='spectral', exclude_inner=1):
    """Tables from one converged state; depth ND (diffusion boundary) excluded.

    force='spectral' is the local comoving force of the export. 'native_local'
    rescales native g_RAD by L(r)/L(ND) and exists only for the input audit.
    """
    s = inputs['state']; r = s['radius']; k = slice(exclude_inner, None)
    if force == 'spectral': g = inputs['moments']['radiative_acceleration']
    elif force == 'native_local':
        L = inputs['moments']['luminosity_comoving']; g = inputs['hydro']['g_rad']*L/L[0]
    else: raise ValueError('Unknown force source')
    x = np.log(r if coordinate == 'radius' else inputs['tau_rosseland'])
    return FrozenSources.tabulated(coordinate, x[k], inputs['sound_speed_squared'][k], (g*r*r/(G*mass))[k],
                                   s['kappa_rosseland'][k], provenance=dict(run_dir=inputs.get('run_dir'), force=force,
                                   exclude_inner=exclude_inner, mass=mass))


# --------------------------------------------------------------- residuals

def momentum_residual(radius, velocity, gas_pressure, density, radiative_acceleration, mass, sigma=None):
    """v dv/dr + rho^-1 dP/dr + GM/r^2 - g_rad on a tabulated state.

    v dv/dr=v^2(sigma+1)/r from the CMFGEN velocity law when sigma is given,
    else PCHIP in ln r; dlnP/dlnr by PCHIP. 'relative' uses the sum of all
    term magnitudes; 'hydro_percent' repeats the native 0-200 normalization.
    """
    r, v, P, rho, g = map(lambda a: np.asarray(a, float), (radius, velocity, gas_pressure, density, radiative_acceleration))
    x = np.log(r)
    vdv = v*v*(np.asarray(sigma)+1)/r if sigma is not None else PchipInterpolator(x, v*v/2).derivative()(x)/r
    pterm = P/rho*PchipInterpolator(x, np.log(P)).derivative()(x)/r
    grav = G*mass/r**2; res = vdv+pterm+grav-g
    return dict(residual=res, relative=res/(abs(vdv)+abs(pterm)+grav+abs(g)),
                hydro_percent=200*res/(abs(vdv)+abs(pterm)+abs(g-grav)),
                vdvdr=vdv, dpdr_on_rho=pterm, gravity=grav, gamma=g/grav)


def critical_coefficients(sources, mass, x, lnv, lntau=None, mdot=None):
    """N and D=v^2-a^2 (velocity^2) at x=ln r; tau mode needs ln tau and mdot."""
    r = np.exp(x); v = np.exp(lnv)
    if sources.coordinate == 'radius':
        a2, da2, gam, k = sources(x); dlt = None; tangent = da2
    else:
        if lntau is None or mdot is None: raise ValueError('Tau mode needs ln tau and mdot')
        a2, da2, gam, k = sources(lntau); dlt = -k*mdot/(4*np.pi*r*v*np.exp(lntau)); tangent = da2*dlt
    grav = G*mass*(1-gam)/r
    return dict(numerator=2*a2-tangent-grav, denominator=v*v-a2, sound_speed_squared=a2,
                temperature_term=tangent, effective_gravity_term=grav, gamma=gam, kappa=k, dlntau_dlnr=dlt)


def _integrate(rhs, span, y0, rtol, events=(), stage_rejection=False, sonic_gap=None):
    """DOP853; a table exit or nonfinite state is a rejected local trial.

    stage_rejection (tau mode, where the table coordinate ln tau is a
    dependent variable): a table exit or nonfinite value at a trial RK stage
    returns NaN, so DOP853 rejects that step (NaN error norm, step x0.2)
    instead of aborting. An accepted step therefore never used an
    out-of-table value; a solution that truly leaves the table ends with
    too small a step and is reported as a table exit. sonic_gap(x,y) =
    |v^2/a^2-1|: a stall below SONIC_STALL is the D->0, N!=0 fold of a
    non-regular sonic approach (the only singularity of the step's ODE),
    a rejected branch; any other stall stays a numerical failure.
    """
    exits = []
    def reject(x):
        exits.append(x)  # bounded: a solver that kept retrying NaN stages would otherwise never return
        if len(exits) > 10000: raise WindBranchError('Trial left the frozen table: 10000 rejected RK stages')
    def wrapped(x, y):
        try:
            with np.errstate(all='ignore'): d = np.asarray(rhs(x, y), float)
        except ValueError:
            if not stage_rejection: raise
            reject(x); return np.full(len(y), np.nan)
        if not np.all(np.isfinite(d)):
            if not stage_rejection: raise WindBranchError('Nonfinite hydro-step derivative')
            reject(x)
        return d
    try:
        s = solve_ivp(wrapped, span, y0, method='DOP853', rtol=rtol, atol=rtol*1e-2, dense_output=True,
                      events=list(events))
    except ValueError as exc:
        raise WindBranchError('Trial left the frozen table or became nonfinite: '+str(exc))
    if s.status < 0:
        try: gap = float(sonic_gap(s.t[-1], s.y[:, -1])) if sonic_gap is not None else np.inf
        except ValueError: gap = np.inf
        if gap < SONIC_STALL:
            raise WindBranchError(f'Branch reaches the sonic surface non-regularly: integration stalls at '
                                  f'|v^2/a^2-1|={gap:.2e} (D->0 with N!=0)')
        if exits: raise WindBranchError(f'Trial left the frozen table or became nonfinite ({len(exits)} rejected stages): '+s.message)
        raise RuntimeError('Hydro-step integration failed; numerical, not physical: '+s.message)
    if not np.all(np.isfinite(s.y)): raise WindBranchError('Nonfinite accepted hydro-step state')
    return s


def _log_shell(lm, sols):
    """ln of int 4 pi r^2 rho dr = mdot int r/v dlnr on 8 dense points per solver step."""
    xs = []; lvs = []
    for sol in sols:
        t = np.sort(sol.t); x = np.unique(np.concatenate([np.linspace(a, b, 9) for a, b in zip(t[:-1], t[1:])]))
        xs.append(x); lvs.append(sol.sol(x)[0])
    x = np.concatenate(xs); lv = np.concatenate(lvs)
    return lm+_log_trapezoid(x-lv, x)


def _sonic(lna2):
    """Terminal v=a crossing event and the |v^2/a^2-1| gap for stall diagnosis."""
    cross = lambda x, y: y[0]-.5*lna2(x, y); cross.terminal = True
    return [cross], lambda x, y: abs(np.expm1(2*y[0]-lna2(x, y)))


def _check_sonic(sol, side):
    if sol.t_events[0].size: raise WindBranchError(f'{side} branch re-encounters the sonic surface')


# -------------------------------------------------------- radius-mode step

def _radius_branch(src, mass, xc, slope, xb, xo, rtol, offset):
    """Two-sided regular branch; tau integral I=int kappa/(4 pi r v) dlnr.

    The inner integral is carried as ln I: exp(-lambda) structure overflows.
    """
    a2c = src(xc)[0]; lvc = .5*np.log(a2c); eps = offset/max(1., abs(slope))
    if not xb < xc-eps or not xc+eps < xo: raise WindBranchError('Critical point too close to the domain edge')
    def parts(x, lv):
        a2, da2, gam, k = src(x)
        return (2*a2-da2-G*mass*(1-gam)/np.exp(x))/(np.exp(2*lv)-a2), np.log(k/(4*np.pi))-x-lv
    def rhs_in(x, y):
        s, lf = parts(x, y[0]); return [s, -np.exp(lf-y[1])]
    def rhs_out(x, y):
        s, lf = parts(x, y[0]); return [s, np.exp(lf)]
    lf = parts_lnf(src, xc, lvc)
    sonic, gap = _sonic(lambda x, y: np.log(src(x)[0]))
    inner = _integrate(rhs_in, [xc-eps, xb], [lvc-slope*eps, lf+np.log(eps)], rtol, sonic, sonic_gap=gap)
    outer = _integrate(rhs_out, [xc+eps, xo], [lvc+slope*eps, np.exp(lf)*eps], rtol, sonic, sonic_gap=gap)
    for sol, side in [(inner, 'inner'), (outer, 'outer')]: _check_sonic(sol, side)
    return inner, outer, eps


def _radius_step(src, mass, *, rb, ro, tau_b, tau_o, rho_b, rtol, offset, step, samples, radii):
    xb, xo = np.log(rb), np.log(ro)
    if not src.domain[0] <= xb < xo <= src.domain[1]: raise ValueError('Base/outer radius outside the frozen radius table')
    nodes = src.nodes[(src.nodes > xb) & (src.nodes < xo)] if src.nodes is not None else np.array([])
    edges = np.unique(np.r_[xb, nodes, xo])
    if src.nodes is None: edges = np.linspace(xb, xo, 257)
    grid = np.unique(np.concatenate([np.linspace(a, b, samples+1) for a, b in zip(edges[:-1], edges[1:])]))
    def Nfun(x):
        a2, da2, gam, k = src(x); return 2*a2-da2-G*mass*(1-gam)/np.exp(x)
    N = Nfun(grid); scaled = N*np.exp(grid)/(G*mass)
    base = dict(coordinate='radius', mass=mass, base_radius=rb, outer_radius=ro, base_tau=tau_b, outer_tau=tau_o,
                base_density=rho_b, sampled_points=len(grid), max_scaled_numerator=float(scaled.max()),
                all_branches_enumerated=True, physical_impossibility_claimed=False, wind_certified=False)
    if scaled[0] >= 0:
        return dict(base, status='base_not_accelerating', reason='N(R_b)>=0: no accelerating subsonic base in this reduction')
    cross = np.flatnonzero(np.sign(scaled[:-1]) != np.sign(scaled[1:]))
    if cross.size == 0:
        if scaled.max() > -1e-6:
            return dict(base, status='marginal_tangency', reason='N touches zero without a resolved sign change')
        return dict(base, status='no_critical_point', **_hydrostatic_bound(src, mass, xb, xo, tau_b, tau_o, rho_b, rtol))
    candidates = []; admissible = []; found = []
    for i in cross:
        xc = brentq(lambda z: float(Nfun(z)), grid[i], grid[i+1], xtol=1e-15, rtol=1e-15)
        if any(abs(xc-z) < 1e-12 for z in found): continue
        found.append(xc)
        a2c = src(xc)[0]; row = dict(sonic_radius=float(np.exp(xc)))
        def coeff(y, xc=xc, a2c=a2c):
            a2, da2, gam, k = src(xc+y[0])
            return np.array([2*a2-da2-G*mass*(1-gam)/np.exp(xc+y[0]), a2c*np.exp(2*y[1])-a2])/a2c
        try:
            slope, roots = regular_sonic_slope(coeff, 0., step)
            inner, outer, eps = _radius_branch(src, mass, xc, slope, xb, xo, rtol, offset)
            row.update(status='admissible', slope=slope, slope_roots=roots)
            admissible.append((xc, slope, roots, inner, outer, eps))
        except WindBranchError as exc:
            row.update(status='rejected_crossing', reason=str(exc))
        except RuntimeError as exc:  # unresolved candidate: enumeration continues, no admissibility claim
            row.update(status='candidate_integration_failure', reason=str(exc))
        candidates.append(row)
    base['critical_points'] = candidates
    if any(c['status'] == 'candidate_integration_failure' for c in candidates):
        return dict(base, status='numerical_failure', reason='critical-point candidate(s) unresolved by integration; '
                    'admissibility undetermined, numerical, not physical')
    if not admissible:
        return dict(base, status='no_admissible_crossing', reason='N changes sign but no admissible X-type crossing')
    if len(admissible) > 1:
        return dict(base, status='multiple_admissible', reason='More than one admissible crossing; continuation required')
    xc, slope, roots, inner, outer, eps = admissible[0]
    a2c = src(xc)[0]; lvc = .5*np.log(a2c)
    lI_in = inner.y[1, -1]; I_out = outer.y[1, -1]
    if tau_b is not None:
        if not tau_b > tau_o: raise ValueError('Base optical depth must exceed the outer value')
        lm = np.log(tau_b-tau_o)-np.logaddexp(lI_in, np.log(I_out))
    else:
        lm = np.log(4*np.pi*rb*rb*rho_b)+inner.y[0, -1]
    x = np.unique(np.r_[np.log(radii) if radii is not None else np.linspace(xb, xo, 401), xc])
    if x[0] < xb-1e-12 or x[-1] > xo+1e-12: raise ValueError('Requested radii outside the hydro domain')
    lv = np.empty_like(x); tau = np.empty_like(x)
    im = x < xc-eps; om = x > xc+eps; near = ~(im | om)
    y = inner.sol(np.clip(x[im], xb, xc-eps)); lv[im] = y[0]
    tau[im] = tau_o+np.exp(lm+np.logaddexp(np.log(I_out), y[1]))  # tau_o+mdot*(I_out(xo)+I_in(x))
    y = outer.sol(np.clip(x[om], xc+eps, xo)); lv[om] = y[0]
    tau[om] = tau_o+np.exp(lm)*np.maximum(I_out-y[1], 0.)        # tau_o+mdot*(I_out(xo)-I_out(x))
    lv[near] = lvc+slope*(x[near]-xc); tau[near] = tau_o+np.exp(lm)*I_out
    a2, da2, gam, k = src(x)
    with np.errstate(all='ignore'): sl = (2*a2-da2-G*mass*(1-gam)/np.exp(x))/(np.exp(2*lv)-a2)
    sl[near] = slope
    sonic_tau = tau_o+np.exp(lm)*I_out
    return dict(base, status='transonic', ln_mdot=float(lm), mdot=float(np.exp(lm)), sonic_radius=float(np.exp(xc)),
                sonic_velocity=float(np.sqrt(a2c)), sonic_tau=float(sonic_tau), sonic_slope=float(slope),
                ln_sonic_tau_excess=float(lm+np.log(I_out)),  # ln(tau_c-tau_o), resolved when << tau_o
                slope_roots=roots, dlntau_dlnr_sonic=float(-np.exp(lm+parts_lnf(src, xc, lvc))/sonic_tau),
                regularity_relative=float(abs(Nfun(xc))*np.exp(xc)/(G*mass)),
                base_log_residual=0., shell_mass=float(np.exp(_log_shell(lm, [inner, outer]))),
                profile=_profile(x, lv, lm, tau, sl, a2))


def parts_lnf(src, x, lv):
    """ln[kappa/(4 pi r v)]: dtau/dlnr = -mdot exp(this)."""
    return np.log(src(x)[3]/(4*np.pi))-x-lv


def _log_trapezoid(lf, x):
    """ln of the trapezoid integral of exp(lf) dx without overflow."""
    o = np.argsort(x); lf, x = np.asarray(lf)[o], np.asarray(x)[o]
    terms = np.logaddexp(lf[1:], lf[:-1])+np.log(.5*np.diff(x)+1e-300)
    return float(np.logaddexp.reduce(terms[np.diff(x) > 0])) if np.any(np.diff(x) > 0) else -np.inf


def _profile(x, lv, lm, tau, slope, a2):
    r = np.exp(x); v = np.exp(lv)
    return dict(radius=r, velocity=v, log_density=lm-np.log(4*np.pi*r*r)-lv, tau=tau,
                sigma=slope-1, mach=v/np.sqrt(a2))


def _hydrostatic_bound(src, mass, xb, xo, tau_b, tau_o, rho_b, rtol):
    """Upper bound when N<0 on the whole domain (plan rule 5).

    dlnP/dlnr <= -GM(1-Gamma)/(r a^2) for every subsonic accelerating flow,
    so P(r_o) <= P_hs(r_o) with the same base normalization (tau or density)
    and mdot < 4 pi r_o^2 P_hs(r_o)/a(r_o). Logs avoid exp(-lambda) underflow.
    """
    def rhs(x, y):
        a2, da2, gam, k = src(x); r = np.exp(x)
        return [-G*mass*(1-gam)/(r*a2), k*r/a2*np.exp(y[0])]
    # y0 = ln(P/P_b); y1 = int kappa (P/P_b)/a^2 dr, bounded because N<0 keeps P/P_b of order one or less
    s = _integrate(rhs, [xb, xo], [0., 0.], rtol)
    if not s.success: raise RuntimeError('Hydrostatic bound integration failed; numerical, not physical')
    lPo_rel = s.y[0, -1]; a2o = src(xo)[0]; a2b = src(xb)[0]
    lPb = np.log(tau_b-tau_o)-np.log(s.y[1, -1]) if tau_b is not None else np.log(rho_b*a2b)
    lb = np.log(4*np.pi)+2*xo+lPb+lPo_rel-.5*np.log(a2o)
    if not np.isfinite(lb): raise RuntimeError('Hydrostatic bound is not finite; numerical, not physical')
    return dict(ln_mdot_bound=float(lb), mdot_bound=float(np.exp(lb)), hydrostatic_base_pressure=float(np.exp(lPb)),
                hydrostatic_outer_pressure=float(np.exp(lPb+lPo_rel)),
                bound_assumption='regular stationary flow subsonic and accelerating on the whole domain (N<0, D<0)')


# ----------------------------------------------------------- tau-mode step

def _ln_log1p_exp(d):
    """ln ln(1+e^d) without overflow or underflow (d -> -inf gives d)."""
    if d > 0: return float(np.log(np.logaddexp(0., d)))
    r = np.exp(d)
    return float(d+(np.log(np.log1p(r)/r) if r > 0 else 0.))


def _tau_branch(src, mass, u, *, rb, ro, tau_b, tau_o, rho_b, rtol, offset, step, floor, final=False):
    """Residual of the frozen-tau problem at u=(ln mdot, ln(tau_c-tau_o)).

    r_c is closed form: N=0 at v=a_c gives
    r_c=[GM(1-Gamma_c)-(da^2/dlntau)_c kappa_c mdot/(4 pi a_c tau_c)]/(2 a_c^2).
    mdot enters only through logs, lq=ln[kappa mdot/(4 pi r v)] (exp(-lambda)
    underflows). The outer branch carries w=ln(tau/tau_o)/z_c, z_c=ln(tau_c/tau_o)
    (ln z_c in closed form), so tau_c-tau_o << tau_o is resolved; the inner
    branch carries ln tau. Newton steps on the event-radius residuals
    (ln r(tau_o)-ln r_o, ln r(tau_b)-ln R_b: smooth and bounded); scale
    converts them to first-order (1+z_c) w(r_o) and ln tau(R_b)-ln tau_b,
    the forms reported and thresholded, on which Newton stops.
    """
    lm, ls = map(float, u); lto = np.log(tau_o); l4p = np.log(4*np.pi)
    lzc = _ln_log1p_exp(ls-lto); zc = np.exp(lzc); ltc = lto+zc  # zc may underflow; ltc is then tau_o to rounding
    def S(lt):
        try: return src(lt)
        except WindBranchError: raise
        except ValueError as exc: raise WindBranchError('Trial left the frozen table: '+str(exc))
    def parts(x, lv, lt):
        """N, D and lq; dlnTau/dlnr=-exp(lq-lt)."""
        a2_, da2_, gam_, k_ = S(lt); lq = np.log(k_)+lm-l4p-x-lv
        return 2*a2_+da2_*np.exp(lq-lt)-G*mass*(1-gam_)/np.exp(x), np.exp(2*lv)-a2_, lq
    a2, da2, gam, k = S(ltc)
    rc = (G*mass*(1-gam)-da2*np.exp(np.log(k)+lm-l4p-.5*np.log(a2)-ltc))/(2*a2)
    if not rb < rc < ro: raise WindBranchError('Critical radius outside the base-to-outer domain')
    xc = np.log(rc); lvc = .5*np.log(a2)
    def coeff(y):
        N, D, _ = parts(xc+y[0], lvc+y[1], ltc+y[2]); return np.array([N, D])/a2
    lqc = parts(xc, lvc, ltc)[2]; g = -np.exp(lqc-ltc); gw = -np.exp(lqc-ltc-lzc)  # dlnTau/dlnr, dw/dlnr at r_c
    slope, roots = regular_sonic_slope(coeff, g, step)
    eps = offset/max(1., abs(slope), abs(g), abs(gw))
    lt_out = lambda w: lto+zc*w
    rhs_out = lambda x, y: (lambda N, D, lq: [N/D, -np.exp(lq-lt_out(y[1])-lzc)])(*parts(x, y[0], lt_out(y[1])))
    rhs_in = lambda x, y: (lambda N, D, lq: [N/D, -np.exp(lq-y[1])])(*parts(x, y[0], y[1]))
    sonic_out, gap_out = _sonic(lambda x, y: np.log(S(lt_out(y[1]))[0]))
    sonic_in, gap_in = _sonic(lambda x, y: np.log(S(y[1])[0]))
    # Sources depend on tau only, so the outward search may pass r_o (up to 100 r_o):
    # ln r(tau_o) - ln r_o is smooth, unlike switching between two residual forms.
    ev_out = lambda x, y: y[1]; ev_out.terminal = True; ev_out.direction = -1
    out = _integrate(rhs_out, [xc+eps, np.log(ro) if final else np.log(100*ro)], [lvc+eps*slope, 1+eps*gw], rtol,
                     sonic_out if final else sonic_out+[ev_out], stage_rejection=True, sonic_gap=gap_out)
    _check_sonic(out, 'outer')
    if final: f1, s1 = (1+zc)*out.y[1, -1], 1.
    elif out.t_events[1].size: f1, s1 = out.t[-1]-np.log(ro), (1+zc)*np.exp(parts(out.t[-1], out.y[0, -1], lto)[2]-lto-lzc)
    else: raise WindBranchError('Outer optical depth not reached within 100 r_o')
    y0 = [lvc-eps*slope, ltc-eps*g]
    if tau_b is not None and not final:
        ev_b = lambda x, y: y[1]-np.log(tau_b); ev_b.terminal = True
        inn = _integrate(rhs_in, [xc-eps, np.log(rb*floor)], y0, rtol, sonic_in+[ev_b], stage_rejection=True, sonic_gap=gap_in)
        _check_sonic(inn, 'inner')
        if not inn.t_events[1].size: raise WindBranchError('Base optical depth not reached above the floor radius')
        f2, s2 = inn.t[-1]-np.log(rb), np.exp(parts(inn.t[-1], inn.y[0, -1], inn.y[1, -1])[2]-inn.y[1, -1])
    else:
        inn = _integrate(rhs_in, [xc-eps, np.log(rb)], y0, rtol, sonic_in, stage_rejection=True, sonic_gap=gap_in)
        _check_sonic(inn, 'inner')
        f2, s2 = (inn.y[1, -1]-np.log(tau_b)) if tau_b is not None else (lm-np.log(4*np.pi*rb*rb)-inn.y[0, -1]-np.log(rho_b)), 1.
    value = np.array([f1, f2])
    if not final: return value, np.array([s1, s2])
    return value, dict(xc=xc, lvc=lvc, ltc=ltc, lto=lto, zc=zc, a2c=a2, slope=slope, roots=roots, g=g, eps=eps,
                       inner=inn, outer=out, parts=parts, outer_excess=float(out.y[1, -1]))


def _newton(fun, u0, tolerance, max_nfev, trace, trust=1., h=1e-6):
    """Admissible damped Newton (spherical_wind.match_base pattern), 2 unknowns.

    fun(u) -> (residual, scale): steps and line search use residual; the
    stop is max|residual*scale| <= tolerance (the reported, thresholded
    form). Every trial, including rejected ones, is appended to trace.
    """
    def F(u, purpose):
        if len(trace) >= max_nfev: raise RuntimeError('Tau-mode evaluation budget exhausted; numerical, not physical')
        row = dict(log_coordinates=[float(z) for z in u], purpose=purpose)
        try:
            value, sc = fun(u)
            if not np.all(np.isfinite(value)) or not np.all(np.isfinite(sc)): raise WindBranchError('Nonfinite tau-mode residual')
            row.update(status='admissible', residual=value.tolist(), scaled_residual=(value*sc).tolist())
        except WindBranchError as exc:
            row.update(status='rejected_local_trial', reason=str(exc)); value = sc = None
        trace.append(row); return value, sc
    u = np.asarray(u0, float); f, sc = F(u, 'initial')
    if f is None: raise RuntimeError('Initial tau-mode trial inadmissible; seed from the radius mode')
    while max(abs(f*sc)) > tolerance:
        cols = []
        for i in range(2):
            hh = h
            for _ in range(8):
                e = np.eye(2)[i]*hh; p = F(u+e, 'jacobian_plus')[0]; m = F(u-e, 'jacobian_minus')[0]
                if p is not None and m is not None: cols.append((p-m)/(2*hh)); break
                hh *= .5
            else: raise RuntimeError('No admissible centered Jacobian stencil; numerical, not physical')
        J = np.column_stack(cols)
        stepv, _, rank, sv = np.linalg.lstsq(J, -f, rcond=1e-12)
        if rank < 2: raise RuntimeError('Tau-mode Jacobian lost rank; continuation required')
        stepv *= min(1., trust/max(abs(stepv)))
        old = np.linalg.norm(f)
        for b in range(20):
            t, tsc = F(u+2.**-b*stepv, 'line_search')
            if t is not None and np.linalg.norm(t) < old*(1-1e-4*2.**-b):
                u = u+2.**-b*stepv; f, sc = t, tsc; break
        else: raise RuntimeError('Tau-mode line search stalled; numerical, not physical')
    return u, f*sc, trace


def _tau_step(src, mass, *, rb, ro, tau_b, tau_o, rho_b, rtol, offset, step, seed, tolerance, max_nfev, floor, radii):
    if seed is None: raise ValueError('Tau mode needs an explicit seed (ln mdot, ln(tau_c-tau_o)), e.g. tau_seed(radius_step)')
    seed = np.asarray(seed, float)
    if seed.shape != (2,) or not np.all(np.isfinite(seed)): raise ValueError('Tau-mode seed must be finite (ln mdot, ln(tau_c-tau_o))')
    if not (src.domain[0] < np.log(tau_o) and (tau_b is None or np.log(tau_b) < src.domain[1])):
        raise ValueError('Outer/base optical depth must lie strictly inside the frozen tau table (margin node)')
    xb, xo = np.log(rb), np.log(ro)
    x = np.log(radii) if radii is not None else np.linspace(xb, xo, 401)
    if np.min(x) < xb-1e-12 or np.max(x) > xo+1e-12: raise ValueError('Requested radii outside the hydro domain')
    kw = dict(rb=rb, ro=ro, tau_b=tau_b, tau_o=tau_o, rho_b=rho_b, rtol=rtol, offset=offset, step=step, floor=floor)
    base = dict(coordinate='tau', mass=mass, base_radius=rb, outer_radius=ro, base_tau=tau_b, outer_tau=tau_o,
                base_density=rho_b, all_branches_enumerated=False, physical_impossibility_claimed=False, wind_certified=False)
    trace = []
    try:
        u, f, trace = _newton(lambda u: _tau_branch(src, mass, u, **kw), seed, tolerance, max_nfev, trace)
        value, b = _tau_branch(src, mass, u, final=True, **kw)
    except (RuntimeError, WindBranchError) as exc:
        return dict(base, status='numerical_failure', reason=str(exc), trials=trace)
    lm, ls = u; xc, eps, slope, ltc = b['xc'], b['eps'], b['slope'], b['ltc']
    x = np.unique(np.r_[x, xc])
    lv = np.empty_like(x); lt = np.empty_like(x)
    for mask, sol, out in [(x < xc-eps, b['inner'], False), (x > xc+eps, b['outer'], True)]:
        lo, hi = sorted(sol.t[[0, -1]]); y = sol.sol(np.clip(x[mask], lo, hi))  # clip <= 1e-12 (checked above)
        lv[mask], lt[mask] = y[0], (b['lto']+b['zc']*y[1] if out else y[1])
    near = abs(x-xc) <= eps; lv[near] = b['lvc']+slope*(x[near]-xc); lt[near] = ltc+b['g']*(x[near]-xc)
    N, D, _ = b['parts'](x, lv, lt); a2 = src(lt)[0]
    sl = np.where(near, slope, N/np.where(near, 1., D))
    return dict(base, status='transonic', ln_mdot=float(lm), mdot=float(np.exp(lm)), sonic_radius=float(np.exp(xc)),
                sonic_velocity=float(np.sqrt(b['a2c'])), sonic_tau=float(np.exp(ltc)), ln_sonic_tau_excess=float(ls),
                sonic_slope=float(slope), slope_roots=b['roots'], dlntau_dlnr_sonic=float(b['g']),
                regularity_relative=float(abs(b['parts'](xc, b['lvc'], ltc)[0])*np.exp(xc)/(G*mass)),
                outer_log_residual=float(b['zc']*b['outer_excess']), outer_excess_residual=b['outer_excess'],
                base_log_residual=float(value[1]), newton_residual=[float(z) for z in f],
                shell_mass=float(np.exp(_log_shell(lm, [b['inner'], b['outer']]))),
                profile=_profile(x, lv, lm, np.exp(lt), sl, a2), trials=trace, method='admissible_damped_newton')


def hydro_step(sources, *, mass, base_radius, outer_radius, outer_tau, base_tau=None, base_density=None,
               seed=None, rtol=1e-10, critical_offset=1e-5, derivative_step=1e-5, samples_per_interval=32,
               tolerance=1e-9, max_nfev=200, floor_ratio=.5, radii=None):
    """One frozen-source hydro step: mdot and the sonic point are outputs.

    Exactly one base condition: tau_R(R_b)=base_tau (primary) or
    rho(R_b)=base_density. Radius mode enumerates every sign change of N on
    the whole domain; tau mode solves (ln mdot, ln(tau_c-tau_o)) by admissible
    damped Newton from an explicit seed. Physical outcomes never raise; invalid
    inputs raise ValueError; numerical failures are returned as such.
    """
    if (base_tau is None) == (base_density is None): raise ValueError('Give exactly one base condition')
    vals = [mass, base_radius, outer_radius, outer_tau, base_tau if base_tau is not None else base_density]
    if not np.all(np.isfinite(vals)) or min(vals) <= 0 or not base_radius < outer_radius:
        raise ValueError('Finite positive mass, radii, optical depths and ordered domain required')
    if base_tau is not None and not base_tau > outer_tau: raise ValueError('Base optical depth must exceed the outer value')
    if radii is not None:
        rr = np.asarray(radii, float)
        if not np.all(np.isfinite(rr)) or rr.min() < base_radius*(1-1e-12) or rr.max() > outer_radius*(1+1e-12):
            raise ValueError('Requested radii outside the hydro domain')
    if not 0 < rtol < 1e-6 or not 0 < critical_offset < 1e-2 or samples_per_interval < 4:
        raise ValueError('Invalid hydro-step controls')
    kw = dict(rb=base_radius, ro=outer_radius, tau_b=base_tau, tau_o=outer_tau, rho_b=base_density,
              rtol=rtol, offset=critical_offset, step=derivative_step, radii=radii)
    if sources.coordinate == 'radius':
        try: return _radius_step(sources, mass, samples=samples_per_interval, **kw)
        except (RuntimeError, WindBranchError) as exc:
            return dict(coordinate='radius', status='numerical_failure', reason=str(exc), physical_impossibility_claimed=False)
    return _tau_step(sources, mass, seed=seed, tolerance=tolerance, max_nfev=max_nfev, floor=floor_ratio, **kw)


def tau_seed(radius_step):
    """Seed (ln mdot, ln(tau_c-tau_o)) for the tau mode from a transonic step.

    The excess is carried in logs by both modes, so tau_c-tau_o << tau_o and
    mdot below the float range keep a finite seed.
    """
    if radius_step.get('status') != 'transonic': raise ValueError('Seed needs a transonic radius-mode step')
    return np.array([radius_step['ln_mdot'], radius_step['ln_sonic_tau_excess']])


# ---------------------------------------------------------- screens/guard

def reduction_screens(inputs, step, *, mass, limits=None):
    """Plan screens at the step's own critical point (rule 9); never clipped.

    State quantities are read where the step froze a^2, Gamma and kappa:
    at ln tau_R = ln sonic_tau in tau mode, at ln r = ln sonic_radius in
    radius mode (before the fixed point the closure's r(tau_c) is not the
    hydro r_c). A point outside the state table raises.
    Thermal relaxation is a gross-exchange PROXY, not a bound: the
    relaxation rate T dLambda_net/dT can be well below gross exchange. It
    uses translational+atomic energy only when the resolved material energy
    is loaded (translational-only cannot pass) and, while the plan
    amendment is pending, must pass with both the registered all-channel
    estimate and the radiative-only one (adiabatic flow-work rows dropped).
    """
    lim = dict(thermal=1e-2, trapping=1e-2, relativistic=1e-2, shell=1e-2, gencool=1e-2)
    lim.update(limits or {})
    s = inputs['state']; mode = step.get('coordinate')
    if mode == 'tau': xg = np.log(inputs['tau_rosseland'])[::-1]; xc = np.log(step['sonic_tau']); o = slice(None, None, -1)
    elif mode == 'radius': xg = np.log(s['radius']); xc = np.log(step['sonic_radius']); o = slice(None)
    else: raise ValueError('Screens need a radius- or tau-mode step')
    if not np.all(np.diff(xg) > 0) or not xg[0] <= xc <= xg[-1]:
        raise ValueError(f'Critical point outside the state {mode} table; screens are never clipped or extrapolated')
    lin = lambda a: float(np.interp(xc, xg, np.asarray(a, float)[o]))
    def at(a):
        with np.errstate(divide='ignore'): return float(np.exp(lin(np.log(np.asarray(a, float)))))
    full = 'thermal_energy_density' in inputs
    energy = at(inputs['thermal_energy_density'] if full else inputs['translational_energy_density'])
    c = inputs['cooling']; closure = lin(c['channel_closure_relative'])
    vc, rc, sc = step['sonic_velocity'], step['sonic_radius'], step['sonic_slope']
    t_flow = rc/abs(2+sc)/vc
    with np.errstate(divide='ignore'):  # zero exchange -> inf -> fails
        ratio = lambda key: float(np.float64(energy)/at(c[key])/t_flow) if key in c else None
        q = dict(thermal_relaxation_over_flow=ratio('gross_exchange_lower_bound'),
                 thermal_relaxation_over_flow_radiative_only=ratio('radiative_gross_exchange'),
                 photon_trapping=vc*step['sonic_tau']/C, max_v_over_c=float(np.max(step['profile']['velocity'])/C),
                 shell_mass_fraction=step['shell_mass']/mass, gencool_channel_closure=closure, atomic_energy_included=full,
                 coordinate=mode, thermal_screen='gross_exchange_proxy')
    thermal = [q[k] for k in ('thermal_relaxation_over_flow', 'thermal_relaxation_over_flow_radiative_only')]
    passed = dict(thermal=full and all(t is not None and t <= lim['thermal'] for t in thermal),
                  trapping=q['photon_trapping'] <= lim['trapping'], relativistic=q['max_v_over_c'] < lim['relativistic'],
                  shell=q['shell_mass_fraction'] <= lim['shell'], gencool=closure <= lim['gencool'])
    return dict(values=q, passed=passed, all_passed=all(passed.values()), limits=lim)


class ConvectiveFluxExit(RuntimeError):
    """Subsonic MLT would carry flux that the radiative CMFGEN structure omits."""
    def __init__(self, diagnostics):
        super().__init__('Convective flux not negligible at a CMFGEN depth; radiative structure invalid there')
        self.diagnostics = diagnostics


def mach_one_guard(inputs, eos, mesa, *, mass, opacity='both', tau_min=2/3, alpha=2., model='Henyey',
                   flux_limit=.01, composition_tolerance=1e-4):
    """Existing Mach-one MLT exclusion at subsonic CMFGEN depths with tau_R>=tau_min.

    LTE FreeEOS (gas+radiation) and native MESA MLT at each CMFGEN (rho,T),
    g=GM/r^2, flux=L_comoving/(4 pi r^2), grad_rad as atmosphere_reference.
    Opacity 'cmfgen' (Rosseland mean), 'mesa' (LTE) or 'both' (either trips).
    Raises AtmosphereDomainExit at Mach>=1, then ConvectiveFluxExit above
    flux_limit. Neither is clipped or a physical impossibility claim.
    """
    if opacity not in ('cmfgen', 'mesa', 'both'): raise ValueError('Unknown guard opacity')
    if eos.hydrogen != mesa.composition.X: raise ValueError('Guard EOS and opacity composition differ')
    if abs(eos.hydrogen-inputs['hydrogen_mass_fraction']) > composition_tolerance:
        raise ValueError('Guard composition differs from the CMFGEN state')
    arad = eos.radiation_constant
    if abs(arad/mesa.radiation_constant-1) > 1e-12: raise ValueError('Guard radiation constants differ')
    s = inputs['state']
    idx = np.flatnonzero((s['velocity']**2 < inputs['sound_speed_squared']) & (inputs['tau_rosseland'] >= tau_min))
    if idx.size == 0: return dict(guarded_depths=0, max_convective_mach=None, max_convective_flux_fraction=None)
    r, rho, T = s['radius'][idx], s['density'][idx], s['temperature'][idx]
    q = eos.evaluate(rho, T); g = G*mass/r**2; F = inputs['moments']['luminosity_comoving'][idx]/(4*np.pi*r*r)
    kap = {}
    if opacity in ('cmfgen', 'both'): kap['cmfgen'] = s['kappa_rosseland'][idx]
    if opacity in ('mesa', 'both'): kap['mesa'] = np.asarray(mesa.opacity_external_eos(rho, T, q['free_e'], q['eta']))
    out = {}
    for name, k in kap.items():
        gradr = 3*k*q['p']*F/(4*arad*C*g*T**4)
        z = mesa.mlt(np.array([q['chiT'], q['chiRho'], q['cp'], g, rho, q['p'], T, k, gradr, q['grad_ad']]), model=model, alpha=alpha)
        mach = z[1]/q['cs']; frac = np.where(gradr > 0, 1-z[0]/gradr, 0.)
        out[name] = dict(radiative_gradient=gradr, gradient=z[0], convective_velocity=z[1], convective_mach=mach,
                         convective_flux_fraction=frac)
        if np.any(mach >= 1):
            i = int(np.argmax(mach))
            raise AtmosphereDomainExit(dict(depth_index=int(idx[i]), radius=float(r[i]), tau=float(inputs['tau_rosseland'][idx[i]]),
                density=float(rho[i]), temperature=float(T[i]), convective_velocity=float(z[1][i]),
                sound_speed=float(q['cs'][i]), convective_mach=float(mach[i]), opacity_source=name,
                model=model, alpha=alpha, physical_impossibility_claimed=False))
    worst = max(float(np.max(v['convective_flux_fraction'])) for v in out.values())
    if worst > flux_limit:
        raise ConvectiveFluxExit(dict(max_convective_flux_fraction=worst, limit=flux_limit,
                                      physical_impossibility_claimed=False))
    return dict(guarded_depths=int(idx.size), depth_index=idx, opacities=out, model=model, alpha=alpha,
                max_convective_mach=max(float(np.max(v['convective_mach'])) for v in out.values()),
                max_convective_flux_fraction=worst)


def guard_status(inputs, eos, mesa, **kw):
    """mach_one_guard as data: 'passed', 'mach_one', 'convective_flux' or 'invalid_input'."""
    try: return dict(status='passed', **mach_one_guard(inputs, eos, mesa, **kw))
    except AtmosphereDomainExit as exc: return dict(status='mach_one', diagnostics=exc.diagnostics)
    except ConvectiveFluxExit as exc: return dict(status='convective_flux', diagnostics=exc.diagnostics)
    except ValueError as exc: return dict(status='invalid_input', reason=str(exc))


def input_checks(inputs, *, sync_tolerance=1e-3):
    """Plan rule-1 checks on one loaded state; guard composition is checked by the guard."""
    sync = inputs.get('spectral_sync') or {}
    return dict(nonrelativistic=bool(np.max(abs(inputs['state']['velocity']))/C < .01),
                spectral_loaded='moments' in inputs,
                spectral_rvtj_sync=bool(sync) and max(sync.values()) <= sync_tolerance)


# ---------------------------------------------------------- classification

def _step_unconverged(step):
    """Plan tolerances a transonic step must meet (missing values fail).

    base_log_residual: plan hydro_step_base_condition_log_residual. The same
    1e-8 on the tau-mode outer ln tau residual is a stricter application
    recorded as pending user approval. regularity: critical_regularity_relative.
    """
    tol = THRESHOLDS['hydro_step_base_condition_log_residual']
    keys = ['base_log_residual']+(['outer_log_residual'] if step.get('coordinate') == 'tau' else [])
    bad = [k for k in keys if not abs(step.get(k, np.nan)) <= tol]
    if not step.get('regularity_relative', np.nan) <= THRESHOLDS['critical_regularity_relative']: bad.append('regularity_relative')
    return bad


def classify(tau_step, radius_step, *, mass, luminosity, base_radius, input_checks=None, guard=None, screens=None):
    """First matching rule of the plan's classification_order.

    Never a physical-impossibility claim; 'stationary_transonic_wind' is a
    component outcome, not a certified wind. A 'passed' guard counts only
    if it guarded >=1 depth with both opacities (plan default while the
    opacity choice is pending); a transonic step above the plan's residual
    or regularity tolerances is a nonconvergent root (rule 7).
    """
    def out(rule, name, reason, **kw):
        return dict(rule=rule, classification=name, reason=reason, physical_impossibility_claimed=False,
                    wind_certified=False, **kw)
    failed = sorted(k for k, v in (input_checks or {}).items() if not v)
    if failed: return out(1, 'reduction_invalid', 'input checks failed: '+', '.join(failed))
    if guard is None or guard['status'] not in ('passed', 'mach_one', 'convective_flux'):
        return out(1, 'reduction_invalid', 'Mach-one MLT guard not evaluated: '+str((guard or {}).get('reason', 'missing')))
    if guard['status'] == 'passed' and (not guard.get('guarded_depths', 0) > 0 or set(guard.get('opacities', ())) != {'cmfgen', 'mesa'}):
        return out(1, 'reduction_invalid', 'Mach-one MLT guard vacuous (no guarded depth) or not run with both cmfgen and mesa opacities')
    if guard['status'] == 'mach_one': return out(2, 'mach_one_mlt_exit', 'convective Mach >= 1', diagnostics=guard['diagnostics'])
    if guard['status'] == 'convective_flux': return out(3, 'reduction_invalid', 'MLT flux fraction above limit', diagnostics=guard['diagnostics'])
    rs = radius_step['status']
    if rs in ('base_not_accelerating', 'no_admissible_crossing', 'marginal_tangency'):
        return out(4, 'reduction_invalid', radius_step.get('reason', rs))
    if rs == 'no_critical_point':
        return out(5, 'hydrostatic_no_wind', 'N<0 on the whole radius-mode domain',
                   ln_mdot_bound=radius_step['ln_mdot_bound'], confirmation_required=True)
    if rs == 'multiple_admissible': return out(6, 'multiple_admissible_branches', radius_step['reason'])
    if rs != 'transonic' or tau_step['status'] != 'transonic':
        return out(7, 'numerical_failure', (radius_step if rs != 'transonic' else tau_step).get('reason', 'failed'))
    bad = {k: v for k, v in (('tau', _step_unconverged(tau_step)), ('radius', _step_unconverged(radius_step))) if v}
    if bad: return out(7, 'numerical_failure', 'hydro-step root above the plan tolerances', unconverged=bad)
    ltiring = tau_step['ln_mdot']+np.log(G*mass/(base_radius*luminosity)); tiring = float(np.exp(ltiring))
    if ltiring >= 0: return out(8, 'photon_tiring_exit', 'mdot GM/(R_b L) >= 1', tiring=tiring)
    if screens is None or not screens['all_passed']:
        return out(9, 'reduction_invalid', 'critical-point screens failed or not evaluated',
                   failed=[k for k, v in (screens or {}).get('passed', {}).items() if not v])
    return out(10, 'stationary_transonic_wind', 'component outcome; outer fixed point pending', tiring=tiring)


# ---------------------------------------------------- outer fixed point

def outer_update(history, *, momentum_passed, damping=.5, trust=np.log(10.), tolerance=3e-4, budget=12):
    """Plan policy for the next ln mdot of the outer fixed point.

    history rows: (x_k, xhat_k, ln r_s,k, ln rhat_s,k), all four required.
    momentum_passed: closure-k momentum-residual criterion of the latest row.
    'converged' only when |G_k|, |ln rhat_s-ln r_s| and momentum all pass;
    mdot and r_s at the fixed point with momentum failing is
    'mdot_converged_pending_checks', never 'converged'. First step damped;
    then secant unless it leaves the latest sign-change bracket; trust
    region ln 10 recorded when it truncates. Deterministic.
    """
    h = [tuple(map(float, row)) for row in history]
    if not h: raise ValueError('Outer history is empty')
    if any(len(r) != 4 for r in h): raise ValueError('Outer history rows must be (x, xhat, ln r_s, ln rhat_s)')
    if not isinstance(momentum_passed, (bool, np.bool_)): raise ValueError('momentum_passed must be a bool for the latest closure')
    x = np.array([r[0] for r in h]); g = np.array([r[1]-r[0] for r in h])
    checks = dict(mdot=bool(abs(g[-1]) <= tolerance), sonic_radius=bool(abs(h[-1][3]-h[-1][2]) <= tolerance),
                  momentum=bool(momentum_passed))
    failing = [k for k, v in checks.items() if not v]
    if not failing: return dict(status='converged', next=float(x[-1]), iterations=len(h), checks=checks)
    if failing == ['momentum']:
        return dict(status='mdot_converged_pending_checks', next=float(x[-1]), iterations=len(h), checks=checks,
                    failing=failing, reason='mdot and sonic radius at the fixed point; closure momentum criterion fails, not converged')
    if len(h) >= budget:
        return dict(status='numerical_failure', reason='outer fixed-point budget exhausted; not a physical statement',
                    iterations=len(h), checks=checks, failing=failing)
    opposite = np.flatnonzero(np.sign(g[:-1]) != np.sign(g[-1]))
    bracket = tuple(sorted((float(x[opposite[-1]]), float(x[-1])))) if opposite.size else None
    rule = 'damped'; stepv = damping*g[-1]
    if len(h) > 1 and g[-1] != g[-2]:
        xs = x[-1]-g[-1]*(x[-1]-x[-2])/(g[-1]-g[-2])
        if bracket is None or bracket[0] < xs < bracket[1]: rule = 'secant'; stepv = xs-x[-1]
    truncated = abs(stepv) > trust
    if truncated: stepv = np.sign(stepv)*trust
    return dict(status='continue', next=float(x[-1]+stepv), rule=rule, truncated=bool(truncated),
                bracket=bracket, iterations=len(h), checks=checks, failing=failing)


# --------------------------------------------------------- ledger and audit

def wind_ledger(inputs, *, mass, depths, specific_internal_energy=None, gauge='neutral_H_He_ground'):
    """Steady total energy flux through static base and outer faces.

    Uses outer_interface.moving_flux with export lab moments. With a static
    point-mass potential the steady total flux is radius independent, so the
    base-to-outer difference is the layer's energy defect. Translational-only
    gas energy (no levels loaded) is flagged and cannot pass.
    """
    s = inputs['state']; m = inputs['moments']; f = s['filling_factor']
    u = inputs.get('specific_internal_energy') if specific_internal_energy is None else specific_internal_energy
    atomic = u is not None
    if u is None: u = 1.5*inputs['sound_speed_squared']
    faces = []
    for i in depths:
        st = InterfaceState(radius=s['radius'][i], enclosed_mass=mass, density=s['density'][i], velocity=s['velocity'][i],
                            gas_pressure=inputs['gas_pressure'][i]*f[i], gas_internal_energy=u[i],
                            radiation_energy=m['E_lab'][i], radiation_flux=m['F_lab'][i], radiation_pressure_rr=m['P_lab'][i],
                            potential=-G*mass/s['radius'][i], energy_gauge=gauge)
        try: faces.append(moving_flux(st))
        except ValueError as exc: return dict(status='unrealizable_interface', reason=str(exc), depth=int(i))
    b, o = faces
    scale = max(sum(abs(v) for v in b['energy_channels'].values()), sum(abs(v) for v in o['energy_channels'].values()))
    i = depths[1]; r = s['radius'][i]; v = s['velocity'][i]
    return dict(status='evaluated', base=b, outer=o, energy_relative=(o['total_energy']-b['total_energy'])/scale,
                mass_relative=(o['mass']-b['mass'])/max(abs(o['mass']), abs(b['mass'])),
                outer_gas_bernoulli=v*v/2+u[i]+inputs['gas_pressure'][i]/s['density_clump'][i]-G*mass/r,
                outer_velocity_over_escape=v/np.sqrt(2*G*mass/r), atomic_energy_included=atomic,
                escaping_radiation_certified=False)


def grad_input_audit(inputs, *, mass, exclude_inner=1):
    """Export g_rad versus native HYDRO g_RAD, and the native density convention.

    Native g_RAD uses L(ND) (hydro_terms_v5.f:214); the definitional
    comparison is g_spec versus g_RAD L(r)/L(ND); since native g_RAD divides
    by POP_ATOM*mu*amu, agreement also tests the export density convention.
    Vsound uses the pass electron density, which the final correction may
    change, so its difference is informational only.
    """
    s = inputs['state']; m = inputs['moments']; h = inputs['hydro']; r = s['radius']
    keep = np.arange(len(r)) >= exclude_inner
    gs, gn, L = m['radiative_acceleration'], h['g_rad'], m['luminosity_comoving']
    raw = gs/gn-1; corr = gs/(gn*L/L[0])-1
    mu_rvtj = s['density_clump']/(inputs['atom_density']*MU)
    chi = m['radiative_force_density']*C/(m['F_comoving']*s['filling_factor'])/s['chi_flux']-1
    mom = momentum_residual(r, s['velocity'], inputs['gas_pressure'], s['density_clump'], gs, mass, sigma=inputs['sigma'])
    mx = lambda a: float(np.max(abs(a[keep])))
    return dict(depths=len(r), excluded_inner_depths=exclude_inner,
                raw_relative_max=mx(raw), raw_relative_inner=float(raw[0]),
                L_corrected_relative_max=mx(corr), L_corrected_relative_inner=float(corr[0]),
                luminosity_ratio_range=[float(np.min(L/L[0])), float(np.max(L/L[0]))],
                rvtj_mean_atomic_weight=dict(median=float(np.median(mu_rvtj)), spread=float(np.ptp(mu_rvtj))),
                flux_mean_opacity_relative_max=mx(chi),
                sound_speed_relative_max=mx(np.sqrt(inputs['sound_speed_squared'])/h['sound_speed']-1),
                momentum_relative_max=mx(mom['relative']), momentum_relative_median=float(np.median(abs(mom['relative'][keep]))),
                native_percent_error_abs_max=float(np.nanmax(abs(h['percent_error'][keep]))),
                native_percent_error_abs_median=float(np.nanmedian(abs(h['percent_error'][keep]))),
                gamma_spectral_range=[float(np.min(mom['gamma'][keep])), float(np.max(mom['gamma'][keep]))],
                tau_rosseland_independent_vs_native=dict(
                    inner=[float(inputs['tau_rosseland'][0]), float(h['tau'][0])],
                    outer=[float(inputs['tau_rosseland'][-1]), float(h['tau'][-1])]),
                spectral_sync=inputs.get('spectral_sync'), hydro_sync_radius=inputs.get('hydro_sync_radius'),
                arrays=dict(radius=r, raw=raw, L_corrected=corr, momentum_relative=mom['relative'],
                            native_percent_error=h['percent_error']))
