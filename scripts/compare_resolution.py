#!/usr/bin/env python3
"""Spatial / frequency resolution comparison of converged cool-closure states (Report 21; Report 23).

State mode (positional coarse, fine base passes): T, n_e, neutral H fraction and He I/II/III fractions at
the shared radii. The 167-radius grid nests the 84-radius grid (84-depth k = 167-depth 2k-1, native
order outer first) and the 333-radius grid the 167-radius one; a doubled-frequency state lives on the same
radii. Threshold: benchmarks/cmfgen_validation_plan.json profile_refinement_relative (1 %), applied to T,
n_e and the dominant-ion fractions; trace fractions are compared in log only where both exceed 1e-6.

Profile mode (--profiles COARSE FINE: the profiles.npz of the two audited verify passes, or their audit
directories): radiation, force, luminosity and heating at the shared radii, with the definitions and
thresholds of benchmarks/cmfgen_cool_frequency_plan.json refinement_profile_criteria (declared before any
result; b = the finer run). Metrics whose definition is pending the user's approval are listed with the
result, so a PASS that relies on them is never reported as unconditional. The plan's per_state_prerequisites
bind: each audit directory's validation.json checks and cool_validation.json additional_checks must all PASS
(per_state_audit), otherwise the pair is 'NOT COMPARABLE' and never a pass; a bare profiles.npz is unaudited.
pass_all needs comparable --profiles (state mode alone never sets it true).
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_trace_record  # noqa: E402

PLAN = ROOT / 'benchmarks/cmfgen_cool_frequency_plan.json'
RELATIVE = ('E_comoving', 'F_comoving', 'P_comoving', 'radiative_acceleration', 'kappa_flux')


def profiles(x):
    H = x[1:31].sum(0) + x[31]
    He = x[32:132].sum(0)
    return dict(T=x[133] * 1e4, ne=x[132], fHI=x[1:31].sum(0) / H, fHeI=x[32:101].sum(0) / He,
                fHeII=x[101:131].sum(0) / He, fHeIII=x[131] / He)


def one_record(path):
    # every reader takes record 0: a trace with a record inherited from its restart is refused (Report 23)
    n = Path(path).stat().st_size
    nt, nd = np.frombuffer(Path(path).read_bytes()[:12], '<i4')[1:]
    if n != 12 + 16 * int(nt) * int(nd):
        raise SystemExit(f'{path}: not a single-record trace ({n} bytes)')
    return read_trace_record(path)


def compare(a, b, nested):
    xa, _ = one_record(Path(a) / 'QSTARS_RESID_TRACE.bin')
    xb, _ = one_record(Path(b) / 'QSTARS_RESID_TRACE.bin')
    if nested:
        xb = xb[:, ::2]
    if xa.shape != xb.shape:
        raise SystemExit(f'shapes differ {xa.shape} {xb.shape}')
    pa, pb = profiles(xa), profiles(xb)
    out = {}
    for k in pa:
        va, vb = pa[k], pb[k]
        if k in ('T', 'ne'):
            rel = np.abs(va - vb) / np.abs(vb)
        else:
            dom = (va > 0.01) & (vb > 0.01)
            rel = np.where(dom, np.abs(va - vb), 0.0)          # absolute difference of a fraction > 1 %
            tr = (va > 1e-6) & (vb > 1e-6) & ~dom
            out[k + '_trace_max_dlog10'] = float(np.max(np.abs(np.log10(va[tr] / vb[tr])))) if tr.any() else 0.0
        out[k + '_max'] = float(np.max(rel))
        out[k + '_argmax_depth'] = int(np.argmax(rel)) + 1
    return out


def load_profiles(p):
    p = Path(p)
    z = np.load(p / 'profiles.npz' if p.is_dir() else p)
    return {k: z[k] for k in z.files}


def per_state_audit(p):
    """Per-state prerequisite (frequency plan per_state_prerequisites): every check of the audit directory's
    validation.json (checks) and cool_validation.json (additional_checks) is PASS. A file path (bare profiles.npz)
    or an audit missing either file is not audited; never comparable."""
    p = Path(p)
    st = dict(path=str(p), checks={}, missing=[])
    if not p.is_dir():
        st.update(status='unaudited', missing=['validation.json', 'cool_validation.json'])
        return st
    for fn, key in (('validation.json', 'checks'), ('cool_validation.json', 'additional_checks')):
        if not (p / fn).exists():
            st['missing'].append(fn)
            continue
        for c in json.loads((p / fn).read_text()).get(key) or []:
            st['checks'][f"{fn.split('.')[0]}:{c.get('name')}"] = c.get('status')
    if len(st['missing']) == 2:
        st['status'] = 'unaudited'
    else:
        st['status'] = 'PASS' if not st['missing'] and st['checks'] and all(v == 'PASS' for v in st['checks'].values()) else 'FAIL'
    return st


def shared_radii(ra, rb, rtol=1e-9):
    """Indices (ia, ib) of the coarse radii in the fine grid (exact nesting or identical grids required)."""
    ib = np.searchsorted(rb, ra)
    ib = np.clip(ib, 0, len(rb) - 1)
    lo = np.clip(ib - 1, 0, len(rb) - 1)
    ib = np.where(np.abs(rb[lo] / ra - 1) < np.abs(rb[ib] / ra - 1), lo, ib)
    if not np.all(np.abs(rb[ib] / ra - 1) <= rtol):
        raise SystemExit('the coarse radii are not a subset of the fine radii (relative 1e-9)')
    return np.arange(len(ra)), ib


def compare_profiles(pa, pb, criteria):
    """Profile metrics of the refinement criteria (definitions in the plan); b is the finer run."""
    oa, ob = np.argsort(pa['radius']), np.argsort(pb['radius'])     # inner radius first
    pa = {k: v[oa] if np.ndim(v) == 1 and len(v) == len(oa) else v for k, v in pa.items()}
    pb = {k: v[ob] if np.ndim(v) == 1 and len(v) == len(ob) else v for k, v in pb.items()}
    ia, ib = shared_radii(pa['radius'], pb['radius'])
    thr = criteria['thresholds']
    out = dict(shared_radii=int(len(ia)), metrics={})

    def put(name, value, where=None, extra=None):
        m = dict(value=float(value), threshold=thr[name], status='PASS' if value <= thr[name] else 'FAIL',
                 definition_status=criteria['status'].get(name, 'existing thresholds'))
        if where is not None:
            m['argmax_radius_cm'] = float(pa['radius'][ia][where])
        m.update(extra or {})
        out['metrics'][name] = m
    for k in RELATIVE:
        a, b = pa[k][ia], pb[k][ib]
        fin = np.isfinite(a) & np.isfinite(b)
        rel = np.where(fin, np.abs(a - b) / np.maximum(np.abs(b), 1e-300), np.inf)
        put(k, float(np.max(rel)), int(np.argmax(rel)), dict(nonfinite_points=int(np.sum(~fin))))
    La, Lb = pa['luminosity_lab'], pb['luminosity_lab']
    Ls = max(abs(Lb[0]), abs(Lb[-1]))
    dL = np.abs(La[ia] - Lb[ib]) / Ls
    put('luminosity_lab_over_Lscale', float(np.max(dL)), int(np.argmax(dL)),
        dict(Lscale=float(Ls), emergent_relative=float(abs(La[-1] - Lb[-1]) / abs(Lb[-1]))))

    def cumulative(p):
        r, h = p['radius'], p['heating_comoving']
        f = r * r * h
        return 4 * np.pi * np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])
    dQ = np.abs(cumulative(pa)[ia] - cumulative(pb)[ib]) / Ls
    hs = np.max(np.abs(pb['heating_comoving']))
    put('cumulative_heating_over_Lscale', float(np.max(dQ)), int(np.argmax(dQ)),
        dict(pointwise_heating_over_max_reported=float(np.max(np.abs(pa['heating_comoving'][ia] - pb['heating_comoving'][ib])) / hs)))
    out['pass'] = all(m['status'] == 'PASS' for m in out['metrics'].values())
    out['pending_user_approval'] = sorted(k for k, m in out['metrics'].items()
                                          if m['definition_status'] == 'pending_user_approval')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('coarse', type=Path, nargs='?', help='converged coarser base pass (state mode)')
    ap.add_argument('fine', type=Path, nargs='?', help='converged finer (333-radius or doubled-frequency) base pass')
    ap.add_argument('--nested', action='store_true', help='fine is the bisected (nested) grid')
    ap.add_argument('--profiles', type=Path, nargs=2, metavar=('COARSE', 'FINE'),
                    help='profiles.npz (or audit directories) of the two audited verify passes')
    ap.add_argument('--plan', type=Path, default=PLAN, help='plan holding refinement_profile_criteria')
    ap.add_argument('--label', default='spatial')
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if (a.coarse is None) != (a.fine is None) or (a.coarse is None and a.profiles is None):
        ap.error('give COARSE FINE base passes, --profiles COARSE FINE, or both')
    res = dict(label=a.label)
    if a.coarse is not None:
        plan = json.loads((ROOT / 'benchmarks/cmfgen_validation_plan.json').read_text())['thresholds']
        thr = plan['profile_refinement_relative']
        st = compare(a.coarse, a.fine, a.nested)
        gate_keys = ['T_max', 'ne_max', 'fHI_max', 'fHeI_max', 'fHeII_max', 'fHeIII_max']
        res.update(st)
        res['threshold'] = thr
        res['pass'] = bool(all(res[k] <= thr for k in gate_keys))
        res['coarse'] = str(a.coarse); res['fine'] = str(a.fine)
        if a.profiles is None:
            res['per_state_audits'] = 'not read in state mode; pass_all needs --profiles with both audit directories'
    if a.profiles is not None:
        crit = json.loads(a.plan.read_text())['refinement_profile_criteria']
        pr = compare_profiles(load_profiles(a.profiles[0]), load_profiles(a.profiles[1]), crit)
        pr.update(coarse=str(a.profiles[0]), fine=str(a.profiles[1]), plan=str(a.plan))
        pr['per_state_audits'] = [per_state_audit(p) for p in a.profiles]
        pr['comparable'] = all(st['status'] == 'PASS' for st in pr['per_state_audits'])
        pr['metrics_pass'] = pr['pass']
        if not pr['comparable']:
            pr['pass'] = False                  # never a pass without both per-state audits
        pr['status'] = 'NOT COMPARABLE' if not pr['comparable'] else 'PASS' if pr['pass'] else 'FAIL'
        res['profiles'] = pr
    res['pass_all'] = bool(a.profiles is not None and res['profiles']['comparable'] and res.get('pass', True)
                           and res['profiles']['pass'])
    a.out.write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
