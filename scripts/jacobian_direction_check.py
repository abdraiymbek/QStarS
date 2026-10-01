#!/usr/bin/env python3
"""Directional check of CMFGEN Newton linearizations against the exact residual.

For a Newton step from x0 with Jacobian J, the fractional correction s satisfies
J.dx = -R0 with dx = -x0*s (POPS <- POPS*(1-SCALE*s)).  If J is accurate along
dx, the exact residual at x0*(1-eps*s) is R0 + eps*J.dx = (1-eps)*R0 to first
order.  This script compares the measured residual change with that prediction.

Layout (built by the jac_check campaign):
  <root>/s_<J>.bin                 float64 NT*ND column-major fractional step
  <root>/R0/QSTARS_RESID_TRACE.bin unperturbed residual (first record used)
  <root>/<J>_e2/, <J>_e3/          perturbed runs with eps = 1e-2 / 1e-3 / max|s|
"""
import argparse, json, struct, sys
from pathlib import Path

import numpy as np

GROUPS = {'H0': (1, 1), 'HI': (2, 31), 'HII': (32, 32), 'HeI': (33, 101),
          'HeII': (102, 131), 'HeIII': (132, 132), 'ne': (133, 133), 'T': (134, 134)}
LEVELS = {'e2': 1e-2, 'e3': 1e-3}


def first_record(path):
    data = path.read_bytes()
    cnt, nt, nd = struct.unpack_from('<3i', data, 0)
    n = nt * nd
    if len(data) < 12 + 16 * n:
        raise SystemExit(f'{path}: truncated first record')
    pops = np.frombuffer(data, '<f8', n, 12).reshape((nt, nd), order='F')
    steq = np.frombuffer(data, '<f8', n, 12 + 8 * n).reshape((nt, nd), order='F')
    return cnt, pops.copy(), steq.copy()


def norm(a):
    return float(np.sqrt(np.sum(a * a)))


def compare(R0, Re, eps, mask):
    dR = (Re - R0)[mask]
    pred = (-eps * R0)[mask]
    err = dR - pred
    out = dict(eps=eps,
               e_glob=norm(err) / norm(pred),
               cos=float(np.dot(dR, pred) / (norm(dR) * norm(pred))) if norm(dR) > 0 else float('nan'),
               merit=norm(Re[mask]) / norm(R0[mask]), merit_pred=1.0 - eps)
    return out


def per_depth(R0, Re, eps, mask):
    nd = R0.shape[1]
    e = np.full(nd, np.nan)
    for d in range(nd):
        m = mask[:, d]
        p = -eps * R0[m, d]
        if norm(p) > 0:
            e[d] = norm((Re[m, d] - R0[m, d]) - p) / norm(p)
    return e


def per_group(R0, Re, eps, mask):
    out = {}
    for g, (a, b) in GROUPS.items():
        m = np.zeros_like(mask)
        m[a - 1:b, :] = mask[a - 1:b, :]
        p = (-eps * R0)[m]
        if norm(p) > 0:
            out[g] = norm((Re - R0)[m] - p) / norm(p)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root', type=Path)
    ap.add_argument('--directions', nargs='+', default=['fac1', 'fac05', 'fac0', 'mom'])
    args = ap.parse_args()
    root = args.root
    _, x0, R0 = first_record(root / 'R0' / 'QSTARS_RESID_TRACE.bin')
    nt, nd = R0.shape
    mask = np.ones_like(R0, dtype=bool)
    mask[np.all(R0 == 0.0, axis=1), :] = False
    print(f'NT={nt} ND={nd}  rows excluded (R0 identically 0): '
          f'{[int(i) + 1 for i in np.where(np.all(R0 == 0.0, axis=1))[0]]}')
    print(f'||R0||={norm(R0[mask]):.4e}  max|R0|={np.abs(R0).max():.4e}')
    summary = {}
    for J in args.directions:
        s = np.fromfile(root / f's_{J}.bin', '<f8')
        if s.size != nt * nd:
            raise SystemExit(f's_{J}.bin size {s.size} != {nt * nd}')
        s = s.reshape((nt, nd), order='F')
        smax = np.abs(s).max()
        res = {}
        dR = {}
        for lev, scale in LEVELS.items():
            d = root / f'{J}_{lev}'
            f = d / 'QSTARS_RESID_TRACE.bin'
            if not f.exists():
                print(f'{J}_{lev}: missing {f}')
                continue
            eps = scale / smax
            _, xp, Re = first_record(f)
            xe = x0 * (1.0 - eps * s)
            state_err = float(np.max(np.abs(xp - xe) / np.abs(xe)))
            r = compare(R0, Re, eps, mask)
            r['state_rel_err'] = state_err
            ed = per_depth(R0, Re, eps, mask)
            r['depth_median'] = float(np.nanmedian(ed))
            worst = np.argsort(np.nan_to_num(ed, nan=-1.0))[::-1][:5]
            r['depth_worst'] = [(int(k) + 1, float(ed[k])) for k in worst]
            r['groups'] = per_group(R0, Re, eps, mask)
            res[lev] = r
            dR[lev] = Re - R0
        if 'e2' in dR and 'e3' in dR:
            ratio = LEVELS['e2'] / LEVELS['e3']
            res['nonlin'] = norm((dR['e2'] - ratio * dR['e3'])[mask]) / norm(dR['e2'][mask])
        summary[J] = dict(max_abs_s=float(smax), **res)
        print(f'\n=== direction {J}: max|s|={smax:.4e}')
        for lev in LEVELS:
            if lev not in res:
                continue
            r = res[lev]
            print(f'  {lev}: eps={r["eps"]:.3e} state_err={r["state_rel_err"]:.1e} '
                  f'e_glob={r["e_glob"]:.3e} cos={r["cos"]:+.4f} '
                  f'merit={r["merit"]:.6f} (pred {r["merit_pred"]:.6f}) '
                  f'depth_median={r["depth_median"]:.3e}')
            print('     worst depths: ' + ', '.join(f'd{k}:{v:.2e}' for k, v in r['depth_worst']))
            print('     groups: ' + ', '.join(f'{g}:{v:.2e}' for g, v in r['groups'].items()))
        if 'nonlin' in res:
            print(f'  nonlinearity ||dR_e2-10 dR_e3||/||dR_e2|| = {res["nonlin"]:.3e}')
    out = root / 'jacobian_direction_check.json'
    out.write_text(json.dumps(summary, indent=1))
    print(f'\nwrote {out}')


if __name__ == '__main__':
    sys.exit(main())
