#!/usr/bin/env python3
"""Complete Jacobian action vs an exact-residual oracle directory (Report 21).

The direct part is built exactly as the Newton-Krylov driver builds it (cool_newton_krylov.Iterate:
native blocks, block-tridiagonal when INC_ADV/INC_AD are on, frozen-J n_e/T/H II columns corrected
for the native off-diagonal part, common POPION term); the radiation part is the analytic tangent
(QSTARS_JV_RAD.bin of the base pass, or of --rad). e = ||J.s - FD|| / ||FD|| over active rows, overall,
for the RE row alone and per species group, for each oracle step size.
Oracle layout (run_tfree.py): base/, col_{ne,hii,T}_{p,m}/, <dir>_{a,b}/ with QSTARS_STEP.bin = eps*s.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_fortran, read_jv_rad, read_trace_record  # noqa: E402

NAMES = ['T', 'ne', 'rand', 'mix']
GROUPS = {'H0': (1, 1), 'HI': (2, 31), 'HII': (32, 32), 'HeI': (33, 101), 'HeII': (102, 131),
          'HeIII': (132, 132), 'ne': (133, 133), 'T': (134, 134)}


def l2(a, b, m):
    return float(np.linalg.norm((a - b)[m]) / max(np.linalg.norm(b[m]), 1e-300))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir', type=Path)
    ap.add_argument('--rad', type=Path, default=None, help='product pass whose QSTARS_JV_RAD.bin to use (default base)')
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    d = a.dir.resolve()
    x, R0 = read_trace_record(d / 'base' / 'QSTARS_RESID_TRACE.bin')
    ck.NT, ck.ND = x.shape
    s_all, rad = read_jv_rad((a.rad or d / 'base') / 'QSTARS_JV_RAD.bin')
    cols = {n: (d / f'col_{n}_p', d / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(d / 'base', x, cols, 1e-4, True)
    m = np.broadcast_to(~np.all(R0 == 0.0, axis=1)[:, None], R0.shape).copy()
    mT = np.zeros_like(m); mT[133] = True
    res = dict(tridiagonal_direct=it.tridiagonal, RE_row_active=bool(np.any(R0[133])))
    print(f'direct part block-tridiagonal: {it.tridiagonal}; RE row active: {res["RE_row_active"]}')
    for k, n in enumerate(NAMES):
        sk = s_all[:, :, k]
        jv = it.Jdir.matvec(sk) + rad[:, :, k]
        fds = {}
        for lab in 'ab':
            p = d / f'{n}_{lab}'
            if (p / 'QSTARS_RESID_TRACE.bin').exists():
                st = read_fortran(p / 'QSTARS_STEP.bin', x.shape)
                nz = np.abs(sk) > 0
                eps = float(np.median(st[nz] / sk[nz]))
                _, Re = read_trace_record(p / 'QSTARS_RESID_TRACE.bin')
                fds[lab] = (Re - R0) / (-eps)
        if not fds:
            continue
        out = {}
        for lab, fd in fds.items():
            g = {}
            for gn, (a0, b0) in GROUPS.items():
                mg = np.zeros_like(m); mg[a0 - 1:b0] = m[a0 - 1:b0]
                if np.any(fd[mg] != 0):
                    g[gn] = l2(jv, fd, mg)
            out[lab] = dict(e=l2(jv, fd, m), e_RE_row=l2(jv, fd, mT), groups=g)
        if len(fds) == 2:
            out['oracle_self_consistency'] = l2(fds['a'], fds['b'], m)
            out['oracle_self_consistency_RE_row'] = l2(fds['a'], fds['b'], mT)
        res[n] = out
        print(f'{n:5s} ' + '  '.join(f'{lab}: e={v["e"]:.2e} RE={v["e_RE_row"]:.2e}' for lab, v in out.items()
                                     if isinstance(v, dict))
              + f'  oracle {out.get("oracle_self_consistency", float("nan")):.2e} / RE '
                f'{out.get("oracle_self_consistency_RE_row", float("nan")):.2e}')
    (a.out or d / 'action_audit.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
