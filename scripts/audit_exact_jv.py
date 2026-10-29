#!/usr/bin/env python3
"""Validate the complete Report-20 Jacobian action against exact-residual oracles.

  J.(x*s) = J_dir.s + dSTEQ_rad(s)

J_dir: native fixed-J blocks (block-diagonal; QS_ALI_LIN, QS_LIN_FAC=0) with the
       n_e and T columns replaced by frozen-J central differences and every ion
       column (H II, He II levels, He III) completed by the common POPION term
       measured on the H II column (per unit absolute population).
rad:   analytic nonlocal formal-transfer tangent through the native dSTEQ/dJ.

Oracles: Report 18 fixed-T runs (results/jac_check, eps per queue.sh) and the
T/n_e/population-direction runs (results/jv/oracle_t; FIX_T=F was requested but the
restart's BAMATPNT kept T fixed, so the RE row is inactive there: T COLUMN validated). e = ||J.s - FD|| / ||FD|| over active
rows; the oracle's own consistency is measured between its two step sizes.
"""
import argparse, json, struct, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_ba_blocks, read_fortran, read_jv_rad, read_trace_record  # noqa: E402

NT, ND = 134, 84
ION_ROWS = [32] + list(range(102, 133))
GROUPS = {'H0': (1, 1), 'HI': (2, 31), 'HII': (32, 32), 'HeI': (33, 101),
          'HeII': (102, 131), 'HeIII': (132, 132), 'ne': (133, 133), 'T': (134, 134)}


def exact_direct(base, cols, h, x, popion=True):
    _, di, _ = read_ba_blocks(Path(base) / 'QSTARS_BA_TRACE.bin', NT, ND)
    di = di.copy()
    fd = {}
    for n in ('ne', 'hii', 'T'):
        _, rp = read_trace_record(Path(cols) / f'{n}_p' / 'QSTARS_RESID_TRACE.bin')
        _, rm = read_trace_record(Path(cols) / f'{n}_m' / 'QSTARS_RESID_TRACE.bin')
        fd[n] = (rp - rm) / (2 * h)
    miss = np.stack([fd['hii'][:, d] - di[d][:, 31] for d in range(ND)], axis=1) / x[31][None, :]
    for d in range(ND):
        di[d][:, 132] = fd['ne'][:, d]
        di[d][:, 133] = fd['T'][:, d]
        if popion:
            for r in ION_ROWS:
                di[d][:, r - 1] += miss[:, d] * x[r - 1, d]
    return di


def apply_bd(di, s):
    return np.einsum('kij,jk->ik', di, s)


def l2(a, b, m):
    return float(np.linalg.norm((a - b)[m]) / np.linalg.norm(b[m]))


def compare(name, jv, fds, R0, rows_exclude=()):
    m = np.broadcast_to(~np.all(R0 == 0.0, axis=1)[:, None], R0.shape).copy()
    for r in rows_exclude:
        m[r - 1, :] = False
    out = {}
    for lab, fd in fds.items():
        out[lab] = dict(e=l2(jv, fd, m),
                        groups={g: l2(jv, fd, _grp(m, a, b)) for g, (a, b) in GROUPS.items()
                                if np.any(fd[_grp(m, a, b)] != 0)})
    labs = list(fds)
    if len(labs) == 2:
        out['oracle_self_consistency'] = l2(fds[labs[0]], fds[labs[1]], m)
    msg = ', '.join(f'{k}: e={v["e"]:.3e}' for k, v in out.items() if isinstance(v, dict))
    print(f'{name:10s} {msg}  oracle self-consistency {out.get("oracle_self_consistency", float("nan")):.3e}')
    return out


def _grp(m, a, b):
    g = np.zeros_like(m)
    g[a - 1:b] = m[a - 1:b]
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', type=Path, default=ROOT / 'results/jv/b21baseT')
    ap.add_argument('--cols', type=Path, default=ROOT / 'results/jv/frozenj_cols')
    ap.add_argument('--h', type=float, default=1e-4)
    ap.add_argument('--out', type=Path, default=ROOT / 'results/jv/exact_jv_audit.json')
    a = ap.parse_args()
    x, Rb = read_trace_record(a.base / 'QSTARS_RESID_TRACE.bin')
    s_all, rad = read_jv_rad(a.base / 'QSTARS_JV_RAD.bin')
    res = {}
    for popion in (False, True):
        di = exact_direct(a.base, a.cols, a.h, x, popion)
        tag = 'with_popion' if popion else 'ne_T_columns_only'
        res[tag] = {}
        print(f'--- direct part: {tag}')
        # T/n_e/population-direction oracle (dirs 0-3; RE row inactive, see docstring)
        o = ROOT / 'results/jv/oracle_t'
        _, R0 = read_trace_record(o / 'R0' / 'QSTARS_RESID_TRACE.bin')
        for k, n in enumerate(['T', 'ne', 'rand', 'mix']):
            fds = {}
            for lev, eps in (('a', None), ('b', None)):
                f = o / f'{n}_{lev}' / 'QSTARS_RESID_TRACE.bin'
                if f.exists():
                    eps = float((o / f'{n}_{lev}' / 'QSTARS_PERTURB_APPLIED').read_text().split()[1])
                    s_file = read_fortran(o / f's_{n}.bin', (NT, ND))
                    if not np.array_equal(s_file, s_all[:, :, k]):
                        raise SystemExit(f'direction {n} differs from the product direction')
                    _, Re = read_trace_record(f)
                    fds[lev] = (Re - R0) / (-eps)
            if fds:
                jv = apply_bd(di, s_all[:, :, k]) + rad[:, :, k]
                res[tag][n] = compare(n, jv, fds, R0)
        # Report-18 fixed-T oracle (dirs 4-7); T row excluded (identity there).
        j = ROOT / 'results/jac_check'
        _, R0f = read_trace_record(j / 'R0' / 'QSTARS_RESID_TRACE.bin')
        eps = {}
        for line in (j / 'queue.sh').read_text().splitlines():
            p = line.replace('SPECS="', '').strip('"').split()
            if len(p) == 3 and p[1].startswith('s_'):
                eps[p[0]] = float(p[2])
        for k, n in enumerate(['fac0', 'fac05', 'fac1', 'mom']):
            fds = {}
            for lev in ('e2', 'e3'):
                _, Re = read_trace_record(j / f'{n}_{lev}' / 'QSTARS_RESID_TRACE.bin')
                fds[lev] = (Re - R0f) / (-eps[f'{n}_{lev}'])
            jv = apply_bd(di, s_all[:, :, 4 + k]) + rad[:, :, 4 + k]
            res[tag][n] = compare(n, jv, fds, R0f, rows_exclude=(134,))
    a.out.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
