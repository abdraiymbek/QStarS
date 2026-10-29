#!/usr/bin/env python3
"""Complete Jacobian action with the radiative-equilibrium row ACTIVE (T released).

Oracle and products from results/jv/tfree/run_tfree.py (every pass released T by
removing the restart's BAMATPNT). J.(x*s) = J_dir.s + dSTEQ_rad(s): J_dir = native
fixed-J blocks with n_e and T columns replaced by frozen-J central differences and the
common POPION term (from the H II column) added to every ion column; dSTEQ_rad from
the analytic tangent (now including the RE row). e = ||J.s - FD|| / ||FD|| over active
rows, overall, per species group and for the RE row alone.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from qstars.nlte_newton import read_jv_rad, read_trace_record  # noqa: E402
import audit_exact_jv as A  # noqa: E402

NAMES = ['T', 'ne', 'rand', 'mix']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', type=Path, default=ROOT / 'results/jv/tfree')
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    d = a.dir
    x, R0 = read_trace_record(d / 'base' / 'QSTARS_RESID_TRACE.bin')
    if not np.any(R0[133]):
        raise SystemExit('RE row residual is identically zero: T was not released')
    s_all, rad = read_jv_rad(d / 'base' / 'QSTARS_JV_RAD.bin')
    cols = d / 'colsdir'
    cols.mkdir(exist_ok=True)
    for n in ('ne', 'hii', 'T'):
        for sg in 'pm':
            link = cols / f'{n}_{sg}'
            if not link.exists():
                link.symlink_to((d / f'col_{n}_{sg}').resolve())
    res = {'RE_row_residual_max': float(np.max(np.abs(R0[133])))}
    for popion in (False, True):
        di = A.exact_direct(d / 'base', cols, 1e-4, x, popion)
        tag = 'with_popion' if popion else 'ne_T_columns_only'
        res[tag] = {}
        print(f'--- direct part: {tag}')
        m = np.broadcast_to(~np.all(R0 == 0.0, axis=1)[:, None], R0.shape).copy()
        mT = np.zeros_like(m); mT[133] = True
        for k, n in enumerate(NAMES):
            jv = A.apply_bd(di, s_all[:, :, k]) + rad[:, :, k]
            jdir = A.apply_bd(di, s_all[:, :, k])
            fds = {}
            for lab in 'ab':
                p = d / f'{n}_{lab}'
                if (p / 'QSTARS_RESID_TRACE.bin').exists():
                    st, sk = read_step(p), s_all[:, :, k]
                    nz = np.abs(sk) > 0
                    eps = float(np.median(st[nz] / sk[nz]))
                    if not np.allclose(st, eps * sk, rtol=1e-12, atol=0):
                        raise SystemExit(f'{p}: step is not eps * product direction')
                    _, Re = read_trace_record(p / 'QSTARS_RESID_TRACE.bin')
                    fds[lab] = (Re - R0) / (-eps)
            if not fds:
                continue
            out = {}
            for lab, fd in fds.items():
                out[lab] = dict(e=A.l2(jv, fd, m), e_RE_row=A.l2(jv, fd, mT),
                                e_RE_row_direct_only=A.l2(jdir, fd, mT),
                                groups={g: A.l2(jv, fd, A._grp(m, a_, b_)) for g, (a_, b_) in A.GROUPS.items()
                                        if np.any(fd[A._grp(m, a_, b_)] != 0)})
            if len(fds) == 2:
                out['oracle_self_consistency'] = A.l2(fds['a'], fds['b'], m)
                out['oracle_self_consistency_RE_row'] = A.l2(fds['a'], fds['b'], mT)
            res[tag][n] = out
            print(f'{n:5s} ' + '  '.join(f'{lab}: e={v["e"]:.2e} RE={v["e_RE_row"]:.2e} (direct-only RE {v["e_RE_row_direct_only"]:.2e})'
                                         for lab, v in out.items() if isinstance(v, dict))
                  + f'  oracle {out.get("oracle_self_consistency", float("nan")):.2e} / RE {out.get("oracle_self_consistency_RE_row", float("nan")):.2e}')
    (a.out or d / 'tfree_jv_audit.json').write_text(json.dumps(res, indent=1))


def read_step(p):
    from qstars.nlte_newton import read_fortran
    return read_fortran(p / 'QSTARS_STEP.bin', (134, 84))


if __name__ == '__main__':
    main()
