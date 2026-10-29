#!/usr/bin/env python3
"""Locate the missing dependencies of the native fixed-J (direct) Jacobian.

The direct Jacobian of the CMFGEN rate equations is block-diagonal in depth
(verified from the dumped blocks), so perturbing ONE variable at ALL depths in a
frozen-J replay pass gives, by central differences, that variable's exact column
at every depth. Compared with the dumped analytic column, the difference is the
dependence the native linearization omits.

Hypothesis tested: the omitted terms enter only through the electron density and
the total ion density POPION (level dissolution / occupation probabilities), so
for any ion level i the omitted column per unit ABSOLUTE population is the same.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_ba_blocks, read_trace_record  # noqa: E402

NT, ND = 134, 84


def column_fd(dirp, dirm, h):
    _, rp = read_trace_record(Path(dirp) / 'QSTARS_RESID_TRACE.bin')
    _, rm = read_trace_record(Path(dirm) / 'QSTARS_RESID_TRACE.bin')
    return (rp - rm) / (2 * h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', type=Path, default=ROOT / 'results/jv/b21baseT')
    ap.add_argument('--cols', type=Path, default=ROOT / 'results/jv/frozenj_cols')
    ap.add_argument('--h', type=float, default=1e-4)
    ap.add_argument('--out', type=Path, default=ROOT / 'results/jv/frozenj_cols/columns_audit.json')
    a = ap.parse_args()
    x, R = read_trace_record(a.base / 'QSTARS_RESID_TRACE.bin')
    lo, di, up = read_ba_blocks(a.base / 'QSTARS_BA_TRACE.bin', NT, ND)
    if np.any(lo) or np.any(up):
        raise SystemExit('direct Jacobian is not block diagonal; column method invalid')
    rows = ~np.all(R == 0.0, axis=1)
    rows[133] = False   # T row: identity under FIX_T, compared separately
    res = {}
    miss = {}
    for name, var in (('ne', 133), ('hii', 32), ('he2g', 102), ('T', 134)):
        # State check: the plus pass must be x*(1+h) in row var only.
        xp, _ = read_trace_record(a.cols / f'{name}_p' / 'QSTARS_RESID_TRACE.bin')
        ratio = xp / x - 1.0
        off = np.delete(ratio, var - 1, axis=0)
        if np.max(np.abs(ratio[var - 1] - a.h)) > 1e-12 or np.max(np.abs(off)) > 1e-14:
            print(name, np.max(np.abs(ratio[var - 1] - a.h)), np.max(np.abs(off)))
            raise SystemExit(f'{name}: perturbed state is not x*(1+h e_{var})')
        fd = column_fd(a.cols / f'{name}_p', a.cols / f'{name}_m', a.h)   # d R / d ln x_var  (all depths)
        ba = np.stack([di[d][:, var - 1] for d in range(ND)], axis=1)       # same, analytic (POPS-scaled)
        m = fd - ba
        miss[name] = m
        den = np.linalg.norm(fd[rows])
        res[name] = dict(rel_missing=float(np.linalg.norm(m[rows]) / den) if den else None,
                         worst_depths=[int(i) + 1 for i in np.argsort(-np.linalg.norm(m[rows], axis=0))[:6]])
        print(f'{name:5s} column: ||fd - BA|| / ||fd|| = {res[name]["rel_missing"]:.3e}  worst depths {res[name]["worst_depths"]}')
    # Per unit absolute population: missing_hii / x_hii vs missing_he2g / x_he2g
    mh = miss['hii'] / x[31][None, :]
    mg = miss['he2g'] / x[101][None, :]
    num = np.linalg.norm((mh - mg)[rows]); den = np.linalg.norm(mh[rows])
    res['popion_hypothesis_rel_diff'] = float(num / den)
    print(f'missing(HII)/x_HII vs missing(HeII g)/x_HeIIg: rel diff {num / den:.3e}')
    a.out.write_text(json.dumps(res, indent=1))
    np.savez(a.out.with_suffix('.npz'), m_ne=miss['ne'], m_ion_abs=mh, m_he2g_abs=mg)


if __name__ == '__main__':
    main()
