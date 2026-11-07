#!/usr/bin/env python3
"""Report 25: the least squares of the front-pieces iterate (it04_from_it04_trial_l1/jv*) in three measures, against the
pre-front attempt's 56 products at the same base. For each direction set and metric (OLD = local-ALI + floor; NEW = with
the exact front columns): linear residual in its own metric, in the fixed measure (it01 reference) and row-scaled, and the
largest step components by region. Output: pieces_check.json."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import front_check as fc  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402


def main():
    np_ = fc.np
    S_old, J_old = fc.load_pairs(sorted(fc.OLDP.glob('jv*')))
    S_new, J_new = fc.load_pairs(sorted(fc.NEWD.glob('jv*')))
    FZ, FJ = fc.load_pairs(sorted(fc.NEWD.glob('front*')))
    xb, _ = read_trace_record(fc.O / fc.B / 'QSTARS_RESID_TRACE.bin')
    win = fc.cnk_trap.front_window(xb, 0.5, 2)
    meta = []
    for zz in FZ:
        r = int(np_.nonzero(np_.abs(zz).sum(1))[0][0]) + 1
        meta.append((r, int(np_.nonzero(zz[r - 1])[0][0]) % 3))
    ref = fc.iterate('it01_trial_l0.5', floor=0.0)
    rows_it = fc.iterate(fc.B, floor=0.0, ali=False)
    regions = dict(outer=(0, win[0]), front=(win[0], win[1] + 1), he_front=(win[1] + 1, 75), mid=(75, 118), deep=(118, 167))
    out = dict(window=[win[0] + 1, win[1] + 1], n_old=len(S_old), n_new=len(S_new), cases={})
    its = dict(OLD=fc.iterate(fc.B), NEW=fc.iterate(fc.B, front=(FZ, FJ, meta, 3, win)))
    for mname, it in its.items():
        for sname, (S, JS) in (('old56', (S_old, J_old)), ('new', (S_new, J_new)), ('both', (S_old + S_new, J_old + J_new))):
            if not S:
                continue
            Wm = lambda v: it.W * it.M.solve(v)
            A = np_.stack([Wm(j).ravel() for j in JS], axis=1); b = Wm(it.R).ravel()
            y = np_.linalg.lstsq(A, b, rcond=1e-12)[0]
            s = sum(yi * si for yi, si in zip(y, S)); Js = sum(yi * ji for yi, ji in zip(y, JS))
            own = float(np_.linalg.norm(A @ y - b) / np_.linalg.norm(b))
            fixed = float(np_.linalg.norm(ref.W * ref.M.solve(it.R - Js)) / np_.linalg.norm(ref.W * ref.M.solve(it.R)))
            rw = lambda r: rows_it.W * rows_it.D * r
            rows = float(np_.linalg.norm(rw(it.R - Js)) / np_.linalg.norm(rw(it.R)))
            rows_reg = {k_: float(np_.linalg.norm(rw(it.R - Js)[:133, r0:r1]) / max(np_.linalg.norm(rw(it.R)[:133, r0:r1]), 1e-300))
                        for k_, (r0, r1) in regions.items()}
            smax = {k_: float(np_.max(np_.abs(s[:, r0:r1][it.sig[:, r0:r1]]))) for k_, (r0, r1) in regions.items()}
            c = dict(n=len(S), lin_res_own=own, lin_res_fixed=fixed, lin_res_rowscale=rows, rowscale_by_region=rows_reg,
                     max_s_by_region=smax, max_abs_coef=float(np_.max(np_.abs(y))))
            out['cases'][f'{mname}/{sname}'] = c
            print(f'{mname:4s} {sname:6s} n={len(S):3d}  lin res own {own:.3f} fixed {fixed:.3f} rowscale {rows:.3f} '
                  f'(by region ' + ' '.join(f'{k_} {v_:.2f}' for k_, v_ in rows_reg.items()) + f')  max|coef| {c["max_abs_coef"]:.3g}')
            print('                 max|s| ' + ' '.join(f'{k_} {v_:.3g}' for k_, v_ in smax.items()))
    Path(__file__).with_name('pieces_check.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
