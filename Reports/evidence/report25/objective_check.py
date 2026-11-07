#!/usr/bin/env python3
"""Report 25: does the span of the exact products at it04_trial_l1 contain a step that lowers BOTH measures? The least
squares solved directly in each measure (row-scaled W D, fixed W_ref M_ref^-1, and the run's metric W M^-1), each step
evaluated in all three by the exact linear model. Output: objective_check.json."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import front_check as fc  # noqa: E402


def main():
    S_old, J_old = fc.load_pairs(sorted(fc.OLDP.glob('jv*')))
    S_new, J_new = fc.load_pairs(sorted(fc.NEWD.glob('jv*')))
    ref = fc.iterate('it01_trial_l0.5', floor=0.0)
    rows_it = fc.iterate(fc.B, floor=0.0, ali=False)
    run = fc.iterate(fc.B)
    R = run.R
    meas = dict(rowscale=lambda r: (rows_it.W * rows_it.D * r)[:133], fixed=lambda r: ref.W * ref.M.solve(r) * ref.active,
                run=lambda r: run.W * run.M.solve(r) * run.active)
    out = {}
    for sname, (S, JS) in (('old56', (S_old, J_old)), ('new', (S_new, J_new)), ('both', (S_old + S_new, J_old + J_new))):
        for oname, f in meas.items():
            A = np.stack([f(j).ravel() for j in JS], axis=1); b = f(R).ravel()
            U, sv, Vt = np.linalg.svd(A, full_matrices=False)
            y = Vt.T @ ((U.T @ b) / sv * (sv > 1e-10 * sv[0]))
            Js = sum(yi * ji for yi, ji in zip(y, JS)); s = sum(yi * si for yi, si in zip(y, S))
            ev = {m: float(np.linalg.norm(g(R - Js)) / np.linalg.norm(g(R))) for m, g in meas.items()}
            smax = float(np.max(np.abs(s[run.sig])))
            out[f'{sname}/{oname}'] = dict(n=len(S), **{f'lin_res_{m}': v for m, v in ev.items()}, max_s=smax,
                                           max_abs_coef=float(np.max(np.abs(y))))
            print(f'{sname:6s} n={len(S):3d} LS in {oname:8s}: lin res rowscale {ev["rowscale"]:.3f} fixed {ev["fixed"]:.3f} '
                  f'run {ev["run"]:.3f}  max|s| {smax:.3g}  max|coef| {np.max(np.abs(y)):.3g}')
    Path(__file__).with_name('objective_check.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
