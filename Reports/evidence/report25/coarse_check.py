#!/usr/bin/env python3
"""Report 25: the two-level (coarse-corrected) preconditioner at it04_trial_l1, offline, on every exact product made at
that base: the pre-front attempt (56), the front-pieces pass (51) and the slow modes Z (coarse* passes).
Fine levels OLD (local ALI + floor) and NEW (+ exact front columns); coarse modes ls and galerkin. For each:
q over the 107 non-coarse pairs, and the least squares in that metric over several direction sets, evaluated by the exact
linear model in the row-scaled and fixed measures. Output: coarse_check.json."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import front_check as fc  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
CUR = fc.NEWD
PIECES = fc.O / 'superseded/it04_front_pieces_2211'


def main():
    S_old, J_old = fc.load_pairs(sorted(fc.OLDP.glob('jv*')))
    S_pc, J_pc = fc.load_pairs(sorted(PIECES.glob('jv00p*')))
    S_z, J_z = fc.load_pairs(sorted(CUR.glob('coarse*')))
    S_cur, J_cur = fc.load_pairs(sorted(CUR.glob('jv*')))
    FZ, FJ = fc.load_pairs(sorted(CUR.glob('front*')))
    xb, _ = read_trace_record(fc.O / fc.B / 'QSTARS_RESID_TRACE.bin')
    win = cnk_trap.front_window(xb, 0.5, 2)
    meta = []
    for zz in FZ:
        r = int(np.nonzero(np.abs(zz).sum(1))[0][0]) + 1
        meta.append((r, int(np.nonzero(zz[r - 1])[0][0]) % 3))
    front = (FZ, FJ, meta, 3, win)
    pre = fc.iterate(fc.B, front=front, floor=0.0)
    Y, Z, info = cnk_trap.slow_modes(pre._pc[1], pre.active, win, 1e-2, rows=np.r_[0:32, 132], min_weight=0.5, exclude=[133])
    assert len(Z) == len(S_z) and all(np.array_equal(a_, b_) for a_, b_ in zip(Z, S_z)), 'coarse directions differ'
    ref = fc.iterate('it01_trial_l0.5', floor=0.0)
    rows_it = fc.iterate(fc.B, floor=0.0, ali=False)
    meas = dict(rowscale=lambda r: (rows_it.W * rows_it.D * r)[:133], fixed=lambda r: ref.W * ref.M.solve(r) * ref.active)
    regions = dict(outer=(0, win[0]), front=(win[0], win[1] + 1), he_front=(win[1] + 1, 75), mid=(75, 118), deep=(118, 167))
    sets = dict(old56=(S_old, J_old), pieces51=(S_pc, J_pc), Z=(S_z, J_z), old56_Z=(S_old + S_z, J_old + J_z),
                pieces51_Z=(S_pc + S_z, J_pc + J_z), all=(S_old + S_pc + S_z, J_old + J_pc + J_z))
    if S_cur:
        sets['current'] = (S_cur, J_cur); sets['current_Z'] = (S_z + S_cur, J_z + J_cur)
    out = dict(n_modes=len(Z), mode_depths=sorted(set(d + 1 for d, _ in info)), variants={})
    for fine in ('OLD', 'NEW'):
        base_it = fc.iterate(fc.B) if fine == 'OLD' else fc.iterate(fc.B, front=front)
        for cmode in (None, 'ls', 'galerkin'):
            it = base_it
            M = base_it.M if cmode is None else cnk_trap.CoarseCorrected(base_it.M, Y, Z, J_z, W=base_it.W * base_it.active,
                                                                         mode=cmode)
            Wm = lambda v: it.W * M.solve(v) * it.active
            q = [float(np.linalg.norm(Wm(j) - it.W * s_) / max(np.linalg.norm(it.W * s_), 1e-300))
                 for s_, j in zip(S_old + S_pc, J_old + J_pc)]
            name = f'{fine}/{cmode or "none"}'
            v = dict(q_median=float(np.median(q)), q_p90=float(np.percentile(q, 90)), q_max=float(np.max(q)),
                     cond=getattr(M, 'cond', None), sets={})
            print(f'{name:14s} q med/p90/max {v["q_median"]:.2f}/{v["q_p90"]:.2f}/{v["q_max"]:.1f}'
                  + (f'  coarse cond {M.cond:.2e}' if cmode else ''))
            for sname, (S, JS) in sets.items():
                A = np.stack([Wm(j).ravel() for j in JS], axis=1); b = Wm(it.R).ravel()
                y = np.linalg.lstsq(A, b, rcond=1e-12)[0]
                Js = sum(yi * ji for yi, ji in zip(y, JS)); s = sum(yi * si for yi, si in zip(y, S))
                ev = {m: float(np.linalg.norm(g(it.R - Js)) / np.linalg.norm(g(it.R))) for m, g in meas.items()}
                own = float(np.linalg.norm(A @ y - b) / np.linalg.norm(b))
                smax = {k_: float(np.max(np.abs(s[:, r0:r1][it.sig[:, r0:r1]]))) for k_, (r0, r1) in regions.items()}
                v['sets'][sname] = dict(n=len(S), own=own, **ev, max_s=smax, max_abs_coef=float(np.max(np.abs(y))))
                print(f'    {sname:10s} n={len(S):3d} own {own:.3f} rowscale {ev["rowscale"]:.3f} fixed {ev["fixed"]:.3f}  '
                      f'max|s| ' + ' '.join(f'{k_} {x_:.3g}' for k_, x_ in smax.items()) + f'  max|coef| {np.max(np.abs(y)):.3g}')
            out['variants'][name] = v
    Path(__file__).with_name('coarse_check.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
