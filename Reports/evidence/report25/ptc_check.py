#!/usr/bin/env python3
"""Report 25: population pseudo-transient continuation at it04_trial_l1, offline, on the exact products at that base.
(J + delta D) s = R with D = sign(J_ii) x (row l1-norm of the direct blocks) on every statistical-equilibrium row (not
the H/He number, charge and T rows): the products of the shifted operator are exact (J z + delta D z), the preconditioner
is M + delta D. Each delta's least-squares step (in its own metric W M_delta^-1) is evaluated by the exact linear model
of the TRUE Jacobian: predicted row-scaled and fixed measures of R - J s, and the step size by region.
Output: ptc_check.json (--front: with the exact front columns, ptc_check_front.json)."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import front_check as fc  # noqa: E402
import cnk_trap  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import BlockTridiag, read_trace_record  # noqa: E402
CONSTRAINT_ROWS = (31, 131, 132, 133)          # H number, He number, charge, T (0-based)


def main():
    S_old, J_old = fc.load_pairs(sorted(fc.OLDP.glob('jv*')))
    S_pc, J_pc = fc.load_pairs(sorted((fc.O / 'superseded/it04_front_pieces_2211').glob('jv00p*')))
    S, JS = S_old + S_pc, J_old + J_pc
    ref = fc.iterate('it01_trial_l0.5', floor=0.0)
    rows_it = fc.iterate(fc.B, floor=0.0, ali=False)
    meas = dict(rowscale=lambda r: (rows_it.W * rows_it.D * r)[:133], fixed=lambda r: ref.W * ref.M.solve(r) * ref.active)
    xb, _ = read_trace_record(fc.O / fc.B / 'QSTARS_RESID_TRACE.bin')
    win = cnk_trap.front_window(xb, 0.5, 2)
    regions = dict(outer=(0, win[0]), front=(win[0], win[1] + 1), he_front=(win[1] + 1, 75), mid=(75, 118), deep=(118, 167))
    front = None
    if '--front' in sys.argv:                  # fine level with the exact front columns (--front-probe)
        FZ, FJ = fc.load_pairs(sorted(fc.NEWD.glob('front*')))
        meta = []
        for zz in FZ:
            r = int(np.nonzero(np.abs(zz).sum(1))[0][0]) + 1
            meta.append((r, int(np.nonzero(zz[r - 1])[0][0]) % 3))
        front = (FZ, FJ, meta, 3, win)
    base = fc.iterate(fc.B, floor=0.0, front=front)
    lo, di, up = base._pc
    jd = base.Jdir
    rsum = 1.0 / np.where(rows_it.D > 0, rows_it.D, np.inf)          # row l1-norm of the direct blocks (NT, ND)
    diag = np.stack([np.diagonal(jd_d) for jd_d in base._pc[1]], axis=1)   # preconditioner diagonal (NT, ND)
    Dm = np.where(diag >= 0, 1.0, -1.0) * rsum * base.active
    Dm[list(CONSTRAINT_ROWS)] = 0.0
    out = {}
    for delta in (0.0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0):
        di2 = di.copy()
        for d in range(di.shape[0]):
            di2[d][np.arange(di.shape[1]), np.arange(di.shape[1])] += delta * Dm[:, d]
        di3, _ = cnk_trap.floor_blocks(di2, 1e-2, base.active)
        M = BlockTridiag(lo, di3, up).factor()
        W = base.W
        Wm = lambda v: W * M.solve(v) * base.active
        Jd = [j + delta * Dm * s_ for s_, j in zip(S, JS)]
        A = np.stack([Wm(j).ravel() for j in Jd], axis=1); b = Wm(base.R).ravel()
        y = np.linalg.lstsq(A, b, rcond=1e-12)[0]
        s = sum(yi * si for yi, si in zip(y, S)); Js = sum(yi * ji for yi, ji in zip(y, JS))
        ev = {m: float(np.linalg.norm(g(base.R - Js)) / np.linalg.norm(g(base.R))) for m, g in meas.items()}
        ev_half = {m: float(np.linalg.norm(g(base.R - 0.5 * Js)) / np.linalg.norm(g(base.R))) for m, g in meas.items()}
        own = float(np.linalg.norm(A @ y - b) / np.linalg.norm(b))
        smax = {k_: float(np.max(np.abs(s[:, r0:r1][base.sig[:, r0:r1]]))) for k_, (r0, r1) in regions.items()}
        out[f'{delta:g}'] = dict(own=own, **ev, half=ev_half, max_s=smax)
        print(f'delta {delta:7.0e}: own {own:.3f}  true-J linear model: rowscale {ev["rowscale"]:.3f} fixed {ev["fixed"]:.3f} '
              f'(half step {ev_half["rowscale"]:.3f}/{ev_half["fixed"]:.3f})  max|s| ' +
              ' '.join(f'{k_} {v_:.3g}' for k_, v_ in smax.items()))
    Path(__file__).with_name('ptc_check_front.json' if front else 'ptc_check.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
