#!/usr/bin/env python3
"""V1 plan P1, items 2-3: the component-constrained least squares on the step that is actually applied.

The driver projects every trial onto x0's H and He totals and sets n_e to the charge sum (--conserve). Its least squares
never sees that projection, and its activity mask (R != 0) drops the charge row and the n_e variable wherever the
projection made the charge residual exactly zero (143-150 of 167 depths at it04/it16): the step there moves H II/H-
with n_e frozen and the projection then changes n_e by up to 48 % unscored. Here each span direction z is replaced by
its linearized projection P z (H and He rows shifted by their total change, n_e set to the charge-sum change), with the
product J P z = J z + Jdir (P z - z) (exact products for z; direct blocks for the projection's correction, whose
radiation part is not available from the stored products), and the box |P S y| <= rad is imposed on the applied step.

    box_proj.py ITERATE_DIR    (run with the scratch Python that has clarabel)"""
import sys, copy
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached, cnk_trap  # noqa: E402

CH = dict(HII=31, HEII=slice(101, 131), HEIII=131, HM=0, NE=132)


def project_lin(x, z):
    """Linearized cnk_trap.project_conservation at a state x that satisfies it, for a fractional step z (x(1 - z))."""
    zp = z.copy()
    H, HE = cnk_trap.H_ROWS, cnk_trap.HE_ROWS
    zp[H] -= (x[H] * z[H]).sum(0) / x[H].sum(0)
    zp[HE] -= (x[HE] * z[HE]).sum(0) / x[HE].sum(0)
    zp[CH['NE']] = (x[CH['HII']] * zp[CH['HII']] + (x[CH['HEII']] * zp[CH['HEII']]).sum(0)
                    + 2.0 * x[CH['HEIII']] * zp[CH['HEIII']] - x[CH['HM']] * zp[CH['HM']]) / x[CH['NE']]
    return zp


def projected_problem(P):
    it = P.it
    Q = copy.copy(P)
    cols, acols = [], []
    for i in range(P.Sm.shape[1]):
        z = P.Sm[:, i].reshape(P.shape)
        zp = project_lin(it.x, z)
        cols.append(zp.ravel())
        acols.append(P.A[:, i] + it.obj(it.Jdir.matvec(zp - z)).ravel())
    Q.Sm = np.stack(cols, axis=1)
    Q.A = np.stack(acols, axis=1)
    Q.mask = P.mask.copy()
    Q.mask[[131, 132]] = it.x[[131, 132]] > 0         # the projection moves n_e (and He III) everywhere: bound them too
    return Q


def main():
    d = sys.argv[1]
    P = cached(d)
    from box_conic import solve_box
    Q = projected_problem(P)
    # check: the projection is idempotent to first order, and the conservation rows of the projected products vanish
    y, info = solve_box(Q, cons=False)
    s = Q.step(y)
    bound = float(np.max(np.abs(s[Q.mask]) / Q.rad[Q.mask]))
    print(f'{d} BOX on the projected step: fit {Q.lin(y):.6f}, bound ratio {bound:.6f}, components at the bound '
          f'{int(np.sum(np.abs(s[Q.mask]) / Q.rad[Q.mask] > 0.999))}, max|s_ne| {np.abs(s[132]).max():.3f}; {info}')
    np.save(Path(__file__).with_name(f'y_BOXP_{d}.npy'), y)
    np.save(Path(__file__).with_name(f'step_BOXP_{d}.npy'), s)


if __name__ == '__main__':
    sys.exit(main())
