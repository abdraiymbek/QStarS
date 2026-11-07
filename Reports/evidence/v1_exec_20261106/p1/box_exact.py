#!/usr/bin/env python3
"""V1 plan P1: the box-constrained least squares restricted to conservation-consistent steps, (P - I) S y = 0 at every
depth (P: the linearized projection of box_proj.py): the applied step equals the planned one to first order, so every
product in the model is an exact stored product (no direct-block approximation of the n_e response).
    box_exact.py ITERATE_DIR    (scratch Python with clarabel)"""
import sys, copy
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached  # noqa: E402
from box_proj import project_lin  # noqa: E402
import box_conic  # noqa: E402

def consistency_rows(P):
    it = P.it
    E = np.stack([(project_lin(it.x, P.Sm[:, i].reshape(P.shape)) - P.Sm[:, i].reshape(P.shape)).ravel()
              for i in range(P.Sm.shape[1])], axis=1)
    rows = np.nonzero(np.any(E != 0, axis=1))[0]
    E = E[rows]
    sc = np.linalg.norm(E, axis=1)
    return E / sc[:, None]


def solve(Pq, E):
    import scipy.sparse as sp, clarabel, time
    A, b, Sm = Pq.A, Pq.b, Pq.Sm
    msk = Pq.mask.ravel()
    cn = np.linalg.norm(A, axis=0); cn[cn == 0] = 1.0
    As = A / cn; Qr, RA = np.linalg.qr(As); c = Qr.T @ b
    n = RA.shape[1]
    Cs = sp.csc_matrix(np.where(np.abs(Sm[msk]) > 0, Sm[msk], 0.0) / Pq.rad.ravel()[msk][:, None] / cn[None, :])
    Es = sp.csc_matrix(E / cn[None, :]); nb = Cs.shape[0]; z1 = lambda m: sp.csc_matrix((m, 1))
    Ac = sp.vstack([sp.hstack([sp.csc_matrix((1, n)), sp.csc_matrix(np.array([[-1.0]]))]), sp.hstack([sp.csc_matrix(-RA), z1(n)]),
                    sp.hstack([Cs, z1(nb)]), sp.hstack([-Cs, z1(nb)]), sp.hstack([Es, z1(Es.shape[0])])]).tocsc()
    bc = np.concatenate([np.zeros(1), -c, np.ones(2 * nb), np.zeros(Es.shape[0])])
    cones = [clarabel.SecondOrderConeT(1 + n), clarabel.NonnegativeConeT(2 * nb), clarabel.ZeroConeT(Es.shape[0])]
    q = np.zeros(n + 1); q[-1] = 1.0
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 400
    sol = clarabel.DefaultSolver(sp.csc_matrix((n + 1, n + 1)), q, Ac, bc, cones, st).solve()
    return np.array(sol.x[:n]) / cn, str(sol.status)


def main():
    d = sys.argv[1]
    P = cached(d); it = P.it
    E = consistency_rows(P)
    print(f'{d}: {E.shape[0]} consistency rows; rank of E {np.linalg.matrix_rank(E, 1e-10)} of {E.shape[1]} directions')
    Q = copy.copy(P)
    Q.mask = P.mask.copy(); Q.mask[[131, 132]] = it.x[[131, 132]] > 0
    y, status = solve(Q, E)
    s = P.step(y)
    print(f'{d} BOX consistent: fit {P.lin(y):.6f}, bound ratio {np.max(np.abs(s[Q.mask]) / Q.rad[Q.mask]):.6f}, '
          f'consistency residual {np.abs(E @ y).max():.2e}, components at the bound {int(np.sum(np.abs(s[Q.mask]) / Q.rad[Q.mask] > 0.999))}; {status}')
    np.save(Path(__file__).with_name(f'y_BOXE_{d}.npy'), y)


if __name__ == '__main__':
    main()
