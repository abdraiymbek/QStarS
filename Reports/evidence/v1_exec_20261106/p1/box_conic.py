#!/usr/bin/env python3
"""V1 plan P1, item 2: the explicit component-constrained least squares, solved as a second-order cone program
(Clarabel interior point; run with a Python that has clarabel, e.g. a scratch venv with --system-site-packages):

    min ||A y - b||  s.t.  |S y| <= rad on every active component (0.5; T 0.03)   [BOX]
                           and A_c y = b_c on the H number, He number and charge rows   [BOX+CONS]

A enters through its QR factor (||A y - b||^2 = ||R_A y - Q^T b||^2 + const, no A^T A), in the column-scaled basis.
Prints the fit ||A y - b|| / ||b||, the bound ratio, the conservation residual and the solver status; saves y.
Read-only.   box_conic.py ITERATE_DIR [--cons] [--rank-cut C]"""
import argparse, sys, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import clarabel
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached  # noqa: E402

CONS_ROWS = (31, 131, 132)


def solve_box(P, cons=False, rank_cut=0.0, log=print):
    A, b, Sm = P.A, P.b, P.Sm
    act = P.it.active.ravel(); msk = P.mask.ravel()
    cn = np.linalg.norm(A, axis=0); cn[cn == 0] = 1.0
    As = A / cn                                       # y = z / cn
    Q, RA = np.linalg.qr(As)
    c = Q.T @ b
    rest = float(np.sqrt(max(np.linalg.norm(b) ** 2 - np.linalg.norm(c) ** 2, 0.0)))
    V = None
    if rank_cut > 0:                                  # optional: work in the leading right singular subspace of A
        U_, s_, Vt = np.linalg.svd(RA)
        r = int(np.sum(s_ > rank_cut * s_[0]))
        V = Vt[:r].T
        RA, n = RA @ V, r
    n = RA.shape[1]
    Cs = Sm[msk] / P.rad.ravel()[msk][:, None] / cn[None, :]
    if V is not None:
        Cs = Cs @ V
    Csp = sp.csc_matrix(np.where(np.abs(Cs) > 0, Cs, 0.0))
    nb = Csp.shape[0]
    # x = [z (n), t]; cones: SOC (t, RA z - c), nonneg (1 - Cz, 1 + Cz), zero (A_c z - b_c)
    blocks, rhs, cones = [], [], []
    soc = sp.hstack([sp.csc_matrix((n + 0, n)), sp.csc_matrix((n, 1))])
    blocks.append(sp.hstack([sp.csc_matrix((1, n)), sp.csc_matrix(np.array([[-1.0]]))]))
    blocks.append(sp.hstack([sp.csc_matrix(-RA), sp.csc_matrix((RA.shape[0], 1))]))
    rhs += [np.zeros(1), -c]
    cones.append(clarabel.SecondOrderConeT(1 + RA.shape[0]))
    blocks.append(sp.hstack([Csp, sp.csc_matrix((nb, 1))])); blocks.append(sp.hstack([-Csp, sp.csc_matrix((nb, 1))]))
    rhs += [np.ones(nb), np.ones(nb)]
    cones.append(clarabel.NonnegativeConeT(2 * nb))
    if cons:
        nt, nd = P.shape
        cidx = np.concatenate([np.ravel_multi_index((np.full(nd, r), np.arange(nd)), (nt, nd)) for r in CONS_ROWS])
        cidx = cidx[act[cidx]]
        Aeq = As[cidx] if V is None else As[cidx] @ V
        blocks.append(sp.hstack([sp.csc_matrix(Aeq), sp.csc_matrix((len(cidx), 1))]))
        rhs.append(b[cidx])
        cones.append(clarabel.ZeroConeT(len(cidx)))
    Ac = sp.vstack(blocks).tocsc()
    bc = np.concatenate(rhs)
    Pm = sp.csc_matrix((n + 1, n + 1))
    q = np.zeros(n + 1); q[-1] = 1.0
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 400
    st.tol_gap_abs = 1e-12; st.tol_gap_rel = 1e-10; st.tol_feas = 1e-10
    t0 = time.time()
    sol = clarabel.DefaultSolver(Pm, q, Ac, bc, cones, st).solve()
    z = np.array(sol.x[:n])
    if V is not None:
        z = V @ z
    y = z / cn
    info = dict(status=str(sol.status), iterations=sol.iterations, seconds=time.time() - t0, rest=rest)
    return y, info


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('itdir'); ap.add_argument('--cons', action='store_true')
    ap.add_argument('--rank-cut', type=float, default=0.0)
    a = ap.parse_args()
    P = cached(a.itdir)
    y, info = solve_box(P, a.cons, a.rank_cut)
    s = P.step(y)
    bound = float(np.max(np.abs(s[P.mask]) / P.rad[P.mask]))
    nt, nd = P.shape
    act = P.it.active.ravel()
    cidx = np.concatenate([np.ravel_multi_index((np.full(nd, r), np.arange(nd)), (nt, nd)) for r in CONS_ROWS])
    cidx = cidx[act[cidx]]
    consr = float(np.linalg.norm(P.A[cidx] @ y - P.b[cidx]))
    tag = f'BOX{"+CONS" if a.cons else ""}{f" rank-cut {a.rank_cut:g}" if a.rank_cut else ""}'
    print(f'{a.itdir} {tag}: fit {P.lin(y):.6f}, bound ratio {bound:.6f}, conservation residual {consr:.2e}, '
          f'max|s| {np.abs(s[P.mask]).max():.3f}, components at the bound {int(np.sum(np.abs(s[P.mask]) / P.rad[P.mask] > 0.999))}; {info}')
    np.save(Path(__file__).with_name(f'y_{tag.replace(" ", "_").replace("+", "_")}_{a.itdir}.npy'), y)


if __name__ == '__main__':
    sys.exit(main())
