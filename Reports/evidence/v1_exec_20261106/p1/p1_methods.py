#!/usr/bin/env python3
"""V1 plan P1, items 1-2: step algebra on identical stored products, identical row weights. Read-only.

For one iterate (cached Problem: A = obj(J S), b = obj(R), span S, bound |s| <= 0.5, T 0.03 on every active component):
  conditioning  singular values of A (raw and with unit-norm columns), numerical rank by cutoff
  LM-NE         the driver: y(mu) from the normal equations (A^T A + mu G) y = A^T b, G = S^T S over active rows,
                bisection in mu until the bound holds
  LM-QR         the same family from the augmented system [A; sqrt(mu) R_s] y ~ [b; 0] (R_s from QR of S): no A^T A
  TSVD          truncated SVD of the column-scaled A at several cutoffs, then the LM-QR bound inside that subspace
  BOX           the explicit component-constrained least squares min ||A y - b|| s.t. |S y| <= rad (ADMM on the
                augmented QR factor, then an active-set polish with a KKT check)
  BOX+CONS      BOX with the linearized H number, He number and charge rows as equalities (A_c y = b_c)
For each: linear residual ratio ||A y - b|| / ||b||, bound ratio max|s|/rad, backward error, distance between methods.

    p1_methods.py ITERATE_DIR [--json OUT]"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached  # noqa: E402

CONS_ROWS = (31, 131, 132)          # H number, He number, charge (0-based)


def lm_qr_family(A, Rs, b):
    """y(mu) = argmin ||A y - b||^2 + mu ||R_s y||^2 through the augmented QR (scipy lstsq on the stacked system)."""
    def solve(mu):
        if mu <= 0:
            return sla.lstsq(A, b, cond=1e-12, lapack_driver='gelsd')[0]
        M = np.vstack([A, np.sqrt(mu) * Rs])
        rhs = np.concatenate([b, np.zeros(Rs.shape[0])])
        Q, R = np.linalg.qr(M)
        return sla.solve_triangular(R, Q.T @ rhs)
    return solve


def bisect(solve, ratio, scale):
    y = solve(0.0)
    if ratio(y) <= 1.0:
        return y, 0.0
    lo, hi = -14.0, 14.0
    while ratio(solve(scale * 10.0 ** hi)) > 1.0 and hi < 40:
        hi += 4.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if ratio(solve(scale * 10.0 ** mid)) > 1.0:
            lo = mid
        else:
            hi = mid
    mu = scale * 10.0 ** hi
    return solve(mu), mu


def box_lsq(A, b, C, Aeq=None, beq=None, rho=None, iters=4000, tol=1e-7, log=print):
    """min 0.5||A y - b||^2  s.t. -1 <= C y <= 1 [, Aeq y = beq]. ADMM with the splitting z = C y (and the equality as
    a heavily weighted extra block of the least squares, then enforced exactly in the polish). Polish: equality-
    constrained least squares with the ADMM active set held at its bounds, re-checked for feasibility and KKT signs."""
    n = A.shape[1]
    if rho is None:
        rho = (np.linalg.norm(A, 2) / max(np.linalg.norm(C, 2), 1e-300)) ** 2
    Ms = [A, np.sqrt(rho) * C]
    if Aeq is not None:
        Ms.append(1e3 * np.sqrt(rho) * Aeq / max(np.linalg.norm(Aeq, 2), 1e-300) * np.linalg.norm(C, 2))
    Q, R = np.linalg.qr(np.vstack(Ms))
    Atb = A.T @ b
    weq = (1e3 * np.sqrt(rho) / max(np.linalg.norm(Aeq, 2), 1e-300) * np.linalg.norm(C, 2)) if Aeq is not None else 0.0
    z = np.zeros(C.shape[0]); u = np.zeros_like(z)
    y = np.zeros(n)
    for k in range(iters):
        rhs = Atb + rho * C.T @ (z - u)
        if Aeq is not None:
            rhs = rhs + weq ** 2 * Aeq.T @ beq
        y = sla.solve_triangular(R, sla.solve_triangular(R, rhs, trans='T'))
        Cy = C @ y
        z_old = z
        z = np.clip(Cy + u, -1.0, 1.0)
        u = u + Cy - z
        rp = np.linalg.norm(Cy - z) / max(np.linalg.norm(Cy), 1.0)
        rd = rho * np.linalg.norm(C.T @ (z - z_old)) / max(np.linalg.norm(Atb), 1e-300)
        if k % 500 == 0:
            log(f'      ADMM {k}: primal {rp:.2e} dual {rd:.2e} fit {np.linalg.norm(A @ y - b) / np.linalg.norm(b):.6f}')
        if rp < tol and rd < tol:
            break
    # polish: active set from the ADMM multipliers, equality-constrained least squares, then check
    act = np.abs(u) > 1e-9
    sgn = np.sign(u[act])
    E = [C[act]]; e = [sgn]
    if Aeq is not None:
        E.append(Aeq); e.append(beq)
    E = np.vstack(E); e = np.concatenate(e)
    y_p = eq_lsq(A, b, E, e)
    feas = float(np.max(np.abs(C @ y_p)))
    return (y_p if feas <= 1.0 + 1e-6 else y), dict(admm_iters=k + 1, primal=rp, dual=rd, active=int(act.sum()),
                                                    polish_feasible=feas <= 1.0 + 1e-6, polish_max=feas)


def eq_lsq(A, b, E, e):
    """min ||A y - b|| s.t. E y = e (null-space method: particular solution + least squares in null(E))."""
    if E.shape[0] == 0:
        return sla.lstsq(A, b, cond=1e-12)[0]
    U, sv, Vt = np.linalg.svd(E, full_matrices=True)
    r = int(np.sum(sv > 1e-12 * sv[0]))
    y_p = Vt[:r].T @ ((U[:, :r].T @ e) / sv[:r])
    N = Vt[r:].T
    t = sla.lstsq(A @ N, b - A @ y_p, cond=1e-12)[0]
    return y_p + N @ t


def backward_error(A, b, y, G=None, mu=0.0):
    """Relative normwise residual of the optimality condition A^T (A y - b) + mu G y = 0."""
    g = A.T @ (A @ y - b) + (mu * (G @ y) if G is not None and mu > 0 else 0.0)
    nA = np.linalg.norm(A, 2)
    return float(np.linalg.norm(g) / (nA ** 2 * np.linalg.norm(y) + nA * np.linalg.norm(b)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('itdir'); ap.add_argument('--json', type=Path, default=None)
    a = ap.parse_args()
    t0 = time.time()
    P = cached(a.itdir)
    A, b, Sm = P.A, P.b, P.Sm
    act = P.it.active.ravel()
    msk = P.mask.ravel()
    radv = P.rad.ravel()
    Sa = Sm[act]                                     # penalty rows (the driver's wvec = active)
    C = Sm[msk] / radv[msk][:, None]                 # bound rows, scaled
    out = dict(iterate=a.itdir, directions=A.shape[1], dense=P.ndense, base_merit=P.beta)
    print(f'{a.itdir}: {A.shape[1]} directions, {C.shape[0]} bounded components; base merit {P.beta:.6e}')
    # conditioning
    sv = np.linalg.svd(A, compute_uv=False)
    cn_ = np.linalg.norm(A, axis=0); cn_[cn_ == 0] = 1.0
    sv_s = np.linalg.svd(A / cn_, compute_uv=False)
    out['conditioning'] = dict(cond_raw=float(sv[0] / sv[-1]), cond_scaled=float(sv_s[0] / sv_s[-1]),
                               rank_raw={f'{c:g}': int(np.sum(sv > c * sv[0])) for c in (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4)},
                               rank_scaled={f'{c:g}': int(np.sum(sv_s > c * sv_s[0])) for c in (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4)},
                               cond_normal_eq=float((sv[0] / sv[-1]) ** 2))
    print(f'  cond(A) {sv[0] / sv[-1]:.3e} (normal equations {(sv[0] / sv[-1]) ** 2:.1e}); column-scaled {sv_s[0] / sv_s[-1]:.3e}; '
          f'rank by cutoff {out["conditioning"]["rank_raw"]}')
    # QR of the penalty rows: ||S y||_active = ||R_s y||
    Rs = np.linalg.qr(Sa, mode='r')
    G = Rs.T @ Rs
    AtA = A.T @ A
    scale = max(np.trace(AtA) / max(np.trace(G), 1e-300), 1e-300)

    def ratio(y):
        return float(np.max(np.abs(C @ y)))
    res = {}
    # unbounded reference
    y_u = sla.lstsq(A, b, cond=1e-12, lapack_driver='gelsd')[0]
    res['unbounded'] = dict(fit=P.lin(y_u), bound_ratio=ratio(y_u))
    # LM-NE (driver)
    from p1_common import ck
    _, mu_ne, y_ne = ck.lm_step_bound((P.S, A, b), P.it.active.astype(float), P.rad, P.mask)
    res['LM-NE'] = dict(mu=mu_ne, fit=P.lin(y_ne), bound_ratio=ratio(y_ne), backward_error=backward_error(A, b, y_ne, G, mu_ne))
    # LM-QR at the driver's mu, and its own bisection
    solve = lm_qr_family(A, Rs, b)
    y_qr_same = solve(mu_ne)
    y_qr, mu_qr = bisect(solve, ratio, scale)
    res['LM-QR'] = dict(mu=mu_qr, fit=P.lin(y_qr), bound_ratio=ratio(y_qr), backward_error=backward_error(A, b, y_qr, G, mu_qr),
                        rel_diff_NE_vs_QR_same_mu=float(np.linalg.norm(Sm @ (y_ne - y_qr_same)) / np.linalg.norm(Sm @ y_qr_same)),
                        fit_same_mu=P.lin(y_qr_same), bound_ratio_same_mu=ratio(y_qr_same))
    print(f'  LM-NE (driver): mu {mu_ne:.3e}, fit {res["LM-NE"]["fit"]:.6f}, bound {res["LM-NE"]["bound_ratio"]:.4f}, '
          f'backward error {res["LM-NE"]["backward_error"]:.2e}')
    print(f'  LM-QR: mu {mu_qr:.3e}, fit {res["LM-QR"]["fit"]:.6f}, bound {res["LM-QR"]["bound_ratio"]:.4f}, backward error '
          f'{res["LM-QR"]["backward_error"]:.2e}; at the driver mu: step differs by {res["LM-QR"]["rel_diff_NE_vs_QR_same_mu"]:.2e} '
          f'(rel.), fit {res["LM-QR"]["fit_same_mu"]:.6f}, bound {res["LM-QR"]["bound_ratio_same_mu"]:.4f}   ({time.time() - t0:.0f} s)')
    # TSVD in the column-scaled basis, then the LM-QR bound inside the subspace
    U, s_, Vt = np.linalg.svd(A / cn_, full_matrices=False)
    res['TSVD'] = {}
    for c in (1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2):
        r = int(np.sum(s_ > c * s_[0]))
        V = (Vt[:r].T / cn_[:, None])                   # y = V t
        AV = A @ V
        y_t = V @ ((U[:, :r].T @ b) / s_[:r])
        Rsv = np.linalg.qr(Sa @ V, mode='r')
        solve_v = lm_qr_family(AV, Rsv, b)
        t_b, mu_t = bisect(solve_v, lambda t: ratio(V @ t), max(np.trace(AV.T @ AV) / np.trace(Rsv.T @ Rsv), 1e-300))
        res['TSVD'][f'{c:g}'] = dict(rank=r, fit_unbounded=P.lin(y_t), bound_ratio_unbounded=ratio(y_t),
                                     fit_bounded=P.lin(V @ t_b), mu=mu_t)
        print(f'  TSVD cutoff {c:g}: rank {r}, unbounded fit {P.lin(y_t):.4f} (bound ratio {ratio(y_t):.1f}), bounded fit '
              f'{P.lin(V @ t_b):.6f}')
    # explicit box constraint
    print(f'  BOX: ADMM on {C.shape[0]} bound rows   ({time.time() - t0:.0f} s)')
    y_b, info = box_lsq(A, b, C, log=print)
    res['BOX'] = dict(fit=P.lin(y_b), bound_ratio=ratio(y_b), **info)
    print(f'  BOX: fit {res["BOX"]["fit"]:.6f}, bound {res["BOX"]["bound_ratio"]:.4f}, {info}   ({time.time() - t0:.0f} s)')
    # with the linearized conservation rows as equalities (rows of A at the constraint rows: obj(J S) there)
    nt, nd = P.shape
    cidx = np.concatenate([np.ravel_multi_index((np.full(nd, r), np.arange(nd)), (nt, nd)) for r in CONS_ROWS])
    cidx = cidx[act[cidx]]
    Aeq, beq = A[cidx], b[cidx]
    y_c, info_c = box_lsq(A, b, C, Aeq=Aeq, beq=beq, log=print)
    res['BOX+CONS'] = dict(fit=P.lin(y_c), bound_ratio=ratio(y_c), cons_residual=float(np.linalg.norm(Aeq @ y_c - beq)),
                           cons_residual_LM_NE=float(np.linalg.norm(Aeq @ y_ne - beq)), base_cons=float(np.linalg.norm(beq)), **info_c)
    print(f'  BOX+CONS: fit {res["BOX+CONS"]["fit"]:.6f}, bound {res["BOX+CONS"]["bound_ratio"]:.4f}, linearized '
          f'conservation residual {res["BOX+CONS"]["cons_residual"]:.2e} (LM-NE {res["BOX+CONS"]["cons_residual_LM_NE"]:.2e}, '
          f'base {res["BOX+CONS"]["base_cons"]:.2e})   ({time.time() - t0:.0f} s)')
    out['methods'] = res
    steps = dict(LM_NE=y_ne, LM_QR=y_qr, BOX=y_b, BOX_CONS=y_c, unbounded=y_u)
    np.savez(Path(__file__).with_name(f'steps_{a.itdir}.npz'), **steps)
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
