#!/usr/bin/env python3
"""Offline GMRES driver for an exact-residual Newton-Krylov step in CMFGEN.

Unknowns are fractional corrections s (POPS <- POPS*(1-s)).  The exact Newton
equation J.(x0*s) = R0 is preconditioned by the frozen-J CMFGEN Jacobian M
(QS_ALI_LIN with QS_LIN_FAC=0):

    A s = b,   A s := M^{-1} J (x0*s),   b := M^{-1} R0 = s_fac0.

Each Krylov run k perturbs x0 -> x0*(1-eps_k*v_k) with eps_k = 1e-2/max|v_k|,
forms rhs = (R(x0*(1-eps_k*v_k)) - R0)/eps_k ~ -J(x0*v_k) in CMFGEN and dumps the
raw solution u_k = M^{-1} rhs, so A v_k = -u_k.

Layout under <root>:
  b.bin                      float64 NT*ND (copy of s_fac0)
  k<NN>/v.bin                Krylov input vector (written by this script)
  k<NN>/eps                  eps used (written by this script)
  k<NN>/QSTARS_SOL_TRACE.bin raw solution record (written by CMFGEN)
  step/s_star.bin            current GMRES Newton step (written by `solve`)

Commands:
  next      ingest every completed run, write the next v (or nothing if done)
  solve     write the GMRES step for all completed runs and its estimate
  evaluate  compare residual traces of line-search runs under step/
"""
import argparse, json, struct, sys
from pathlib import Path

import numpy as np

NT, ND = 134, 84
EPS_LEVEL = 1e-2


def read_first(path):
    data = path.read_bytes()
    cnt, nt, nd = struct.unpack_from('<3i', data, 0)
    n = nt * nd
    if (nt, nd) != (NT, ND) or len(data) < 12 + 16 * n:
        raise SystemExit(f'{path}: bad record (nt={nt} nd={nd} bytes={len(data)})')
    pops = np.frombuffer(data, '<f8', n, 12).copy()
    arr = np.frombuffer(data, '<f8', n, 12 + 8 * n).copy()
    return pops, arr


def load(path):
    a = np.fromfile(path, '<f8')
    if a.size != NT * ND:
        raise SystemExit(f'{path}: size {a.size}')
    return a


def save(path, a):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.asarray(a, '<f8').tofile(path)


def arnoldi(root):
    b = load(root / 'b.bin')
    beta = np.linalg.norm(b)
    V = [b / beta]
    H = []
    k = 0
    while True:
        d = root / f'k{k + 1:02d}'
        if not (d / 'v.bin').exists():
            break
        v = load(d / 'v.bin')
        if np.max(np.abs(v - V[k])) > 1e-12 * np.max(np.abs(V[k])):
            raise SystemExit(f'{d}/v.bin does not match the recomputed Arnoldi vector')
        sol = d / 'QSTARS_SOL_TRACE.bin'
        if not sol.exists():
            return b, beta, V, H, k, False
        _, u = read_first(sol)
        w = -u
        h = np.zeros(k + 2)
        for i in range(k + 1):
            h[i] = np.dot(w, V[i])
            w = w - h[i] * V[i]
        for i in range(k + 1):
            c = np.dot(w, V[i])
            h[i] += c
            w = w - c * V[i]
        h[k + 1] = np.linalg.norm(w)
        H.append(h)
        V.append(w / h[k + 1])
        k += 1
    return b, beta, V, H, k, True


def hess(H, k):
    Hm = np.zeros((k + 1, k))
    for j in range(k):
        Hm[:j + 2, j] = H[j][:j + 2]
    return Hm


def gmres_solution(beta, V, H, k):
    Hm = hess(H, k)
    e1 = np.zeros(k + 1)
    e1[0] = beta
    y, *_ = np.linalg.lstsq(Hm, e1, rcond=None)
    rel = np.linalg.norm(Hm @ y - e1) / beta
    s = np.zeros(NT * ND)
    for j in range(k):
        s += y[j] * V[j]
    return s, y, rel


def cmd_next(root):
    b, beta, V, H, k, complete = arnoldi(root)
    if not complete:
        print(f'run k{k + 1:02d} has no QSTARS_SOL_TRACE.bin yet')
        return
    if k:
        _, _, rel = gmres_solution(beta, V, H, k)
        print(f'after {k} Krylov vectors: preconditioned GMRES residual {rel:.4e}; '
              f'h(k+1,k)={H[-1][-1]:.3e}')
    v = V[k]
    eps = EPS_LEVEL / np.max(np.abs(v))
    d = root / f'k{k + 1:02d}'
    save(d / 'v.bin', v)
    (d / 'eps').write_text(f'{eps:.17e}\n')
    print(f'wrote {d}/v.bin  eps={eps:.6e}  max|eps*v|={eps * np.max(np.abs(v)):.3e}')


def cmd_solve(root):
    b, beta, V, H, k, _ = arnoldi(root)
    if k == 0:
        raise SystemExit('no completed Krylov runs')
    s, y, rel = gmres_solution(beta, V, H, k)
    save(root / 'step' / 's_star.bin', s)
    S = s.reshape((NT, ND), order='F')
    info = dict(k=k, gmres_rel=rel, y=y.tolist(), max_s=float(S.max()), min_s=float(S.min()),
                n_ge_1=int(np.sum(S >= 1.0)))
    (root / 'step' / 'solve.json').write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))
    for lam in (1.0, 0.5, 0.25, 0.1):
        print(f'lambda={lam}: min(1-lambda*s)={np.min(1 - lam * S):+.4e}')


def norms(R0, R):
    fl = np.maximum(np.abs(R0), 1e-3 * np.abs(R0).max(axis=0, keepdims=True))
    raw = np.linalg.norm(R) / np.linalg.norm(R0)
    sc = np.linalg.norm(R / fl) / np.linalg.norm(R0 / fl)
    dep = np.linalg.norm(R, axis=0) / np.linalg.norm(R0, axis=0)
    return raw, sc, dep


def cmd_evaluate(root, r0_trace):
    x0, R0 = read_first(r0_trace)
    R0 = R0.reshape((NT, ND), order='F')
    R0[NT - 1, :] = 0.0
    rows = []
    for d in sorted((root / 'step').glob('*/')):
        f = d / 'QSTARS_RESID_TRACE.bin'
        if not f.exists():
            continue
        _, R = read_first(f)
        R = R.reshape((NT, ND), order='F')
        R[NT - 1, :] = 0.0
        raw, sc, dep = norms(R0, R)
        rows.append(dict(run=d.name, raw=raw, scaled=sc, depth_median=float(np.median(dep)),
                         depth_max=float(dep.max()), worst_depth=int(np.argmax(dep)) + 1))
        print(f'{d.name:>16}: ||R||/||R0|| raw={raw:.4e} row-scaled={sc:.4e} '
              f'per-depth median={np.median(dep):.4e} max={dep.max():.3e}@d{np.argmax(dep) + 1}')
    (root / 'step' / 'evaluate.json').write_text(json.dumps(rows, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('command', choices=['next', 'solve', 'evaluate'])
    ap.add_argument('root', type=Path)
    ap.add_argument('--r0', type=Path, help='R0 QSTARS_RESID_TRACE.bin (evaluate)')
    a = ap.parse_args()
    if a.command == 'next':
        cmd_next(a.root)
    elif a.command == 'solve':
        cmd_solve(a.root)
    else:
        cmd_evaluate(a.root, a.r0)


if __name__ == '__main__':
    sys.exit(main())
