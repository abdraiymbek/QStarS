#!/usr/bin/env python3
"""Offline (no native runs) re-solve of the cool_newton_krylov augmented least squares on the ARCHIVED exact
products of results/nk167_hot, with the driver's objective (ultra weight 1e-8) and with --ultra-weight 1e-10.
Products/directions are the archived exact pairs (pairs.npz); only the least-squares weights differ.
Usage: .venv/bin/python offline_krylov_eval.py [k ...]   (default 12 13 14)"""
import sys, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad')
from nkoff import setup_globals, load_iterate, load_pairs
ks = [int(a) for a in sys.argv[1:]] or [12, 13, 14]
x0 = setup_globals()
for k in ks:
    it, rec = load_iterate(k, x0)
    Z, JZ, pas, _ = load_pairs(k, it.dir.name)
    act = it.active; sig = (it.frac > 1e-10) & act; ultra = (it.frac <= 1e-20) & act
    Mb = it.M.solve(it.R); Mi = [it.M.solve(j) for j in JZ]
    print(f'it{k:02d} (logged driver residuals {[round(g, 4) for g in rec["gmres"][1:]]}, max|s_sig| {rec["max_step_significant"]:.4f})')
    for name, uw in (('ultra-weight 1e-8 (driver)', 1e-8), ('ultra-weight 1e-10', 1e-10)):
        Wt = np.where(ultra, uw, it.W) * act
        for p in range(pas.max() + 1):
            idx = np.nonzero(pas <= p)[0]
            A = np.stack([(Wt * Mi[i]).ravel() for i in idx], 1); bb = (Wt * Mb).ravel()
            y = np.linalg.lstsq(A, bb, rcond=1e-12)[0]
            s = sum(yi * Z[i] for yi, i in zip(y, idx)); rM = Mb - sum(yi * Mi[i] for yi, i in zip(y, idx))
            print(f'  {name:27s} passes {p + 1} ({len(idx)} products): residual {np.linalg.norm(bb - A @ y) / np.linalg.norm(bb):.4f} '
                  f'significant-row {np.linalg.norm(rM[sig]) / np.linalg.norm(Mb[sig]):.4f} max|s_sig| {np.abs(s[sig]).max():.4f} '
                  f's_T d31-38 {np.round(s[133, 31:39], 3)}')
