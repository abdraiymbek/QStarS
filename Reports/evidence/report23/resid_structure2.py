"""Residual structure of a recent nk167_hot iterate with the driver's CURRENT weights (ultra 1e-10)."""
import sys, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad')
import nkoff
from nkoff import load_iterate, load_pairs
import cool_newton_krylov as ck
k = int(sys.argv[1])
x0 = nkoff.setup_globals(); ck.ULTRA_WEIGHT = 1e-10
it, rec = load_iterate(k, x0)
Z, JZ, pas, _ = load_pairs(k, it.dir.name)
act = it.active
Mb = it.M.solve(it.R); Mi = [it.M.solve(j) for j in JZ]
W = it.W
A = np.stack([(W * m).ravel() for m in Mi], 1); b = (W * Mb).ravel()
y = np.linalg.lstsq(A, b, rcond=1e-12)[0]
r = b - A @ y; r = r.reshape(Mb.shape); b2 = b.reshape(Mb.shape)
tot = np.sum(r**2); tot0 = np.sum(b2**2)
print(f'it{k}: LS residual {np.sqrt(tot/tot0):.4f} (logged {rec["gmres"][-1]:.4f}); products {len(Z)}')
ultra = (it.frac <= 1e-20) & act; trace = (it.frac > 1e-20) & (it.frac <= 1e-10) & act; sig = (it.frac > 1e-10) & act
rows = np.zeros_like(act); rows[132] = True; rT = np.zeros_like(act); rT[133] = True
for name, m in [('ultra', ultra & ~rows & ~rT), ('trace', trace & ~rows & ~rT), ('sig pops', sig & ~rows & ~rT), ('n_e', rows & act), ('T', rT & act)]:
    print(f'  {name:9s} final share {np.sum(r[m]**2)/tot:6.3f} initial share {np.sum(b2[m]**2)/tot0:6.3f} reduced to {np.sqrt(np.sum(r[m]**2)/max(np.sum(b2[m]**2),1e-300)):.3f}')
groups = {'H (1-32)': slice(0, 32), 'He I (33-101)': slice(32, 101), 'He II (102-131)': slice(101, 131), 'He III (132)': slice(131, 132), 'ne': slice(132,133), 'T': slice(133,134)}
for g, sl in groups.items():
    m = np.zeros_like(act); m[sl] = True
    print(f'  {g:15s} final share {np.sum(r[m & act]**2)/tot:6.3f}')
ND = r.shape[1]
prof = np.array([np.sum(r[:, d][act[:, d]]**2) for d in range(ND)]) / tot
top = np.argsort(-prof)[:10]
print('  top depths (1-based, share):', [(int(d) + 1, round(float(prof[d]), 3)) for d in top])
s = sum(yi * z for yi, z in zip(y, Z))
print('  step s_T d28-48:', np.round(s[133, 27:48], 3))
# T row alone
WT = np.zeros_like(W); WT[133, 27:42] = 1.0
A2 = np.stack([(WT * m).ravel() for m in Mi], 1); b3 = (WT * Mb).ravel()
y2 = np.linalg.lstsq(A2, b3, rcond=1e-12)[0]
print(f'  front T row (d28-42) alone reducible to {np.linalg.norm(b3 - A2 @ y2)/np.linalg.norm(b3):.3f}')
