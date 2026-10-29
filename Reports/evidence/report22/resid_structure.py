"""Where does the augmented least-squares residual of nk167_hot iterate 14 live (74 exact products, driver weights)?
Split ||W M^-1 (R - J s)||^2 by row class: ultra-trace (frac<=1e-20, weight 1e-8), trace (1e-20<frac<=1e-10, 1e-3),
significant populations, n_e row, T (RE) row; and by depth region. Also: residual of the T row alone if the T/front
subproblem were solved exactly (projection of the T-row residual onto the span is what matters)."""
import sys, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad')
from nkoff import setup_globals, load_iterate, load_pairs
x0 = setup_globals()
it, rec = load_iterate(14, x0)
Z, JZ, pas, _ = load_pairs(14, it.dir.name)
act = it.active
b = it.obj(it.R)
A = np.stack([it.obj(j).ravel() for j in JZ], 1)
y = np.linalg.lstsq(A, b.ravel(), rcond=1e-12)[0]
r = (b.ravel() - A @ y).reshape(b.shape)
tot = np.sum(r**2); tot0 = np.sum(b**2)
ultra = (it.frac <= 1e-20) & act; trace = (it.frac > 1e-20) & (it.frac <= 1e-10) & act
sig = (it.frac > 1e-10) & act
rows = np.zeros_like(act); rows[132] = True; rT = np.zeros_like(act); rT[133] = True
sigpop = sig & ~rows & ~rT
print(f'relative LS residual {np.sqrt(tot/tot0):.4f} (logged {rec["gmres"][-1]:.4f})')
for name, m in [('ultra-trace rows (w 1e-8)', ultra & ~rows & ~rT), ('trace rows (w 1e-3)', trace & ~rows & ~rT),
                ('significant population rows', sigpop), ('n_e row', rows & act), ('T (RE) row', rT & act)]:
    print(f'  {name:30s}: share of final residual^2 {np.sum(r[m]**2)/tot:7.4f}   share of initial {np.sum(b[m]**2)/tot0:7.4f}'
          f'   reduced to {np.sqrt(np.sum(r[m]**2)/max(np.sum(b[m]**2),1e-300)):.3f}')
# species split of ultra-trace residual
groups = {'H-/H I/H II (1-32)': slice(0, 32), 'He I (33-101)': slice(32, 101), 'He II (102-131)': slice(101, 131), 'He III (132)': slice(131, 132)}
for g, sl in groups.items():
    m = np.zeros_like(act); m[sl] = True
    print(f'  {g:20s}: final residual^2 share {np.sum(r[m & act]**2)/tot:7.4f} (ultra part {np.sum(r[m & ultra]**2)/tot:7.4f})')
ND = r.shape[1]
for lo, hi in [(0, 19), (19, 48), (48, 100), (100, ND)]:
    print(f'  depths {lo+1:3d}-{hi:3d}: final residual^2 share {np.sum(r[:, lo:hi][act[:, lo:hi]]**2)/tot:7.4f}')
# Is the T-row stall a missing subspace or a weighting conflict? Same 74 exact products, different LS row weights.
Mb = it.M.solve(it.R); Mi = [it.M.solve(j) for j in JZ]
def solve(Wt):
    A_ = np.stack([(Wt * m).ravel() for m in Mi], 1); b_ = (Wt * Mb).ravel()
    y_ = np.linalg.lstsq(A_, b_, rcond=1e-12)[0]
    rM = Mb - sum(yi * m for yi, m in zip(y_, Mi))
    s_ = sum(yi * z for yi, z in zip(y_, Z))
    return rM, s_
for name, uw in [('driver (ultra 1e-8)', 1e-8), ('ultra 1e-10', 1e-10), ('ultra rows excluded', 0.0)]:
    Wt = np.where(ultra, uw, it.W) * act
    rM, s_ = solve(Wt)
    tr = np.linalg.norm((it.W * rM)[133, 19:48]) / np.linalg.norm((it.W * Mb)[133, 19:48])
    ur = np.linalg.norm(rM[ultra]) / np.linalg.norm(Mb[ultra])
    print(f'{name:22s}: front T-row residual reduced to {tr:.3f}; ultra-row M^-1 residual to {ur:.3f}; '
          f'max|s_T| front {np.abs(s_[133, 19:48]).max():.4f}; max|s| ultra {np.abs(s_[ultra]).max():.3f}')
# T-row only objective: the best the span can do for the front RE row alone
WT = np.zeros_like(it.W); WT[133, 19:48] = 1.0
rM, s_ = solve(WT)
print(f'T row (front) alone     : reduced to {np.linalg.norm(rM[133,19:48])/np.linalg.norm(Mb[133,19:48]):.3f}')
