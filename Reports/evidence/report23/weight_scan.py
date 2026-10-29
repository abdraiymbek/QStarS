import sys, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad')
import nkoff
from nkoff import load_iterate, load_pairs
import cool_newton_krylov as ck
k = int(sys.argv[1])
x0 = nkoff.setup_globals(); ck.ULTRA_WEIGHT = 1e-10
it, rec = load_iterate(k, x0)
Z, JZ, pas, _ = load_pairs(k, it.dir.name)
act = it.active; ultra = (it.frac <= 1e-20) & act; sig = (it.frac > 1e-10) & act
Mb = it.M.solve(it.R); Mi = [it.M.solve(j) for j in JZ]
print(f'it{k}: max|M^-1 R| ultra {np.abs(Mb[ultra]).max():.2e}, significant {np.abs(Mb[sig]).max():.2e}')
for uw in [1e-10, 1e-12, 1e-14, 0.0]:
    Wt = np.where(ultra, uw, it.W) * act
    A = np.stack([(Wt * m).ravel() for m in Mi], 1); b = (Wt * Mb).ravel()
    y = np.linalg.lstsq(A, b, rcond=1e-12)[0]
    rM = Mb - sum(yi * m for yi, m in zip(y, Mi)); s = sum(yi * z for yi, z in zip(y, Z))
    Wn = np.where(ultra, 0.0, it.W) * act          # objective restricted to non-ultra rows
    rel_nonultra = np.linalg.norm(Wn * rM) / np.linalg.norm(Wn * Mb)
    tr = np.linalg.norm(rM[133, 31:39]) / np.linalg.norm(Mb[133, 31:39])
    print(f'  ultra weight {uw:7.0e}: non-ultra LS residual {rel_nonultra:.3f}; front T row (d32-39) to {tr:.3f}; '
          f'max|s_sig| {np.abs(s[sig]).max():.3f}; front max|s_T| {np.abs(s[133,31:39]).max():.3f}; max|s_ultra| {np.abs(s[ultra]).max():.2e}')
