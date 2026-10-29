import sys, struct, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record
Q = '/Users/abraiymbek/Projects/QStarS/'
gd = Q + sys.argv[1]
path = gd + '/QSTARS_SPECTRAL.bin'
f = open(path, 'rb'); hdr = f.read(8); nd, ncf = struct.unpack('<2i', f.read(8)); f.close()
rec = np.memmap(path, '<f8', 'r', 16 + 5 * nd * 8, (ncf, 2 + 6 * nd))
x, R = read_trace_record(gd + '/QSTARS_RESID_TRACE.bin')
nu = np.array(rec[:, 0]) / 1e15; fqw = np.array(rec[:, 1])
print('hdr', hdr, nd, ncf, 'nu range', nu.min(), nu.max(), 'sum fqw', fqw.sum())
TWOHCSQ = 0.0147452575; HDKT = 4.7994145; hP = 6.62607e-27; h15 = hP * 1e15
E = 3.288087 / np.arange(1, 31) ** 2; nu0 = E[0] - E[1]; A2s = 8.2249
def AY(y):
    u = y * (1 - y); fu = 4 * u
    return 24.56 * A2s * (u * (1 - fu ** 0.8) + 0.88 * u ** 1.53 * fu ** 0.8)
m = (nu < nu0)
y = nu[m] / nu0; ay = AY(y)
print('norm check: sum 0.5*A(y)dy*0.25 =', (0.5 * ay / (nu0 * 1e15) * fqw[m]).sum() * 0.25)
order = np.argsort(nu)
for d in [int(v) for v in sys.argv[2:]]:
    k = d - 1
    T = x[133, k]; n1 = x[1, k]; n2 = x[2, k]
    chi = np.array(rec[:, 2 + k]); chs = np.array(rec[:, 2 + nd + k]); eta = np.array(rec[:, 2 + 2 * nd + k]); J = np.array(rec[:, 2 + 3 * nd + k])
    RE_direct = 4 * np.pi * (fqw * ((chi - chs) * J - eta)).sum()
    YJ = J[m] / (TWOHCSQ * nu[m] ** 3)
    nub = nu0 - nu[m]
    Jb = np.interp(nub, nu[order], J[order])
    Yb = Jb / (TWOHCSQ * np.maximum(nub, 1e-30) ** 3)
    Yb = np.where(nub < nu.min(), 0.0, Yb)   # partner below grid: dropped (A(y)->0 there)
    Yw = np.exp(-HDKT * nub / T)
    w = 0.5 * ay / (nu0 * 1e15) * fqw[m]
    Nup_old = n1 * (w * YJ * Yw).sum(); Nup_ex = n1 * (w * YJ * Yb).sum()
    Eabs_old = n1 * (2 * w * h15 * nu[m] * YJ * Yw).sum(); Eabs_ex = n1 * (2 * w * h15 * nu[m] * YJ * Yb).sum()
    ybar = Eabs_old / (2 * Nup_old * h15 * nu0)
    D = Eabs_old - h15 * nu0 * Nup_old
    # two-photon emission (spontaneous) energy
    Eem = n2 * 0.25 * (2 * w * h15 * nu[m]).sum()
    RTc = 4 * np.pi * 1e-10 * R[133, k]
    # T-derivative of OLD term at fixed J, pops: d/dlnT of Eabs_old (in STEQ units)
    dlnT = (2 * w * h15 * nu[m] * YJ * Yw * (HDKT * nub / T)).sum() * n1 / (4 * np.pi * 1e-10)
    print('d%3d T=%.0fK n1=%.3e n2=%.3e | RE direct %.4e trace %.4e | 2g abs OLD %.4e (%.1f%% of RE) exact %.4e ratio %.1f | Nup OLD %.2f (per atom %.4e) exact %.3f ratio %.1f | <y>=%.3f D=%.3e | 2g emis %.2e | d(2g OLD)/dlnT [STEQ] %.3f' % (
        d, T * 1e4, n1, n2, RE_direct, RTc, Eabs_old, 100 * Eabs_old / RTc, Eabs_ex, Eabs_old / Eabs_ex, Nup_old, Nup_old / n1, Nup_ex, Nup_old / Nup_ex, ybar, D, Eem, dlnT))
    # brightness temperatures at a few wavelengths
    for lam in (1215.67 * 0.995, 1300, 1665, 2500, 3600, 4500, 6000, 10000):
        f0 = 2.99792458e18 / lam / 1e15; i = np.argmin(np.abs(nu - f0)); Yv = J[i] / (TWOHCSQ * nu[i] ** 3)
        print('     lam %7.1f A  J=%.3e  occupation %.3e  T_b=%.0f K  Wien(Te) occ %.3e' % (lam, J[i], Yv, HDKT * nu[i] / np.log(1 + 1 / Yv) * 1e4, np.exp(-HDKT * nu[i] / T)))
