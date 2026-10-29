import sys, struct, numpy as np
sys.path.insert(0, '/Users/abraiymbek/Projects/QStarS/src')
from qstars.nlte_newton import read_trace_record
gd = '/Users/abraiymbek/Projects/QStarS/results/gencool_167/pin100_it05'
path = gd + '/QSTARS_SPECTRAL.bin'
f = open(path, 'rb'); f.read(8); nd, ncf = struct.unpack('<2i', f.read(8)); f.close()
rec = np.memmap(path, '<f8', 'r', 16 + 5 * nd * 8, (ncf, 2 + 6 * nd))
x, R = read_trace_record(gd + '/QSTARS_RESID_TRACE.bin')
nu = np.array(rec[:, 0]) / 1e15; fqw = np.array(rec[:, 1])
TWOHCSQ = 0.0147452575; HDKT = 4.7994145; hP = 6.62607e-27; h15 = hP * 1e15
E = 3.288087 / np.arange(1, 31) ** 2; nu0 = E[0] - E[1]; A2s = 8.2249
def AY(y):
    u = y * (1 - y); fu = 4 * u
    return 24.56 * A2s * (u * (1 - fu ** 0.8) + 0.88 * u ** 1.53 * fu ** 0.8)
for d in (35,):
    k = d - 1; T = x[133, k]; n1 = x[1, k]; n2 = x[2, k]
    chi = np.array(rec[:, 2 + k]); chs = np.array(rec[:, 2 + nd + k]); eta = np.array(rec[:, 2 + 2 * nd + k]); J = np.array(rec[:, 2 + 3 * nd + k])
    term = 4 * np.pi * fqw * ((chi - chs) * J - eta)
    m = nu < nu0; y = nu[m] / nu0
    t2 = np.zeros(ncf); t2[m] = 4 * np.pi * fqw[m] * (1.0 / (TWOHCSQ * nu[m] ** 3)) * (hP / (4 * np.pi)) * y * AY(y) * n1 * np.exp(-HDKT * (nu0 - nu[m]) / T) * J[m]
    rest = term - t2
    dop = nu0 * 50 / 2.998e5
    sel = np.abs(nu - nu0) < 8 * dop
    print('d%d: total RE %.4e ; two-photon OLD %.4e ; non-2g %.4e' % (d, term.sum(), t2.sum(), rest.sum()))
    print('  Ly-alpha zone (|x|<8 Doppler widths of 50 km/s): net non-2g %.4e erg = %.1f cm-3 s-1 at 10.2 eV ; npts %d' % (rest[sel].sum(), rest[sel].sum() / (h15 * nu0), sel.sum()))
    idx = np.nonzero(sel)[0]
    for i in idx:
        xD = (nu[i] - nu0) / dop
        print('     x=%+6.2f  nu=%.6f  J=%.4e  chi-chs=%.3e  eta/(chi-chs)=%.4e  net=%+.3e' % (xD, nu[i], J[i], chi[i] - chs[i], eta[i] / (chi[i] - chs[i]), rest[i]))
    oth = 0.0
    for n in range(3, 31):
        f0 = E[0] - E[n - 1]; s2 = np.abs(nu - f0) < 8 * f0 * 50 / 2.998e5
        oth += rest[s2].sum()
    print('  other Lyman lines n=3..30: net %.4e erg' % oth)
    print('  remaining (continua + subordinate lines): %.4e' % (rest.sum() - rest[sel].sum() - oth))
