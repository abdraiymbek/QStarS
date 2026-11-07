"""Read-only quantification of the outward density increase in the Ball 2011 row-0 models.

Reads the preserved relaxed row-0 profiles of the 2026-10-14 Ball-track campaign
(benchmarks/results/2026-10-14-ball-track/*_relax_C*_n*/profile_row0_<seed>rsun.npz)
and the smoke-run history/final profile, and writes row0_inversion.json next to this
script. Nothing outside this directory is written.

Quantities per profile (cell-centred density, outward cell order):
  adjacent_max     max_i (rho_{i+1}-rho_i)/rho_i       (the registered Q2 metric)
  resolved_rise    max_{i<j} rho_j/rho_i - 1           (largest outward rise over any interval)
  rise on fixed physical windows: same max restricted to cells with
      T in [20, 40] kK, and with exterior mass M*-m in [1e-3, 1e-1] Msun
  location of the trough/peak of the resolved rise in r/R*, exterior mass, T, tau
  dlnrho/dlnr min/max in the rising run
  local radiative Eddington factor kappa*L_r/(4 pi c G m) at the trough/peak
  rho(T) on a common T grid (26-30 kK) for cross-resolution comparison
"""
import glob, json, math, os, sys
import numpy as np

ROOT = '/Users/abraiymbek/Projects/QStarS'
RES = os.path.join(ROOT, 'benchmarks/results/2026-10-14-ball-track')
OUT = os.path.dirname(os.path.abspath(__file__))
MSUN = 1.98847e33   # value irrelevant at the quoted precision; ratios only
RSUN = 6.957e10
G = 6.674e-8
C = 2.99792458e10

RUNS = ['20261014T045753.248504Z_relax_C1_n200',
        '20261014T045958.492497Z_relax_C0_n400',
        '20261014T050203.014930Z_relax_C2_n800',
        '20261014T050616.157669Z_relax_C5_n400',
        '20261014T050739.087855Z_relax_C6_n400',
        '20261014T050913.937107Z_relax_C7_n400']


def max_rise(rho, mask=None):
    """max over i<j (both in mask) of rho_j/rho_i - 1, with (i, j)."""
    idx = np.arange(len(rho)) if mask is None else np.flatnonzero(mask)
    best, bi, bj = 0.0, None, None
    run_min, run_i = math.inf, None
    for j in idx:
        if run_min < math.inf and rho[j] / run_min - 1 > best:
            best, bi, bj = rho[j] / run_min - 1, run_i, j
        if rho[j] < run_min:
            run_min, run_i = rho[j], j
    return best, bi, bj


def analyse(path):
    z = np.load(path)
    rho = z['density_cgs']; T = z['temperature_K']; kap = z['opacity_cgs']
    r = z['r_edges_cm']; m = z['enclosed_mass_edges_g']; L = z['luminosity_edges_erg_s']
    n = len(rho); R = r[-1]; M = m[-1]
    rc = 0.5 * (r[1:] + r[:-1]); mc = 0.5 * (m[1:] + m[:-1]); Lc = 0.5 * (L[1:] + L[:-1])
    dr = np.diff(r)
    # optical depth at cell centres measured inward from the outer face, plus 2/3
    # (approximate: assumes the outer face is the tau = 2/3 photosphere of the grey Ball surface)
    dtau = kap * rho * dr
    tau_face_out = 2 / 3 + np.r_[np.cumsum(dtau[::-1])[::-1][1:], 0.0]
    tau_c = tau_face_out + 0.5 * dtau
    q_ext = (M - mc) / MSUN
    gam = kap * Lc / (4 * math.pi * C * G * mc)
    inc = np.diff(rho) / rho[:-1]
    ia = int(np.argmax(inc))
    rise, i, j = max_rise(rho)
    rT, iT, jT = max_rise(rho, (T >= 2.0e4) & (T <= 4.0e4))
    rq, iq, jq = max_rise(rho, (q_ext >= 1e-3) & (q_ext <= 1e-1))
    loc = lambda k: dict(cell=int(k), r_over_R=float(rc[k] / R), exterior_mass_msun=float(q_ext[k]),
                         T_K=float(T[k]), tau_approx=float(tau_c[k]), rho_cgs=float(rho[k]),
                         kappa_cgs=float(kap[k]), gamma_rad=float(gam[k]))
    dlnrho_dlnr = np.diff(np.log(rho)) / np.diff(np.log(rc))
    seg = slice(i, j) if i is not None else slice(0, 0)
    # rho on a common temperature grid through the inversion (T decreases outward there)
    Tg = np.linspace(26.0e3, 30.0e3, 9)
    sel = (T > 2.0e4) & (T < 4.0e4)
    o = np.argsort(T[sel])
    rho_T = np.interp(Tg, T[sel][o], np.log(rho[sel][o]))
    return dict(path=os.path.relpath(path, ROOT), zones=n, R_star_rsun=float(R / RSUN),
                Teff_proxy_outer_cell_T_K=float(T[-1]),
                adjacent_max=float(inc[ia]), adjacent_max_location=loc(ia + 1),
                n_positive_adjacent_increments=int(np.sum(inc > 0)),
                positive_increment_cells=[int(k) for k in np.flatnonzero(inc > 0)],
                resolved_rise=float(rise), trough=loc(i) if i is not None else None, peak=loc(j) if j is not None else None,
                rise_T_20_40kK=float(rT), rise_extmass_1em3_1em1=float(rq),
                rising_run_width_r_over_R=float((rc[j] - rc[i]) / R) if i is not None else 0.0,
                rising_run_width_extmass_msun=float(q_ext[i] - q_ext[j]) if i is not None else 0.0,
                dlnrho_dlnr_max_in_run=float(np.max(dlnrho_dlnr[seg])) if i is not None else None,
                max_gamma_rad_profile=float(np.max(gam)), max_gamma_rad_location=loc(int(np.argmax(gam))),
                ln_rho_on_T_grid=dict(T_K=Tg.tolist(), ln_rho=rho_T.tolist()))


def smoke():
    import csv
    d = sorted(glob.glob(os.path.join(RES, '*_smoke_C1_n200')))[-1]
    rows = list(csv.DictReader(open(os.path.join(d, 'history.csv'))))
    f = lambda x: float(x) if x not in ('', 'None', 'nan') else math.nan
    acc = [r for r in rows if r['event'] in ('step', 'landing', 'initial', 'start')] or rows
    teff = np.array([f(r['effective_temperature_K']) for r in acc])
    mbh = np.array([f(r['bh_mass_msun']) for r in acc])
    inv = np.array([f(r['max_density_inversion']) for r in acc])
    marks = {}
    for thr in (1e-2, 0.1, 1.0, 3.0):
        k = np.flatnonzero(inv > thr)
        marks[f'first_inversion_gt_{thr:g}'] = None if not len(k) else dict(bh_mass_msun=float(mbh[k[0]]), Teff_K=float(teff[k[0]]), inversion=float(inv[k[0]]))
    for tt in (12000., 10000., 9000., 8000., 7000., 6000.):
        k = np.flatnonzero(teff < tt)
        marks[f'first_Teff_below_{tt:g}'] = None if not len(k) else dict(bh_mass_msun=float(mbh[k[0]]), Teff_K=float(teff[k[0]]), inversion=float(inv[k[0]]))
    sample = [dict(bh_mass_msun=float(mbh[k]), Teff_K=float(teff[k]), max_density_inversion=float(inv[k]))
              for k in np.unique(np.r_[0, np.linspace(0, len(acc) - 1, 25).astype(int), len(acc) - 1])]
    final = analyse(os.path.join(d, 'profile_final.npz'))
    return dict(run=os.path.relpath(d, ROOT), n_records=len(acc), events=sorted({r['event'] for r in rows}),
                marks=marks, sample=sample, final_profile=final)


if __name__ == '__main__':
    out = dict(row0={}, smoke=None)
    for run in RUNS:
        for seed in (6000, 3464):
            p = os.path.join(RES, run, f'profile_row0_{seed}rsun.npz')
            out['row0'][f'{run.split("_relax_")[1]}_seed{seed}'] = analyse(p)
    out['smoke'] = smoke()
    json.dump(out, open(os.path.join(OUT, 'row0_inversion.json'), 'w'), indent=1)
    hdr = f"{'run':16s} {'adjmax':>9s} {'npos':>4s} {'rise':>8s} {'rise_T':>8s} {'rise_q':>8s} {'r/R tr-pk':>15s} {'q_ext tr-pk Msun':>19s} {'T tr-pk kK':>13s} {'tau tr-pk':>13s} {'Gam pk':>6s} {'dlnr/dlnr':>9s}"
    print(hdr)
    for k, a in out['row0'].items():
        t, p = a['trough'], a['peak']
        print(f"{k:16s} {a['adjacent_max']:9.2e} {a['n_positive_adjacent_increments']:4d} {a['resolved_rise']:8.4f} {a['rise_T_20_40kK']:8.4f} {a['rise_extmass_1em3_1em1']:8.4f} "
              f"{t['r_over_R']:.4f}-{p['r_over_R']:.4f} {t['exterior_mass_msun']:.4f}-{p['exterior_mass_msun']:.4f} {t['T_K']/1e3:5.2f}-{p['T_K']/1e3:5.2f} "
              f"{t['tau_approx']:5.0f}-{p['tau_approx']:5.0f} {p['gamma_rad']:6.3f} {a['dlnrho_dlnr_max_in_run']:9.3f}  R*={a['R_star_rsun']:.1f}")
    for k, a in out['row0'].items():
        if k.endswith('seed6000'):
            print(k, ' '.join(f'{x:.5f}' for x in a['ln_rho_on_T_grid']['ln_rho']))
    s = out['smoke']
    print('smoke', s['run'], s['n_records'], s['events'])
    print(json.dumps(s['marks'], indent=1))
    for r in s['sample']:
        print(r)
    fp = s['final_profile']
    print('smoke final: adj', fp['adjacent_max'], 'rise', fp['resolved_rise'], 'trough', fp['trough'], 'peak', fp['peak'])
