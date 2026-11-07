#!/usr/bin/env python3
"""V1 plan P0: the frozen reference norm of a cool-closure campaign, next to the local norm of each iterate.

The driver builds and accepts each step in a LOCAL measure (--merit rowscale --accept-rowscale): ||w_k D_k R|| with the
row weights w_k and gross scales D_k of the current base k. Those change from base to base, so the driver's merits of
different iterates are not one convergence curve. Here the weights and scales of one reference iterate are frozen in a
file (`freeze`), and any pass is scored with them (`score`):

  frozen      ||w_ref D_ref R||, by depth region, by species group, the largest cells
  local       the same with the pass's own weights/scales (when its blocks are present), from the raw native blocks:
              without the driver's finite-difference n_e/ion/T column completion, so ~0.02 % off the driver's merit
  unfloored   R / (row l1 scale without --eq-floor): rows the floor keeps out of the measure stay visible (largest)
  new rows    active significant rows of the pass that weigh 0 in the frozen measure, scored with the pass's own
              scales: freezing must not hide an equation that has become significant
  physics     H and He number and charge conservation against x0, RE net/gross heating

    campaign_norms.py freeze ITERATE_DIR OUT.npz [--x0-base PASS]
    campaign_norms.py score FROZEN.npz PASS [PASS ...] [--json OUT.json]

ITERATE_DIR is a driver iterate directory itNN_from_BASE (its col_* passes complete the direct blocks); PASS a pass
directory (its <PASS>_tp two-photon replay is scored when present, as the driver does for a base). The measure
settings are those of the stage-2 runs (Report 25/26): sig-frac 1e-10, trace weight 1e-3, ultra-frac 1e-20, ultra
weight 0, weight ramp 2, eq-floor 1e-15, pc-t-floor 1e-3, RE row included, H I n=2 significant everywhere."""
from __future__ import annotations

import argparse, hashlib, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_ba_blocks, read_trace_record  # noqa: E402

SETTINGS = dict(sig_frac=1e-10, trace_weight=1e-3, ultra_frac=1e-20, ultra_weight=0.0, weight_ramp=2.0,
                merit='rowscale', eq_floor=1e-15, pc_t_floor=1e-3, sig_rows=[3], nrow=134)
X0_BASE = ROOT / 'results/nk167_trap/fixedT_b35/it00_base'
GROUPS = {'H-': [0], 'H I n=1/bal': [1], 'H I n=2': [2], 'H I n>=3': list(range(3, 31)), 'H num': [31],
          'He I': list(range(32, 101)), 'He II': list(range(101, 131)), 'He num': [131], 'charge': [132], 'RE': [133]}
REGIONS = [(1, 9), (10, 18), (19, 32), (33, 47), (48, 54), (55, 61), (62, 100), (101, 167)]


def configure(x0_base=X0_BASE):
    ck.SIG_FRAC, ck.TRACE_WEIGHT = SETTINGS['sig_frac'], SETTINGS['trace_weight']
    ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = SETTINGS['ultra_frac'], SETTINGS['ultra_weight'], SETTINGS['weight_ramp']
    ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = SETTINGS['merit'], SETTINGS['eq_floor'], SETTINGS['pc_t_floor']
    x0, _ = read_trace_record(Path(x0_base) / 'QSTARS_RESID_TRACE.bin')
    ck.NT, ck.ND = x0.shape
    return x0


def residual_dir(p):
    p = Path(p)
    tp = p.with_name(p.name + '_tp')
    return tp if (tp / 'QSTARS_RESID_TRACE.bin').exists() else p


def iterate_at(itdir, x0):
    """The driver's Iterate of an iterate directory itNN_from_BASE (no preconditioner changes: the measure only)."""
    itdir = Path(itdir)
    base = itdir.parent / itdir.name.split('_from_', 1)[1]
    cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(residual_dir(base), x0, cols, 1e-4, True, None, rj=base / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=SETTINGS['sig_rows'])
    return it


def gross_scales(lo, di, up):
    return (np.abs(lo) + np.abs(di) + np.abs(up)).sum(axis=2).T


def freeze(itdir, out, x0_base=X0_BASE):
    x0 = configure(x0_base)
    it = iterate_at(itdir, x0)
    nr = SETTINGS['nrow']
    w = (it.row_weights() * it.D * it.active)[:nr]
    np.savez(out, w=w, rs_gross=it.rs_gross[:nr], active=it.active[:nr], sig=it.sig[:nr], x_ref=it.x,
             reference=str(Path(itdir).resolve()), settings=json.dumps(SETTINGS), x0_base=str(Path(x0_base).resolve()))
    h = hashlib.sha256(Path(out).read_bytes()).hexdigest()
    print(f'frozen measure of {Path(itdir).name} -> {out} (sha256 {h}); reference merit {np.linalg.norm(w * it.R[:nr]):.6e}')
    return h


def decompose(r2):
    tot = float(r2.sum())
    dep = r2.sum(0) / tot
    cells = np.argsort(r2.ravel())[::-1][:12]
    nd = r2.shape[1]
    return dict(
        regions={f'd{a}-{b}': float(dep[a - 1:b].sum()) for a, b in REGIONS},
        groups={g: float(r2[i].sum() / tot) for g, i in GROUPS.items()},
        region_groups={f'd{a}-{b}': {g: float(r2[i, a - 1:b].sum() / tot) for g, i in GROUPS.items()
                                     if r2[i, a - 1:b].sum() / tot > 0.005} for a, b in REGIONS},
        largest_cells=[dict(row=int(c // nd) + 1, depth=int(c % nd) + 1, share=float(r2.ravel()[c] / tot)) for c in cells])


def score(frozen, p, x0):
    z = np.load(frozen)
    nr = SETTINGS['nrow']
    w = z['w']
    rd = residual_dir(p)
    x, R = read_trace_record(rd / 'QSTARS_RESID_TRACE.bin')
    R = R[:nr]
    r2 = (w * R) ** 2
    out = dict(pass_=str(Path(p).resolve()), residual_from=rd.name, frozen=float(np.sqrt(r2.sum())), **decompose(r2))
    # unfloored: every active row relative to the reference's own l1 scale (no eq-floor); the largest values only
    rs = z['rs_gross']
    act = (R != 0.0) & (rs > 0)
    rel = np.where(act, np.abs(R) / np.where(rs > 0, rs, 1.0), 0.0)
    top = np.argsort(rel.ravel())[::-1][:8]
    out['unfloored_largest'] = [dict(row=int(c // ck.ND) + 1, depth=int(c % ck.ND) + 1, rel=float(rel.ravel()[c]),
                                     frozen_weight=float(w.ravel()[c] * rs.ravel()[c])) for c in top]
    # the pass's own (local) measure and the rows it finds significant that the frozen measure does not weigh
    if (rd / 'QSTARS_BA_TRACE.bin').exists():
        lo, di, up = read_ba_blocks(rd / 'QSTARS_BA_TRACE.bin', ck.NT, ck.ND)
        rsl = gross_scales(lo, di, up)[:nr]
        rsf = np.maximum(rsl, SETTINGS['eq_floor'] * rsl.max(axis=0, keepdims=True))
        D = np.where(rsf > 0, 1.0 / np.where(rsf > 0, rsf, 1.0), 0.0)
        xH = x[1:32].sum(0); xHe = x[32:132].sum(0)
        frac = np.ones_like(x); frac[0:32] = x[0:32] / xH; frac[32:132] = x[32:132] / xHe
        active = (R != 0.0)
        Wl = ck.significance_weights(frac)[:nr] * active
        Wl[[31, 131, 132]] = active[[31, 131, 132]]
        for r_ in SETTINGS['sig_rows']:
            Wl[r_ - 1] = np.where(active[r_ - 1], 1.0, Wl[r_ - 1])
        loc2 = (Wl * D * R) ** 2
        out['local'] = float(np.sqrt(loc2.sum()))
        new = (Wl >= 1.0) & (w == 0.0) & active
        out['new_significant_rows'] = dict(count=int(new.sum()), local_share=float(loc2[new].sum() / max(loc2.sum(), 1e-300)),
                                           cells=[dict(row=int(i) + 1, depth=int(d) + 1) for i, d in zip(*np.nonzero(new))][:20])
    else:
        out['local'] = None
    # physics: exact conservation rows against x0 and the RE imbalance relative to the reference gross heating
    # (the largest |factor - 1| of the projection onto x0's H and He totals and the charge sum: 0 when exact)
    out['conservation_max_rel'] = cnk_trap.project_conservation(x, cnk_trap.conservation_totals(x0))[1]
    out['re_net_over_gross_max'] = float(np.max(np.abs(R[133]) / np.maximum(rs[133], 1e-300)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('freeze'); f.add_argument('itdir', type=Path); f.add_argument('out', type=Path)
    f.add_argument('--x0-base', type=Path, default=X0_BASE)
    s = sub.add_parser('score'); s.add_argument('frozen', type=Path); s.add_argument('passes', type=Path, nargs='+')
    s.add_argument('--json', type=Path, default=None)
    a = ap.parse_args()
    if a.cmd == 'freeze':
        freeze(a.itdir, a.out, a.x0_base)
        return 0
    z = np.load(a.frozen)
    x0 = configure(str(z['x0_base']))
    res = [score(a.frozen, p, x0) for p in a.passes]
    for r in res:
        print(f'{Path(r["pass_"]).name} ({r["residual_from"]}): frozen {r["frozen"]:.6e}'
              + (f', local {r["local"]:.6e}' if r['local'] is not None else '')
              + f'; regions ' + ', '.join(f'{k} {v:.3f}' for k, v in r['regions'].items())
              + f'; conservation {", ".join(f"{k} {v:.1e}" for k, v in r["conservation_max_rel"].items())}'
              + f'; RE net/gross max {r["re_net_over_gross_max"]:.2e}')
        if r.get('new_significant_rows', {}).get('count'):
            print(f'   {r["new_significant_rows"]["count"]} significant rows unweighted in the frozen measure '
                  f'({r["new_significant_rows"]["local_share"]:.3f} of the local measure^2)')
    if a.json:
        a.json.write_text(json.dumps(dict(frozen=str(a.frozen.resolve()),
                                          frozen_sha256=hashlib.sha256(a.frozen.read_bytes()).hexdigest(), scores=res), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
