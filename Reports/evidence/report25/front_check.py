#!/usr/bin/env python3
"""Report 25: the exact-front-column preconditioner (--front-probe) tested offline at base it04_trial_l1 of fixedT_b35
against the 56 exact products of the pre-front iterate-4 attempt at the same base (superseded/it04_from_it04_trial_l1_pre_front).

Preconditioners: OLD = fixed-J + local-ALI increment + SVD floor 1e-2 (the Report 24 run); NEW = the same with exact
front columns (front* probe passes of it04_from_it04_trial_l1) installed before the floor; NEW_NF = NEW without the
floor in the front window (--front-no-floor). For each:
  q per pair            ||W (M^-1 J z - z)|| / ||W z|| (0 perfect, 1 useless);
  front leakage         share of W M^-1 J z inside the front window for directions supported outside it;
  least squares         over the 56 pairs in that metric: linear residual in its own metric, in the fixed measure (it01
                        reference, Report 24) and row-scaled (no inverse); coefficients; max|s| by region.
Output: front_check.json (and a printed table)."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.EQ_FLOOR, ck.PC_T_FLOOR = 0.0, 1e-3
O = ROOT / 'results/nk167_trap/fixedT_b35'
B = 'it04_trial_l1'
OLDP = O / 'superseded/it04_from_it04_trial_l1_pre_front'
NEWD = O / f'it04_from_{B}'
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
z = np.load(O / 'homotopy_frontpin.npz'); hom = dict(Tfix=z['Tfix'], c=z['c'], w=z['w'], R0=None, mode='pin', lam=1.0)


def iterate(b_, front=None, floor=1e-2, no_floor_front=False, ali=True):
    ck.MERIT_MODE = 'precond'
    cols = {n: (O / f'cols_{b_}' / f'col_{n}_p', O / f'cols_{b_}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(O / (b_ + '_tp'), x0, cols, 1e-4, True, hom, rj=O / b_ / 'QSTARS_RJ.bin')
    if ali:
        it.install_ali_increment(O / (b_ + '_ali'))
    if front is not None:
        dirs, jzs, meta, stride, win = front
        it.install_exact_columns(dirs, jzs, meta, stride)
    if floor:
        it.floor_preconditioner(floor, skip=(front[4] if front is not None and no_floor_front else None))
    it.extend_significance(0.0, rows=[3])
    return it


def load_pairs(dirs):
    S, JS = [], []
    for d in dirs:
        if not (d / "pairs.npz").exists():
            continue
        Q = np.load(d / "pairs.npz")
        S += [Q['z'][:, :, i] for i in range(Q['z'].shape[2])]; JS += [Q['jz'][:, :, i] for i in range(Q['jz'].shape[2])]
    return S, JS


def main():
    S, JS = load_pairs(sorted(OLDP.glob('jv*')))
    fdirs = sorted(NEWD.glob('front*'))
    FZ, FJ = load_pairs(fdirs)
    rec = json.loads((O / 'history.json').read_text())
    xb, _ = read_trace_record(O / B / 'QSTARS_RESID_TRACE.bin')
    win = cnk_trap.front_window(xb, 0.5, 2)
    stride = 3
    meta = []
    for zz in FZ:                                     # (row, colour) of each probe direction
        r = int(np.nonzero(np.abs(zz).sum(1))[0][0]) + 1
        meta.append((r, int(np.nonzero(zz[r - 1])[0][0]) % stride))
    front = (FZ, FJ, meta, stride, win)
    ref = iterate('it01_trial_l0.5', floor=0.0)       # the fixed measure of Report 24 (it01 reference, ALI, no floor)
    rows_it = iterate(B, floor=0.0, ali=False)
    rowscale = lambda r: rows_it.W * rows_it.D * r    # row-scaled measure (Iterate.obj, merit rowscale)
    out = dict(base=B, window=[win[0] + 1, win[1] + 1], probe_dirs=len(FZ), products=len(S), variants={})
    regions = dict(outer=(0, win[0]), front=(win[0], win[1] + 1), he_front=(win[1] + 1, 75), mid=(75, 118), deep=(118, 167))
    for name, kw in (('OLD', {}), ('NEW', dict(front=front)), ('NEW_NF', dict(front=front, no_floor_front=True)),
                     ('NEW_NOFLOOR', dict(front=front, floor=0.0))):
        it = iterate(B, **kw)
        cond = cnk_trap.block_conditioning(it._pc[1], it.active, win)
        Wm = lambda v: it.W * it.M.solve(v)
        q = [float(np.linalg.norm(Wm(j) - it.W * s_) / max(np.linalg.norm(it.W * s_), 1e-300)) for s_, j in zip(S, JS)]
        A = np.stack([Wm(j).ravel() for j in JS], axis=1); b = Wm(it.R).ravel()
        y = np.linalg.lstsq(A, b, rcond=1e-12)[0]; s = sum(yi * si for yi, si in zip(y, S))
        Js = sum(yi * ji for yi, ji in zip(y, JS))
        own = float(np.linalg.norm(A @ y - b) / np.linalg.norm(b))
        fixed = float(np.linalg.norm(ref.W * ref.M.solve(it.R - Js)) / np.linalg.norm(ref.W * ref.M.solve(it.R)))
        rows = float(np.linalg.norm(rowscale(it.R - Js)) / np.linalg.norm(rowscale(it.R)))
        leak = []
        for i, s_ in enumerate(S):
            sup = np.abs(s_).sum(0) > 0
            if sup[win[0]:win[1] + 1].any():
                continue
            a = A[:, i].reshape(s_.shape)
            leak.append((i, float(np.linalg.norm(a[:, win[0]:win[1] + 1]) / max(np.linalg.norm(a), 1e-300)), float(y[i])))
        smax = {k_: float(np.max(np.abs(s[:, r0:r1][it.sig[:, r0:r1]]))) for k_, (r0, r1) in regions.items()}
        big = sorted(range(len(y)), key=lambda i: -abs(y[i]))[:5]
        v = dict(q_median=float(np.median(q)), q_p90=float(np.percentile(q, 90)), q_max=float(np.max(q)),
                 front_cond_median=float(np.median(cond)), front_cond_min=float(np.min(cond)),
                 lin_res_own=own, lin_res_fixed_measure=fixed, lin_res_rowscale=rows,
                 leak_median=float(np.median([l_[1] for l_ in leak])), leak_max=float(np.max([l_[1] for l_ in leak])),
                 leak_weighted=float(sum(abs(l_[2]) * l_[1] for l_ in leak) / sum(abs(l_[2]) for l_ in leak)),
                 max_s_by_region=smax, max_abs_coef=float(np.max(np.abs(y))),
                 largest_coefs=[(int(i), float(y[i]), int(np.nonzero(np.abs(S[i]).sum(0))[0][0]) + 1,
                                 int(np.nonzero(np.abs(S[i]).sum(0))[0][-1]) + 1) for i in big])
        out['variants'][name] = v
        print(f"{name:12s} q med/p90/max {v['q_median']:.2f}/{v['q_p90']:.2f}/{v['q_max']:.1f}  front cond med/min "
              f"{v['front_cond_median']:.1e}/{v['front_cond_min']:.1e}  lin res own/fixed/rowscale {own:.3f}/{fixed:.3f}/{rows:.3f}  "
              f"leak into front med/max/coef-weighted {v['leak_median']:.1%}/{v['leak_max']:.1%}/{v['leak_weighted']:.1%}  "
              f"max|coef| {v['max_abs_coef']:.3g}")
        print('             max|s| ' + ' '.join(f'{k_} {x_:.3g}' for k_, x_ in smax.items()))
    (Path(__file__).with_name('front_check.json')).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
