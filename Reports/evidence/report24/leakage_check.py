#!/usr/bin/env python3
"""Report 24: far-field leakage of the augmented least-squares directions (fixedT_b35).
For iterate 1 (raw local-ALI preconditioner, the run's 18 exact products) and iterate 4 (SVD floor 1e-2, 46 products):
the least-squares coefficient of every direction, the share of its preconditioned exact action W M^-1 J z that lies
outside the depths where z is nonzero, and where that action is largest. Also iterate 1 with the floor applied to the
same 18 directions. Output: leakage_check.json."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'precond', 0.0, 1e-3
O = ROOT / 'results/nk167_trap/fixedT_b35'; B = O / 'branch_unconstrained_it2plus'
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
z = np.load(O / 'homotopy_frontpin.npz'); hom = dict(Tfix=z['Tfix'], c=z['c'], w=z['w'], R0=None, mode='pin', lam=1.0)


def iterate(base_dir, b_, floor):
    cols = {n: (base_dir / f'cols_{b_}' / f'col_{n}_p', base_dir / f'cols_{b_}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(base_dir / (b_ + '_tp'), x0, cols, 1e-4, True, hom, rj=base_dir / b_ / 'QSTARS_RJ.bin')
    it.install_ali_increment(base_dir / (b_ + '_ali'))
    if floor:
        it.floor_preconditioner(floor)
    it.extend_significance(0.0, rows=[3])
    return it


def pairs(dirs):
    S, JS = [], []
    for d in dirs:
        for j in sorted(d.glob('jv*')):
            Q = np.load(j / 'pairs.npz')
            S += [Q['z'][:, :, i] for i in range(Q['z'].shape[2])]; JS += [Q['jz'][:, :, i] for i in range(Q['jz'].shape[2])]
    return S, JS


def analyse(it, S, JS, label):
    A = np.stack([it.obj(j).ravel() for j in JS], axis=1); b = it.obj(it.R).ravel()
    y = np.linalg.lstsq(A, b, rcond=1e-12)[0]; s = sum(yi * si for yi, si in zip(y, S))
    rows = []
    for i in np.argsort(np.abs(y))[::-1]:
        zz = S[i]; a = A[:, i].reshape(x0.shape); sup = np.abs(zz).sum(0) > 0; dep = np.nonzero(sup)[0]
        rr = np.nonzero(np.abs(zz).sum(1))[0]
        rows.append(dict(direction=int(i), coefficient=float(y[i]), support_depths=[int(dep[0]) + 1, int(dep[-1]) + 1],
                         support_rows=[int(rr[0]) + 1, int(rr[-1]) + 1],
                         action_outside_support=float(np.linalg.norm(a[:, ~sup]) / max(np.linalg.norm(a), 1e-300)),
                         largest_action_depth=int(np.argmax(np.abs(a).max(0))) + 1))
    ev = np.arange(117, 166, 2)
    out = dict(label=label, n_dirs=len(S), linear_residual=float(np.linalg.norm(A @ y - b) / np.linalg.norm(b)),
               max_abs_s_significant=float(np.max(np.abs(s[it.sig]))),
               max_abs_s_deep_even=float(np.max(np.abs(s[:, ev][it.sig[:, ev]]))), directions=rows)
    print(f"{label}: {len(S)} dirs, lin res {out['linear_residual']:.3f}, max|s| sig {out['max_abs_s_significant']:.3g}, "
          f"deep even {out['max_abs_s_deep_even']:.3g}")
    for r in rows[:4]:
        print(f"   dir {r['direction']:2d} y {r['coefficient']:9.3g} support d{r['support_depths'][0]}-{r['support_depths'][1]} rows "
              f"{r['support_rows'][0]}-{r['support_rows'][1]}: {r['action_outside_support']:.1%} outside, largest at d{r['largest_action_depth']}")
    return out


S1, J1 = pairs([O / 'it01f0_from_it01_trial_l0.5'])
res = [analyse(iterate(O, 'it01_trial_l0.5', 0.0), S1, J1, 'iterate 1, raw local ALI'),
       analyse(iterate(O, 'it01_trial_l0.5', 1e-2), S1, J1, 'iterate 1 directions, floor 1e-2')]
S4, J4 = pairs([B / 'it04pre_from_it04_trial_l0.25', B / 'it04f2_from_it04_trial_l0.25'])
res.append(analyse(iterate(B, 'it04_trial_l0.25', 1e-2), S4, J4, 'iterate 4, floor 1e-2'))
(Path(__file__).parent / 'leakage_check.json').write_text(json.dumps(res, indent=1) + '\n')
