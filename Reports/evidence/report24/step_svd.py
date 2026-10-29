#!/usr/bin/env python3
"""Report 24: singular-value anatomy of the augmented least-squares Newton step of fixedT_b35 iterate 4 (base
it04_trial_l0.25, SVD-floored local-ALI preconditioner, the 46 exact products of the run). Solves min ||A y - b|| with
A = [W M^-1 J z_i], b = W M^-1 R as the driver does (rcond 1e-12), checks the result against the driver's saved step,
and compares truncated-SVD solutions: linear residual, max|s| over significant variables (by group), ||W s||."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402
ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'precond', 0.0, 1e-3
O = ROOT / 'results/nk167_trap/fixedT_b35'; b_ = 'it04_trial_l0.25'
B = O / 'branch_unconstrained_it2plus'           # iterates 2-5 were moved here at 2026-10-22 01:15
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
z = np.load(O / 'homotopy_frontpin.npz'); hom = dict(Tfix=z['Tfix'], c=z['c'], w=z['w'], R0=None, mode='pin', lam=1.0)
cols = {n: (B / f'cols_{b_}' / f'col_{n}_p', B / f'cols_{b_}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
it = ck.Iterate(B / (b_ + '_tp'), x0, cols, 1e-4, True, hom, rj=B / b_ / 'QSTARS_RJ.bin')
it.install_ali_increment(B / (b_ + '_ali'))
it.floor_preconditioner(1e-2)
it.extend_significance(0.0, rows=[3])
S, JS = [], []
for d in ('it04pre_from_' + b_, 'it04f2_from_' + b_):         # it04_from_ renamed it04f2_from_ at the restart
    for j in sorted((B / d).glob('jv*')):
        Q = np.load(j / 'pairs.npz')
        S += [Q['z'][:, :, i] for i in range(Q['z'].shape[2])]; JS += [Q['jz'][:, :, i] for i in range(Q['jz'].shape[2])]
A = np.stack([it.obj(j).ravel() for j in JS], axis=1); bo = it.obj(it.R).ravel()
Sm = np.stack([s.ravel() for s in S], axis=1)
U, sv, Vt = np.linalg.svd(A, full_matrices=False)
sig = it.sig; W = it.W * it.active
groups = [('H I', 1, 31), ('H II', 31, 32), ('He I', 32, 101), ('He II', 101, 131), ('He III', 131, 132), ('n_e', 132, 133)]
saved = np.load(B / ('it04f2_from_' + b_) / 'newton_step.npy')
out = dict(n_dirs=len(S), singular_values=(sv / sv[0]).tolist(), rows=[])
for rc in (1e-12, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
    y = np.linalg.lstsq(A, bo, rcond=rc)[0]; s = (Sm @ y).reshape(x0.shape)
    row = dict(rcond=rc, kept=int(np.sum(sv > rc * sv[0])), lin_res=float(np.linalg.norm(A @ y - bo) / np.linalg.norm(bo)),
               max_s_sig=float(np.max(np.abs(s[sig]))), W_s=float(np.linalg.norm(W * s)),
               by_group={g: float(np.max(np.abs(s[r0:r1][sig[r0:r1]]), initial=0.0)) for g, r0, r1 in groups},
               n_sig_over_0p5=int(np.sum(np.abs(s[sig]) > 0.5)))
    if rc == 1e-12:
        row['matches_driver_step'] = float(np.max(np.abs(s - saved)) / np.max(np.abs(saved)))
    out['rows'].append(row)
    print(f"rcond {rc:7.0e} kept {row['kept']:2d}/{len(sv)} res {row['lin_res']:.3f} max|s|sig {row['max_s_sig']:9.3e} "
          f"||Ws|| {row['W_s']:9.3e} n|s|>0.5 {row['n_sig_over_0p5']:5d} " + ' '.join(f"{g} {v:.2g}" for g, v in row['by_group'].items()))
print('singular values / max:', np.array2string(sv / sv[0], precision=2, max_line_width=200))
print('driver step match (rcond 1e-12):', out['rows'][0].get('matches_driver_step'))
Path(ROOT / 'Reports/evidence/report24/step_svd_it04.json').write_text(json.dumps(out, indent=1) + '\n')
