#!/usr/bin/env python3
"""Report 24: line inversions created by a trial step at the He front (fixedT_b35 iterate-3 trial at lambda = 0.5).
Prints n_u g_l / (n_l g_u) of H I lines at depths 55-70 for the base and the trial state, and the trial/base ratio of
H I n = 8 and n = 12 (rows n+1, 1-based). Values > 1 are inverted lines."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_trace_record  # noqa: E402
O = ROOT / 'results/nk167_trap/fixedT_b35'
B = O / 'branch_unconstrained_it2plus'           # iterates 2-5 moved here at 2026-10-22 01:15
xb, _ = read_trace_record(B / 'it04_trial_l0.25/QSTARS_RESID_TRACE.bin')   # base of the iterate-4 trials
xp, _ = read_trace_record(B / 'it03_trial_l0.25/QSTARS_RESID_TRACE.bin')   # base of the iterate-3 trials (it04_trial_*)
x0, _ = read_trace_record(O / 'it00_base/QSTARS_RESID_TRACE.bin')
st = np.fromfile(B / 'it04_trial_l0.5/QSTARS_STEP.bin', dtype='<f8')[-x0.size:].reshape(x0.shape, order='F')
xt = x0 * (1 - st)
d = np.arange(54, 70)
out = dict(depths=(d + 1).tolist(), base='it03_trial_l0.25', trial='it04_trial_l0.5', lines={})
for u, l in [(12, 8), (10, 9), (11, 10), (16, 9)]:
    q = lambda x: (x[u] / (2 * u * u)) / (x[l] / (2 * l * l))
    out['lines'][f'{u}-{l}'] = dict(base=np.round(q(xp)[d], 4).tolist(), trial=np.round(q(xt)[d], 4).tolist())
for n in (8, 12):
    out[f'trial_over_base_n{n}'] = np.round((xt[n] / xp[n])[d], 4).tolist()
print(json.dumps(out, indent=1))
