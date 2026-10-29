#!/usr/bin/env python3
"""Report 22: numerical agreement of one Newton-Krylov iterate re-run on a faster path with the archived iterate.

Both runs start from the same base pass (same state, residual, BA, stored J). Compared: the frozen-J columns, every
exact product pair (directions z and products J z of each Krylov pass), the least-squares history, the Newton step,
the trial merits and the accepted trial state/residual.

    compare_iterate.py ARCHIVED_RUN NEW_RUN K [--out JSON]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_trace_record  # noqa: E402


def rel(a, b):
    """Worst over rows (first axis: variable/equation) of max|a-b| / max|a| within the row, so that every equation
    is judged on its own scale (rows span many decades)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    ax = tuple(range(1, a.ndim))
    s = np.max(np.abs(a), axis=ax); d = np.max(np.abs(a - b), axis=ax)
    ok = s > 0
    return float(np.max(d[ok] / s[ok])) if ok.any() else float(np.max(d))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('archived', type=Path); ap.add_argument('new', type=Path); ap.add_argument('k', type=int)
    ap.add_argument('--out', type=Path)
    a = ap.parse_args()
    k = a.k
    ha = json.loads((a.archived / 'history.json').read_text()); hb = json.loads((a.new / 'history.json').read_text())
    ra = next(r for r in ha if r['it'] == k); rb = next(r for r in hb if r['it'] == k)
    base = Path(ra['base']).name
    out = dict(iterate=k, base=base)
    # frozen-J columns
    cols = {}
    for c in sorted((a.archived / f'cols_{base}').glob('col_*')):
        cb = a.new / f'cols_{base}' / c.name
        if (cb / 'QSTARS_RESID_TRACE.bin').exists():
            _, Ra = read_trace_record(c / 'QSTARS_RESID_TRACE.bin'); _, Rb = read_trace_record(cb / 'QSTARS_RESID_TRACE.bin')
            nz = Ra != 0
            cols[c.name] = dict(max_rel=float(np.max(np.abs(Ra[nz] - Rb[nz]) / np.abs(Ra[nz]))),
                                rel_l2=float(np.linalg.norm(Ra - Rb) / np.linalg.norm(Ra)))
    out['columns'] = cols
    # product pairs
    passes = []
    da, db = a.archived / f'it{k:02d}_from_{base}', a.new / f'it{k:02d}_from_{base}'
    for pa in sorted(da.glob('jv*')):
        pb = db / pa.name
        if not (pb / 'pairs.npz').exists():
            passes.append(dict(name=pa.name, missing=True)); continue
        A = np.load(pa / 'pairs.npz'); B = np.load(pb / 'pairs.npz')
        n = min(A['z'].shape[2], B['z'].shape[2])
        passes.append(dict(name=pa.name, ndirs=[int(A['z'].shape[2]), int(B['z'].shape[2])],
                           z_max_rel=max(rel(A['z'][..., i], B['z'][..., i]) for i in range(n)),
                           jz_max_rel=max(rel(A['jz'][..., i], B['jz'][..., i]) for i in range(n))))
    out['product_passes'] = passes
    out['gmres'] = dict(archived=ra.get('gmres'), new=rb.get('gmres'))
    out['products'] = dict(archived=ra.get('products'), new=rb.get('products'), reused_new=rb.get('products_reused'))
    sa, sb = np.load(da / 'newton_step.npy'), np.load(db / 'newton_step.npy')
    out['newton_step_max_rel'] = rel(sa, sb)
    out['max_step_significant'] = dict(archived=ra.get('max_step_significant'), new=rb.get('max_step_significant'))
    out['merit'] = dict(archived=ra.get('merit'), new=rb.get('merit'))
    out['trials'] = dict(archived=ra.get('trials'), new=rb.get('trials'))
    out['accepted'] = dict(archived=ra.get('accepted'), new=rb.get('accepted'))
    lab = ra.get('accepted')
    ta, tb = a.archived / f'it{k + 1:02d}_trial_{lab}', a.new / f'it{k + 1:02d}_trial_{lab}'
    if (tb / 'QSTARS_RESID_TRACE.bin').exists():
        xa, Ra = read_trace_record(ta / 'QSTARS_RESID_TRACE.bin'); xb, Rb = read_trace_record(tb / 'QSTARS_RESID_TRACE.bin')
        out['accepted_state_max_rel'] = float(np.max(np.abs(xa - xb) / np.abs(xa)))
        out['accepted_residual_rel_l2'] = float(np.linalg.norm(Ra - Rb) / np.linalg.norm(Ra))
    txt = json.dumps(out, indent=1)
    if a.out:
        a.out.write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
