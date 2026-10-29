#!/usr/bin/env python3
"""Summarize a cool-closure Newton-Krylov run: residual, frozen-J correction and
exact Newton step per iterate, split by species group and depth region."""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_trace_record  # noqa: E402

GROUPS = {'H0': (1, 1), 'HI': (2, 31), 'HII': (32, 32), 'HeI': (33, 101),
          'HeII': (102, 131), 'HeIII': (132, 132), 'ne': (133, 133), 'T': (134, 134)}
REGIONS = {'outer d1-9': (1, 9), 'cool d10-30': (10, 30), 'mid d31-60': (31, 60), 'inner d61-84': (61, 84)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path)
    a = ap.parse_args()
    hist = json.loads((a.run / 'history.json').read_text())
    rows = []
    for h in hist:
        it = h['it']
        base = Path(h['base'])
        x, R = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
        stepf = a.run / f'it{it:02d}_from_{base.name}' / 'newton_step.npy'
        s = np.load(stepf) if stepf.exists() else None
        rec = dict(it=it, resid_norm=float(np.linalg.norm(R)), merit_own_M=h['merit'],
                   max_frozen=h['max_frozen_correction'], gmres=h['gmres'], products=h['products'],
                   lam=h.get('accepted_lam'), merit_after=h.get('merit_new'))
        if s is not None:
            rec['newton_max'] = float(np.max(np.abs(s)))
            rec['newton_by_group'] = {g: float(np.max(np.abs(s[a_ - 1:b]))) for g, (a_, b) in GROUPS.items()}
            rec['newton_by_region'] = {r: float(np.max(np.abs(s[:, a_ - 1:b]))) for r, (a_, b) in REGIONS.items()}
            major = np.r_[np.arange(1, 31), np.arange(32, 101), [132]]   # HI, HeI levels, n_e (0-based rows)
            rec['newton_major_max'] = float(np.max(np.abs(s[major])))
        rows.append(rec)
        print(f"it{it:02d} ||R||={rec['resid_norm']:.3e} merit={rec['merit_own_M']:.3e} "
              f"maxNewton={rec.get('newton_max', float('nan')):.3e} major={rec.get('newton_major_max', float('nan')):.3e} "
              f"gmres={[round(g, 3) for g in h['gmres']]} lam={rec['lam']} merit_after={rec['merit_after']}")
        if s is not None:
            print('      by region: ' + ', '.join(f'{k}:{v:.2e}' for k, v in rec['newton_by_region'].items()))
            print('      by group : ' + ', '.join(f'{k}:{v:.2e}' for k, v in rec['newton_by_group'].items()))
    (a.run / 'progress.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
