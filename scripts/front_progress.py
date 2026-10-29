#!/usr/bin/env python3
"""Report 23: progress of a released cool closure at the H I front, per accepted base.

For each base in history order (plus the newest base not yet in history): T and the RE-row imbalance
net/gross (R_RE over the RE row's l1 scale, as the PTC switch uses) over a depth window, with the iterate's
max significant Newton step and PTC delta. A coupled finite-T root is being approached when the front T keeps
rising while its net heating falls; net heating that stays flat or grows while T rises means the front is
still on the thermally unstable stretch.

    front_progress.py RUN [--depths 30 40] [--json OUT]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_ba_blocks, read_trace_record  # noqa: E402


def front(base, d0, d1):
    x, R = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
    nt, nd = x.shape
    lo, di, up = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', nt, nd)
    c = (np.abs(lo[:, nt - 1, :]) + np.abs(di[:, nt - 1, :]) + np.abs(up[:, nt - 1, :])).sum(axis=1)
    r = R[nt - 1] / np.where(c > 0, c, 1.0)
    w = slice(d0 - 1, d1)
    T = x[nt - 1] * 1e4
    return dict(T_front=[round(float(v), 1) for v in T[w]], net_gross_front=[float(v) for v in r[w]],
                max_net_gross_front=float(np.max(np.abs(r[w]))), net_gross_l2_all=float(np.linalg.norm(r)),
                max_net_gross_outside=float(np.max(np.abs(np.delete(r, np.arange(d0 - 1, d1))))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path)
    ap.add_argument('--depths', type=int, nargs=2, default=[30, 40])
    ap.add_argument('--json', type=Path)
    a = ap.parse_args()
    hist = json.loads((a.run / 'history.json').read_text())
    rows, seen = [], set()
    for h in hist:
        b = a.run / Path(h['base']).name
        if b.name in seen or not (b / 'QSTARS_BA_TRACE.bin').exists():
            continue
        seen.add(b.name)
        rows.append(dict(it=h['it'], base=b.name, max_step_significant=h.get('max_step_significant'),
                         ptc_delta=h.get('ptc_delta'), accepted=h.get('accepted'), **front(b, *a.depths)))
    log = (a.run / 'stdout.log').read_text().splitlines()
    newest = [ln.split(' base ')[1].split(':')[0] for ln in log if ' base it' in ln]
    if newest and newest[-1] not in seen and (a.run / newest[-1] / 'QSTARS_BA_TRACE.bin').exists():
        rows.append(dict(it=None, base=newest[-1], **front(a.run / newest[-1], *a.depths)))
    d0, d1 = a.depths
    print(f'{"base":20s} {"step":>7s} {"delta":>7s} {"Tmax(K)":>8s} {"Tmean":>7s} {"max net/gross":>13s} {"outside":>9s}')
    for r in rows:
        print(f'{r["base"]:20s} {r.get("max_step_significant") or float("nan"):7.3f} {r.get("ptc_delta") or float("nan"):7.4f} '
              f'{max(r["T_front"]):8.0f} {np.mean(r["T_front"]):7.0f} {r["max_net_gross_front"]:13.3e} {r["max_net_gross_outside"]:9.1e}')
    if a.json:
        a.json.write_text(json.dumps(dict(run=str(a.run), depths=a.depths, rows=rows), indent=1))


if __name__ == '__main__':
    main()
