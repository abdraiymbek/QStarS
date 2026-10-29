#!/usr/bin/env python3
"""Report 22: wall/CPU profile of the archived Newton-Krylov passes of a run (no new computation).

Classifies every pass directory of a driver --out run:
  base   : start_base and it*_trial_* (full native pass: formal-ALI J solve, BA dump, stored J)
  column : cols_*/col_* (frozen-J replay, residual only)
  product: it*_from_*/jv* (replay + exact tangent for QSTARS_DIRS.bin directions)
and reads, per pass, /usr/bin/time -l (real, user, sys, max RSS) from native.log, the thread count and
extra_env from run_manifest.json, the number of tangent directions, and the formal-ALI statistics
(QSTARS_ALI_STATS: formal solutions, Ng steps, max iterations, unconverged frequencies). Timings of
passes that ran concurrently with other work are marked by the run's own concurrency (threads) only;
machine load is not recorded and must be stated when comparing.
"""
import argparse, json, re
from pathlib import Path

import numpy as np


def time_l(native):
    if not native.exists():
        return {}
    t = native.read_text(errors='replace')
    m = re.search(r'([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys', t)
    out = {}
    if m:
        out.update(real=float(m.group(1)), user=float(m.group(2)), sys=float(m.group(3)))
    m = re.search(r'(\d+)\s+maximum resident set size', t)
    if m:
        out['rss_gb'] = int(m.group(1)) / 1e9
    return out


def ali_stats(d):
    p = d / 'QSTARS_ALI_STATS'
    if not p.exists() or p.stat().st_mtime < (d / 'run_manifest.json').stat().st_mtime - 1:
        return {}
    rows = []
    for line in p.read_text().splitlines():
        m = re.search(r'max_it\s+(\d+)\s+unconverged\s+(\d+).*?max_res\s+(\S+)\s+formal_solutions\s+(\d+)\s+ng_steps\s+(\d+)', line)
        if m:
            rows.append(dict(max_it=int(m.group(1)), unconverged=int(m.group(2)), max_res=float(m.group(3)),
                             formal=int(m.group(4)), ng=int(m.group(5))))
    if not rows:
        return {}
    # a pass directory starts as a copy of its source, so earlier lines are inherited: only the LAST line
    # belongs to this pass (one formal-ALI J solve per native pass)
    r = rows[-1]
    return dict(formal_this_pass=r['formal'], ng_this_pass=r['ng'], max_it=r['max_it'], unconverged=r['unconverged'],
                max_res=r['max_res'])


def ndirs(d):
    p = d / 'QSTARS_DIRS.bin'
    if not p.exists():
        return 0
    return int(np.frombuffer(p.read_bytes()[:4], '<i4')[0])


def classify(d, run):
    rel = d.relative_to(run).as_posix()
    if rel.startswith('cols_'):
        return 'column'
    if re.match(r'it\d+_from_.*/jv\d+$', rel):
        return 'product'
    if rel == 'start_base' or re.match(r'it\d+_trial_', rel):
        return 'base'
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+', type=Path)
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    rows = []
    for run in a.runs:
        run = run.resolve()
        for man in run.rglob('run_manifest.json'):
            d = man.parent
            if 'superseded' in d.parts:
                continue
            kind = classify(d, run)
            if kind is None:
                continue
            m = json.loads(man.read_text())
            r = dict(run=run.name, pass_=d.relative_to(run).as_posix(), kind=kind, threads=m.get('threads'),
                     ndirs=ndirs(d), warm='QSTARS_ALI_WARM' in m.get('extra_env', {}),
                     exe=Path(m.get('executable', '')).parts[-4] if m.get('executable') else '')
            r.update(time_l(d / 'native.log'))
            r.update(ali_stats(d))
            rows.append(r)
    for kind in ('base', 'column', 'product'):
        k = [r for r in rows if r['kind'] == kind and 'real' in r]
        if not k:
            continue
        print(f'== {kind}: {len(k)} passes')
        for thr in sorted({r["threads"] for r in k}):
            kk = [r for r in k if r['threads'] == thr]
            real = np.array([r['real'] for r in kk]); user = np.array([r['user'] for r in kk])
            line = f'   threads {thr}: n={len(kk)} real median {np.median(real):.0f}s [{real.min():.0f},{real.max():.0f}]  user median {np.median(user):.0f}s'
            if kind == 'product':
                nd = np.array([r['ndirs'] for r in kk])
                line += f'  dirs median {np.median(nd):.0f} [{nd.min()},{nd.max()}]  user/dir median {np.median(user / np.maximum(nd, 1)):.0f}s'
            if kind == 'base':
                fs = [r['formal_this_pass'] for r in kk if 'formal_this_pass' in r]
                if fs:
                    line += f'  formal solutions/pass median {np.median(fs):.3g}'
                    unc = sum(1 for r in kk if r.get('unconverged', 0) > 0)
                    line += f'  passes with an unconverged frequency: {unc}'
            print(line)
    if a.out:
        a.out.write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
