#!/usr/bin/env python3
"""Report 22: seconds of system sleep (macOS `pmset -g log` Sleep -> next Wake/DarkWake) inside a wall-clock window,
so that benchmark wall times can be stated net of sleep.

    sleep_overlap.py START END        (times as 'YYYY-MM-DD HH:MM:SS', local)
    sleep_overlap.py --pass DIR       (window from DIR/run_result.json started/finished, else file mtimes)
"""
import argparse, datetime as dt, json, re, subprocess
from pathlib import Path

FMT = '%Y-%m-%d %H:%M:%S'


def sleep_intervals():
    log = subprocess.run(['pmset', '-g', 'log'], capture_output=True, text=True).stdout.splitlines()
    out, t_sleep = [], None
    for line in log:
        m = re.match(r'^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) [+-]\d{4} (Sleep|Wake|DarkWake)\s', line)
        if not m:
            continue
        t = dt.datetime.strptime(m.group(1), FMT)
        if m.group(2) == 'Sleep' and 'Entering Sleep' in line:
            t_sleep = t
        elif m.group(2) in ('Wake', 'DarkWake') and ' from ' in line and t_sleep is not None:
            out.append((t_sleep, t)); t_sleep = None
    return out


def overlap(t0, t1):
    tot = 0.0
    for a, b in sleep_intervals():
        lo, hi = max(a, t0), min(b, t1)
        if hi > lo:
            tot += (hi - lo).total_seconds()
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('start', nargs='?'); ap.add_argument('end', nargs='?'); ap.add_argument('--pass', dest='pdir', type=Path)
    a = ap.parse_args()
    if a.pdir:
        rr = json.loads((a.pdir / 'run_result.json').read_text()) if (a.pdir / 'run_result.json').exists() else {}
        t0 = dt.datetime.fromtimestamp(min(p.stat().st_mtime for p in a.pdir.iterdir() if p.is_file() and not p.is_symlink()))
        t1 = dt.datetime.fromtimestamp(max(p.stat().st_mtime for p in a.pdir.iterdir() if p.is_file() and not p.is_symlink()))
        if 'started' in rr:
            t0 = dt.datetime.fromisoformat(rr['started'])
    else:
        t0, t1 = dt.datetime.strptime(a.start, FMT), dt.datetime.strptime(a.end, FMT)
    print(json.dumps(dict(start=t0.strftime(FMT), end=t1.strftime(FMT), wall_s=(t1 - t0).total_seconds(),
                          sleep_s=overlap(t0, t1))))


if __name__ == '__main__':
    main()
