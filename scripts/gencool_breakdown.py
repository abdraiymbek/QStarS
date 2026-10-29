#!/usr/bin/env python3
"""Report 21: per-process gas heating/cooling at chosen depths from CMFGEN's GENCOOL.

GENCOOL lists, per block of depths, the thermal-balance terms of the RE equation: bound-free
(per level, summed here), collisional and free-free cooling for each species, plus the adiabatic
ratios and the net cooling. Positive = net cooling of the gas, negative = net heating. This is a
diagnostic of the state that produced the file (an unconverged state is labelled as such).
"""
import argparse, json, re
from pathlib import Path

import numpy as np


def parse(path):
    lines = Path(path).read_text().splitlines()
    data, depths, cur, block = {}, [], None, []
    num = re.compile(r'^\s*[-+]?\d')
    i = 0
    while i < len(lines):
        s = lines[i]
        if s.strip().startswith('Depth') and i + 1 < len(lines):
            block = [int(v) for v in lines[i + 1].split()]
            depths += block
            i += 2
            continue
        if s.strip() and not num.match(s):
            cur = s.strip()
            i += 1
            continue
        if s.strip() and num.match(s) and cur is not None and block:
            vals = [float(v.replace('D', 'E')) for v in s.split()]
            if len(vals) == len(block):
                arr = data.setdefault(cur, {})
                for d, v in zip(block, vals):
                    arr[d] = arr.get(d, 0.0) + v      # per-level rows are summed
        i += 1
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gencool', type=Path)
    ap.add_argument('--depths', type=int, nargs=2, required=True, metavar=('FIRST', 'LAST'))
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    data = parse(a.gencool)
    ds = list(range(a.depths[0], a.depths[1] + 1))
    rows = {}
    for name, arr in data.items():
        if not any(d in arr for d in ds):
            continue
        rows[name] = [arr.get(d, float('nan')) for d in ds]
    keys = [k for k in rows if 'Cooling' in k or 'Net' in k or 'heating' in k or 'Temperature' in k or 'Electron' in k]
    print('depth'.ljust(46) + ''.join(f'{d:>11d}' for d in ds))
    for k in keys:
        print(k[:45].ljust(46) + ''.join(f'{v:11.3e}' for v in rows[k]))
    if a.out:
        a.out.write_text(json.dumps(dict(file=str(a.gencool), depths=ds, rows=rows), indent=1))


if __name__ == '__main__':
    main()
