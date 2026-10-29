#!/usr/bin/env python3
"""Report 22: aggregate macOS `sample` call-graph files into self/inclusive time per function (all threads)."""
import re, sys
from collections import Counter
from pathlib import Path

def parse(path):
    incl, selfc = Counter(), Counter()
    total = 0
    stack = []
    for line in Path(path).read_text(errors='replace').splitlines():
        m = re.match(r'^(\s*)[+!:|\s]*(\d+)\s+(.+?)\s+\(in ([^)]+)\)', line)
        if not m:
            continue
        depth = len(line) - len(line.lstrip(' +!:|'))
        n = int(m.group(2)); fn = m.group(3).split(' + ')[0].strip(); lib = m.group(4)
        key = f'{fn} [{lib}]'
        while stack and stack[-1][0] >= depth:
            stack.pop()
        if not stack:
            total += n
        if key not in [s[1] for s in stack]:
            incl[key] += n
        if stack:
            selfc[stack[-1][1]] -= n
        selfc[key] += n
        stack.append((depth, key))
    return incl, selfc, total

def main():
    files = sys.argv[1:]
    I, S, T = Counter(), Counter(), 0
    for f in files:
        if 'Call graph' not in Path(f).read_text(errors='replace')[:20000]:
            continue
        i, s, t = parse(f); I.update(i); S.update(s); T += t
    print(f'total thread-samples {T}')
    print('-- self time (top 25)')
    for k, v in S.most_common(25):
        print(f'{100*v/max(T,1):6.1f}%  {k}')
    print('-- inclusive (top 25)')
    for k, v in I.most_common(25):
        print(f'{100*v/max(T,1):6.1f}%  {k}')

if __name__ == '__main__':
    main()
