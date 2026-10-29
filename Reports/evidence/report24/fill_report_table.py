#!/usr/bin/env python3
"""Report 24: write the constrained-restart trial table of Reports/Report_24.md from restart.json (between the markers)."""
import json
from pathlib import Path
import numpy as np
E = Path(__file__).resolve().parent; R = E.parents[1] / 'Report_24.md'
rs = json.loads((E / 'restart.json').read_text())
lines = ['| iterate | trial | native pass | driver natural merit | fixed measure | H fronts (23–45) | deep (119–167) | conservation |',
         '|---|---|---|---|---|---|---|---|']
for q in rs['rows']:
    deep = float(np.linalg.norm(q['by_region'][4:])) if 'by_region' in q else None
    lines.append(f"| {q['iterate']} | λ = {q['label'][1:]}{' (accepted)' if q['accepted'] else ''} | {q['native']} | "
                 f"{q['natural']:.2f} | {q['total']:.1f} | {q['by_region'][1]:.1f} | {deep:.1f} | {q['conservation_worst']:.0e} |"
                 if 'total' in q else f"| {q['iterate']} | λ = {q['label'][1:]} | {q['native']} | — | — | — | — | — |")
tab = '\n'.join(lines)
s = R.read_text()
a, b = '<!-- constrained-table:start -->', '<!-- constrained-table:end -->'
if 'CONSTRAINED_TABLE' in s:
    s = s.replace('CONSTRAINED_TABLE', f'{a}\n{tab}\n{b}')
else:
    s = s[:s.index(a)] + f'{a}\n{tab}\n{b}' + s[s.index(b) + len(b):]
R.write_text(s)
print(tab)
