#!/usr/bin/env python3
"""Report 23: exact Jacobian columns of one level at a range of depths (one batched tangent pass).

Directions are unit fractional perturbations of row ROW (1-based) at each depth D0..D1 of BASE. The pass replays the
base's stored J; products are the exact radiation parts (QSTARS_JV_RAD.bin), saved with the directions in pairs-like
form (z, drad) for offline use with the driver's direct blocks.

    probe_level_columns.py BASE OUT --row 3 --depths 26 46 --x0 X0 --exe EXE [--threads 6]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_fortran, read_jv_rad, read_trace_record  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('base', type=Path); ap.add_argument('out', type=Path)
ap.add_argument('--row', type=int, required=True); ap.add_argument('--depths', type=int, nargs=2, required=True)
ap.add_argument('--x0', type=Path, required=True); ap.add_argument('--exe', type=Path, required=True)
ap.add_argument('--threads', type=int, default=6)
ap.add_argument('--vadat', nargs=2, action='append', default=[])
a = ap.parse_args()
base = a.base.resolve(); out = a.out.resolve(); out.parent.mkdir(parents=True, exist_ok=True)
man = json.loads((base / 'run_manifest.json').read_text())['extra_env']
vad = [(p[-1][1:-1], p[0]) for p in map(str.split, (base / 'VADAT').read_text().splitlines())
       if len(p) >= 2 and p[-1] in ('[INC_AD]', '[INC_ADV]')] + [tuple(v) for v in a.vadat]
xb, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
nt, nd = xb.shape; ck.NT, ck.ND = nt, nd
d0, d1 = a.depths
dirs = np.zeros((nt, nd, d1 - d0 + 1))
for i, d in enumerate(range(d0, d1 + 1)):
    dirs[a.row - 1, d - 1, i] = 1.0
R = ck.Runner(a.x0, a.exe, Path(man['QSTARS_REPLACE_FORCE']), False, a.threads, out.parent / (out.name + '.log'), vad)
R.slim = False
R.run(out, read_fortran(base / 'QSTARS_STEP.bin', xb.shape), replay=base / 'QSTARS_RJ.bin', dirs=dirs, ba=False, store_j=False, threads=a.threads)
_, drad = read_jv_rad(out / 'QSTARS_JV_RAD.bin')
np.savez(out / 'probe.npz', z=dirs, drad=drad, row=a.row, depths=np.arange(d0, d1 + 1))
print('saved', out / 'probe.npz', drad.shape)
