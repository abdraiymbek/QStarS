#!/usr/bin/env python3
"""Report 21 reproducibility test of a tangent product pass (after two transient SIGSEGVs).

Reruns an existing product pass (same state QSTARS_STEP.bin, same replay field, same QSTARS_DIRS.bin,
same executable and VADAT) into a new directory and compares QSTARS_JV_RAD.bin. With the dynamic
ray schedule and the moment reduction, only the summation order can differ (relative ~1e-15);
anything larger would indicate a race in the threaded tangent.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_fortran, read_jv_rad  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pass_dir', type=Path)
    ap.add_argument('dest', type=Path)
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--exe', type=Path, default=None, help='executable for the repeat (default: the source pass executable)')
    a = ap.parse_args()
    src = a.pass_dir.resolve()
    man = json.loads((src / 'run_manifest.json').read_text())
    env = man['extra_env']
    raw = (src / 'QSTARS_DIRS.bin').read_bytes()
    nv, nt, nd = np.frombuffer(raw[:12], "<i4")
    dirs = np.frombuffer(raw[12:], '<f8').reshape((nt, nd, nv), order='F')
    ck.NT, ck.ND = int(nt), int(nd)
    s_cum = read_fortran(src / 'QSTARS_STEP.bin', (nt, nd))
    vad = []
    for line in (src / 'VADAT').read_text().splitlines():
        for key in ('INC_AD', 'INC_ADV'):
            if line.split() and line.split()[-1] == f'[{key}]':
                vad.append((key, line.split()[0]))
    R = ck.Runner(a.x0, a.exe or Path(man['executable']), a.replace, False, a.threads, a.dest.parent / 'repeat.log', vad)
    R.run(a.dest, s_cum, replay=env['QSTARS_JV_REPLAY'], dirs=dirs, ba=False, store_j=False)
    _, r1 = read_jv_rad(src / 'QSTARS_JV_RAD.bin')
    _, r2 = read_jv_rad(a.dest / 'QSTARS_JV_RAD.bin')
    out = dict(source=str(src), repeat=str(a.dest), directions=int(nv),
               source_exe=man['executable'], repeat_exe=str(a.exe or man['executable']),
               max_rel_diff=[float(np.max(np.abs(r1[..., k] - r2[..., k])) / max(np.max(np.abs(r1[..., k])), 1e-300))
                             for k in range(int(nv))],
               bitwise_identical=bool(np.array_equal(r1, r2)))
    (a.dest.parent / f'{a.dest.name}_compare.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
