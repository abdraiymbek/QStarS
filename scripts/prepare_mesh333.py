#!/usr/bin/env python3
"""Report 21: 333-radius restart nesting the 167- and 84-radius cool grids (spatial resolution check).

Same physical inputs as the 167 chain (benchmarks/cmfgen_cool_validation_plan.json, case T9, same
source model and atomic data): backends/cmfgen/prepare_cool.py with mesh factor 4 (every 84-radius
log-r interval split in four, so 333[::2] = 167 and 333[::4] = 84), then ONE native iteration
(the bootstrap stage of backends/cmfgen/recover_cool.py). The restart's own populations only
initialize the files: the Newton-Krylov passes set the state through the perturbation hook.
The numerical controls (VADAT, IN_ITS) are then copied from the 167-radius restart
results/cmfgen_ali_lambda167, so that 167 and 333 passes differ only in ND/NP and the radius grid.
Nesting of radius and velocity is verified against the 167 and 84 restarts (relative 1e-12).
"""
import argparse, json, shutil, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'backends' / 'cmfgen')); sys.path.insert(0, str(ROOT / 'src'))
from prepare_cool import prepare  # noqa: E402
from run_component import run  # noqa: E402
from qstars.cmfgen_backend import read_rvtj, sha256  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', type=Path, default=ROOT / 'sources/cache/cmfgen_models_2023/sample_grid/zpup_like_model')
    ap.add_argument('--atomic', type=Path, default=ROOT / 'sources/cache/cmfgen_atomic_2023/atomic')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_28/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--ref167', type=Path, default=ROOT / 'results/cmfgen_ali_lambda167')
    ap.add_argument('--ref84', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--boot', type=Path, default=ROOT / 'results/cmfgen_recovery_T9_mesh333_bootstrap')
    ap.add_argument('--out', type=Path, default=ROOT / 'results/cmfgen_ali_lambda333')
    ap.add_argument('--threads', type=int, default=2)
    a = ap.parse_args()
    if not a.boot.exists():
        prepare(a.source, a.boot, a.atomic, case='T9', iterations=1, factor=4)
    if not (a.boot / 'run_result.json').exists():
        r = run(a.boot, a.exe, 7200, a.threads)
        print(json.dumps(dict(status=r['status'], error=r['error'], seconds=r['seconds'])), flush=True)
        if not r['status'].startswith('completed'):
            raise SystemExit('bootstrap failed (preserved)')
    if a.out.exists():
        raise SystemExit(f'{a.out} exists')
    # the bootstrap iteration's QSTARS_*TRACE* outputs describe the pre-bootstrap state: not part of the restart
    shutil.copytree(a.boot, a.out, symlinks=True, copy_function=shutil.copy2,
                    ignore=shutil.ignore_patterns('QSTARS_*TRACE*'))
    for name in ('VADAT', 'IN_ITS'):
        shutil.copy2(a.ref167 / name, a.out / name)
    for name in ('run_manifest.json', 'run_result.json'):
        (a.out / name).rename(a.out / f'bootstrap_{name}')
    r333 = read_rvtj(a.out / 'RVTJ'); r167 = read_rvtj(a.ref167 / 'RVTJ'); r84 = read_rvtj(a.ref84 / 'RVTJ')
    chk = {}
    for key in ('radius', 'velocity'):
        x3 = np.asarray(r333[key]); x1 = np.asarray(r167[key]); x8 = np.asarray(r84[key])
        chk[key] = dict(nd=len(x3), vs167=float(np.max(np.abs(x3[::2] / x1 - 1))),
                        vs84=float(np.max(np.abs(x3[::4] / x8 - 1))))
    ok = all(v['vs167'] < 1e-12 and v['vs84'] < 1e-12 and v['nd'] == 333 for v in chk.values())
    (a.out / 'mesh333_control.json').write_text(json.dumps(dict(
        purpose='333-radius restart nesting the 167/84 cool grids; Report 21 spatial resolution check',
        bootstrap=str(a.boot), controls_from=str(a.ref167), physical_inputs_changed=False,
        plan_sha256=sha256(ROOT / 'benchmarks/cmfgen_cool_validation_plan.json'),
        vadat_sha256=sha256(a.out / 'VADAT'), nesting=chk, nested=ok), indent=1) + '\n')
    print(json.dumps(chk, indent=1), 'nested' if ok else 'NOT NESTED')
    if not ok:
        raise SystemExit('grids do not nest')


if __name__ == '__main__':
    main()
