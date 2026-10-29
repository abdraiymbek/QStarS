#!/usr/bin/env python3
"""Report 22: deterministic check that QSTARS_NO_OBS_FORMAL changes nothing the solver reads. Two single-thread
replay passes (frozen stored J of a base pass, so no ALI iteration and no thread-order effects) with the BA dump,
observer-frame formal solution on and off; residual, BA and nothing else compared bit for bit.

    check_noobs_ba.py BASE_PASS OUTDIR --x0 X0 --exe EXE
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_ba_blocks, read_fortran, read_trace_record  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('base', type=Path); ap.add_argument('out', type=Path)
ap.add_argument('--x0', type=Path, required=True); ap.add_argument('--exe', type=Path, required=True)
ap.add_argument('--base-mode', action='store_true',
                help='full base passes (ALI-iterated J, warm start as archived) instead of frozen-J replay passes')
a = ap.parse_args()
base = a.base.resolve(); out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
man = json.loads((base / 'run_manifest.json').read_text())['extra_env']
vad = [(p[-1][1:-1], p[0]) for p in map(str.split, (base / 'VADAT').read_text().splitlines())
       if len(p) >= 2 and p[-1] in ('[INC_AD]', '[INC_ADV]')]
x, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
ck.NT, ck.ND = x.shape
s_cum = read_fortran(base / 'QSTARS_STEP.bin', x.shape)
res = {}
for tag, no_obs in (('obs_on', False), ('obs_off', True)):
    R = ck.Runner(a.x0, a.exe, Path(man['QSTARS_REPLACE_FORCE']), False, 1, out / 'check.log', vad)
    R.slim = False; R.no_obs = no_obs
    if a.base_mode:
        R.run(out / tag, s_cum, ba=True, store_j=True, threads=1, warm=man.get('QSTARS_ALI_WARM'))
    else:
        R.run(out / tag, s_cum, replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False, threads=1)
xa, Ra = read_trace_record(out / 'obs_on' / 'QSTARS_RESID_TRACE.bin'); xb, Rb = read_trace_record(out / 'obs_off' / 'QSTARS_RESID_TRACE.bin')
nt, nd = xa.shape
A = read_ba_blocks(out / 'obs_on' / 'QSTARS_BA_TRACE.bin', nt, nd); B = read_ba_blocks(out / 'obs_off' / 'QSTARS_BA_TRACE.bin', nt, nd)
res = dict(state_identical=bool(np.array_equal(xa, xb)), residual_bitwise_identical=bool(np.array_equal(Ra, Rb)),
           ba_bitwise_identical=all(np.array_equal(p, q) for p, q in zip(A, B)),
           residual_max_abs_diff=float(np.max(np.abs(Ra - Rb))),
           ba_max_abs_diff=float(max(np.max(np.abs(p - q)) for p, q in zip(A, B))))
if a.base_mode:
    ja, jb = (out / 'obs_on' / 'QSTARS_RJ.bin').read_bytes(), (out / 'obs_off' / 'QSTARS_RJ.bin').read_bytes()
    res['stored_J_bitwise_identical'] = ja == jb
    Ja = np.frombuffer(ja[12:], np.uint8).reshape(-1, 12 + 8 * nd)[:, 12:].copy().view('<f8')
    Jb = np.frombuffer(jb[12:], np.uint8).reshape(-1, 12 + 8 * nd)[:, 12:].copy().view('<f8')
    res['stored_J_max_rel_diff'] = float(np.max(np.abs(Ja - Jb) / np.maximum(np.abs(Ja), 1e-300)))
(out / 'check.json').write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
