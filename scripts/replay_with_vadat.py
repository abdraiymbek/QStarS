#!/usr/bin/env python3
"""Report 23: one frozen-J replay pass of a base state with extra VADAT settings (diagnostic), compared with the base.

    replay_with_vadat.py BASE OUT --x0 X0 --exe EXE [--vadat KEY VALUE ...] [--lin-fac F] [--depths D0 D1] [--two-phot-j RJ]
Reports the residual rows of interest (RE row NT and the H I n=2 row 3) at the chosen depths for the base and the
replay; the replay shares the base's state and stored J, so differences come only from the VADAT settings.
--two-phot-j RJ (build >= 32): keep the two-photon method USE_RAD with the conjugate-frequency J from the stored-J
file RJ (QSTARS_TWO_PHOT_J; 'base' = the base's own QSTARS_RJ.bin) instead of the OLD_DEFAULT substitution.
--lin-fac F: QS_LIN_FAC of the pass (default 0, the fixed-J blocks; 1 = the local-ALI blocks, fast with build >= 33).
The runner takes QS_LIN_FAC only through Runner.run(lin_fac=...); a '--vadat QS_LIN_FAC F' pair is read as --lin-fac F.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_fortran, read_trace_record  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('base', type=Path); ap.add_argument('out', type=Path)
ap.add_argument('--x0', type=Path, required=True); ap.add_argument('--exe', type=Path, required=True)
ap.add_argument('--vadat', nargs=2, action='append', default=[]); ap.add_argument('--depths', type=int, nargs=2, default=[30, 41])
ap.add_argument('--threads', type=int, default=1); ap.add_argument('--two-phot-j', default=None)
ap.add_argument('--lin-fac', type=float, default=None)
a = ap.parse_args()
lfs = {float(v.replace('D', 'E').replace('d', 'e')) for k, v in a.vadat if k == 'QS_LIN_FAC'}
if a.lin_fac is not None:
    lfs.add(a.lin_fac)
if len(lfs) > 1:
    ap.error(f'conflicting QS_LIN_FAC values {sorted(lfs)}')
lin_fac = lfs.pop() if lfs else 0.0
base = a.base.resolve(); out = a.out.resolve(); out.parent.mkdir(parents=True, exist_ok=True)
man = json.loads((base / 'run_manifest.json').read_text())['extra_env']
vad = [(p[-1][1:-1], p[0]) for p in map(str.split, (base / 'VADAT').read_text().splitlines())
       if len(p) >= 2 and p[-1] in ('[INC_AD]', '[INC_ADV]')] + [tuple(v) for v in a.vadat if v[0] != 'QS_LIN_FAC']
xb, Rb = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
ck.NT, ck.ND = xb.shape
R = ck.Runner(a.x0, a.exe, Path(man['QSTARS_REPLACE_FORCE']), False, a.threads, out.parent / (out.name + '.log'), vad)
R.slim = False
tpj = None if a.two_phot_j is None else (base / 'QSTARS_RJ.bin' if a.two_phot_j == 'base' else Path(a.two_phot_j))
R.run(out, read_fortran(base / 'QSTARS_STEP.bin', xb.shape), replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False, threads=a.threads,
      two_phot_j=tpj, lin_fac=lin_fac)
xr, Rr = read_trace_record(out / 'QSTARS_RESID_TRACE.bin')
nt = xb.shape[0]; d0, d1 = a.depths
pth = out / 'QSTARS_ALI_LIN_PATH'            # build >= 33, QS_LIN_FAC /= 0: which local-ALI path ran
ali_path = pth.read_text().strip() if pth.is_file() else None
rows = dict(state_identical=bool(np.array_equal(xb, xr)), vadat=a.vadat, lin_fac=lin_fac, ali_lin_path=ali_path,
            two_phot_j=None if tpj is None else str(tpj),
            depths=list(range(d0, d1 + 1)),
            RE_base=Rb[nt - 1, d0 - 1:d1].tolist(), RE_replay=Rr[nt - 1, d0 - 1:d1].tolist(),
            n2_base=Rb[2, d0 - 1:d1].tolist(), n2_replay=Rr[2, d0 - 1:d1].tolist())
(out.parent / (out.name + '.json')).write_text(json.dumps(rows, indent=1))
print('two-photon line in OUTGEN:', [ln.strip() for ln in (out / 'OUTGEN').read_text(errors='replace').splitlines() if 'two-photon' in ln.lower() or 'two photon' in ln.lower()][-2:])
print(f'QS_LIN_FAC {lin_fac:g}; local-ALI path: {ali_path}')
print('depth   RE base   RE replay   ratio    n2 base  n2 replay')
for i, d in enumerate(rows['depths']):
    print(f'{d:5d} {rows["RE_base"][i]:9.3f} {rows["RE_replay"][i]:10.3f} {rows["RE_replay"][i]/rows["RE_base"][i]:7.3f} {rows["n2_base"][i]:9.1f} {rows["n2_replay"][i]:9.1f}')
