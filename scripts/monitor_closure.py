#!/usr/bin/env python3
"""Report 23: per-iterate monitor of a cool_newton_krylov.py run (read-only).

    monitor_closure.py OUT [--x0 X0]

For every history record: base, merit / natural measure, accepted trial and its merit, largest significant step,
and for the base state: normalized luminosity check L(r)/L(ND) from the photosphere inward (min, max; CMFGEN's
warning band is 0.8-1.25, its stop band 0.5-2), T and H ionization fraction at depths 25/35/45/60, and the H I ground
departure coefficient b1 at those depths."""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from qstars.nlte_newton import read_trace_record  # noqa: E402

kB, h, me, eV = 1.380649e-16, 6.62607015e-27, 9.1093837e-28, 1.602176634e-12


def lum(d):
    p = d / 'OBSFLUX'
    if not p.exists():
        return None
    t = p.read_text().splitlines()
    i = [k for k, l in enumerate(t) if 'Normalized luminosity check' in l]
    if not i:
        return None
    v = []
    for l in t[i[0] + 2:i[0] + 40]:
        if not l.strip():
            if v:
                break
            continue
        v += [float(x) for x in l.split()]
    return np.array(v)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out', type=Path); a = ap.parse_args()
    out = a.out
    hist = json.loads((out / 'history.json').read_text()) if (out / 'history.json').exists() else []
    print(f'{"it":>3} {"base":24s} {"merit":>10s} {"accepted":>10s} {"new":>10s} {"max|s|sig":>9s} {"L/L(ND) min-max":>16s} '
          'T(kK) d25/35/45/60          xHII d25/35/45/60              b1 d25/35/45/60')
    for r in hist:
        b = Path(r['base']); bd = b if b.is_absolute() else ROOT / b
        L = lum(bd)
        x, _ = read_trace_record(bd / 'QSTARS_RESID_TRACE.bin')
        T = x[133] * 1e4; ne = x[132]; npr = x[31]; n1 = x[1]; xh = npr / x[1:32].sum(0)
        b1 = n1 / (ne * npr * (h * h / (2 * np.pi * me * kB * T)) ** 1.5 * np.exp(13.598434 * eV / (kB * T)))
        ds = [24, 34, 44, 59]
        print(f'{r["it"]:3d} {bd.name:24s} {r["merit"]:10.3e} {str(r.get("accepted")):>10s} {r.get("merit_new", float("nan")):10.3e} '
              f'{r.get("max_step_significant", float("nan")):9.2e} '
              f'{"%.2f-%.2f" % (L.min(), L.max()) if L is not None else "n/a":>16s} '
              + ' '.join(f'{T[d] / 1e3:5.2f}' for d in ds) + '   ' + ' '.join(f'{xh[d]:7.1e}' for d in ds)
              + '   ' + ' '.join(f'{b1[d]:8.1e}' for d in ds))


if __name__ == '__main__':
    main()
