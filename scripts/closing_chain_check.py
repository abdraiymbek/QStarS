#!/usr/bin/env python3
"""V1 plan P3 (Report 26 review): refuse a closing chain whose physics differs from the converged solve's.

    closing_chain_check.py STATE --x0 X0 --exe EXE --replace REPLACE [--two-phot-j] [--vadat KEY VALUE ...] [--free-t]

STATE is the accepted base pass the chain starts from (scripts/post_convergence.sh BASE). Its run_manifest.json and
VADAT record what the solve ran with. The planned chain must use the same executable (sha256), the same equation table
(QSTARS_REPLACE_FORCE), the same two-photon treatment (stored J iff the state's pass had QSTARS_TWO_PHOT_J), FIX_T as
the solve (--free-t iff FIX_T=F), and every VADAT setting the solve changed from the restart X0's VADAT (other than
the per-pass QS_* linearization switches and TWO_METH, which --two-phot-j sets) with the same value. Exit 0 and print
the checks, or exit 2 with every mismatch."""
from __future__ import annotations

import argparse, hashlib, json, re, sys
from pathlib import Path

PASS_KEYS = {'QS_ALI_LIN', 'QS_LIN_FAC', 'TWO_METH', 'FIX_T'}     # set per pass by the tools themselves


def vadat_settings(path):
    """{KEY: value} of a CMFGEN VADAT file ('value [KEY]' lines)."""
    out = {}
    for ln in Path(path).read_text().splitlines():
        m = re.match(r'\s*(\S+)\s+\[([A-Za-z0-9_]+)\]', ln)
        if m:
            out[m.group(2)] = m.group(1)
    return out


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check(state, x0, exe, replace, two_phot_j, vadat, free_t):
    state = Path(state)
    man = json.loads((state / 'run_manifest.json').read_text())
    env = man.get('extra_env', {})
    bad, ok = [], []

    def cmp(name, want, got):
        (ok if want == got else bad).append(f'{name}: solve {want!r}, chain {got!r}')
    cmp('executable sha256', man.get('executable_sha256'), sha(exe))
    cmp('equation table', str(Path(env.get('QSTARS_REPLACE_FORCE', '')).resolve()), str(Path(replace).resolve()))
    cmp('two-photon stored J', bool(env.get('QSTARS_TWO_PHOT_J')), bool(two_phot_j))
    sv, xv = vadat_settings(state / 'VADAT'), vadat_settings(Path(x0) / 'VADAT')
    cmp('FIX_T', sv.get('FIX_T'), 'F' if free_t else 'T')
    planned = dict(vadat)
    for key, val in sorted(sv.items()):
        if key in PASS_KEYS or xv.get(key) == val:
            continue
        cmp(f'VADAT {key} (set by the solve)', val, planned.get(key))
    for key, val in sorted(planned.items()):
        if key not in sv:
            bad.append(f'VADAT {key}: not in the solve\'s VADAT, chain {val!r}')
        elif sv[key] != val and xv.get(key) == sv[key]:
            bad.append(f'VADAT {key}: solve {sv[key]!r} (restart default), chain {val!r}')
    return ok, bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('state', type=Path)
    ap.add_argument('--x0', type=Path, required=True)
    ap.add_argument('--exe', type=Path, required=True)
    ap.add_argument('--replace', type=Path, required=True)
    ap.add_argument('--two-phot-j', action='store_true')
    ap.add_argument('--free-t', action='store_true')
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'))
    a = ap.parse_args()
    ok, bad = check(a.state, a.x0, a.exe, a.replace, a.two_phot_j, [tuple(v) for v in a.vadat], a.free_t)
    for line in ok:
        print('ok      ' + line)
    for line in bad:
        print('MISMATCH ' + line)
    return 2 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
