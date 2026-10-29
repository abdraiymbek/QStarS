#!/usr/bin/env python3
"""Exact radiation response of selected columns by depth colouring (Report 21 preconditioner).

The fixed-J direct Jacobian M lacks the radiation response, and in the cool layer the fixed-J RE
row is nearly insensitive to T: the coupled T sensitivity comes almost entirely from the thick
radiation field (Report 20 §4.5). This measures the EXACT radiation part of chosen columns (T, n_e,
...) with one batched tangent replay pass. For column v and colour c in 0..p-1, the direction is
z = e_v at every depth d with d % p == c (fractional variables s). Its radiation response at row
depth k is attributed to the column depth of the same colour within +-(p-1)/2 of k. The offsets -1, 0
and +1 enter the block-tridiagonal preconditioner, and the response beyond +-1 is dropped. Coupling at
offsets that are multiples of p aliases onto offset 0; this is measured and reported (alias check,
second colouring), not hidden.

Only the PRECONDITIONER changes. Products, residual and convergence test stay exact, so the solution
cannot depend on this. The probed blocks may be reused at nearby states (stale preconditioner).

Output: <dir>/probe.npz with rlo, rdi, rup of shape (ND, NT, len(cols)) (radiation entries of the
probed columns at offsets -1/0/+1) plus column indices and colour period.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import BlockTridiag, read_jv_rad, read_trace_record  # noqa: E402

COLS = {'T': 133, 'ne': 132, 'HII': 31, 'Hm': 0, 'HI1': 1, 'HI2': 2}


def directions(cols, p, nt, nd):
    zs, tags = [], []
    for v in cols:
        for c in range(p):
            z = np.zeros((nt, nd)); z[COLS[v], c::p] = 1.0
            zs.append(z); tags.append((v, c))
    return zs, tags


def assemble(rad, tags, cols, p, nt, nd):
    """rad[:, :, i] response of direction i -> radiation entries of each probed column at offsets -1/0/+1."""
    rlo = np.zeros((nd, nt, len(cols))); rdi = np.zeros_like(rlo); rup = np.zeros_like(rlo)
    for i, (v, c) in enumerate(tags):
        j = cols.index(v)
        for k in range(nd):
            off = ((k - c) % p)                       # k - d for the nearest same-colour column depth d <= k
            if off > p // 2:
                off -= p                              # nearest column depth may lie deeper (d > k)
            d = k - off
            if not (0 <= d < nd):
                continue
            if off == 0:
                rdi[k, :, j] = rad[:, k, i]
            elif off == 1:                            # column one depth above (outer): lower block of row k
                rlo[k, :, j] = rad[:, k, i]
            elif off == -1:
                rup[k, :, j] = rad[:, k, i]
    return rlo, rdi, rup


def add_to(it, probe, pc_floor=0.0):
    """Preconditioner M + probed radiation columns (block tridiagonal), optional T floor as in the driver."""
    lo, di, up = it.Jdir.lo.copy(), it.Jdir.di.copy(), it.Jdir.up.copy()
    for j, v in enumerate(probe['cols']):
        c = COLS[str(v)]
        di[:, :, c] += probe['rdi'][:, :, j]
        lo[:, :, c] += probe['rlo'][:, :, j]
        up[:, :, c] += probe['rup'][:, :, j]
    if pc_floor > 0:
        for d in range(lo.shape[0]):
            rT = np.abs(lo[d][133]).sum() + np.abs(di[d][133]).sum() + np.abs(up[d][133]).sum()
            sgn = 1.0 if di[d][133, 133] >= 0 else -1.0
            di[d][133, 133] += sgn * pc_floor * rT
    return BlockTridiag(lo, di, up).factor()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path, help='NK run directory (driver --out)')
    ap.add_argument('base', help='accepted base pass name inside run')
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_27/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--cols', nargs='+', default=['T', 'ne'])
    ap.add_argument('--period', type=int, default=5)
    ap.add_argument('--vadat', nargs=2, action='append', default=[])
    ap.add_argument('--evaluate', type=Path, nargs='*', default=[],
                    help='iterate dirs holding jv*/pairs.npz at this base: report q = ||W(P^-1 J z - z)||/||W z||')
    ap.add_argument('--pc-t-floor', type=float, default=1e-3)
    ap.add_argument('--hom', type=Path, default=None, help='homotopy npz (shift mode) to build the Iterate as the driver')
    a = ap.parse_args()
    run = a.run.resolve(); base = run / a.base
    x0, _ = read_trace_record(run / 'it00_base' / 'QSTARS_RESID_TRACE.bin')
    ck.NT, ck.ND = x0.shape
    nt, nd = x0.shape
    pdir = run / f'probe_{a.base}'
    tagfile = pdir / 'probe.npz'
    if not tagfile.exists():
        R = ck.Runner(a.x0, a.exe, a.replace, False, a.threads, run / 'probe.log', [tuple(v) for v in a.vadat])
        xb, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
        s_cum = 1.0 - xb / x0
        zs, tags = directions(a.cols, a.period, nt, nd)
        R.log(f'probe {a.base}: {len(zs)} directions ({a.cols}, period {a.period})')
        R.run(pdir / 'pass', s_cum, replay=base / 'QSTARS_RJ.bin', dirs=np.stack(zs, axis=2), ba=False, store_j=False)
        _, rad = read_jv_rad(pdir / 'pass' / 'QSTARS_JV_RAD.bin')
        rlo, rdi, rup = assemble(rad, tags, a.cols, a.period, nt, nd)
        # far-field share: response at |offset| >= 2 (dropped) relative to the kept band, per column
        far = {}
        for j, v in enumerate(a.cols):
            kept = np.sum(rdi[:, :, j] ** 2) + np.sum(rlo[:, :, j] ** 2) + np.sum(rup[:, :, j] ** 2)
            tot = sum(np.sum(rad[:, :, i] ** 2) for i, t in enumerate(tags) if t[0] == v)
            far[v] = float(np.sqrt(max(tot - kept, 0.0) / max(tot, 1e-300)))
        np.savez(tagfile, rlo=rlo, rdi=rdi, rup=rup, cols=np.array(a.cols), period=a.period, base=str(base))
        (pdir / 'probe.json').write_text(json.dumps(dict(base=str(base), cols=a.cols, period=a.period,
                                                         far_field_share=far), indent=1))
        R.log(f'probe {a.base}: far-field share (|offset| >= 2, dropped) {far}')
    probe = dict(np.load(tagfile))
    probe['cols'] = [str(c) for c in probe['cols']]
    if not a.evaluate:
        return
    cols = {n: (run / f'cols_{a.base}' / f'col_{n}_p', run / f'cols_{a.base}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    hom = None
    if a.hom:
        z = np.load(a.hom)
        hom = dict(Tfix=z['Tfix'], c=z['c'], R0=z['R0'], mode='shift', lam=0.0)
    ck.ULTRA_FRAC = 1e-20; ck.ULTRA_WEIGHT = 1e-8
    res = {}
    for fl in (0.0, a.pc_t_floor):
        ck.PC_T_FLOOR = fl
        it = ck.Iterate(base, x0, cols, 1e-4, True, hom)
        Ms = {f'M floor {fl:g}': it.M, f'M+probe floor {fl:g}': add_to(it, probe, fl)}
        for name, M in Ms.items():
            qs = []
            for d in a.evaluate:
                for j in sorted(Path(d).glob('jv*')):
                    P = np.load(j / 'pairs.npz')
                    for i in range(P['z'].shape[2]):
                        zz, jz = P['z'][:, :, i], P['jz'][:, :, i]
                        e = (M.solve(jz) - zz) * it.W
                        qs.append(float(np.linalg.norm(e) / max(np.linalg.norm(zz * it.W), 1e-300)))
            res[name] = qs
            print(f'{name:22s} q = ' + ' '.join(f'{q:.3f}' for q in qs))
    (pdir / 'evaluation.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
