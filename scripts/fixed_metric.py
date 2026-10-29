#!/usr/bin/env python3
"""Report 23: progress of a cool_newton_krylov.py run in ONE fixed measure (read-only).

    fixed_metric.py OUT REF_BASE [--ali] [--extra DIR ...]

The driver's merit and natural measure use each base's own preconditioner (and the natural measure its own Krylov
subspace), so they are not comparable between iterates. Here every accepted base's residual (its <base>_tp replay when
present) is measured as ||W_ref M_ref^-1 R_k|| with the preconditioner and weights of REF_BASE (local-ALI increment with
--ali), overall and by species group (RE row: without the PTC term). --extra: also measure these pass directories
(names under OUT, e.g. line-search trials), as they are (no _tp replay)."""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out', type=Path); ap.add_argument('ref'); ap.add_argument('--ali', action='store_true')
    ap.add_argument('--extra', nargs='+', default=[])
    a = ap.parse_args(); out = a.out
    x0, _ = read_trace_record(out / 'it00_base' / 'QSTARS_RESID_TRACE.bin'); ck.NT, ck.ND = x0.shape
    ck.SIG_FRAC, ck.TRACE_WEIGHT, ck.ULTRA_FRAC, ck.ULTRA_WEIGHT, ck.WEIGHT_RAMP = 1e-10, 1e-3, 1e-20, 0.0, 2.0
    ck.MERIT_MODE, ck.EQ_FLOOR, ck.PC_T_FLOOR = 'precond', 0.0, 1e-3
    ref = out / a.ref
    cols = {n: (out / f'cols_{ref.name}' / f'col_{n}_p', out / f'cols_{ref.name}' / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    resdir = out / (ref.name + '_tp') if (out / (ref.name + '_tp')).exists() else ref
    it = ck.Iterate(resdir, x0, cols, 1e-4, True, None, rj=ref / 'QSTARS_RJ.bin')
    if a.ali:
        it.install_ali_increment(out / (ref.name + '_ali'))
    it.extend_significance(0.0, rows=[3])
    hist = json.loads((out / 'history.json').read_text())
    bases = []
    for r in hist:
        for b in (Path(r['base']).name, f"it{r['it'] + 1:02d}_trial_{r['accepted']}" if r.get('accepted') and
                  r.get('merit_new', np.inf) < r.get('merit', 0) else None):
            if b and b not in bases and (out / b / 'QSTARS_RESID_TRACE.bin').exists():
                bases.append(b)
    groups = [('H I', 1, 31), ('H II', 31, 32), ('He', 32, 132), ('n_e', 132, 133), ('RE', 133, 134)]
    print(f'{"base":24s} {"||W M_ref^-1 R||":>16s} ' + ' '.join(f'{g:>9s}' for g, _, _ in groups) + '   max|M_ref^-1 R| sig')
    for b in bases + [e for e in a.extra if e not in bases]:
        d = out / (b + '_tp') if b in bases and (out / (b + '_tp') / 'QSTARS_RESID_TRACE.bin').exists() else out / b
        if not (d / 'QSTARS_RESID_TRACE.bin').exists():
            print(f'{b:24s} (no residual)'); continue
        _, R = read_trace_record(d / 'QSTARS_RESID_TRACE.bin')
        y = it.M.solve(R); w = it.W * y * it.active
        print(f'{b:24s} {np.linalg.norm(w):16.4e} ' + ' '.join(f'{np.linalg.norm(w[r0:r1]):9.2e}' for _, r0, r1 in groups)
              + f'   {np.max(np.abs(y[it.sig])):.3e}')


if __name__ == '__main__':
    main()
