#!/usr/bin/env python3
"""Council block: evaluate one finished iterate of results/nk167_trap/v1_p2 against the agreed rules. Read-only.

    block_eval.py IT [--radius R]     (IT: the iterate number, e.g. 20; R: the radius the iterate ran with)

rho = (1 - (phi_new/phi)^2) / (1 - pred^2) for the accepted trial (pred: the linear-model merit at that lambda, history
'predicted_rel_by_lam'); next radius: rho >= 0.7 keep, else max(0.125, lambda_acc x radius). Frozen it16 measure of the
base and the accepted trial (absolute region squares d33-47, d62-100), RE net/gross max, significant max|s|, wall
time from the driver log, frozen progress per hour. Appends to block_eval.json."""
import argparse, json, re, sys, time
from datetime import datetime
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402

RUN = ROOT / 'results/nk167_trap/v1_p2'
FROZEN = ROOT / 'Reports/evidence/v1_exec_20261106/frozen_norm_it16.npz'

ap = argparse.ArgumentParser(); ap.add_argument('it', type=int); ap.add_argument('--radius', type=float, required=True)
a = ap.parse_args()
rec = next(r for r in reversed(json.loads((RUN / 'history.json').read_text())) if r['it'] == a.it)
x0 = cn.configure()
z = np.load(FROZEN)
base = Path(rec['base'])
out = dict(it=a.it, radius=a.radius, base=base.name, status=rec.get('status'), accepted=rec.get('accepted'),
           merit=rec['merit'], trials=rec.get('trials'), predicted=rec.get('predicted_rel_by_lam'),
           max_step_significant=rec.get('max_step_significant'), lm_bound=rec.get('lm_bound'),
           lm_bound_linear=rec.get('lm_bound_linear'))
sb = cn.score(FROZEN, base, x0)
out['frozen_base'] = sb['frozen']
if rec.get('accepted'):
    lab = rec['accepted']; lam = float(lab[1:].split('_')[0]) / (16.0 if '_d16' in lab else 1.0)
    pred = rec['predicted_rel_by_lam'].get(f'l{float(lab[1:].split("_")[0]):.4g}') if rec.get('predicted_rel_by_lam') else None
    ratio = rec['merit_new'] / rec['merit']
    out.update(lam_accepted=lam, ratio=ratio)
    if pred is not None:
        if '_d16' in lab:
            pred = None
        else:
            out['rho'] = (1 - ratio ** 2) / max(1 - pred ** 2, 1e-300)
    rho = out.get('rho')
    out['next_radius'] = a.radius if (rho is not None and rho >= 0.7) else max(0.125, lam * a.radius)
    st = cn.score(FROZEN, RUN / f'it{a.it + 1:02d}_trial_{lab}', x0)
    out['frozen_new'] = st['frozen']
    out['frozen_drop'] = 1 - st['frozen'] / sb['frozen']
    w = z['w']
    for name, (d0, d1) in (('d33-47', (33, 47)), ('d62-100', (62, 100))):
        for tag_, sc in (('base', sb), ('new', st)):
            _, Rr = __import__('qstars.nlte_newton', fromlist=['read_trace_record']).read_trace_record(
                cn.residual_dir(base if tag_ == 'base' else RUN / f'it{a.it + 1:02d}_trial_{lab}') / 'QSTARS_RESID_TRACE.bin')
            out[f'abs2_{name}_{tag_}'] = float(((w * Rr[:134])[:, d0 - 1:d1] ** 2).sum())
    out['re_net_gross_base'], out['re_net_gross_new'] = sb['re_net_over_gross_max'], st['re_net_over_gross_max']
    out['conservation_new'] = st['conservation_max_rel']
# wall time of the iterate: from its first log line to the last trial merit line
t = [datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S') for ln in (RUN / 'driver.log').read_text().splitlines()
     if (m := re.match(r'(\S+ \S+) it%02d ' % a.it, ln))]
if t:
    out['wall_h'] = (t[-1] - t[0]).total_seconds() / 3600.0
    if 'frozen_drop' in out:
        out['progress_pct_per_h'] = 100.0 * out['frozen_drop'] / out['wall_h']
print(json.dumps(out, indent=1))
f = Path(__file__).with_name('block_eval.json')
prev = json.loads(f.read_text()) if f.exists() else []
f.write_text(json.dumps([p for p in prev if p['it'] != a.it] + [out], indent=1))
