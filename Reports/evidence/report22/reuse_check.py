"""Offline check of the same-base reuse identity on archived nk167_hot iterate 14: for every archived exact pair,
Jdir(delta).z + drad (drad read from the pass's QSTARS_JV_RAD.bin, as the reuse code does) must equal the stored jz;
with delta -> 0 the difference is exactly the PTC diagonal term delta*c*z_T."""
import sys, json, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad')
from nkoff import setup_globals, load_iterate, RUN
from qstars.nlte_newton import read_jv_rad
x0 = setup_globals()
it, rec = load_iterate(14, x0)
it0, _ = load_iterate(14, x0, delta=0.0)
itd = RUN / 'it14_from_it14_trial_l1'
worst = worst0 = 0.0; n = 0
for jvd in sorted(itd.glob('jv*')):
    P = np.load(jvd / 'pairs.npz'); _, drad = read_jv_rad(jvd / 'QSTARS_JV_RAD.bin')
    assert (jvd / 'QSTARS_STEP.bin').read_bytes() == np.asarray(it.s_cum, '<f8').tobytes(order='F')
    man = json.loads((jvd / 'run_manifest.json').read_text())['extra_env']
    assert man['QSTARS_JV_REPLAY'] == str(it.rj.resolve())
    for i in range(P['z'].shape[2]):
        z, jz = P['z'][:, :, i], P['jz'][:, :, i]
        r = it.Jdir.matvec(z) + drad[:, :, i]
        worst = max(worst, np.max(np.abs(r - jz)) / np.max(np.abs(jz)))
        r0 = it0.Jdir.matvec(z) + drad[:, :, i]
        expect = jz.copy(); expect[133] += it.hom['delta'] * it.ptc_c * z[133]
        worst0 = max(worst0, np.max(np.abs(r0 - expect)) / np.max(np.abs(jz)))
        n += 1
print(f'{n} archived pairs; reuse reconstruction vs stored jz (same delta): {worst:.2e}; '
      f'delta->0 reconstruction vs jz + delta c z_T: {worst0:.2e}')
