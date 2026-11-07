#!/usr/bin/env python3
"""V1 plan P2: the least squares over the v1_p2 Krylov span with bounded coefficient amplification. The directions are
normalized and the span is truncated at relative singular value tau (state space, SVD of [z_i / |z_i|]); the products of
the orthonormal basis are J U = J Z_n V S^-1, whose error is at most the products' own error / tau (the full span had
coefficients up to 4e14 cancelling to |s| = 0.4: the model was product noise). For each tau: bounded LM-QR and box
steps, their predicted fit, and the amplification ||y * |z_i| || / ||s||. With --verify, one exact product pass of the
candidate steps gives each step's TRUE linear fit (fresh product). Read-only apart from that pass (RUN/tsvd_verify).

    p2_tsvd.py [--passes N] [--taus 1e-2 1e-4 1e-6 1e-8] [--verify]"""
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'p1'))
import campaign_norms as cn  # noqa: E402
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_jv_rad  # noqa: E402
from p1_native import EXE, REPLACE, X0DIR  # noqa: E402

RUN, NAME, BASE = ROOT / 'results/nk167_trap/v1_p2', 'it16_from_it16_trial_l0.5', ROOT / 'results/nk167_trap/free_r25/it16_trial_l0.5'


def load(npass):
    x0b = cn.configure(); ck.STRUCT_ACTIVE = True
    itdir = RUN / NAME
    cols = {n: (itdir / f'col_{n}_p', itdir / f'col_{n}_m') for n in ('ne', 'hii', 'T')}
    it = ck.Iterate(cn.residual_dir(BASE), x0b, cols, 1e-4, True, None, rj=BASE / 'QSTARS_RJ.bin')
    it.extend_significance(0.0, rows=[3])
    passes = sorted({p.name[:4] for p in itdir.glob('jv[0-9][0-9]p[0-9]') if (p / 'pairs.npz').exists()})
    if npass:
        passes = passes[:npass]
    S, JS = [], []
    for pn in passes:
        for q in sorted(itdir.glob(f'{pn}p[0-9]')):
            z = np.load(q / 'pairs.npz')
            S += [z['z'][:, :, i] for i in range(z['z'].shape[2])]; JS += [z['jz'][:, :, i] for i in range(z['z'].shape[2])]
    return it, passes, S, JS


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--passes', type=int, default=0)
    ap.add_argument('--taus', type=float, nargs='+', default=[1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8])
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    it, passes, S, JS = load(a.passes)
    b = it.obj(it.R).ravel(); beta = float(np.linalg.norm(b))
    rad = np.full(it.x.shape, 0.5); rad[133] = 0.03
    Zm = np.stack([z.ravel() for z in S], 1); nz = np.linalg.norm(Zm, axis=0)
    Jm = np.stack([j.ravel() for j in JS], 1)
    U, sv, Vt = np.linalg.svd(Zm / nz, full_matrices=False)
    print(f'{len(passes)} passes, {len(S)} directions; singular values of the normalized span: max {sv[0]:.3e}, '
          f'min {sv[-1]:.3e}; count above 1e-2/1e-4/1e-6/1e-8: ' + '/'.join(str(int(np.sum(sv > t * sv[0]))) for t in (1e-2, 1e-4, 1e-6, 1e-8)))
    cands, rows = [], []
    for tau in a.taus:
        k = int(np.sum(sv > tau * sv[0]))
        M = (Vt[:k].T / sv[:k]) / nz[:, None]           # y = M c: coefficients of the orthonormal basis U_k
        Sb = [U[:, i].reshape(it.x.shape) for i in range(k)]
        Jb = Jm @ M
        Ab = np.stack([it.obj(Jb[:, i].reshape(it.x.shape)).ravel() for i in range(k)], 1)
        for sv_ in ('qr', 'box'):
            s, _, c = ck.bounded_step((Sb, Ab, b), it.active, rad, it.active, sv_)
            y = M @ c
            r = dict(tau=tau, basis=k, solver=sv_, fit=float(np.linalg.norm(Ab @ c - b) / beta),
                     amplification=float(np.linalg.norm(y * nz) / max(np.linalg.norm(s), 1e-300)),
                     max_s=float(np.max(np.abs(s[it.active]))))
            rows.append(r); cands.append(s)
            print(f'  tau {tau:g} ({k} basis vectors) {sv_}: predicted fit {r["fit"]:.4f}, amplification {r["amplification"]:.2e}', flush=True)
    out = Path(__file__).with_name(f'tsvd_{len(passes)}passes.json')
    np.save(Path(__file__).with_name(f'tsvd_{len(passes)}passes_steps.npy'), np.stack(cands, 2))
    if a.verify:
        R = ck.Runner(X0DIR, EXE, REPLACE, False, 4, RUN / 'tsvd_verify.log', [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')])
        R.no_obs = True
        d = RUN / f'tsvd_verify_{len(passes)}passes'
        R.run(d, it.s_cum, replay=it.rj, dirs=np.stack(cands, 2), ba=False, store_j=False, two_phot_j=BASE / 'QSTARS_RJ.bin', threads=4)
        _, drad = read_jv_rad(d / 'QSTARS_JV_RAD.bin')
        for i, (r, s) in enumerate(zip(rows, cands)):
            Js = it.Jdir.matvec(s) + drad[:, :, i]
            r['true_fit'] = float(np.linalg.norm((b.reshape(it.x.shape) - it.obj(Js))[:134]) / beta)
            print(f'  tau {r["tau"]:g} {r["solver"]}: predicted {r["fit"]:.4f}, TRUE (fresh product) {r["true_fit"]:.4f}', flush=True)
    out.write_text(json.dumps(dict(passes=len(passes), directions=len(S), base_merit=beta, rows=rows), indent=1))


if __name__ == '__main__':
    sys.exit(main())
