#!/usr/bin/env python3
"""Audit the build-18 exact directional Jacobian against the Report 18 oracle.

Three layers, all on the actual 84-radius cool state x0 and in native depth order
(index 0 = outer boundary):

  spectra   dCHI, dETA, dCHI_SCAT and dJ per frequency from QSTARS_JV_SPECTRAL.bin
            versus finite differences of the oracle QSTARS_SPECTRAL.bin exports
            (perturbed run x0*(1-eps*s) minus R0, divided by -eps).
  residual  J.(x0*s) = BA_direct.s + dSTEQ_rad  versus the exact residual change
            (R(x0*(1-eps*s)) - R0)/(-eps), using the same equation selection.

The oracle runs are never rerun or edited here; this reads their archived output.
"""
import argparse, json, struct, sys
from pathlib import Path

import numpy as np

GROUPS = {'H0': (1, 1), 'HI': (2, 31), 'HII': (32, 32), 'HeI': (33, 101),
          'HeII': (102, 131), 'HeIII': (132, 132), 'ne': (133, 133), 'T': (134, 134)}


def read_native_spectral(path):
    """QSCMF001 in native depth order. Returns nu, fqw and (nf, 6, nd) arrays."""
    path = Path(path)
    with path.open('rb') as f:
        if f.read(8) != b'QSCMF001':
            raise ValueError(f'{path}: not a QSCMF001 stream')
        nd, nf = struct.unpack('<ii', f.read(8))
    a = np.memmap(path, dtype='<f8', mode='r', offset=16)
    if a.size != 5 * nd + nf * (2 + 6 * nd):
        raise ValueError(f'{path}: truncated')
    data = a[5 * nd:].reshape(nf, 2 + 6 * nd)
    return data[:, 0] * 1e-15, data[:, 1], data[:, 2:].reshape(nf, 6, nd)


def read_jv_spectral(path):
    """QSJV0001: per frequency fl,fqw,rj(nd),replay(nd) then nv*[dchi,deta,dscat,dj](nd)."""
    path = Path(path)
    with path.open('rb') as f:
        if f.read(8) != b'QSJV0001':
            raise ValueError(f'{path}: not a QSJV0001 stream')
        nv, nd = struct.unpack('<ii', f.read(8))
    rec = 2 + 2 * nd + nv * 4 * nd
    a = np.memmap(path, dtype='<f8', mode='r', offset=16)
    if a.size % rec:
        raise ValueError(f'{path}: partial record')
    data = a.reshape(-1, rec)
    fl, fqw = data[:, 0], data[:, 1]
    rj = data[:, 2:2 + nd]
    replay = data[:, 2 + nd:2 + 2 * nd]
    d = data[:, 2 + 2 * nd:].reshape(-1, nv, 4, nd)
    return fl, fqw, rj, replay, d


def first_record(path):
    data = Path(path).read_bytes()
    cnt, nt, nd = struct.unpack_from('<3i', data, 0)
    n = nt * nd
    pops = np.frombuffer(data, '<f8', n, 12).reshape((nt, nd), order='F').copy()
    steq = np.frombuffer(data, '<f8', n, 12 + 8 * n).reshape((nt, nd), order='F').copy()
    return pops, steq


def read_jv_rad(path):
    data = Path(path).read_bytes()
    nv, nt, nd = struct.unpack_from('<3i', data, 0)
    n = nt * nd * nv
    s = np.frombuffer(data, '<f8', n, 12).reshape((nt, nd, nv), order='F')
    d = np.frombuffer(data, '<f8', n, 12 + 8 * n).reshape((nt, nd, nv), order='F')
    return s.copy(), d.copy()


def read_ba(path, nt, nd):
    """Full blocks keyed by (depth, band); depth/band are 1-based native indices."""
    data = Path(path).read_bytes()
    rec = 12 + 8 * nt * nt
    if len(data) % rec:
        raise ValueError(f'{path}: partial BA record')
    blocks = {}
    for off in range(0, len(data), rec):
        k, b, n = struct.unpack_from('<3i', data, off)
        if n != nt:
            raise ValueError('BA dimension mismatch')
        blocks[(k, b)] = np.frombuffer(data, '<f8', nt * nt, off + 12).reshape((nt, nt), order='F').copy()
    return blocks


def ba_apply(blocks, s, nd, diag=2):
    """y(:,K) = sum_band C(K,band) s(:,K+band-diag); blocks already scaled by POPS."""
    nt = s.shape[0]
    y = np.zeros_like(s)
    for (k, b), c in blocks.items():
        l = k + b - diag
        if 1 <= l <= nd:
            y[:, k - 1] += c @ s[:, l - 1]
    return y


def rel(err, ref):
    den = np.max(np.abs(ref))
    return float(np.max(np.abs(err)) / den) if den > 0 else float('nan')


def l2rel(err, ref):
    den = np.sqrt(np.sum(ref * ref))
    return float(np.sqrt(np.sum(err * err)) / den) if den > 0 else float('nan')


def audit_spectra(jv_dir, oracle, names, eps_table, out):
    fl, fqw, rj, replay, d = read_jv_spectral(jv_dir / 'QSTARS_JV_SPECTRAL.bin')
    nu0, w0, s0 = read_native_spectral(oracle / 'R0' / 'QSTARS_SPECTRAL.bin')
    nf = len(fl)
    if nf != len(nu0) or np.max(np.abs(fl / nu0 - 1)) > 1e-12:
        raise SystemExit('frequency grids differ between jv run and oracle')
    res = {'nfreq': nf,
           'primal_replay_vs_native_J': float(np.max(np.abs(replay - rj) / np.maximum(np.abs(rj), 1e-300))),
           'jv_J_vs_oracle_R0_J_rel': float(np.max(np.abs(rj - s0[:, 3, :]) / np.maximum(np.abs(s0[:, 3, :]), 1e-300)))}
    for k, name in enumerate(names):
        for lev in ('e3', 'e2'):
            run = oracle / f'{name}_{lev}'
            if not (run / 'QSTARS_SPECTRAL.bin').exists():
                continue
            eps = eps_table[f'{name}_{lev}']
            _, _, sp = read_native_spectral(run / 'QSTARS_SPECTRAL.bin')
            # Exported CHI/CHI_SCAT/ETA are scaled by 1e-10; native derivatives are not.
            fd = (np.asarray(sp, dtype=float) - np.asarray(s0, dtype=float)) / (-eps)
            fd_chi, fd_sca, fd_eta, fd_j = fd[:, 0] * 1e10, fd[:, 1] * 1e10, fd[:, 2] * 1e10, fd[:, 3]
            an_chi, an_eta, an_sca, an_j = d[:, k, 0], d[:, k, 1], d[:, k, 2], d[:, k, 3]
            # Per-frequency response-norm errors, then global summaries.
            def perfreq(an, fdv):
                num = np.max(np.abs(an - fdv), axis=1)
                den = np.max(np.abs(fdv), axis=1)
                return num / np.where(den > 0, den, np.nan), den
            out_k = {}
            for lab, an, fdv in (('chi', an_chi, fd_chi), ('eta', an_eta, fd_eta), ('scat', an_sca, fd_sca), ('J', an_j, fd_j)):
                pf, den = perfreq(an, fdv)
                # Weight frequencies by their FD response so negligible responses do not dominate.
                big = den >= 1e-6 * np.nanmax(den)
                out_k[lab] = dict(global_l2=l2rel(an - fdv, fdv),
                                  global_max=rel(an - fdv, fdv),
                                  median_perfreq=float(np.nanmedian(pf[big])),
                                  p95_perfreq=float(np.nanpercentile(pf[big], 95)),
                                  worst_freqs=[(float(fl[i]), float(pf[i])) for i in np.argsort(np.nan_to_num(pf * big, nan=-1))[::-1][:5]])
                # Relative J error weighted by local J (matters for rates).
            jrel = np.abs(an_j - fd_j) / np.maximum(np.abs(rj), 1e-300)
            out_k['J_err_over_J_max'] = float(np.max(jrel))
            out_k['eps'] = eps
            res[f'{name}_{lev}'] = out_k
            print(f'{name}_{lev}: ' + ', '.join(f'{lab} L2={out_k[lab]["global_l2"]:.2e} max={out_k[lab]["global_max"]:.2e}'
                                               for lab in ('chi', 'eta', 'scat', 'J')))
    (out / 'spectra_audit.json').write_text(json.dumps(res, indent=1))
    return res


def audit_residual(jv_dir, oracle, names, eps_table, out, r0_dir=None):
    pops0, R0 = first_record((r0_dir or oracle / 'R0') / 'QSTARS_RESID_TRACE.bin')
    s, drad = read_jv_rad(jv_dir / 'QSTARS_JV_RAD.bin')
    nt, nd, nv = s.shape
    blocks = read_ba(jv_dir / 'QSTARS_BA_TRACE.bin', nt, nd)
    pj, Rj = first_record(jv_dir / 'QSTARS_RESID_TRACE.bin')
    res = {'state_rel_diff_jv_vs_oracle': float(np.max(np.abs(pj - pops0) / np.abs(pops0))),
           'R_jv_vs_R0_l2rel': l2rel(Rj - R0, R0)}
    print(f"state diff {res['state_rel_diff_jv_vs_oracle']:.2e}  ||R_jv-R0||/||R0|| {res['R_jv_vs_R0_l2rel']:.2e}")
    mask = np.broadcast_to(~np.all(R0 == 0.0, axis=1)[:, None], R0.shape)
    for k, name in enumerate(names):
        direct = ba_apply(blocks, s[:, :, k], nd)
        jv = direct + drad[:, :, k]
        for lev in ('e3', 'e2'):
            f = oracle / f'{name}_{lev}' / 'QSTARS_RESID_TRACE.bin'
            if not f.exists():
                continue
            eps = eps_table[f'{name}_{lev}']
            _, Re = first_record(f)
            fd = (Re - R0) / (-eps)            # exact action of J on x0*s (to O(eps))
            m = mask
            e_full = l2rel((jv - fd)[m], fd[m])
            e_direct = l2rel((direct - fd)[m], fd[m])
            cos = float(np.dot(jv[m], fd[m]) / (np.linalg.norm(jv[m]) * np.linalg.norm(fd[m])))
            # Report-18 metric: predicted merit change along the native step.
            pred = R0 - eps * jv
            merit_pred = float(np.linalg.norm(pred[m]) / np.linalg.norm(R0[m]))
            merit_meas = float(np.linalg.norm(Re[m]) / np.linalg.norm(R0[m]))
            groups = {}
            for g, (a, b) in GROUPS.items():
                mm = np.zeros_like(m)
                mm[a - 1:b] = m[a - 1:b]
                if np.any(fd[mm] != 0):
                    groups[g] = l2rel((jv - fd)[mm], fd[mm])
            depth = []
            for dd in range(nd):
                md = m[:, dd]
                depth.append(l2rel((jv - fd)[md, dd], fd[md, dd]) if np.any(fd[md, dd] != 0) else float('nan'))
            depth = np.array(depth)
            worst = np.argsort(np.nan_to_num(depth, nan=-1))[::-1][:6]
            r = dict(eps=eps, e_full=e_full, e_direct_only=e_direct, cos=cos,
                     merit_measured=merit_meas, merit_predicted=merit_pred,
                     radiation_to_direct_norm=float(np.linalg.norm(drad[:, :, k][m]) / np.linalg.norm(direct[m])),
                     groups=groups, depth_median=float(np.nanmedian(depth)),
                     depth_worst=[(int(i) + 1, float(depth[i])) for i in worst])
            res[f'{name}_{lev}'] = r
            print(f'{name}_{lev}: e_full={e_full:.3e} e_direct={e_direct:.3e} cos={cos:+.4f} '
                  f'merit meas {merit_meas:.6f} pred {merit_pred:.6f}')
            print('   groups: ' + ', '.join(f'{g}:{v:.2e}' for g, v in groups.items()))
            print('   worst depths: ' + ', '.join(f'd{i}:{v:.2e}' for i, v in r['depth_worst']))
    (out / 'residual_audit.json').write_text(json.dumps(res, indent=1))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('jv_dir', type=Path)
    ap.add_argument('--oracle', type=Path, default=Path('results/jac_check'))
    ap.add_argument('--names', nargs='+', default=['fac0', 'fac05', 'fac1', 'mom'])
    ap.add_argument('--layers', nargs='+', default=['spectra', 'residual'])
    ap.add_argument('--fd-dir', type=Path, default=Path('results/jv/frozenj'))
    ap.add_argument('--r0', type=Path, default=None, help='alternative base residual run')
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    out = a.out or a.jv_dir
    out.mkdir(parents=True, exist_ok=True)
    eps_table = {}
    for line in (a.oracle / 'queue.sh').read_text().splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[1].startswith('s_') and parts[1].endswith('.bin'):
            eps_table[parts[0].replace('SPECS=', '').strip('"')] = float(parts[2].strip('"'))
    if 'spectra' in a.layers:
        audit_spectra(a.jv_dir, a.oracle, a.names, eps_table, out)
    if 'residual' in a.layers:
        audit_residual(a.jv_dir, a.oracle, a.names, eps_table, out, a.r0)
    if 'direct' in a.layers:
        audit_direct(a.jv_dir, a.fd_dir, a.names, out)



def audit_direct(base_dir, fd_dir, names, out):
    """Fixed-J (direct) Jacobian part: dumped BA.s versus frozen-J central differences
    (R(x(1+h s); J0) - R(x(1-h s); J0)) / 2h from replay passes that reread J0."""
    import struct as _st
    pops, R = first_record(Path(base_dir) / 'QSTARS_RESID_TRACE.bin')
    nt, nd = R.shape
    blocks = read_ba(Path(base_dir) / 'QSTARS_BA_TRACE.bin', nt, nd)
    mask = np.broadcast_to(~np.all(R == 0.0, axis=1)[:, None], R.shape)
    res = {}
    for name in names:
        s = np.fromfile(Path('results/jac_check') / f's_{name}.bin', '<f8').reshape((nt, nd), order='F')
        direct = ba_apply(blocks, s, nd)
        for lev in (3, 4):
            p = Path(fd_dir) / f'{name}_h{lev}p' / 'QSTARS_RESID_TRACE.bin'
            m = Path(fd_dir) / f'{name}_h{lev}m' / 'QSTARS_RESID_TRACE.bin'
            if not (p.exists() and m.exists()):
                continue
            h = 10.0 ** -lev / np.max(np.abs(s))
            _, Rp = first_record(p)
            _, Rm = first_record(m)
            fd = (Rp - Rm) / (2 * h)
            e = l2rel((direct - fd)[mask], fd[mask])
            groups = {}
            for g, (a, b) in GROUPS.items():
                mm = np.zeros_like(mask); mm[a - 1:b] = mask[a - 1:b]
                if np.any(fd[mm] != 0):
                    groups[g] = l2rel((direct - fd)[mm], fd[mm])
            depth = np.array([l2rel((direct - fd)[mask[:, d], d], fd[mask[:, d], d]) if np.any(fd[mask[:, d], d] != 0) else np.nan
                              for d in range(nd)])
            worst = np.argsort(np.nan_to_num(depth, nan=-1))[::-1][:6]
            res[f'{name}_h{lev}'] = dict(h=h, e_direct=e, groups=groups, depth_worst=[(int(i) + 1, float(depth[i])) for i in worst])
            print(f'direct {name}_h{lev}: e={e:.3e}  ' + ', '.join(f'{g}:{v:.1e}' for g, v in groups.items()))
            print('   worst depths: ' + ', '.join(f'd{i + 1}:{depth[i]:.2e}' for i in worst))
    (Path(out) / 'direct_audit.json').write_text(json.dumps(res, indent=1))
    return res


if __name__ == '__main__':
    sys.exit(main())
