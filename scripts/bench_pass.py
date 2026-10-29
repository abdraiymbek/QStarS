#!/usr/bin/env python3
"""Report 22: re-run an archived Newton-Krylov pass with identical inputs, time it, sample it, compare it.

The pass is reconstructed exactly as cool_newton_krylov.Runner made it: same restart (--x0), equation table,
VADAT settings, state (QSTARS_STEP.bin), replay field, warm-start field and tangent directions, read from the
archived pass's run_manifest.json. Only the executable, thread count and (optionally) a subset of the
directions may change. While the native process runs, macOS `sample` records call stacks (--sample-every
seconds apart, --sample-secs long), giving a statistical per-routine profile. Afterwards the outputs are
compared with the archived pass: residual trace (x and R), radiation products (QSTARS_JV_RAD.bin) and BA
blocks when both exist. With identical arithmetic only summation order may differ.

    bench_pass.py ARCHIVED_PASS DEST --x0 X0 [--exe EXE] [--threads N] [--dirs-subset K] [--sample-every S]
"""
import argparse, json, subprocess, sys, threading, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
from qstars.nlte_newton import read_ba_blocks, read_fortran, read_jv_rad, read_trace_record  # noqa: E402


def find_native(dest, timeout=600):
    t0 = time.time()
    while time.time() - t0 < timeout:
        out = subprocess.run(['pgrep', '-f', 'cmfgen_dev.exe'], capture_output=True, text=True).stdout.split()
        for pid in out:
            cwd = subprocess.run(['lsof', '-a', '-p', pid, '-d', 'cwd', '-Fn'], capture_output=True, text=True).stdout
            if str(dest) in cwd and 'time' not in subprocess.run(['ps', '-o', 'comm=', '-p', pid],
                                                                  capture_output=True, text=True).stdout:
                return pid
        time.sleep(2)
    return None


def sampler(dest, every, secs, outdir, stop):
    pid = find_native(dest)
    if pid is None:
        return
    k = 0
    time.sleep(min(every, 30))
    while not stop.is_set():
        f = outdir / f'sample_{k:02d}.txt'
        r = subprocess.run(['sample', pid, str(secs), '-file', str(f)], capture_output=True, text=True)
        if r.returncode != 0:
            break
        k += 1
        if stop.wait(every):
            break


def rowwise(ra, rb):
    """Worst over rows (variable index) of max_depth|ra-rb| / max_depth|ra|: every equation judged on its own scale
    (rows span many decades; a globally normalized difference would hide errors in the small rows)."""
    sa = np.max(np.abs(ra), axis=1); d = np.max(np.abs(ra - rb), axis=1)
    ok = sa > 0
    return float(np.max(d[ok] / sa[ok])) if ok.any() else 0.0


def compare_products(a, b, k=None):
    _, ra = read_jv_rad(a / 'QSTARS_JV_RAD.bin'); _, rb = read_jv_rad(b / 'QSTARS_JV_RAD.bin')
    k = rb.shape[2] if k is None else k
    return dict(products_max_rel_per_dir=[float(np.max(np.abs(ra[..., i] - rb[..., i])) /
                                                max(np.max(np.abs(ra[..., i])), 1e-300)) for i in range(k)],
                products_rowwise_max_rel_per_dir=[rowwise(ra[..., i], rb[..., i]) for i in range(k)])


def compare(a, b):
    out = {}
    xa, Ra = read_trace_record(a / 'QSTARS_RESID_TRACE.bin'); xb, Rb = read_trace_record(b / 'QSTARS_RESID_TRACE.bin')
    out['state_identical'] = bool(np.array_equal(xa, xb))
    nz = Ra != 0
    out['residual_max_rel'] = float(np.max(np.abs(Ra[nz] - Rb[nz]) / np.abs(Ra[nz]))) if nz.any() else 0.0
    out['residual_rel_l2'] = float(np.linalg.norm(Ra - Rb) / max(np.linalg.norm(Ra), 1e-300))
    out['residual_rowwise_max_rel'] = rowwise(Ra, Rb)
    if (a / 'QSTARS_JV_RAD.bin').exists() and (b / 'QSTARS_JV_RAD.bin').exists():
        out.update(compare_products(a, b))
    if (a / 'QSTARS_BA_TRACE.bin').exists() and (b / 'QSTARS_BA_TRACE.bin').exists():
        nt, nd = xa.shape
        A = read_ba_blocks(a / 'QSTARS_BA_TRACE.bin', nt, nd); B = read_ba_blocks(b / 'QSTARS_BA_TRACE.bin', nt, nd)
        out['ba_rel_l2'] = float(sum(np.linalg.norm(p - q) for p, q in zip(A, B)) / max(sum(np.linalg.norm(p) for p in A), 1e-300))
        # per (depth, equation): max element difference over the row's l1 scale (all three blocks)
        scale = sum(np.abs(p).sum(axis=2) for p in A)
        diff = np.maximum.reduce([np.abs(p - q).max(axis=2) for p, q in zip(A, B)])
        ok = scale > 0
        out['ba_max_rel_rowscale'] = float(np.max(diff[ok] / scale[ok])) if ok.any() else 0.0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('archived', type=Path)
    ap.add_argument('dest', type=Path)
    ap.add_argument('--x0', type=Path, required=True)
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_28/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--dirs-subset', type=int, default=0, help='use only the first K directions (0: all)')
    ap.add_argument('--sample-every', type=float, default=120.0)
    ap.add_argument('--sample-secs', type=float, default=10.0)
    ap.add_argument('--no-compare', action='store_true')
    ap.add_argument('--obs-formal', action='store_true', help='keep the observer-frame formal solution (build >= 29 skips it by default)')
    ap.add_argument('--compare-only', action='store_true', help='recompute the agreement of an existing DEST and exit')
    ap.add_argument('--extra-env', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'))
    a = ap.parse_args()
    src = a.archived.resolve(); dest = a.dest.resolve()
    if a.compare_only:
        bj = dest.parent / (dest.name + '_bench.json')
        res = json.loads(bj.read_text())
        if (dest / 'QSTARS_DIRS.bin').exists() and (src / 'QSTARS_DIRS.bin').exists() and \
                (dest / 'QSTARS_DIRS.bin').stat().st_size != (src / 'QSTARS_DIRS.bin').stat().st_size:
            res['agreement'] = compare_products(src, dest)
        else:
            res['agreement'] = compare(src, dest)
        bj.write_text(json.dumps(res, indent=1)); print(json.dumps(res['agreement'], indent=1))
        return
    man = json.loads((src / 'run_manifest.json').read_text())
    env = man.get('extra_env', {})
    vad = []
    for line in (src / 'VADAT').read_text().splitlines():
        p = line.split()
        if len(p) >= 2 and p[-1] in ('[INC_AD]', '[INC_ADV]'):
            vad.append((p[-1][1:-1], p[0]))
    x0, _ = read_trace_record(ROOT / 'results' / 'nk167_adv' / 'it00_base' / 'QSTARS_RESID_TRACE.bin') if False else (None, None)
    s_cum = None
    raw_step = (src / 'QSTARS_STEP.bin')
    dirs = None
    if (src / 'QSTARS_DIRS.bin').exists():
        raw = (src / 'QSTARS_DIRS.bin').read_bytes()
        nv, nt, nd = np.frombuffer(raw[:12], '<i4')
        dirs = np.frombuffer(raw[12:], '<f8').reshape((nt, nd, nv), order='F')
        if a.dirs_subset:
            dirs = dirs[:, :, :a.dirs_subset]
    nt = nd = None
    xs, _ = read_trace_record(src / 'QSTARS_RESID_TRACE.bin') if (src / 'QSTARS_RESID_TRACE.bin').exists() else (None, None)
    nt, nd = xs.shape
    ck.NT, ck.ND = nt, nd
    s_cum = read_fortran(raw_step, (nt, nd))
    fix_t = any(len(p) >= 2 and p[-1] == '[FIX_T]' and p[0] == 'T' for p in map(str.split, (src / 'VADAT').read_text().splitlines()))
    R = ck.Runner(a.x0, a.exe, Path(env['QSTARS_REPLACE_FORCE']), fix_t, a.threads, dest.parent / 'bench.log', vad)
    R.slim = False
    R.no_obs = not a.obs_formal
    dest.parent.mkdir(parents=True, exist_ok=True)
    stop = threading.Event(); sdir = dest.parent / (dest.name + '_samples'); sdir.mkdir(exist_ok=True)
    th = threading.Thread(target=sampler, args=(dest, a.sample_every, a.sample_secs, sdir, stop), daemon=True)
    th.start()
    t0 = time.time()
    R.run(dest, s_cum, replay=env.get('QSTARS_JV_REPLAY'), dirs=dirs, ba='QSTARS_BA_DUMP' in env,
          store_j='QSTARS_RJ_STORE' in env, warm=env.get('QSTARS_ALI_WARM'), extra_env=dict(a.extra_env) or None)
    wall = time.time() - t0
    stop.set(); th.join(timeout=5)
    res = dict(archived=str(src), dest=str(dest), exe=str(a.exe), threads=a.threads,
               ndirs=0 if dirs is None else int(dirs.shape[2]), no_obs_formal=R.no_obs,
               extra_env=dict(a.extra_env), wall_s=wall)
    nat = (dest / 'native.log').read_text(errors='replace')
    import re
    m = re.search(r'([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys', nat)
    if m:
        res.update(real=float(m.group(1)), user=float(m.group(2)), sys=float(m.group(3)))
    if (dest / 'QSTARS_ALI_STATS').exists():
        last = (dest / 'QSTARS_ALI_STATS').read_text().splitlines()[-1]
        mm = re.search(r'formal_solutions\s+(\d+)', last)
        if mm:
            res['formal_solutions'] = int(mm.group(1))
    if not a.no_compare:
        res['agreement'] = compare(src, dest) if not a.dirs_subset else {}
        if a.dirs_subset and (src / 'QSTARS_JV_RAD.bin').exists():
            res['agreement'] = compare_products(src, dest)
    (dest.parent / (dest.name + '_bench.json')).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
