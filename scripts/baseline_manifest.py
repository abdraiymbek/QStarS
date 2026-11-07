#!/usr/bin/env python3
"""V1 plan P0: the restartable baseline manifest of a cool-closure run, and a native replay check of its accepted state.

    baseline_manifest.py RUN_DIR --frozen FROZEN.npz --out MANIFEST.json [--replay SCRATCH_DIR]

Records, for the accepted state named by RUN_DIR/accepted_state.json (cnk_trap.write_accepted_pointer):
  * the source: git HEAD, the uncommitted diff (saved next to the manifest, with its sha256), sha256 of the driver files
  * the native pass: executable and capability-record sha256, the sha256 of every input file (from the pass's
    run_manifest.json, digested), the equation map (QSTARS_REPLACE_FORCE file), the two-photon J of the residual
    definition, the stored J and step of the state, VADAT settings
  * the launch: the launch script (sha256) and the last recorded driver command line (RUN_DIR/README.txt)
  * the recorded residual in the frozen reference measure and the driver's local measure (campaign_norms.py)
--replay re-creates, in SCRATCH_DIR, (a) the frozen-J two-photon replay that defines the iterate's residual (<base>_tp)
and (b) the full native pass of the state from x0, with the recorded physics, and compares their residuals with the
recorded ones (relative max and l2 per row group, and in the frozen measure)."""
from __future__ import annotations

import argparse, hashlib, json, re, subprocess, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import campaign_norms as cn  # noqa: E402
import cnk_trap  # noqa: E402
from qstars.nlte_newton import read_trace_record  # noqa: E402

DRIVER_FILES = ['scripts/cool_newton_krylov.py', 'scripts/cnk_trap.py', 'scripts/campaign_norms.py',
                'src/qstars/nlte_newton.py', 'backends/cmfgen/run_component.py', 'backends/cmfgen/formal_ali_trial.py']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def pass_record(d):
    d = Path(d)
    man = json.loads((d / 'run_manifest.json').read_text())
    res = json.loads((d / 'run_result.json').read_text()) if (d / 'run_result.json').exists() else {}
    inputs = man.get('inputs', {})
    vad = [ln.strip() for ln in (d / 'VADAT').read_text().splitlines()] if (d / 'VADAT').exists() else []
    env = man.get('extra_env', {})
    return dict(dir=str(d.resolve()), executable=man.get('executable'), executable_sha256=man.get('executable_sha256'),
                inputs_count=len(inputs),
                inputs_digest=hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest(),
                extra_env=env, threads=man.get('threads'), seconds=res.get('seconds'), native_status=res.get('status'),
                native_error=res.get('error'),
                vadat_qstars=[ln for ln in vad if re.search(r'\[(QS_|FIX_T|TWO_METH|INC_AD|STOP_IF_MW|COMP_F)', ln)],
                files={f: sha(d / f) for f in ('QSTARS_RESID_TRACE.bin', 'QSTARS_STEP.bin', 'QSTARS_RJ.bin',
                                               'QSTARS_BA_TRACE.bin') if (d / f).exists()})


def compare(rec_dir, new_dir, frozen, x0):
    xa, Ra = read_trace_record(Path(rec_dir) / 'QSTARS_RESID_TRACE.bin')
    xb, Rb = read_trace_record(Path(new_dir) / 'QSTARS_RESID_TRACE.bin')
    z = np.load(frozen)
    w = z['w']; nr = w.shape[0]
    out = dict(state_identical=bool(np.array_equal(xa, xb)), residual_identical=bool(np.array_equal(Ra, Rb)),
               frozen_recorded=float(np.linalg.norm(w * Ra[:nr])), frozen_replayed=float(np.linalg.norm(w * Rb[:nr])),
               frozen_of_difference=float(np.linalg.norm(w * (Ra - Rb)[:nr])), groups={})
    for g, rows in cn.GROUPS.items():
        a_, b_ = Ra[rows], Rb[rows]
        n_ = np.linalg.norm(a_)
        out['groups'][g] = dict(rel_l2=float(np.linalg.norm(a_ - b_) / n_) if n_ > 0 else 0.0,
                                frozen_diff=float(np.linalg.norm((w[rows] * (a_ - b_)))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run', type=Path)
    ap.add_argument('--frozen', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--replay', type=Path, default=None, help='scratch directory for the native replay check')
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda167')
    a = ap.parse_args()
    run = a.run.resolve()
    ptr = cnk_trap.read_accepted_pointer(run)
    if ptr is None:
        raise SystemExit(f'{run}: no accepted_state.json')
    base = Path(ptr['base'])
    why = cnk_trap.check_restart_base(base)
    if why is not None:
        raise SystemExit(f'the pointer names an unaccepted state: {why}')
    rd = cn.residual_dir(base)
    m = dict(created=time.strftime('%Y-%m-%d %H:%M:%S'), run=str(run), accepted_pointer=ptr,
             restart_check='admissible (cnk_trap.check_restart_base)')
    # source
    diff = git('diff', 'HEAD')
    dp = a.out.with_name(a.out.stem + '_source.diff')
    dp.write_text(diff)
    untracked = [l for l in git('ls-files', '--others', '--exclude-standard').splitlines() if l.endswith(('.py', '.sh'))]
    m['source'] = dict(head=git('rev-parse', 'HEAD').strip(), dirty_files=git('diff', '--stat', 'HEAD').strip().splitlines(),
                       diff_file=str(dp), diff_sha256=sha(dp),
                       files={f: sha(ROOT / f) for f in DRIVER_FILES if (ROOT / f).exists()},
                       untracked_scripts={f: sha(ROOT / f) for f in untracked})
    # native passes of the accepted state
    m['state_pass'] = pass_record(base)
    m['residual_pass'] = pass_record(rd)
    env = m['residual_pass']['extra_env']
    rep = Path(env['QSTARS_REPLACE_FORCE'])
    exe = Path(m['state_pass']['executable'])
    cap = exe.with_suffix('.capabilities.json')
    m['physics'] = dict(
        executable_sha256_now=sha(exe), executable_matches=sha(exe) == m['state_pass']['executable_sha256'],
        capabilities=str(cap), capabilities_sha256=sha(cap) if cap.exists() else None,
        equation_map=str(rep), equation_map_sha256=sha(rep),
        two_photon_j_of_residual=env.get('QSTARS_TWO_PHOT_J'),
        two_photon_j_sha256=sha(env['QSTARS_TWO_PHOT_J']) if env.get('QSTARS_TWO_PHOT_J') else None,
        x0=str(a.x0.resolve()), x0_restart_sha256={f: sha(a.x0 / f) for f in ('MODEL', 'SCRATCH_REC') if (a.x0 / f).exists()},
        x0_base_pass=str(cn.X0_BASE))
    # launch
    lines = [l for l in (run / 'README.txt').read_text().splitlines() if 'driver started' in l]
    launch = re.search(r'\((launch_\S+\.sh)', lines[-1]).group(1) if lines and re.search(r'\((launch_\S+\.sh)', lines[-1]) else None
    lp = run.parent / launch if launch else None
    m['launch'] = dict(last_driver_start=lines[-1] if lines else None, script=str(lp) if lp else None,
                       script_sha256=sha(lp) if lp and lp.exists() else None,
                       script_text=lp.read_text() if lp and lp.exists() else None,
                       resume=f'{lp} <method flags> with --resume replacing BASE K' if lp else None)
    # recorded residual, frozen and local
    z = np.load(a.frozen)
    x0 = cn.configure(str(z['x0_base']))
    m['frozen_measure'] = dict(file=str(a.frozen.resolve()), sha256=sha(a.frozen), reference=str(z['reference']),
                               settings=json.loads(str(z['settings'])))
    m['recorded_residual'] = cn.score(a.frozen, base, x0)
    hist = json.loads((run / 'history.json').read_text())
    last = hist[-1]
    m['last_iterate'] = dict(it=last['it'], base=last['base'], driver_merit=last['merit'],
                             best_trial=last.get('best_trial', last.get('accepted')),
                             accepted=cnk_trap.history_accepted(last),
                             status=last.get('status', 'stopped_no_decrease' if cnk_trap.history_accepted(last) is None
                                             else 'accepted (legacy record)'))
    if a.replay is not None:
        import cool_newton_krylov as ck
        scratch = a.replay.resolve(); scratch.mkdir(parents=True, exist_ok=True)
        vad = [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')]
        R = ck.Runner(a.x0, exe, rep, False, 4, scratch / 'replay.log', vad)
        R.no_obs = True
        xs, _ = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
        s_cum = 1.0 - xs / np.asarray(x0)
        recorded_step = np.frombuffer((base / 'QSTARS_STEP.bin').read_bytes(), dtype='<f8')
        m['replay'] = dict(scratch=str(scratch))
        t0 = time.time()
        tp = R.run(scratch / 'residual_tp', s_cum, replay=base / 'QSTARS_RJ.bin', ba=True, store_j=False, threads=1,
                   two_phot_j=Path(env['QSTARS_TWO_PHOT_J']))
        m['replay']['residual_tp'] = dict(seconds=time.time() - t0, **compare(rd, tp, a.frozen, x0))
        senv = m['state_pass']['extra_env']
        t0 = time.time()
        st_dest = scratch / 'state_pass'
        # the recorded step of the state, byte for byte (x0 * (1 - s) reproduces the state's populations)
        st = recorded_step.reshape(s_cum.shape, order='F')
        sp = R.run(st_dest, st, threads=m['state_pass']['threads'] or 2,
                   warm=senv.get('QSTARS_ALI_WARM'), two_phot_j=senv.get('QSTARS_TWO_PHOT_J'))
        m['replay']['state_pass'] = dict(seconds=time.time() - t0, **compare(base, sp, a.frozen, x0))
    a.out.write_text(json.dumps(m, indent=1, default=str))
    print(json.dumps({k: m[k] for k in ('accepted_pointer', 'last_iterate')}, indent=1))
    print(f'frozen measure of the recorded residual: {m["recorded_residual"]["frozen"]:.6e}')
    if 'replay' in m:
        for k in ('residual_tp', 'state_pass'):
            r = m['replay'][k]
            print(f'replay {k}: {r["seconds"]:.0f} s, state identical {r["state_identical"]}, residual identical '
                  f'{r["residual_identical"]}, frozen {r["frozen_recorded"]:.6e} -> {r["frozen_replayed"]:.6e} '
                  f'(difference {r["frozen_of_difference"]:.2e})')
    print(f'manifest: {a.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
