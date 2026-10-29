"""Report 23 tooling for Gates 1-3: closing chain control flow, trace pairing, profile refinement metrics,
restart copying and the pre-declared plans. Synthetic data only (no native runs)."""
import hashlib
import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    sys.path.insert(0, str(ROOT / 'scripts'))
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------- post_convergence.sh control flow
STUB_CLEANUP = '''import json, sys, os
from pathlib import Path
a = sys.argv[1:]
open(os.environ['STUB_CALLS'], 'a').write('cleanup ' + ' '.join(a) + '\\n')
run, start, tag = Path(a[0]), a[1], a[a.index('--tag') + 1]
ok = os.environ.get('STUB_CONVERGED', '1') == '1'
hist = [dict(sweep=0, base=start, tol=1e-4, max_trace_step=1e-6 if ok else 0.5, converged=ok)]
(run / f'cleanup_{tag}_history.json').write_text(json.dumps(hist))
'''
STUB_VERIFY = '''import sys, os
from pathlib import Path
open(os.environ['STUB_CALLS'], 'a').write('verify ' + ' '.join(sys.argv[1:]) + '\\n')
d = Path(sys.argv[2]); d.mkdir(); d.with_name(d.name + '_audit').mkdir()
'''


def _chain(tmp_path, env_extra):
    (tmp_path / 'scripts').mkdir(exist_ok=True)
    shutil.copy2(ROOT / 'scripts/post_convergence.sh', tmp_path / 'scripts/post_convergence.sh')
    (tmp_path / 'scripts/trace_cleanup.py').write_text(STUB_CLEANUP)
    (tmp_path / 'scripts/verify_cool_state.py').write_text(STUB_VERIFY)
    (tmp_path / '.venv/bin').mkdir(parents=True, exist_ok=True)
    if not (tmp_path / '.venv/bin/python').exists():
        os.symlink(sys.executable, tmp_path / '.venv/bin/python')
    run = tmp_path / 'run'; (run / 'it30_trial_l1').mkdir(parents=True, exist_ok=True)
    env = {**os.environ, 'STUB_CALLS': str(tmp_path / 'calls.txt'), **env_extra}
    return subprocess.run(['bash', str(tmp_path / 'scripts/post_convergence.sh'), str(run), 'it30_trial_l1', 'X0', 'TABLE',
                           'FRAC_SP', '1.04'], env=env, capture_output=True, text=True)


def test_post_convergence_chain_guards_and_per_state_destinations(tmp_path):
    p = _chain(tmp_path, {})
    assert p.returncode == 0, p.stderr
    calls = (tmp_path / 'calls.txt').read_text().splitlines()
    run = tmp_path / 'run'
    # cleanup converged at sweep 0: no cleanup_* pass exists, the base itself is verified (critic B1)
    assert '--threads 2 --pairing structural --tag it30_trial_l1_structural_build_31 ' in calls[0]
    assert 'build_31/cur_cmf/exe/cmfgen_dev.exe' in calls[0] and '--vadat FRAC_SP 1.04' in calls[0]
    assert calls[1].startswith(f'verify {run}/it30_trial_l1 {run}/verify_it30_trial_l1 ') and '--threads 6' in calls[1]
    # a rerun never collides with the earlier verification
    p = _chain(tmp_path, {'EXE': 'other.exe', 'VERIFY_THREADS': '3'})
    assert p.returncode == 0, p.stderr
    calls = (tmp_path / 'calls.txt').read_text().splitlines()
    last = calls[-1]
    assert f'{run}/verify_it30_trial_l1_r2 ' in last and '--exe other.exe' in last and '--threads 3' in last
    # another executable never reuses the cleanup passes of the previous one (tag carries the build)
    assert '--tag it30_trial_l1_structural_other ' in calls[-2] and (run / 'cleanup_it30_trial_l1_structural_other_history.json').exists()
    # an unconverged cleanup stops the chain before verification
    n = len((tmp_path / 'calls.txt').read_text().splitlines())
    p = _chain(tmp_path, {'STUB_CONVERGED': '0'})
    assert p.returncode != 0 and 'NOT converged' in p.stderr
    calls = (tmp_path / 'calls.txt').read_text().splitlines()
    assert len(calls) == n + 1 and calls[-1].startswith('cleanup')
    assert not (run / 'verify_it30_trial_l1_r3').exists()


# ---------------------------------------------------------------- trace classification and pairing
def test_l1_split_is_order_independent():
    ncc = _script('native_correction_check')
    nt, nd = 4, 1
    frac = np.array([[1.0], [1e-12], [1e-12], [1e-12]])
    active = np.ones((nt, nd), bool)
    di = np.zeros((nd, nt, nt))
    # row 1 is significant-dominated; rows 2, 3 couple to row 1's variable. With the old in-loop view
    # mutation, demoting row 1 first turned column 1 'significant' for rows 2 and 3.
    di[0, 1] = [5.0, 1.0, 0.0, 0.0]
    di[0, 2] = [0.0, 3.0, 2.0, 0.0]
    di[0, 3] = [0.0, 3.0, 0.0, 2.0]
    kept, demoted = ncc.l1_split(frac, active, di, 1e-10)
    assert kept[:, 0].tolist() == [False, False, True, True] and demoted[:, 0].tolist() == [False, True, False, False]


def _hhe_frac(trace_rows):
    frac = np.full((134, 1), 0.5)
    frac[trace_rows, 0] = 1e-30
    return frac


def test_structural_pairing_follows_cmfgen_equation_selection():
    ncc = _script('native_correction_check')
    act = np.ones((134, 1), bool)
    rep = np.zeros((1, 6), bool); rep[0, 4] = True                       # He II ground row = C_ION(He II)
    # He III trace, He II ground significant: He III solved with C_ION(He II) (row 102), never row 132
    row_of, kind, unp = ncc.structural_pairing(_hhe_frac([131, 120, 0]), act, rep, 1e-10)
    assert row_of[131, 0] == 101 and kind[131, 0] == 2 and row_of[120, 0] == 120 and row_of[0, 0] == 0
    assert not unp.any() and not np.isin(row_of, [31, 131, 132, 133]).any()
    # no C_ION(He II) row at this depth: He III has no equation of its own -> unpaired, held fixed
    row_of, kind, unp = ncc.structural_pairing(_hhe_frac([131]), act, np.zeros((1, 6), bool), 1e-10)
    assert unp[131, 0] and row_of[131, 0] == -1
    # He II ground (replaced) and He III trace: chain anchored at He I through C_ION(He I) (row 33)
    rep2 = rep.copy(); rep2[0, 3] = True
    row_of, kind, unp = ncc.structural_pairing(_hhe_frac([101, 131]), act, rep2, 1e-10)
    assert row_of[101, 0] == 32 and row_of[131, 0] == 101 and not unp.any()
    # H- with its row replaced: C_ION(H-) determines it, anchored at the significant H I (downward)
    rep3 = np.zeros((1, 6), bool); rep3[0, 0] = True
    row_of, kind, unp = ncc.structural_pairing(_hhe_frac([0]), act, rep3, 1e-10)
    assert row_of[0, 0] == 0 and kind[0, 0] == 2


def test_structural_step_zeroes_paired_rows_exactly():
    ncc = _script('native_correction_check')
    from qstars.nlte_newton import BlockTridiag
    rng = np.random.default_rng(5)
    nt, nd = 134, 3
    lo = rng.normal(size=(nd, nt, nt)) * 1e-2; up = rng.normal(size=(nd, nt, nt)) * 1e-2
    di = rng.normal(size=(nd, nt, nt)) + 20 * np.eye(nt)[None]
    scale = 10.0 ** rng.uniform(-30, 0, size=(nd, nt))                # rows spanning 30 decades
    lo *= scale[:, :, None]; di *= scale[:, :, None]; up *= scale[:, :, None]
    lo[0] = 0; up[-1] = 0
    x = np.full((nt, nd), 0.5); x[[0, 120, 131]] = 1e-30
    R = rng.normal(size=(nt, nd)) * scale.T
    rep = np.zeros((nd, 6), bool); rep[:, 4] = True
    s, solved, unpaired, bw = ncc.structural_step(x, R, lo, di, up, rep, 1e-10)
    assert solved.sum() == 3 * nd and not unpaired.any() and bw < 1e-12
    row_of, _, _ = ncc.structural_pairing(ncc.fractions(x), R != 0, rep, 1e-10)
    Rn = R - BlockTridiag(lo, di, up).matvec(np.where(solved, s, 0.0))  # linear residual after the update
    for v, d in zip(*np.nonzero(solved)):
        assert abs(Rn[row_of[v, d], d]) <= 1e-10 * abs(R[row_of[v, d], d])
    # significant components of the native correction are unchanged by a trace-only update
    M = BlockTridiag(lo, di, up).factor()
    ya, yb = M.solve(R), M.solve(Rn)
    assert np.max(np.abs(ya - yb)[~solved]) <= 1e-8 * np.max(np.abs(ya[~solved]))
    assert np.allclose(ya[solved] - yb[solved], s[solved], rtol=1e-6)


# ---------------------------------------------------------------- native correction check and trace cleanup wiring
def _write_pass(d, x, R, blocks, rep, manifest=None):
    """A synthetic base pass in the native trace formats (one RESID record, 3 nd - 2 BA blocks, equation table)."""
    lo, di, up = blocks
    nt, nd = x.shape
    d.mkdir(parents=True, exist_ok=True)
    (d / 'QSTARS_RESID_TRACE.bin').write_bytes(struct.pack('<3i', 1, nt, nd) + x.tobytes(order='F') + R.tobytes(order='F'))
    with open(d / 'QSTARS_BA_TRACE.bin', 'wb') as f:
        for k in range(1, nd + 1):
            for b, blk in ((1, lo), (2, di), (3, up)):
                if not ((b == 1 and k == 1) or (b == 3 and k == nd)):
                    f.write(struct.pack('<3i', k, b, nt) + np.ascontiguousarray(blk[k - 1]).tobytes(order='F'))
    (d / 'QSTARS_REPLACE_TRACE').write_text(''.join(f'{k + 1} ' + ' '.join(str(int(v)) for v in r) + '\n'
                                                    for k, r in enumerate(rep)))
    if manifest is not None:
        (d / 'run_manifest.json').write_text(json.dumps(manifest))
    return d


def _hhe_system(nd=3, rscale=0.1, seed=7):
    """134 x nd H/He system: rows over 30 decades, trace H-, He II level 120 and He III; R small enough that
    every trace step |s| < 1 (the cleanup update x(1 - s) is then the exact linear solution)."""
    rng = np.random.default_rng(seed)
    nt = 134
    lo = rng.normal(size=(nd, nt, nt)) * 1e-2; up = rng.normal(size=(nd, nt, nt)) * 1e-2
    di = rng.normal(size=(nd, nt, nt)) + 20 * np.eye(nt)[None]
    scale = 10.0 ** rng.uniform(-30, 0, size=(nd, nt))
    lo *= scale[:, :, None]; di *= scale[:, :, None]; up *= scale[:, :, None]
    lo[0] = 0; up[-1] = 0
    x = np.full((nt, nd), 0.5); x[[0, 120, 131]] = 1e-30
    R = rng.normal(size=(nt, nd)) * scale.T * rscale
    rep = np.zeros((nd, 6), bool); rep[:, 4] = True
    return x, R, (lo, di, up), rep


def test_audit_normalized_matches_native_audit_parser():
    ncc = _script('native_correction_check')
    from qstars.cmfgen_backend import nonlinear_corrections
    for y in (np.array([-3.0, 0.1, 2e-4]), np.array([-2e-3, 0.5, 0.0]), np.array([-1e-2, 0.9999995]), np.array([-4e-4, 7e-4])):
        inc, dec = 100 * max(0.0, -y.min()), 100 * max(0.0, y.max())   # POPS -> POPS(1 - y): y > 0 is a decrease
        txt = (' Current great iteration count is   1\n'
               f' Maximum % increase at depth   34 is  {inc:.12E}  (LAMBDA)          --- iteration    1\n'
               f' Maximum % decrease at depth   20 is  {dec:.12E}  (LAMBDA)          --- iteration    1\n'
               ' Maximm changes as returned by SOLVEBA_V13 is   1.0E-30\n')
        assert np.isclose(nonlinear_corrections(txt)[-1]['all_level_percent'] / 100, ncc.audit_normalized(y).max(), rtol=1e-10)


def test_summarize_gate_location_and_fail_at_floor():
    ncc = _script('native_correction_check')
    y = np.array([[0.6], [-1.0], [2e-3]])
    fl = np.array([[0.0], [2e-3], [0.9995e-3]])           # raw 0.9995e-3 < 1e-3, audit-normalized 1.0005e-3 > 1e-3
    frac = np.array([[1.0], [0.5], [1e-20]])
    cls = dict(significant=np.array([[True], [True], [False]]), all=np.ones((3, 1), bool), empty=np.zeros((3, 1), bool))
    out = ncc.summarize(y, cls, frac, fl, 1e-3, model=np.full((3, 1), 7.0))
    sig = out['significant']
    # gate = decrease 0.6 -> 0.6/0.4 = 1.5 at row 1; the largest |y| (1.0) is another variable (row 2)
    assert np.isclose(sig['audit_normalized'], 1.5) and (sig['gate_row'], sig['gate_depth'], sig['gate_y']) == (1, 1, 0.6)
    assert (sig['max_abs'], sig['max_abs_row'], sig['max_abs_frac']) == (1.0, 2, 0.5)
    a = out['all']
    assert a['n_above_threshold'] == 3 and a['n_fail_at_floor'] == 2 and a['n_fail_above_floor'] == 1
    assert [r[:2] for r in a['fail_at_floor']] == [[2, 1], [3, 1]] and sig['n_fail_at_floor'] == 1
    assert a['model_floor_max'] == 7.0 and out['empty'] == dict(n=0)
    assert ncc.summarize(y, cls, frac)['all'].keys().isdisjoint({'n_fail_at_floor', 'model_floor_max'})


def test_check_cleanup_prediction_and_measured_floor(tmp_path):
    ncc = _script('native_correction_check')
    from qstars.nlte_newton import BlockTridiag
    x, R, blocks, rep = _hhe_system()
    man = dict(executable_sha256='a' * 64, threads=2)
    base = _write_pass(tmp_path / 'base', x, R, blocks, rep, man)
    out, y = ncc.check(base, cleanup=True)
    M = BlockTridiag(*blocks)
    Mf = BlockTridiag(*blocks).factor()
    assert np.allclose(y, Mf.solve(R), rtol=1e-12, atol=0)
    s, solved, _, _ = ncc.structural_step(x, R, *blocks, rep, 1e-10)
    yb = Mf.solve(R - M.matvec(s))                       # native correction after the trace update x(1 - s)
    frac, cls = ncc.classes(x, R, blocks[1], rep, 1e-10)
    want = ncc.summarize(yb, cls, frac)
    got = out['after_structural_cleanup']['classes']
    assert out['after_structural_cleanup']['trace_solved'] == int(solved.sum()) == 9
    for name, w in want.items():
        g = got[name]
        assert g['n'] == w['n']
        if w['n']:
            assert np.isclose(g['audit_normalized'], w['audit_normalized'], rtol=1e-8)
            assert np.isclose(g['max_abs'], w['max_abs'], rtol=1e-8) and g['n_above_threshold'] == w['n_above_threshold']
            assert (g['gate_row'], g['gate_depth']) == (w['gate_row'], w['gate_depth'])
    # the significant part of the native correction is unchanged by the trace-only update (first order, exact here)
    assert np.allclose(yb[~solved], y[~solved], rtol=1e-8, atol=0)
    # measured floor: identical-state passes (same build other threads; exact repeat) through BASE's own M
    rng = np.random.default_rng(3)
    d1 = rng.normal(size=R.shape) * np.abs(R) * 1e-9; d2 = rng.normal(size=R.shape) * np.abs(R) * 1e-12
    p1 = _write_pass(tmp_path / 'p1', x, R + d1, blocks, rep, dict(man, threads=8))
    p2 = _write_pass(tmp_path / 'p2', x, R + d2, blocks, rep, man)
    out, _ = ncc.check(base, noise_pairs=[p1, p2], noise_eps=1e-12, cleanup=True)
    fl = np.maximum(np.abs(Mf.solve(d1)), np.abs(Mf.solve(d2)))
    f = out['floor']
    assert f['composition'] == dict(n_pairs=2, same_build_other_threads=True, exact_repeat=True, as_declared=True)
    assert f['use'].startswith('declared floor') and [p['same_threads'] for p in f['pairs']] == [False, True]
    assert np.isclose(out['classes']['all']['floor_max'], fl.max(), rtol=1e-8)
    assert np.isclose(f['pairs'][0]['class_max']['all'], np.abs(Mf.solve(d1)).max(), rtol=1e-8)
    ya = ncc.audit_normalized(y); fa = ncc.audit_normalized(fl); act = R != 0
    assert out['classes']['all']['n_fail_at_floor'] == int(np.sum(act & (ya > 1e-3) & (fa > 1e-3)))
    assert 'model_floor_max' in out['classes']['all'] and out['noise_model']['use'].startswith('report only')
    # one pair only: diagnostic, never the declared floor
    out, _ = ncc.check(base, noise_pairs=[p1])
    assert not out['floor']['composition']['as_declared'] and out['floor']['use'].startswith('diagnostic only')
    # a pass at another state is refused
    bad = _write_pass(tmp_path / 'bad', x * (1 + 1e-15), R, blocks, rep, man)
    with pytest.raises(SystemExit):
        ncc.check(base, noise_pairs=[bad])


def test_trace_cleanup_history_parsed_by_post_convergence(tmp_path, monkeypatch):
    tc = _script('trace_cleanup')
    from qstars.nlte_newton import BlockTridiag, read_trace_record
    x, R, blocks, rep = _hhe_system(rscale=0.05)
    M = BlockTridiag(*blocks)
    run = tmp_path / 'run'
    x0 = x * 1.25
    _write_pass(run / 'it00_base', x0, R, blocks, rep)
    _write_pass(run / 'base', x, R, blocks, rep)
    table = tmp_path / 'table.txt'; table.write_text((run / 'base/QSTARS_REPLACE_TRACE').read_text())
    exe = tmp_path / 'cmfgen_dev.exe'; exe.write_bytes(b'stub executable')

    class StubRunner:
        """Runner.run at the linear trace subsystem: R(x(1 - s)) = R - M s (exact for the synthetic model)."""
        def __init__(self, *a):
            pass

        def log(self, msg):
            pass

        def run(self, dest, s_cum, warm=None, **kw):
            xp, Rp = read_trace_record(Path(warm).parent / 'QSTARS_RESID_TRACE.bin')
            xn = x0 * (1.0 - s_cum)
            return _write_pass(Path(dest), xn, Rp - M.matvec(1.0 - xn / xp), blocks, rep)
    monkeypatch.setattr(tc.ck, 'Runner', StubRunner)
    argv = ['trace_cleanup.py', str(run), str(run / 'base'), '--x0', str(tmp_path), '--exe', str(exe), '--replace', str(table),
            '--free-t', '--pairing', 'structural']
    monkeypatch.setattr(sys, 'argv', argv + ['--tag', 'T1'])
    tc.main()
    h = json.loads((run / 'cleanup_T1_history.json').read_text())
    assert [r['sweep'] for r in h] == [0, 1] and not h[0]['converged'] and h[0]['max_trace_step'] > 1e-4
    assert h[-1]['converged'] and h[-1]['max_trace_step'] < h[-1]['tol'] == 1e-4 and h[-1]['base'] == str(run / 'cleanup_T1_01')
    assert h[0]['trace_variables'] == 9 and h[0]['held_fixed_unpaired'] == 0 and 'T' in h[-1]['native_correction_audit_normalized']
    # the parser of post_convergence.sh selects the cleaned state from this history
    sh = (ROOT / 'scripts/post_convergence.sh').read_text()
    parser = sh.split("<<'EOF'\n", 1)[1].split('\nEOF', 1)[0]
    p = subprocess.run([sys.executable, '-', str(run / 'cleanup_T1_history.json')], input=parser, capture_output=True, text=True)
    assert p.returncode == 0 and p.stdout.strip() == str(run / 'cleanup_T1_01')
    # a cleanup pass made by another executable is never reused
    (run / 'cleanup_T2_01').mkdir(); (run / 'cleanup_T2_01/run_manifest.json').write_text(json.dumps(dict(executable_sha256='0' * 64)))
    monkeypatch.setattr(sys, 'argv', argv + ['--tag', 'T2'])
    with pytest.raises(SystemExit, match='another executable'):
        tc.main()
    # an x0 reference holding a record inherited from its restart is refused
    rec = (run / 'it00_base/QSTARS_RESID_TRACE.bin').read_bytes()
    (run / 'it00_base/QSTARS_RESID_TRACE.bin').write_bytes(rec + rec)
    monkeypatch.setattr(sys, 'argv', argv + ['--tag', 'T3'])
    with pytest.raises(RuntimeError, match='not one'):
        tc.main()


# ---------------------------------------------------------------- profile refinement metrics
def _profiles(r, scale=1.0, heat=1.0):
    L = 3e39 * (1 + 0.1 * np.log(r / r[0]))
    F = L / (4 * np.pi * r ** 2)
    return dict(radius=r, E_comoving=1e3 / r * scale, F_comoving=F * scale, P_comoving=3e2 / r,
                radiative_acceleration=1e2 * r / r[0], kappa_flux=0.34 + 0 * r, luminosity_lab=L * scale,
                heating_comoving=heat * 3e36 / (4 * np.pi * r[-1]) * np.sin(7 * r / r[-1]) / r ** 2)


def test_profile_metrics_on_nested_and_identical_grids():
    cr = _script('compare_resolution')
    crit = json.loads((ROOT / 'benchmarks/cmfgen_cool_frequency_plan.json').read_text())['refinement_profile_criteria']
    rf = np.geomspace(1e12, 1e14, 41)
    fine = _profiles(rf)
    same = cr.compare_profiles(_profiles(rf), fine, crit)
    assert same['shared_radii'] == 41 and same['pass'] and all(m['value'] == 0 for m in same['metrics'].values())
    coarse = _profiles(rf[::2], scale=1.004)                           # nested coarse grid, 0.4 % in E/F/L
    coarse = {k: v[::-1] for k, v in coarse.items()}                  # native order (outer first) is handled
    out = cr.compare_profiles(coarse, fine, crit)
    assert out['shared_radii'] == 21
    assert np.isclose(out['metrics']['E_comoving']['value'], 0.004) and out['metrics']['P_comoving']['value'] < 1e-14
    assert np.isclose(out['metrics']['luminosity_lab_over_Lscale']['emergent_relative'], 0.004)
    assert out['pass'] and out['pending_user_approval'] == ['cumulative_heating_over_Lscale', 'luminosity_lab_over_Lscale']
    bad = _profiles(rf[::2]); bad['kappa_flux'] = bad['kappa_flux'] * 1.02; bad['kappa_flux'][3] = np.nan
    out = cr.compare_profiles(bad, fine, crit)
    assert out['metrics']['kappa_flux']['status'] == 'FAIL' and out['metrics']['kappa_flux']['nonfinite_points'] == 1
    # heating: the cumulative 4 pi r^2 q dr difference, normalized by the luminosity scale
    out = cr.compare_profiles(_profiles(rf, heat=3.0), fine, crit)
    from scipy.integrate import cumulative_trapezoid
    Q = lambda p: 4 * np.pi * cumulative_trapezoid(p['radius'] ** 2 * p['heating_comoving'], p['radius'], initial=0)
    want = np.max(np.abs(Q(_profiles(rf, heat=3.0)) - Q(fine))) / max(fine['luminosity_lab'][[0, -1]])
    assert 1e-5 < want < 1e-2 and np.isclose(out['metrics']['cumulative_heating_over_Lscale']['value'], want, rtol=1e-10, atol=0)
    with pytest.raises(SystemExit):
        cr.compare_profiles(_profiles(rf * 1.01), fine, crit)            # not nested


def test_profile_mode_reproduces_archived_hot_frequency_pair():
    cr = _script('compare_resolution')
    e = ROOT / 'Reports/evidence/backend_20260928/final_audit'
    if not (e / 'transport_frequency/profiles.npz').exists():
        pytest.skip('archived hot frequency pair not present')
    crit = json.loads((ROOT / 'benchmarks/cmfgen_cool_frequency_plan.json').read_text())['refinement_profile_criteria']
    out = cr.compare_profiles(cr.load_profiles(e / 'transport_rel'), cr.load_profiles(e / 'transport_frequency'), crit)
    ref = next(c['values'] for c in json.loads((e / 'resolution.json').read_text())['comparisons'] if c['kind'] == 'frequency')
    for k in ('E_comoving', 'F_comoving', 'P_comoving', 'radiative_acceleration', 'kappa_flux'):
        assert np.isclose(out['metrics'][k]['value'], ref[k], rtol=1e-12)


def _audit_dir(d, prof, failing=(), cool=True):
    d.mkdir(parents=True)
    np.savez(d / 'profiles.npz', **prof)
    checks = [dict(name=n, status='FAIL' if n in failing else 'PASS') for n in ('electron_balance', 'last_nonlinear_correction')]
    (d / 'validation.json').write_text(json.dumps(dict(checks=checks)))
    if cool:
        (d / 'cool_validation.json').write_text(json.dumps(dict(additional_checks=[
            dict(name=n, status='FAIL' if n in failing else 'PASS') for n in ('cumulative_radiation_energy',)])))
    return d


def test_profile_mode_requires_both_per_state_audits(tmp_path, monkeypatch):
    cr = _script('compare_resolution')
    rf = np.geomspace(1e12, 1e14, 41)

    def run(coarse, fine, name):
        monkeypatch.setattr(sys, 'argv', ['compare_resolution.py', '--profiles', str(coarse), str(fine), '--out', str(tmp_path / name)])
        cr.main()
        return json.loads((tmp_path / name).read_text())
    good = _audit_dir(tmp_path / 'fine_audit', _profiles(rf))
    ok = run(_audit_dir(tmp_path / 'ok_audit', _profiles(rf[::2])), good, 'ok.json')
    assert ok['profiles']['comparable'] and ok['profiles']['status'] == 'PASS' and ok['pass_all']
    # identical profiles (metrics pass) but one failing per-state check: not comparable, never a pass
    bad = run(_audit_dir(tmp_path / 'bad_audit', _profiles(rf[::2]), failing=('last_nonlinear_correction',)), good, 'bad.json')
    pr = bad['profiles']
    assert pr['metrics_pass'] and not pr['pass'] and pr['status'] == 'NOT COMPARABLE' and not pr['comparable'] and not bad['pass_all']
    assert pr['per_state_audits'][0]['checks']['validation:last_nonlinear_correction'] == 'FAIL'
    assert pr['per_state_audits'][0]['status'] == 'FAIL' and pr['per_state_audits'][1]['status'] == 'PASS'
    # missing cool audit, bare profiles.npz: not comparable
    nocool = run(_audit_dir(tmp_path / 'nocool_audit', _profiles(rf[::2]), cool=False), good, 'nocool.json')
    assert nocool['profiles']['per_state_audits'][0]['missing'] == ['cool_validation.json'] and not nocool['pass_all']
    bare = run(tmp_path / 'ok_audit/profiles.npz', good, 'bare.json')
    assert bare['profiles']['per_state_audits'][0]['status'] == 'unaudited' and bare['profiles']['status'] == 'NOT COMPARABLE'
    assert not bare['profiles']['pass'] and not bare['pass_all']


# ---------------------------------------------------------------- restart copying and declarations
def test_continuation_never_copies_trace_records(tmp_path):
    sys.path.insert(0, str(ROOT / 'backends/cmfgen'))
    from continue_component import prepare_continuation
    src = tmp_path / 'restart'; src.mkdir()
    for n in ('MOD_SUM', 'IN_ITS', 'VADAT', 'QSTARS_ALI_STATS', 'QSTARS_RESID_TRACE.bin', 'QSTARS_BA_TRACE.bin',
              'QSTARS_REPLACE_TRACE', 'QSTARS_SCALE_TRACE', 'QSTARS_STEQ_TRACE.bin'):
        (src / n).write_text('x\n')
    (src / 'bootstrap_traces').mkdir(); (src / 'bootstrap_traces/QSTARS_SOL_TRACE.bin').write_text('x')
    dst = prepare_continuation(src, tmp_path / 'pass', 1)
    names = {p.name for p in dst.iterdir()}
    assert not any('TRACE' in n for n in names) and 'bootstrap_traces' not in names
    assert {'MOD_SUM', 'QSTARS_ALI_STATS', 'VADAT'} <= names


def _sha(block):
    b = {k: v for k, v in block.items() if k != 'declaration_sha256'}
    return hashlib.sha256(json.dumps(b, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def test_predeclared_plans_are_intact():
    cool = json.loads((ROOT / 'benchmarks/cmfgen_cool_validation_plan.json').read_text())
    cp = cool['closing_phase']
    assert cp['declaration_sha256'] == _sha(cp) and cp['floating_point_floor']['exemption'] == 'pending_user_approval'
    # the floor is measured from identical-state pairs; the random-sign model is report-only, with no fallback EPS
    fp = cp['floating_point_floor']
    assert '--noise-pair' in fp['definition'] and 'REPORT ONLY' in fp['noise_model'] and 'eps' not in fp
    assert cp['revisions'][0]['superseded_declaration_sha256'] == '8c9b07958d84928021c7902f15daef0cf34348a57e9b1bbf555b97f01df65bdd'
    assert cool['thresholds']['cumulative_radiation_energy_relative'] == 0.001          # gates unchanged
    freq = json.loads((ROOT / 'benchmarks/cmfgen_cool_frequency_plan.json').read_text())
    assert freq['declaration_sha256'] == _sha(freq)
    assert freq['settings'] == json.loads((ROOT / 'benchmarks/cmfgen_frequency_plan.json').read_text())['settings']
    assert len(freq['settings']) == 9 and freq['frequencies']['expected_NCF'] == 85359
    st = freq['refinement_profile_criteria']['status']
    assert st['cumulative_heating_over_Lscale'] == st['luminosity_lab_over_Lscale'] == 'pending_user_approval'
