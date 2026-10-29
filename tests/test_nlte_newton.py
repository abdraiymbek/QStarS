"""Offline linear algebra and file layouts used by the Report-20 Newton-Krylov driver."""
import struct
import subprocess
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from qstars.nlte_newton import (BlockTridiag, gmres, read_ba_blocks, read_fortran, read_jv_rad,
                                read_trace_record, write_directions, write_fortran)

ROOT = Path(__file__).resolve().parents[1]


def _random_blocks(rng, nd, nt):
    lo = rng.normal(size=(nd, nt, nt)) * 0.1
    up = rng.normal(size=(nd, nt, nt)) * 0.1
    di = rng.normal(size=(nd, nt, nt)) + 5 * np.eye(nt)[None]
    # widely varying row scales, as in the native rate equations
    scale = 10.0 ** rng.uniform(-20, 20, size=(nd, nt))
    lo *= scale[:, :, None]; di *= scale[:, :, None]; up *= scale[:, :, None]
    lo[0] = 0; up[-1] = 0
    return lo, di, up


def _dense(lo, di, up):
    nd, nt, _ = di.shape
    A = np.zeros((nd * nt, nd * nt))
    for k in range(nd):
        A[k*nt:(k+1)*nt, k*nt:(k+1)*nt] = di[k]
        if k > 0:
            A[k*nt:(k+1)*nt, (k-1)*nt:k*nt] = lo[k]
        if k < nd - 1:
            A[k*nt:(k+1)*nt, (k+1)*nt:(k+2)*nt] = up[k]
    return A


def test_block_tridiagonal_matches_dense():
    rng = np.random.default_rng(0)
    nd, nt = 9, 6
    lo, di, up = _random_blocks(rng, nd, nt)
    M = BlockTridiag(lo, di, up).factor()
    A = _dense(lo, di, up)
    x_true = rng.normal(size=(nt, nd))
    b = (A @ x_true.ravel(order='F')).reshape((nt, nd), order='F')
    assert np.allclose(M.matvec(x_true), b, rtol=1e-12, atol=0)
    # Row scales span 40 decades; the per-pivot equilibrated block LU must still
    # recover a well-scaled solution.
    x = M.solve(b)
    assert np.linalg.norm(x - x_true) / np.linalg.norm(x_true) < 1e-9


def test_fortran_layouts_round_trip(tmp_path):
    a = np.arange(134 * 84, dtype=float).reshape(134, 84)
    write_fortran(tmp_path / 'a.bin', a)
    # ndarray.tofile would write C order; the native reader is column-major.
    assert np.array_equal(read_fortran(tmp_path / 'a.bin', (134, 84)), a)
    raw = np.frombuffer((tmp_path / 'a.bin').read_bytes(), '<f8')
    assert raw[1] == a[1, 0]
    d = np.random.default_rng(1).normal(size=(5, 4, 3))
    write_directions(tmp_path / 'd.bin', d)
    nv, nt, nd = struct.unpack('<3i', (tmp_path / 'd.bin').read_bytes()[:12])
    assert (nv, nt, nd) == (3, 5, 4)
    assert np.array_equal(np.frombuffer((tmp_path / 'd.bin').read_bytes(), '<f8', offset=12)
                          .reshape((5, 4, 3), order='F'), d)


def test_trace_and_block_readers(tmp_path):
    nt, nd = 4, 3
    pops = np.arange(nt * nd, dtype=float).reshape(nt, nd) + 1
    vec = -pops
    with open(tmp_path / 'tr.bin', 'wb') as f:
        f.write(struct.pack('<3i', 1, nt, nd) + pops.tobytes(order='F') + vec.tobytes(order='F'))
    p, v = read_trace_record(tmp_path / 'tr.bin')
    assert np.array_equal(p, pops) and np.array_equal(v, vec)
    rng = np.random.default_rng(2)
    with open(tmp_path / 'ba.bin', 'wb') as f:
        blocks = {}
        for k in range(1, nd + 1):
            for b in (2, 1, 3):
                if (k == 1 and b == 1) or (k == nd and b == 3):
                    continue
                c = rng.normal(size=(nt, nt)); blocks[(k, b)] = c
                f.write(struct.pack('<3i', k, b, nt) + c.tobytes(order='F'))
    lo, di, up = read_ba_blocks(tmp_path / 'ba.bin', nt, nd)
    assert np.array_equal(di[1], blocks[(2, 2)]) and np.array_equal(lo[2], blocks[(3, 1)])
    assert np.array_equal(up[0], blocks[(1, 3)]) and not np.any(lo[0]) and not np.any(up[-1])
    s = rng.normal(size=(nt, nd, 2)); d = rng.normal(size=(nt, nd, 2))
    with open(tmp_path / 'rad.bin', 'wb') as f:
        f.write(struct.pack('<3i', 2, nt, nd) + s.tobytes(order='F') + d.tobytes(order='F'))
    s2, d2 = read_jv_rad(tmp_path / 'rad.bin')
    assert np.array_equal(s2, s) and np.array_equal(d2, d)


def test_gmres_recovers_solution():
    rng = np.random.default_rng(3)
    n = 40
    A = np.eye(n) + 0.3 * rng.normal(size=(n, n)) / np.sqrt(n)
    x = rng.normal(size=n)
    sol, hist = gmres(lambda v: A @ v, A @ x, n, tol=1e-12)
    assert np.linalg.norm(sol - x) / np.linalg.norm(x) < 1e-8
    assert hist[-1] < 1e-10


def test_native_linearized_tangent_matches_dual(tmp_path):
    native = Path(os.environ.get('QSTARS_CMFGEN_ROOT', ROOT / 'backends/cmfgen/build_22/cur_cmf'))
    if not shutil.which('gfortran') or not (native / 'lib/libsubs.a').exists():
        pytest.skip('Initialized compiler and recorded CMFGEN formal-ALI libraries required')
    src = ROOT / 'backends/cmfgen'
    flags = ['-cpp', '-DNDV=20', '-O2', '-fopenmp', '-ffree-line-length-none', f'-I{native / "mod"}', f'-I{tmp_path}']
    libs = [str(native / 'lib' / n) for n in ('libsubs.a', 'libmod.a', 'libmy_lapack.a', 'libmy_blas.a')]
    for f in ('formal_tangent.f90', 'formal_tangent_lin.f90'):
        subprocess.run(['gfortran', *flags, '-c', str(src / f), '-o', str(tmp_path / (f + '.o'))],
                       cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(['gfortran', *flags, str(src / 'formal_tangent_lin_probe.f90'),
                    str(tmp_path / 'formal_tangent.f90.o'), str(tmp_path / 'formal_tangent_lin.f90.o'),
                    *libs, '-framework', 'Accelerate', '-o', str(tmp_path / 'probe')],
                   cwd=tmp_path, check=True, capture_output=True)
    worst = 0.0
    for case in ('1000 100 100 1e6 1e4 T', '1000 100 20 0.001 0.0001 T', '1000 10000 0 2 0.03 F'):
        for mode in (1, 3, 4, 5, 6):
            for threads in ('1', '4'):
                r = subprocess.run([str(tmp_path / 'probe')], input=f'{case} {mode} 1e-5 DIFFUSION\n', cwd=tmp_path,
                                   env={**os.environ, 'OMP_NUM_THREADS': threads}, capture_output=True, text=True,
                                   timeout=120, check=True)
                vals = [float(l.split()[2]) for l in r.stdout.splitlines() if l.startswith('RESULT lin_vs_dual')]
                assert len(vals) == 5 and 'MISMATCH' not in r.stdout
                worst = max(worst, max(vals))
    assert worst < 1e-11


def test_split_pieces_partition():
    import importlib.util
    spec = importlib.util.spec_from_file_location('cnk', ROOT / 'scripts' / 'cool_newton_krylov.py')
    cnk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cnk)
    rng = np.random.default_rng(3)
    for nd, free_t in ((84, False), (167, True)):
        regions = cnk.split_regions(nd)
        assert regions[0][0] == 0 and regions[-1][1] == nd
        assert all(a[1] == b[0] for a, b in zip(regions[:-1], regions[1:]))
        v = rng.normal(size=(134, nd)) * 10.0 ** rng.integers(-8, 2, size=(134, nd))
        pcs = cnk.split_pieces(v, regions, cnk.split_groups(free_t), 1e-3)
        assert np.array_equal(sum(pcs), v)
        # pieces have disjoint support
        support = sum((p != 0).astype(int) for p in pcs)
        assert support.max() <= 1
    # 84-radius regions nest into the 167-radius ones and those into the 333-radius ones (depth k -> 2k-1)
    r84, r167, r333 = cnk.split_regions(84), cnk.split_regions(167), cnk.split_regions(333)
    assert [(max(2 * a - 1, 0), 2 * b - 1) for a, b in r84] == r167
    assert [(max(2 * a - 1, 0), 2 * b - 1) for a, b in r167] == r333
    v = rng.normal(size=(134, 333))
    assert np.array_equal(sum(cnk.split_pieces(v, r333, cnk.split_groups(True), 1e-3)), v)


def _driver():
    import importlib.util
    spec = importlib.util.spec_from_file_location('cnk', ROOT / 'scripts' / 'cool_newton_krylov.py')
    cnk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cnk)
    return cnk


def test_lm_step_newton_limit_and_trust_radius():
    cnk = _driver()
    rng = np.random.default_rng(7)
    shape = (134, 6)
    S = [rng.normal(size=shape) for _ in range(5)]
    A = rng.normal(size=(np.prod(shape), 5))
    b = rng.normal(size=np.prod(shape))
    w = np.ones(shape)
    # radius larger than the Newton T step: mu = 0 and the plain least-squares step
    y0 = np.linalg.lstsq(A, b, rcond=None)[0]
    s_newton = sum(yi * si for yi, si in zip(y0, S))
    s, mu, pred = cnk.lm_step((S, A, b), w, 10.0 * np.max(np.abs(s_newton[133])))
    assert mu == 0.0 and np.allclose(s, s_newton)
    assert np.isclose(pred, np.linalg.norm(A @ y0 - b))
    # binding radius: |s_T| <= radius (to bisection accuracy) and the predicted residual grows
    r = 0.2 * np.max(np.abs(s_newton[133]))
    s2, mu2, pred2 = cnk.lm_step((S, A, b), w, r)
    assert mu2 > 0 and np.max(np.abs(s2[133])) <= r * (1 + 1e-6) and pred2 >= pred


def test_release_t_removes_inherited_ba_pointer(tmp_path):
    cnk = _driver()
    (tmp_path / 'BAMATPNT').write_text('          T        BA structure successfully output\n'
                                       '          T        Temperature was held fixed\n')
    (tmp_path / 'BAMAT').write_bytes(b'\0' * 16)
    cnk.release_t(tmp_path)
    assert not (tmp_path / 'BAMATPNT').exists() and not (tmp_path / 'BAMAT').exists()
    assert (tmp_path / 'QSTARS_T_RELEASED').exists()


def test_t_homotopy_residual_and_jacobian_term_consistent():
    """R_T(lam) = lam*R_RE + (1-lam)*c*(T-Tfix)/Tfix; in fractional variables (x <- x(1-s)) the
    driver adds (1-lam)*c*T/Tfix to the T-T entry. Check that against a finite difference."""
    cnk = _driver()
    nd = 5
    rng = np.random.default_rng(11)
    it = cnk.Iterate.__new__(cnk.Iterate)
    it.hom = dict(lam=0.3, Tfix=rng.uniform(0.3, 2.0, nd), c=rng.uniform(1.0, 10.0, nd))
    x = rng.uniform(0.5, 1.5, (134, nd)); x[133] = it.hom['Tfix'] * rng.uniform(0.8, 1.2, nd)
    R = rng.normal(size=(134, nd))
    h = 1e-7
    for d in range(nd):
        xp = x.copy(); xp[133, d] *= (1 - h)           # s_T = +h
        xm = x.copy(); xm[133, d] *= (1 + h)
        fd = (it.hom_resid(R, xm)[133, d] - it.hom_resid(R, xp)[133, d]) / (2 * h)   # J (x o s) convention
        term = (1 - it.hom['lam']) * it.hom['c'][d] * x[133, d] / it.hom['Tfix'][d]
        assert abs(fd - term) <= 1e-6 * abs(term)
    # lam = 1 leaves the residual unchanged
    it.hom['lam'] = 1.0
    assert np.array_equal(it.hom_resid(R, x), R)


def test_ptc_residual_jacobian_and_heating_direction():
    """Pseudo-transient continuation on T: G_T = R_RE - delta*c*(T - T_k)/T_k. In fractional variables the
    driver adds -delta*c (x T/T_k) to the T-T entry; check it against a finite difference, that the term
    vanishes at the base (fixed points are roots of the unmodified residual) and that, with the term
    dominant, net heating (R_T > 0) raises T."""
    cnk = _driver()
    nd = 5
    rng = np.random.default_rng(12)
    it = cnk.Iterate.__new__(cnk.Iterate)
    it.hom = dict(mode='ptc', lam=1.0, delta=0.2)
    it.ptc_c = rng.uniform(1.0, 10.0, nd)
    it.ptc_Tref = rng.uniform(0.3, 2.0, nd)
    x = rng.uniform(0.5, 1.5, (134, nd)); x[133] = it.ptc_Tref * rng.uniform(0.8, 1.2, nd)
    R = rng.normal(size=(134, nd))
    h = 1e-7
    for d in range(nd):
        xp = x.copy(); xp[133, d] *= (1 - h)           # s_T = +h
        xm = x.copy(); xm[133, d] *= (1 + h)
        fd = (it.hom_resid(R, xm)[133, d] - it.hom_resid(R, xp)[133, d]) / (2 * h)
        term = -it.hom['delta'] * it.ptc_c[d] * x[133, d] / it.ptc_Tref[d]
        assert abs(fd - term) <= 1e-6 * abs(term)
    xb = x.copy(); xb[133] = it.ptc_Tref
    assert np.array_equal(it.hom_resid(R, xb), R)
    # dominant term: s_T = R_T / (-delta*c) < 0 for net heating, so T(1 - s_T) > T
    s_T = 1.0 / (-it.hom['delta'] * it.ptc_c[0])
    assert xb[133, 0] * (1 - s_T) > xb[133, 0]


def test_pass_trace_single_record_and_ali_stats(tmp_path):
    """Report 23: a residual trace with a record inherited from the restart is refused; ALI records written by
    the pass are those after the inherited lines, and an unconverged frequency is logged, never fatal."""
    cnk = _driver()
    nt, nd = 3, 2
    rec = struct.pack('<3i', 1, nt, nd) + np.ones(2 * nt * nd).tobytes()
    (tmp_path / 'one.bin').write_bytes(rec)
    cnk.check_one_record(tmp_path / 'one.bin')
    (tmp_path / 'two.bin').write_bytes(rec + rec)
    with pytest.raises(RuntimeError):
        cnk.check_one_record(tmp_path / 'two.bin')
    line = ('NCF   43314 tol  1.0000E-10 max_it   {} unconverged       {} plain_steps         0 max_res  {} '
            'formal_solutions   1146855 ng_steps    238987\n')
    (tmp_path / 'QSTARS_ALI_STATS').write_text(line.format(298, 0, '1.0000E-10') + line.format(2000, 1, '1.5624E-09'))
    st = cnk.ali_stats(tmp_path / 'QSTARS_ALI_STATS', 1)
    assert len(st) == 1 and st[0]['unconverged'] == 1 and st[0]['max_res'] == 1.5624e-09 and st[0]['NCF'] == 43314
    R = cnk.Runner(tmp_path, tmp_path, tmp_path, False, 1, tmp_path / 'log.txt')
    R.ali_check(tmp_path, 1)
    R.ali_check(tmp_path, 0)
    R.ali_check(tmp_path / 'missing', 0)          # unreadable: warning only
    log = (tmp_path / 'log.txt').read_text()
    assert log.count('ALI unconverged at 1 of 43314') == 2 and 'not parsed' in log
    (tmp_path / 'ok').mkdir()
    (tmp_path / 'ok' / 'QSTARS_ALI_STATS').write_text(line.format(298, 0, '1.0000E-10'))
    n = len(log.splitlines())
    R.ali_check(tmp_path / 'ok', 0)
    assert len((tmp_path / 'log.txt').read_text().splitlines()) == n
    # an earlier unconverged J computation is reported even when the last one converged
    (tmp_path / 'mixed').mkdir()
    (tmp_path / 'mixed' / 'QSTARS_ALI_STATS').write_text(line.format(2000, 3, '2.0000E-09') + line.format(2000, 1, '1.5624E-09')
                                                         + line.format(298, 0, '1.0000E-10'))
    R.ali_check(tmp_path / 'mixed', 0)
    last = (tmp_path / 'log.txt').read_text().splitlines()[-1]
    assert 'ALI unconverged at 3 of 43314' in last and 'max_res 2.000e-09' in last and '2 unconverged J computations out of 3' in last


def test_runner_two_phot_j_plumbing(tmp_path, monkeypatch):
    """Build 32: Runner.run(two_phot_j=...) sets env QSTARS_TWO_PHOT_J and VADAT TWO_METH=USE_RAD for that pass
    only; the default adds neither; a pass is reused only with the same stored two-photon J. The variable given
    through extra_env is the same request (same guards, no override), and the native QSTARS_TWO_PHOT_J_FAILURE
    file is reported as the cause of a stopped pass."""
    import json, sys, types
    cnk = _driver()
    calls = []
    fa = types.ModuleType('formal_ali_trial'); rc = types.ModuleType('run_component')

    def prepare(src, dest, its, tol, max_it, sets):
        Path(dest).mkdir(parents=True); calls.append(dict(sets=list(sets)))

    def run(dest, exe, timeout, threads, extra_env=None):
        calls[-1]['env'] = dict(extra_env)
        (Path(dest) / 'run_manifest.json').write_text(json.dumps(dict(extra_env=extra_env)))
        if native_stop:                                # build 32 stopping on the stored-J file: no residual
            (Path(dest) / 'QSTARS_TWO_PHOT_J_FAILURE').write_text(native_stop[0])
            return dict(error='Native run did not produce fresh RVTJ')
        (Path(dest) / 'QSTARS_RESID_TRACE.bin').write_bytes(struct.pack('<3i', 1, 3, 2) + np.ones(12).tobytes())
        return dict(error=None)
    native_stop = []
    fa.prepare = prepare; rc.run = run
    monkeypatch.setitem(sys.modules, 'formal_ali_trial', fa); monkeypatch.setitem(sys.modules, 'run_component', rc)
    R = cnk.Runner(tmp_path, tmp_path / 'cmfgen_dev.exe', tmp_path, True, 1, tmp_path / 'log.txt'); R.slim = False
    s = np.zeros((3, 2)); rj = tmp_path / 'QSTARS_RJ.bin'; rj.write_bytes(b'QSRJ0001')
    R.run(tmp_path / 'default', s)
    assert 'QSTARS_TWO_PHOT_J' not in calls[-1]['env'] and all(k != 'TWO_METH' for k, _ in calls[-1]['sets'])
    with pytest.raises(RuntimeError):                  # executable without the recorded capability (build < 32)
        R.run(tmp_path / 'on', s, two_phot_j=rj)
    (tmp_path / 'cmfgen_dev.capabilities.json').write_text(json.dumps(dict(qstars_two_phot_stored_j=True)))
    R.run(tmp_path / 'on', s, two_phot_j=rj)
    assert calls[-1]['env']['QSTARS_TWO_PHOT_J'] == str(rj.resolve()) and ('TWO_METH', 'USE_RAD') in calls[-1]['sets']
    n = len(calls)
    R.run(tmp_path / 'on', s, two_phot_j=rj); R.run(tmp_path / 'default', s)       # reuse, no new native pass
    assert len(calls) == n
    with pytest.raises(RuntimeError):
        R.run(tmp_path / 'on', s)                      # made with a stored two-photon J: not the default residual
    with pytest.raises(RuntimeError):
        R.run(tmp_path / 'default', s, two_phot_j=rj)
    with pytest.raises(FileNotFoundError):
        R.run(tmp_path / 'new', s, two_phot_j=tmp_path / 'missing.bin')
    R.vadat = [('TWO_METH', 'NOSTIM')]
    with pytest.raises(ValueError):
        R.run(tmp_path / 'new', s, two_phot_j=rj)
    with pytest.raises(ValueError):                    # the extra_env route does not skip the TWO_METH guard
        R.run(tmp_path / 'new', s, extra_env={'QSTARS_TWO_PHOT_J': str(rj)})
    R.vadat = []
    assert not (tmp_path / 'new').exists() and len(calls) == n

    # QSTARS_TWO_PHOT_J through extra_env (per call or runner-wide) is the same request as two_phot_j
    cap = tmp_path / 'cmfgen_dev.capabilities.json'; cap_text = cap.read_text(); cap.unlink()
    with pytest.raises(RuntimeError, match='qstars_two_phot_stored_j'):     # capability guard not bypassed
        R.run(tmp_path / 'ee', s, extra_env={'QSTARS_TWO_PHOT_J': str(rj)})
    cap.write_text(cap_text)
    with pytest.raises(FileNotFoundError):
        R.run(tmp_path / 'ee', s, extra_env={'QSTARS_TWO_PHOT_J': str(tmp_path / 'missing.bin')})
    R.run(tmp_path / 'ee', s, extra_env={'QSTARS_TWO_PHOT_J': str(rj), 'OTHER': '1'})
    assert calls[-1]['env']['QSTARS_TWO_PHOT_J'] == str(rj.resolve()) and ('TWO_METH', 'USE_RAD') in calls[-1]['sets']
    assert calls[-1]['env']['OTHER'] == '1'
    R.extra_env = {'QSTARS_TWO_PHOT_J': str(rj)}
    R.run(tmp_path / 'ew', s)
    assert calls[-1]['env']['QSTARS_TWO_PHOT_J'] == str(rj.resolve()) and ('TWO_METH', 'USE_RAD') in calls[-1]['sets']
    n = len(calls)
    R.run(tmp_path / 'ee', s, two_phot_j=rj); R.run(tmp_path / 'on', s)      # all three routes: one residual definition
    assert len(calls) == n
    with pytest.raises(RuntimeError):
        R.run(tmp_path / 'default', s)                 # runner-wide variable: not the default residual either
    # the requested file cannot be replaced through extra_env
    other = tmp_path / 'OTHER_RJ.bin'; other.write_bytes(b'QSRJ0001')
    with pytest.raises(ValueError, match='extra_env sets QSTARS_TWO_PHOT_J'):
        R.run(tmp_path / 'new', s, two_phot_j=other)
    R.extra_env = {}
    with pytest.raises(ValueError, match='extra_env sets QSTARS_TWO_PHOT_J'):
        R.run(tmp_path / 'new', s, two_phot_j=rj, extra_env={'QSTARS_TWO_PHOT_J': str(other)})
    assert not (tmp_path / 'new').exists() and len(calls) == n
    # an empty value is 'unset' (as in the native code) and is not handed on
    R.run(tmp_path / 'empty', s, extra_env={'QSTARS_TWO_PHOT_J': ''})
    assert 'QSTARS_TWO_PHOT_J' not in calls[-1]['env'] and all(k != 'TWO_METH' for k, _ in calls[-1]['sets'])

    # a pass stopped by the native stored-J check names QSTARS_TWO_PHOT_J_FAILURE, not just 'no residual'
    native_stop.append('stored-J frequency grid differs at ML=    1000\nfile=' + str(rj) + '\n')
    with pytest.raises(RuntimeError, match='stored two-photon J rejected: stored-J frequency grid differs at ML=    1000'):
        R.run(tmp_path / 'stopped', s, two_phot_j=rj)
    assert not (tmp_path / 'stopped' / 'QSTARS_RESID_TRACE.bin').exists()
    with pytest.raises(FileExistsError):               # the stopped pass is kept and never silently rerun
        R.run(tmp_path / 'stopped', s, two_phot_j=rj)


def test_frontpin_residual_and_jacobian_per_depth():
    """--hom-mode frontpin: pinned depths (w = 0) carry only c*(T - Tfix)/Tfix, all other depths the unmodified
    RE row; the T-T Jacobian term (1 - w)*c*T/Tfix matches a finite difference, and pin_lam scales products."""
    cnk = _driver()
    nd = 6
    rng = np.random.default_rng(13)
    it = cnk.Iterate.__new__(cnk.Iterate)
    w = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
    it.hom = dict(mode='pin', lam=1.0, w=w, Tfix=rng.uniform(0.3, 2.0, nd), c=rng.uniform(1.0, 10.0, nd))
    x = rng.uniform(0.5, 1.5, (134, nd)); x[133] = it.hom['Tfix'] * rng.uniform(0.8, 1.2, nd)
    R = rng.normal(size=(134, nd))
    G = it.hom_resid(R, x)
    assert np.array_equal(G[:133], R[:133])
    assert np.array_equal(G[133, w == 1], R[133, w == 1])
    pinned = w == 0
    assert np.allclose(G[133, pinned], it.hom['c'][pinned] * (x[133, pinned] - it.hom['Tfix'][pinned]) / it.hom['Tfix'][pinned])
    h = 1e-7
    for d in range(nd):
        xp = x.copy(); xp[133, d] *= (1 - h)
        xm = x.copy(); xm[133, d] *= (1 + h)
        fd = (it.hom_resid(R, xm)[133, d] - it.hom_resid(R, xp)[133, d]) / (2 * h)
        term = (1 - w[d]) * it.hom['c'][d] * x[133, d] / it.hom['Tfix'][d]
        assert abs(fd - term) <= 1e-6 * max(abs(term), 1e-12)
    assert np.array_equal(cnk.pin_lam(it.hom), w)
    assert cnk.pin_lam(dict(lam=0.3)) == 0.3


def test_significance_weight_ramp_is_continuous_and_default_unchanged():
    cnk = _driver()
    cnk.SIG_FRAC, cnk.TRACE_WEIGHT, cnk.ULTRA_FRAC, cnk.ULTRA_WEIGHT = 1e-10, 1e-3, 1e-20, 0.0
    frac = np.array([1e-30, 1e-20, 1e-15, 1e-12, 1e-11, 0.999e-10, 1e-10, 1.001e-10, 1e-5])
    cnk.WEIGHT_RAMP = 0.0
    w0 = cnk.significance_weights(frac)
    assert np.array_equal(w0, np.array([0.0, 0.0, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1.0, 1.0]))
    cnk.WEIGHT_RAMP = 2.0
    w = cnk.significance_weights(frac)
    assert w[0] == 0.0 and w[2] == 1e-3 and w[3] == 1e-3 and w[-1] == 1.0
    assert abs(w[4] - 10 ** -1.5) < 1e-12                       # halfway in log: 1e-3^(1/2)
    assert abs(w[6] - 1.0) < 1e-12                                # continuous at SIG_FRAC
    assert abs(w[5] - 1e-3 ** ((np.log10(1e-10) - np.log10(0.999e-10)) / 2.0)) < 1e-12
    f = np.logspace(-14, -8, 2001); ww = cnk.significance_weights(f)
    assert np.all(np.diff(ww) >= 0) and np.max(np.abs(np.diff(np.log10(ww)))) < 0.01
    cnk.WEIGHT_RAMP = 0.0
