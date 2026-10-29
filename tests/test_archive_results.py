"""Device-transfer tooling: lossless packing of results/ units (archive_results.py) and the retention rules of long
runs (retain_runs.py). Synthetic trees only; the packing tests need zstd on PATH."""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
needs_zstd = pytest.mark.skipif(shutil.which('zstd') is None, reason='zstd not on PATH')


def _script(name):
    sys.path.insert(0, str(ROOT / 'scripts'))
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture
def ar(tmp_path, monkeypatch):
    m = _script('archive_results')
    monkeypatch.setattr(m, 'RESULTS', tmp_path / 'results')
    monkeypatch.setattr(m, 'ARCHIVE', tmp_path / 'results' / '_archive')
    return m


def make_run(root):
    d = root / 'run_a'
    (d / 'pass1' / 'empty').mkdir(parents=True)
    (d / 'pass1' / 'OUTGEN').write_text('log line\n' * 50)
    (d / 'pass1' / 'run_manifest.json').write_text('{"ok": true}\n')
    (d / 'pass1' / 'QSTARS_BA_TRACE.bin').write_bytes(os.urandom(3 << 20))
    (d / 'pass1' / 'BCSCRATCH0001').write_bytes(b'\0' * 4096)
    (d / 'pass1' / 'HI_IN').write_text('1.0 2.0\n')
    (d / 'pass1' / 'star').write_bytes(b'\xcf\xfa' * 10)
    (d / 'pass1' / 'tool').write_bytes(os.urandom(2 << 20))
    os.chmod(d / 'pass1' / 'tool', 0o755)
    os.symlink('/nonexistent/atomic/col_data', d / 'pass1' / 'HI_COL_DATA')
    os.link(d / 'pass1' / 'OUTGEN', d / 'pass1' / 'OUTGEN_copy')
    return d


def tree(d):
    out = {}
    for p in sorted(d.rglob('*')):
        st = p.lstat()
        key = str(p.relative_to(d))
        if p.is_symlink():
            out[key] = ('l', os.readlink(p))
        elif p.is_dir():
            out[key] = ('d', st.st_mode & 0o7777, st.st_mtime)
        else:
            out[key] = ('f', p.read_bytes(), st.st_mode & 0o7777, st.st_mtime)
    return out


@needs_zstd
def test_pack_delete_restore_is_byte_exact(ar, tmp_path):
    src = make_run(ar.RESULTS)
    before = tree(src)
    e = ar.pack('run_a', 2, True, None)
    assert not src.exists() and e['source_removed'] and e['symlinks'] == 1
    ar.check(ar.paths('run_a')[0], ar.read_manifest('run_a'), False)
    z = ar.zstd_in(ar.paths('run_a')[0])
    subprocess.run(['tar', '-xpf', '-', '-C', str(ar.RESULTS)], stdin=z.stdout, check=True)
    assert z.wait() == 0
    assert tree(src) == before
    assert os.stat(src / 'pass1' / 'OUTGEN').st_ino == os.stat(src / 'pass1' / 'OUTGEN_copy').st_ino


@needs_zstd
def test_record_holds_only_the_record_rule(ar, tmp_path):
    make_run(ar.RESULTS)
    ar.pack('run_a', 2, False, None)
    dest = tmp_path / 'rec'
    dest.mkdir()
    z = ar.zstd_in(ar.paths('run_a')[1])
    subprocess.run(['tar', '-xpf', '-', '-C', str(dest)], stdin=z.stdout, check=True)
    names = {p.name for p in (dest / 'run_a' / 'pass1').iterdir()}
    assert names == {'OUTGEN', 'OUTGEN_copy', 'run_manifest.json', 'empty'}


@needs_zstd
def test_corrupt_archive_is_refused(ar):
    make_run(ar.RESULTS)
    ar.pack('run_a', 2, False, None)
    m = ar.read_manifest('run_a')
    m[[e['p'] for e in m].index('run_a/pass1/OUTGEN')]['h'] = '0' * 64
    with pytest.raises(ar.ArchiveError, match='content differs'):
        ar.check(ar.paths('run_a')[0], m, False)


def test_units_skip_the_working_set(ar):
    for rel in ('old_run/p', 'nk167_trap/fixedT_b35/it01_trial_l1', 'nk167_trap/fixedT_b35/branch_x/p',
                'nk167_trap/free_ali/p', 'nk167_trap/exe_b36', 'nk167_adv/transfer', 'nk167_adv/it01_trial_l1',
                'cmfgen_ali_lambda167'):
        (ar.RESULTS / rel).mkdir(parents=True)
    assert ar.units() == ['nk167_adv/it01_trial_l1', 'nk167_trap/fixedT_b35/branch_x', 'nk167_trap/free_ali',
                          'old_run']


def test_retention_protects_what_the_driver_reads(tmp_path):
    rr = _script('retain_runs')
    out = tmp_path / 'run'
    names = ['it00_base', 'start_base', 'it01_trial_l0.5', 'cols_it01_trial_l0.5', 'it01_from_it01_trial_l0.5',
             'it02_trial_l1', 'it02_trial_l0.5', 'cols_it02_trial_l1', 'it02_from_it02_trial_l1', 'it03_trial_l1',
             'it03_trial_l0', 'cols_it03_trial_l1', 'it03_from_it03_trial_l1', 'it04_trial_l1', 'it04_trial_l0.25',
             'it04_from_it04_trial_l1', 'it05_trial_l0.25', 'it05_trial_l1', 'branch_old', 'superseded']
    for n in names:
        (out / n).mkdir(parents=True)
    hist = [dict(it=k, base=str(out / b), accepted=acc) for k, b, acc in
            ((1, 'it01_trial_l0.5', 'l1'), (2, 'it02_trial_l1', 'l1'), (3, 'it03_trial_l1', 'l1'),
             (4, 'it04_trial_l1', 'l0.25'))]
    (out / 'history.json').write_text(json.dumps(hist))
    keep = rr.protected(out, 2, [])
    packable = sorted(n for n in names if not keep(n))
    # kept: the reference it01 base, the bases of the last 2 entries and the next one, iterates 4-5
    assert packable == ['branch_old', 'cols_it02_trial_l1', 'it02_from_it02_trial_l1', 'it02_trial_l0.5',
                        'it02_trial_l1', 'it03_trial_l0', 'superseded']
