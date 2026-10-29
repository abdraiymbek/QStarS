"""Build 33 (Report 23): fast local-ALI linearization pass. Record consistency only; the numerical identity with
the full path is a native test (results/fast_ali_b33)."""
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

CM = Path(__file__).resolve().parents[1] / 'backends' / 'cmfgen'
MAN = CM / 'fast_ali_patch_manifest.json'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_manifest_matches_patch_and_record_script():
    man = json.loads(MAN.read_text())
    assert man['build'] == 'build_33' and man['base_release'] == 'build_32'
    assert sha(CM / man['patch_file']) == man['patch_sha256']
    assert set(man['files']) == set(man['source_hashes']) == {'plane/var_jrel_v5.f', 'new_main/mod_subs/do_var_cont.f'}
    patch = (CM / man['patch_file']).read_text()
    for f in man['files']:
        assert f'--- a/cur_cmf/{f}' in patch and f'+++ b/cur_cmf/{f}' in patch
    # the patch layers on the files of the earlier QS_ALI_LIN / fixed-J patches (verified chain of record_build.py)
    assert man['source_hashes']['plane/var_jrel_v5.f']['base'] == \
        json.loads((CM / 'ali_lin_patch_manifest.json').read_text())['source_hashes']['plane/var_jrel_v5.f']['patched']
    assert man['source_hashes']['new_main/mod_subs/do_var_cont.f']['base'] == \
        json.loads((CM / 'jv_fast_patch_manifest.json').read_text())['source_hashes']['new_main/mod_subs/do_var_cont.f']['patched']
    rec = (CM / 'record_build.py').read_text()
    assert rec.count("'fast_ali_patch_manifest.json'") == 2 and "'qstars_fast_local_ali'" in rec
    # no arithmetic statement of the full path is removed: every deleted line is a declaration or control line
    removed = [ln[1:].strip() for ln in patch.splitlines() if ln.startswith('-') and not ln.startswith('---')]
    assert all(('=' not in ln) or ln.startswith(('USE ', 'IF(', '1 ')) for ln in removed), removed


def test_patch_reproduces_build_33_sources(tmp_path):
    man = json.loads(MAN.read_text())
    base = CM / 'build_32' / 'cur_cmf'
    if not all((base / f).is_file() for f in man['files']) or shutil.which('patch') is None:
        pytest.skip('local build_32 sources (ignored build directory) or patch(1) not available')
    for f in man['files']:
        assert sha(base / f) == man['source_hashes'][f]['base']
        d = tmp_path / 'cur_cmf' / f
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(base / f, d)
    subprocess.run(['patch', '-p1', '-i', str(CM / man['patch_file'])], cwd=tmp_path, check=True, capture_output=True)
    for f in man['files']:
        assert sha(tmp_path / 'cur_cmf' / f) == man['source_hashes'][f]['patched']
        built = CM / 'build_33' / 'cur_cmf' / f
        if built.is_file():
            assert sha(built) == man['source_hashes'][f]['patched']


def test_build_33_capability_record():
    exe = CM / 'build_33' / 'cur_cmf' / 'exe' / 'cmfgen_dev.exe'
    cap = exe.with_suffix('.capabilities.json')
    if not (exe.is_file() and cap.is_file()):
        pytest.skip('local build_33 executable not present')
    rec = json.loads(cap.read_text())
    assert rec['executable_sha256'] == sha(exe)
    assert rec['qstars_fast_local_ali'] is True
    # the layers below it are still recognised through the hash chain
    assert rec['qstars_ali_lin'] is True and rec['qstars_jv_fast'] is True and rec['qstars_two_phot_stored_j'] is True


def test_build_33_path_independent_fingerprint():
    # the executable sha256 depends on the build directory (ld debug map); sections and changed objects do not
    import sys
    sys.path.insert(0, str(CM))
    from macho_sections import member_sha256, section_sha256
    rep = json.loads(MAN.read_text())['reproducibility']
    assert set(rep['object_sha256']) == {'var_jrel_v5.o', 'do_var_cont.o'} and '__TEXT,__text' in rep['section_sha256']
    root = CM / 'build_33' / 'cur_cmf'
    exe, lib = root / 'exe' / 'cmfgen_dev.exe', root / 'lib' / 'libsubs.a'
    if not (exe.is_file() and lib.is_file()) or shutil.which('ar') is None:
        pytest.skip('local build_33 executable/library or ar(1) not available')
    assert sha(exe) == rep['executable_sha256']
    assert section_sha256(exe) == rep['section_sha256']
    assert {m: member_sha256(lib, m) for m in rep['object_sha256']} == rep['object_sha256']
