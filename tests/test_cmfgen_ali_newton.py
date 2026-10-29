"""Build 36: Newton fallback for the coherent-scattering fixed point of the QS_FRM_ALI residual pass (comp_j_blank.f).
Record consistency and the shape of the patch; the native validation is in results/nk167_trap/b36_tests (Report 24)."""
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

CM = Path(__file__).resolve().parents[1] / 'backends' / 'cmfgen'
MAN = CM / 'ali_newton_patch_manifest.json'
F = 'new_main/mod_subs/comp_j_blank.f'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _diff(patch):
    add, rem = [], []
    for ln in patch.splitlines():
        if ln.startswith('+++ ') or ln.startswith('--- '):
            continue
        if ln.startswith('+'):
            add.append(ln[1:])
        elif ln.startswith('-'):
            rem.append(ln[1:])
    return add, rem


def test_manifest_layers_on_the_last_comp_j_blank_layer_and_record_script():
    man = json.loads(MAN.read_text())
    assert man['build'] == 'build_36' and man['base_release'] == 'build_35'
    assert sha(CM / man['patch_file']) == man['patch_sha256']
    assert man['files'] == [F] and set(man['source_hashes']) == {F}
    # the base is the patched hash of an earlier layer of this file (no unrecorded edit in between)
    rows = [json.loads(p.read_text()).get('source_hashes', {}) for p in CM.glob('*_patch_manifest.json') if p != MAN]
    patched = {r[F]['patched'] for r in rows if F in r}
    assert man['source_hashes'][F]['base'] in patched
    rec = (CM / 'record_build.py').read_text()
    assert rec.count("'ali_newton_patch_manifest.json'") == 2 and "'qstars_ali_newton_fallback'" in rec


def test_patch_keeps_the_unconverged_report_and_adds_a_guarded_newton_block():
    add, rem = _diff((CM / json.loads(MAN.read_text())['patch_file']).read_text())
    # the only removed statements are the unconverged count and message, re-added (inside IF(.NOT. QSN_OK))
    assert [r.strip() for r in rem] == ['QS_N_UNCONV=QS_N_UNCONV+1',
                                         "WRITE(LUER,*)'QS_FRM_ALI: scattering not converged',FL,QS_RES"]
    body = [a.strip() for a in add]
    assert 'IF(.NOT. QSN_OK)THEN' in body and all(r.strip() in body for r in rem)
    # accepted only below the ALI tolerance, with finite nonnegative formal J and positive iterates
    assert 'IF(QSN_RES .LT. QS_ALI_TOL)THEN' in body
    assert 'IF(.NOT. ALL(TC(1:ND) .GE. 0.0_LDP .AND. TC(1:ND) .LE. HUGE(T1)))THEN' in body
    assert any(b.startswith('DO WHILE(ANY(QSN_J(1:ND)+QSN_LAM*QSN_DJ(1:ND) .LE. 0.0_LDP)') for b in body)
    # on failure TC and the ray store are recomputed from the last ALI iterate RJ (the build-35 result)
    assert body.count('TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)') == 1
    # opacities, emissivities and the negative-opacity treatment are not touched
    assert not any(b.startswith(('CHI_CLUMP', 'ETA_CLUMP', 'CHI_SCAT_CLUMP', 'CHI(', 'ETA(')) for b in body)
    assert "CALL GET_ENVIRONMENT_VARIABLE('QSTARS_ALI_NO_NEWTON',LENGTH=QSN_LEN)" in body
    # fixed form: statements start with a tab, comments with '!'; no line starts with 'D' (-fd-lines-as-comments)
    assert all(a.startswith('\t') or a.startswith('!') for a in add)


def test_patch_reproduces_build_36_source(tmp_path):
    man = json.loads(MAN.read_text())
    base = CM / 'build_35' / 'cur_cmf' / F
    if not base.is_file() or shutil.which('patch') is None:
        pytest.skip('local build_35 sources (ignored build directory) or patch(1) not available')
    assert sha(base) == man['source_hashes'][F]['base']
    d = tmp_path / F
    d.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(base, d)
    subprocess.run(['patch', '-p1', '-i', str(CM / man['patch_file'])], cwd=tmp_path, check=True, capture_output=True)
    assert sha(d) == man['source_hashes'][F]['patched']
    built = CM / 'build_36' / 'cur_cmf' / F
    if built.is_file():
        assert sha(built) == man['source_hashes'][F]['patched']


def test_build_36_capability_record():
    exe = CM / 'build_36' / 'cur_cmf' / 'exe' / 'cmfgen_dev.exe'
    cap = exe.with_suffix('.capabilities.json')
    if not (exe.is_file() and cap.is_file()):
        pytest.skip('local build_36 executable not present')
    rec = json.loads(cap.read_text())
    assert rec['executable_sha256'] == sha(exe) == json.loads(MAN.read_text())['reproducibility']['executable_sha256']
    flags = {k: v for k, v in rec.items() if k.startswith('qstars_') or k.startswith('hminus_')}
    assert rec['qstars_ali_newton_fallback'] is True and all(flags.values()), flags
