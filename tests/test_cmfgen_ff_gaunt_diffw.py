"""Build 35: the Gaunt-factor temperature term of DIFFW(NT) (diffusion inner boundary) in the Jacobian-vector tangent, and
the free-free resonance guard. Record consistency and the shape of the patch; the native validation is in
results/ir_dT_b34 (README section 8, r2_ana.py)."""
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

CM = Path(__file__).resolve().parents[1] / 'backends' / 'cmfgen'
MAN = CM / 'ff_gaunt_dT_diffw_patch_manifest.json'
FILES = ['new_main/cmfgen_sub.f', 'new_main/qs_jv.f90']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _added(patch, fname):
    """Lines added to FNAME by the unified diff PATCH, and lines removed from it."""
    add, rem, cur = [], [], None
    for ln in patch.splitlines():
        if ln.startswith('+++ '):
            cur = ln[4:].split('\t')[0]
        elif ln.startswith('--- '):
            continue
        elif cur and cur.endswith(fname):
            if ln.startswith('+'):
                add.append(ln[1:])
            elif ln.startswith('-'):
                rem.append(ln[1:])
    return add, rem


def test_manifest_layers_on_build_34_and_record_script():
    man = json.loads(MAN.read_text())
    assert man['build'] == 'build_35' and man['base_release'] == 'build_34'
    assert sha(CM / man['patch_file']) == man['patch_sha256']
    assert man['files'] == FILES and set(man['source_hashes']) == set(FILES)
    b34 = json.loads((CM / 'ff_gaunt_dT_patch_manifest.json').read_text())['source_hashes']
    assert man['source_hashes']['new_main/qs_jv.f90']['base'] == b34['new_main/qs_jv.f90']['patched']
    # cmfgen_sub.f: the base is the last patched hash of the earlier layers (build 32's two-photon patch)
    tp = json.loads((CM / 'two_phot_j_patch_manifest.json').read_text())['source_hashes']
    assert man['source_hashes']['new_main/cmfgen_sub.f']['base'] == tp['new_main/cmfgen_sub.f']['patched']
    rec = (CM / 'record_build.py').read_text()
    assert rec.count("'ff_gaunt_dT_diffw_patch_manifest.json'") == 2 and "'qstars_ff_gaunt_dT_diffw'" in rec


def test_patch_only_adds_and_cmfgen_sub_gets_three_read_only_calls():
    patch = (CM / json.loads(MAN.read_text())['patch_file']).read_text()
    add, rem = _added(patch, 'cmfgen_sub.f')
    assert rem == []
    calls = [ln.strip() for ln in add if not ln.startswith('!')]
    assert calls == ['CALL QS_JV_DTDR_GFF(0,0.0_LDP,0.0_LDP,0.0_LDP,1.0_LDP,0.0_LDP,ND)',
                     'CALL QS_JV_DTDR_GFF(1,(FL),(CONT_FREQ),(T3),(CHI(ND)),0.0_LDP,ND)',
                     'CALL QS_JV_DTDR_GFF(2,0.0_LDP,0.0_LDP,0.0_LDP,1.0_LDP,(T1),ND)']
    # fixed form: statements start with a tab; no line starts with 'D' (-fd-lines-as-comments)
    assert all(ln.startswith('\t') or ln.startswith('!') for ln in add)
    add, rem = _added(patch, 'qs_jv.f90')
    assert rem == []
    # the tangent's dDTDR gets the term only for directions with dT(ND) /= 0, in both sweep paths
    assert sum('ddtdr=ddtdr+dwg*dpop(nt,nd,k)' in ln for ln in add) == 1
    assert sum('ddtdr_=ddtdr_+dwg*dpop(nt,nd,kk)' in ln for ln in add) == 1
    # the hook writes only qs_jv_state variables (no native array)
    body = '\n'.join(add)
    hook = body[body.index('subroutine qs_jv_dtdr_gff'):body.index('end subroutine qs_jv_dtdr_gff')]
    assert re.search(r'use qs_jv_state, only: dtg_sum,dtg_sum_all,diffw_gff,diffw_gff_all,dtg_nfreq,dtg_nres', hook)
    assigned = set(re.findall(r'^\s*(?:case\(\d\)\s*)?([a-z_]+)(?:\([^)]*\))?\s*=(?!=)', hook, re.M))
    assert assigned <= {'dtg_sum', 'dtg_sum_all', 'diffw_gff', 'diffw_gff_all', 'dtg_nfreq', 'dtg_nres', 'cf_last', 'dc', 'de',
                        'z', 'eqion', 'pop_sum', 'emh', 'al', 't4', 'dcf'}, assigned


def test_patch_reproduces_build_35_sources(tmp_path):
    man = json.loads(MAN.read_text())
    base = CM / 'build_34' / 'cur_cmf'
    if not all((base / f).is_file() for f in FILES) or shutil.which('patch') is None:
        pytest.skip('local build_34 sources (ignored build directory) or patch(1) not available')
    for f in FILES:
        assert sha(base / f) == man['source_hashes'][f]['base']
        d = tmp_path / 'cur_cmf' / f
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(base / f, d)
    subprocess.run(['patch', '-p1', '-i', str(CM / man['patch_file'])], cwd=tmp_path, check=True, capture_output=True)
    for f in FILES:
        assert sha(tmp_path / 'cur_cmf' / f) == man['source_hashes'][f]['patched']
        built = CM / 'build_35' / 'cur_cmf' / f
        if built.is_file():
            assert sha(built) == man['source_hashes'][f]['patched']


def test_build_35_capability_record_and_fingerprint():
    root = CM / 'build_35' / 'cur_cmf'
    exe = root / 'exe' / 'cmfgen_dev.exe'
    cap = exe.with_suffix('.capabilities.json')
    if not (exe.is_file() and cap.is_file()):
        pytest.skip('local build_35 executable not present')
    rec = json.loads(cap.read_text())
    assert rec['executable_sha256'] == sha(exe)
    flags = {k: v for k, v in rec.items() if k.startswith('qstars_') or k.startswith('hminus_')}
    assert rec['qstars_ff_gaunt_dT_diffw'] is True and rec['qstars_ff_gaunt_dT'] is True and all(flags.values()), flags
    sys.path.insert(0, str(CM))
    from macho_sections import section_sha256
    rep = json.loads(MAN.read_text())['reproducibility']
    assert sha(exe) == rep['executable_sha256']
    assert section_sha256(exe) == rep['section_sha256']
    for o, h in rep['object_sha256'].items():
        assert sha(root / o) == h
