"""Build 34: temperature derivative of the ion free-free Gaunt factor in the Jacobian-vector tangent (qs_jv.f90).
Record consistency, the copied Gaunt table, and the derivative formula; the native validation is in results/ir_dT_b34."""
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

CM = Path(__file__).resolve().parents[1] / 'backends' / 'cmfgen'
MAN = CM / 'ff_gaunt_dT_patch_manifest.json'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_manifest_matches_patch_and_record_script():
    man = json.loads(MAN.read_text())
    assert man['build'] == 'build_34' and man['base_release'] == 'build_33'
    assert sha(CM / man['patch_file']) == man['patch_sha256']
    assert man['files'] == ['new_main/qs_jv.f90'] and set(man['source_hashes']) == {'new_main/qs_jv.f90'}
    # layered on the build-32 qs_jv.f90 (two-photon patch; build 33 did not change it)
    assert man['source_hashes']['new_main/qs_jv.f90']['base'] == \
        json.loads((CM / 'two_phot_j_patch_manifest.json').read_text())['source_hashes']['new_main/qs_jv.f90']['patched']
    rec = (CM / 'record_build.py').read_text()
    assert rec.count("'ff_gaunt_dT_patch_manifest.json'") == 2 and "'qstars_ff_gaunt_dT'" in rec
    # tangent only: the patch removes no line of the build-33 file
    patch = (CM / man['patch_file']).read_text()
    removed = [ln for ln in patch.splitlines() if ln.startswith('-') and not ln.startswith('---')]
    assert removed == []


def test_patch_reproduces_build_34_source(tmp_path):
    man = json.loads(MAN.read_text())
    base = CM / 'build_33' / 'cur_cmf'
    f = man['files'][0]
    if not (base / f).is_file() or shutil.which('patch') is None:
        pytest.skip('local build_33 sources (ignored build directory) or patch(1) not available')
    assert sha(base / f) == man['source_hashes'][f]['base']
    d = tmp_path / 'cur_cmf' / f
    d.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(base / f, d)
    subprocess.run(['patch', '-p1', '-i', str(CM / man['patch_file'])], cwd=tmp_path, check=True, capture_output=True)
    assert sha(d) == man['source_hashes'][f]['patched']
    built = CM / 'build_34' / 'cur_cmf' / f
    if built.is_file():
        assert sha(built) == man['source_hashes'][f]['patched']


def _table(text, start):
    """The 12x11 DATA table that follows `start` in a Fortran source (numbers in source order)."""
    i = text.index(start)
    body = text[i:text.index('/', text.index('/', i) + 1)]
    return [float(v) for v in re.findall(r'(\d+\.\d+)_(?:LDP|ldp)', body)]


def test_gaunt_table_is_gff_vec_table():
    gff = CM / 'build_34' / 'cur_cmf' / 'subs' / 'gff.f'
    if not gff.is_file():
        pytest.skip('local build_34 sources not present')
    ref = gff.read_text()
    ref_vec = _table(ref[ref.index('SUBROUTINE GFF_VEC'):], 'DATA A')
    patch = (CM / json.loads(MAN.read_text())['patch_file']).read_text()
    new = '\n'.join(ln[1:] for ln in patch.splitlines() if ln.startswith('+') and not ln.startswith('+++'))
    ours = _table(new, 'data a')
    assert len(ref_vec) == 132 and ours == ref_vec



def _gff_and_slope(nu, T, Z, table):
    """Python transcription of GFF_VEC and of qs_gff_dt_vec (bilinear cell, clamps, slope inside the cell)."""
    a1 = np.log10(4.7994 * nu); a2 = np.log10(15.789 * Z * Z); lt = np.log10(T)
    u0 = a1 - lt; g0 = a2 - lt
    u = min(max(u0, -3.99999999), 1.499999); gb = min(max(g0, -2.99999999), 1.999999)
    iu = int(u * 2 + 9.0); ig = int(gb * 2 + 7.0)
    du = 2 * u + 9.0 - iu; dg = 2 * gb + 7.0 - ig
    a = lambda i, j: table[(i - 1) + 12 * (j - 1)]
    g = ((1 - du) * a(iu, ig) + du * a(iu + 1, ig)) * (1 - dg) + ((1 - du) * a(iu, ig + 1) + du * a(iu + 1, ig + 1)) * dg
    gdu = (a(iu + 1, ig) - a(iu, ig)) * (1 - dg) + (a(iu + 1, ig + 1) - a(iu, ig + 1)) * dg
    gdg = ((1 - du) * a(iu, ig + 1) + du * a(iu + 1, ig + 1)) - ((1 - du) * a(iu, ig) + du * a(iu + 1, ig))
    cu = 1.0 if -3.99999999 < u0 < 1.499999 else 0.0
    cg = 1.0 if -2.99999999 < g0 < 1.999999 else 0.0
    return g, -2.0 * (cu * gdu + cg * gdg) / (T * np.log(10.0)), (iu, ig, cu, cg)


def test_gaunt_slope_matches_finite_difference():
    gff = CM / 'build_34' / 'cur_cmf' / 'subs' / 'gff.f'
    if not gff.is_file():
        pytest.skip('local build_34 sources not present')
    ref = gff.read_text()
    table = _table(ref[ref.index('SUBROUTINE GFF_VEC'):], 'DATA A')
    rng = np.random.default_rng(34)
    n = 0
    for _ in range(4000):
        nu = 10 ** rng.uniform(-3, 2); T = 10 ** rng.uniform(-1, 1); Z = float(rng.integers(1, 3))
        h = 1e-7 * T
        g, s, c0 = _gff_and_slope(nu, T, Z, table)
        gp, _, cp = _gff_and_slope(nu, T + h, Z, table); gm, _, cm = _gff_and_slope(nu, T - h, Z, table)
        fd = (gp - gm) / (2 * h)
        # skip the measure-zero neighbourhoods of cell boundaries/clamps, where the table is not differentiable
        if not (cp == c0 == cm):
            continue
        assert abs(fd - s) <= 1e-6 * abs(s) + 1e-7 * g / T, (nu, T, Z, fd, s)   # FD round-off ~ 1e-16 g / h
        n += 1
    assert n > 3900


def test_build_34_capability_record_and_fingerprint():
    root = CM / 'build_34' / 'cur_cmf'
    exe = root / 'exe' / 'cmfgen_dev.exe'
    cap = exe.with_suffix('.capabilities.json')
    if not (exe.is_file() and cap.is_file()):
        pytest.skip('local build_34 executable not present')
    rec = json.loads(cap.read_text())
    assert rec['executable_sha256'] == sha(exe)
    assert rec['qstars_ff_gaunt_dT'] is True
    assert rec['qstars_fast_local_ali'] is True and rec['qstars_two_phot_stored_j'] is True and rec['qstars_jv_dissolution'] is True
    import sys
    sys.path.insert(0, str(CM))
    from macho_sections import section_sha256
    rep = json.loads(MAN.read_text())['reproducibility']
    assert sha(exe) == rep['executable_sha256']
    assert section_sha256(exe) == rep['section_sha256']
    assert sha(root / 'new_main' / 'qs_jv.o') == rep['object_sha256']['new_main/qs_jv.o']
