"""V1 plan P3 (Report 26 review): the closing chain refuses physics that differs from the converged solve's."""
import hashlib, json, struct, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import closing_chain_check as ccc  # noqa: E402
import verify_cool_state as vcs  # noqa: E402


def _pass(tmp, exe, two_phot):
    st = tmp / 'state'; st.mkdir()
    env = {'QSTARS_REPLACE_FORCE': str(tmp / 'replace.txt')}
    if two_phot:
        env['QSTARS_TWO_PHOT_J'] = str(st / 'QSTARS_RJ.bin')
    (st / 'run_manifest.json').write_text(json.dumps(dict(executable=str(exe),
                                                          executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                                                          extra_env=env)))
    (st / 'VADAT').write_text('T [INC_AD]\nT [INC_ADV]\nF [FIX_T]\nF [STOP_IF_MW]\nUSE_RAD [TWO_METH]\n0 [QS_LIN_FAC]\n')
    x0 = tmp / 'x0'; x0.mkdir()
    (x0 / 'VADAT').write_text('F [INC_AD]\nF [INC_ADV]\nT [FIX_T]\nT [STOP_IF_MW]\n')
    (tmp / 'replace.txt').write_text('1 0 1\n')
    return st, x0


def test_closing_chain_check_accepts_the_solve_physics_and_names_every_mismatch(tmp_path):
    exe = tmp_path / 'b36.exe'; exe.write_bytes(b'build 36')
    old = tmp_path / 'b31.exe'; old.write_bytes(b'build 31')
    st, x0 = _pass(tmp_path, exe, True)
    full = [('INC_AD', 'T'), ('INC_ADV', 'T'), ('STOP_IF_MW', 'F')]
    ok, bad = ccc.check(st, x0, exe, tmp_path / 'replace.txt', True, full, True)
    assert not bad and len(ok) == 7
    # the pre-Report-26 default chain: build 31, no stored two-photon J, STOP_IF_MW not carried, T fixed
    ok, bad = ccc.check(st, x0, old, tmp_path / 'replace.txt', False, full[:2], False)
    msg = '\n'.join(bad)
    assert len(bad) == 4 and 'executable' in msg and 'two-photon' in msg and 'STOP_IF_MW' in msg and 'FIX_T' in msg


def test_two_photon_selfconsistency_reports_the_relative_change_of_the_stored_j(tmp_path):
    def write(p, js):
        with open(p, 'wb') as f:
            f.write(b'QSRJ0001' + struct.pack('<i', 2))
            for ml, j in enumerate(js, 1):
                f.write(struct.pack('<id', ml, 10.0 - ml) + struct.pack('<2d', *j))
    write(tmp_path / 'a.bin', [(1.0, 2.0), (4.0, 8.0)])
    write(tmp_path / 'b.bin', [(1.0, 2.0), (4.0, 8.8)])
    r = vcs.two_photon_selfconsistency(tmp_path / 'a.bin', tmp_path / 'b.bin')
    assert r['same_grid'] and np.isclose(r['max_rel'], 0.1) and r['median_rel'] == 0.0
