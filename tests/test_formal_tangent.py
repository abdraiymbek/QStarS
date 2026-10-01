"""Native oracle test: actual CMF formal transfer, not a mocked Jacobian."""
import importlib.util
import os
from pathlib import Path
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_native_formal_tangent_against_transfer_reruns(tmp_path):
    native = Path(os.environ.get('QSTARS_CMFGEN_ROOT',
                                ROOT / 'backends/cmfgen/build_17/cur_cmf'))
    if not shutil.which('gfortran') or not (native / 'lib/libsubs.a').exists():
        pytest.skip('Initialized compiler and recorded CMFGEN formal-ALI libraries required')
    spec = importlib.util.spec_from_file_location(
        'audit_formal_tangent', ROOT / 'scripts/audit_formal_tangent.py')
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    result = driver.audit(native, tmp_path)
    assert result['pass'], result['gates']
    assert len(result['cases']) == 7
    assert result['finite_difference_products_in_tangent'] == 0
    # The source's blue derivative must influence the red part of the profile.
    pulse = result['cases']['inserted_front']['blue_emissivity'][-1]
    assert pulse['dropped_history_error_1'] > 1e-3
    assert pulse['ray_points_inserted']
    # Stress-case face FD must use the declared branch-preserving step, and
    # the audit must retain the larger-step failure instead of hiding it.
    mixed = result['cases']['thick_absorption']['mixed_profile']
    assert mixed[0]['error_moment_5'] > 1e-2
    gate = next(r for r in mixed if r['gate_evaluation'])
    assert gate['epsilon'] == 1e-8
    assert gate['error_moment_5'] < 1e-6
