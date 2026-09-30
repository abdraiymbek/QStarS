import struct
from pathlib import Path
import numpy as np
import pytest
from qstars.cmfgen_backend import read_spectral,vector_section,hhe_material,production_outer
from qstars.constants import C


def snapshot(path,truncate=False):
    n=3;nf=2;r=np.array([3.,2.,1.]);v=np.array([1e6,2e6,3e6]);T=np.ones(n)*1e4
    rho=np.array([2.,3.,4.]);fill=np.array([.1,.5,1.])
    J=np.array([[2.,4.,6.],[3.,5.,7.]])
    H=J*.1;K=J/3;chi=np.ones_like(J)*.3;scat=chi*.6
    eta=(chi-scat)*J # coherent scattering contributes no heat
    data=np.stack([chi,scat,eta,J,H,K],axis=1).reshape(nf,6*n)
    data=np.column_stack(([2e15,1e15],[3.,5.],data))
    content=b'QSCMF001'+struct.pack('<ii',n,nf)+np.array([r,v,T,rho,fill],dtype='<f8').tobytes()+data.astype('<f8').tobytes()
    path.write_bytes(content[:-8] if truncate else content)
    return J,H,K,rho,fill


def test_spectral_units_quadrature_and_scattering_exchange(tmp_path):
    p=tmp_path/'spectral';J,H,K,rho,fill=snapshot(p)
    s=read_spectral(p);q=s.moments(coherent_scattering=True)
    integrated=np.array([3.,5.])@J[:,::-1]
    np.testing.assert_allclose(q['E_comoving'],4*np.pi*integrated/C)
    np.testing.assert_allclose(q['heating_comoving'],0,atol=1e-15)
    np.testing.assert_allclose(q['radiative_acceleration'],.3/rho[::-1]*q['F_comoving']/C)
    # Lab-frame radiation work remains present despite zero comoving heating.
    np.testing.assert_allclose(q['energy_to_matter_lab'],q['velocity']*q['radiative_force_density']/np.sqrt(1-(q['velocity']/C)**2))
    np.testing.assert_allclose(q['density'],(rho*fill)[::-1])
    assert not q['wind_eigenvalue_solved']


def test_truncated_native_snapshot_rejected(tmp_path):
    p=tmp_path/'bad';snapshot(p,truncate=True)
    with pytest.raises(ValueError,match='Incomplete'):read_spectral(p)
    p.write_bytes(b'QSCMF001'+struct.pack('<ii',-1,2))
    with pytest.raises(ValueError,match='dimensions'):read_spectral(p)


def test_labeled_vector_cannot_consume_next_section():
    with pytest.raises(ValueError,match='Malformed'):vector_section(['foo','1.0','bar','2.0'],'foo',2)
    with pytest.raises(ValueError):vector_section(['foo','nan'],'foo',1)


def test_population_atomic_energy_keeps_ionization_separate():
    p={'HI':np.array([[2.],[3.]]),'HII':np.array([5.]),'HeI':np.array([[7.]]),
       'HeII':np.array([[11.]]),'HeIII':np.array([13.])}
    s={'temperature':np.array([1e4]),'density_clump':np.array([1e-10]),'hydrogen_density':np.array([10.]),
       'helium_density':np.array([31.]),'electron_density':np.array([42.])}
    l={'HI':{'wavenumber':np.array([0.,10.]),'ionization_wavenumber':20.},
       'HeI':{'wavenumber':np.array([0.]),'ionization_wavenumber':30.},
       'HeII':{'wavenumber':np.array([0.]),'ionization_wavenumber':40.}}
    q=hhe_material(s,p,l)
    expected=6.62607015e-27*C*(3*10+5*20+11*30+13*70)
    np.testing.assert_allclose(q['atomic_energy_density_clump'],expected)
    assert max(q['conservation_relative'].values())==0
    assert q['energy_gauge']=='neutral_H_He_ground'
    assert not q['equilibrium_derivatives_available']


def test_production_wind_remains_gated():
    with pytest.raises(RuntimeError,match='production wind remains gated'):production_outer()


def test_inelastic_scattering_is_not_silently_coherent(tmp_path):
    p=tmp_path/'spectral';snapshot(p)
    with pytest.raises(NotImplementedError,match='Inelastic'):read_spectral(p).moments(coherent_scattering=False)

def test_fortran_exponent_without_e_is_supported():
    np.testing.assert_allclose(vector_section(['foo','1.25-100 2.0D+02'],'foo',2),[1.25e-100,200],atol=0)


def runner_module():
    import importlib.util
    path=Path(__file__).resolve().parents[1]/'backends/cmfgen/run_component.py'
    spec=importlib.util.spec_from_file_location('cmfgen_runner',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_native_zero_exit_is_not_success(tmp_path):
    executable=tmp_path/'stop';executable.write_text('#!/bin/sh\nprintf "Major Error - cant read BA File\\n" > OUTGEN\nexit 0\n')
    executable.chmod(0o755);run=tmp_path/'run';run.mkdir()
    result=runner_module().run(run,executable,timeout=5)
    assert result['exit_code']==0
    assert result['status']=='numerical_failure'
    assert 'fatal error' in result['error']
    assert (run/'run_result.json').is_file()


def test_stale_native_output_rejected(tmp_path):
    import time
    (tmp_path/'RVTJ').write_text('stale')
    with pytest.raises(RuntimeError,match='fresh RVTJ'):
        runner_module().check_completion(tmp_path,time.time_ns(),False)


def test_nan_luminosity_cannot_pass_zero_correction(tmp_path):
    (tmp_path/'OUTGEN').write_text('Luminosity of star is: NaN NaN\nMaximm changes as returned by SOLVEBA_V13 is 0.0\n')
    with pytest.raises(RuntimeError,match='Nonfinite native diagnostic'):
        runner_module().check_completion(tmp_path,0,False)


def test_native_stopping_cannot_discard_largest_levels():
    from qstars.cmfgen_backend import nonlinear_corrections
    text='''Current great iteration count is 10
Maximum % increase at depth 81 is 40.3 (LAMBDA)
Maximum % decrease at depth 84 is 22.1 (LAMBDA)
Maximm changes as returned by SOLVEBA_V13 is 4.8E-10
'''
    q=nonlinear_corrections(text)[0]
    assert q['iteration']==10
    assert q['all_level_percent']==40.3
    assert q['reported_percent']==4.8e-10
    assert nonlinear_corrections('Maximm changes as returned by SOLVEBA_V13 is 0')[0]['all_level_percent'] is None


def test_raw_decrease_normalization_and_three_digit_exponent():
    from qstars.cmfgen_backend import nonlinear_corrections
    q=nonlinear_corrections('''Maximum % increase at depth 1 is 1.0-100
Maximum % decrease at depth 2 is 50.0
Maximm changes as returned by SOLVEBA_V13 is 1.0E-10''')[0]
    assert q['all_level_percent']==100.


def test_native_moment_failure_rejected_despite_finite_output(tmp_path):
    for message in ['Error in MOM_J_REL_V9: excessive iteration count.\n',
                    'Mean intensity blowing up\n',
                    'You have poor flux accuracy in the inner region of the model\n']:
        (tmp_path/'OUTGEN').write_text(message)
        with pytest.raises(RuntimeError,match='radiation iteration failed'):
            runner_module().check_completion(tmp_path,0,False)


def test_long_lived_runner_retains_imported_source_hashes(tmp_path,monkeypatch):
    mod=runner_module()
    runner_hash=mod.RUNNER_SOURCE_SHA256;adapter_hash=mod.LOADED_SOURCE_SHA256
    exe=tmp_path/'native';exe.write_text('unexecuted capacity failure');run=tmp_path/'run';run.mkdir()
    (run/'MODEL_SPEC').write_text('3000 [NP]\n')
    # Simulate changed on-disk source hashes after Python imported the code.
    monkeypatch.setattr(mod,'sha256',lambda path:'changed-on-disk')
    result=mod.run(run,exe,timeout=1)
    assert result['status']=='numerical_failure'
    assert result['runner_sha256']==runner_hash
    assert result['adapter_sha256']==adapter_hash


def test_live_nonfinite_native_run_is_stopped_and_retained(tmp_path):
    exe=tmp_path/'nan_run'
    exe.write_text('#!/bin/sh\nprintf "Luminosity of star is: NaN NaN\\n" > OUTGEN\nsleep 30\n')
    exe.chmod(0o755);run=tmp_path/'run';run.mkdir()
    result=runner_module().run(run,exe,timeout=15)
    assert result['status']=='numerical_failure'
    assert 'Nonfinite native diagnostic' in result['error']
    assert result['seconds']<10
    assert (run/'OUTGEN').read_text().endswith('NaN NaN\n')


def test_nonphysical_frequency_rejected(tmp_path):
    p=tmp_path/'spectral';snapshot(p)
    a=bytearray(p.read_bytes());struct.pack_into('<d',a,16+8*(5*3+2+6*3),-1.)
    p.write_bytes(a)
    with pytest.raises(ValueError,match='grid'):read_spectral(p)


def test_unverified_native_ray_capacity_cannot_overrun(tmp_path):
    directory=tmp_path/'run';directory.mkdir()
    (directory/'MODEL_SPEC').write_text('680 [NP]\n')
    executable=tmp_path/'native';executable.write_text('dummy')
    with pytest.raises(RuntimeError,match='verified capacity 500'):
        runner_module().check_grid_capacity(directory,executable)
    import json
    from qstars.cmfgen_backend import sha256
    executable.with_suffix('.capabilities.json').write_text(json.dumps(dict(executable_sha256=sha256(executable),explicit_native_capacity_guard=True,max_rays=2048)))
    assert runner_module().check_grid_capacity(directory,executable)==2048
    executable.write_text('changed')
    with pytest.raises(RuntimeError,match='does not match'):
        runner_module().check_grid_capacity(directory,executable)


def test_known_hminus_rayleigh_bug_is_gated(tmp_path):
    run=tmp_path/'run';run.mkdir();(run/'MODEL_SPEC').write_text('99 [NP]\n1,1,1 [H0_ISF]\n')
    (run/'VADAT').write_text('T [INC_RAY]\n');exe=tmp_path/'native';exe.write_text('dummy')
    with pytest.raises(RuntimeError,match='HI selection/stride/Jacobian'):
        runner_module().check_grid_capacity(run,exe)
    import json
    from qstars.cmfgen_backend import sha256
    exe.with_suffix('.capabilities.json').write_text(json.dumps(dict(executable_sha256=sha256(exe),hminus_rayleigh_fixed=True)))
    assert runner_module().check_grid_capacity(run,exe)==500


def test_hminus_population_labels_and_charge_energy(tmp_path):
    from qstars.cmfgen_backend import read_hhe_populations
    p=tmp_path/'RVTJ'
    text='ND: 2\n'
    for ion,label,data in [('H0','Hminus','0.1 0.2'),('HI','Hydrogen','2 3'),
                            ('HeI','Helium I','4 5'),('He2','Helium II','6 7')]:
        text+=f'{ion} populations\nNumber of {label} levels: 1\n{ion} oscillator date: test\n{data}\n'
        if ion=='HI':text+='DHYD population\n8 9\n'
        if ion=='He2':text+='DHe2 population\n10 11\n'
    p.write_text(text);q=read_hhe_populations(p)
    np.testing.assert_array_equal(q['Hminus'],[[.2,.1]])
    np.testing.assert_array_equal(q['HeIII'],[11,10])
    # Upstream fixed indices can mislabel H0 as hydrogen. Fail closed.
    p.write_text(text.replace('HI populations','H0 populations'))
    with pytest.raises(ValueError,match='ion label mismatch'):read_hhe_populations(p)


def test_continuation_guard_requires_attested_inherited_prefix(tmp_path):
    import json
    from qstars.cmfgen_backend import sha256
    mod=runner_module();initial=tmp_path/'initial_inputs';initial.mkdir()
    old='Mean intensity blowing up\nNaN\n'
    (initial/'OUTGEN').write_text(old)
    (tmp_path/'run_manifest.json').write_text(json.dumps({'inputs':{'OUTGEN':sha256(initial/'OUTGEN')}}))
    (tmp_path/'OUTGEN').write_text(old+'New finite iteration\n')
    assert mod.current_native_log(tmp_path)=='New finite iteration\n'
    mod.check_native_finiteness(tmp_path)
    (tmp_path/'OUTGEN').write_text(old+'New Mean intensity blowing up\n')
    with pytest.raises(RuntimeError,match='radiation iteration failed'):
        mod.check_completion(tmp_path,0,False)
    (initial/'OUTGEN').write_text(old+'New ')
    # An unverified prefix cannot hide a failure, even if its bytes match.
    with pytest.raises(RuntimeError,match='Nonfinite native diagnostic'):
        mod.check_native_finiteness(tmp_path)


def test_experimental_moment_failures_cannot_pass_finite_outputs(tmp_path):
    for message in ['QSTARS_NODE_MOMENT negative J 86.13',
                    'QSTARS_RAW_NH reversed frequency characteristic',
                    'QSTARS_RAW_NH singular formal H']:
        (tmp_path/'OUTGEN').write_text(message+'\n')
        with pytest.raises(RuntimeError,match='radiation iteration failed'):
            runner_module().check_completion(tmp_path,0,False)
