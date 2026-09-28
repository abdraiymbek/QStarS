"""Opt-in compiled MESA integration tests; no mocked microphysics."""
import os
import numpy as np
import pytest
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.config import resolve
from qstars.evolution import run

@pytest.fixture(scope='module')
def mesa():
    library=os.environ.get('QSTARS_MESA_LIBRARY');directory=os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Set QSTARS_MESA_LIBRARY and MESA_DIR for compiled MESA tests')
    return MesaMicrophysics(Composition(),library,directory)


def test_mesa_ionization_and_thermodynamic_consistency(mesa):
    rho=np.array([1e-9,1e-9,1e-7]);T=np.array([5000.,1e4,1e6])
    q=mesa.evaluate(rho,T)
    assert q['free_e'][0]<.001 and .1<q['free_e'][1]<.8 and q['free_e'][2]>.8
    assert max(mesa.consistency(rho,T))<.001
    for key in ('u','p'):
        np.testing.assert_allclose(mesa.temperature(rho,q[key],key),T,rtol=1e-9)


def test_mesa_blend_counterexample_preserved(mesa):
    assert mesa.consistency(1e-4,6811.544788207988)>.1


def test_mesa_derivatives_drive_complex_jacobian(mesa):
    rho=np.array([1e-8,1e-7]);T=np.array([2e6,1e6]);h=1e-4
    for key in ('p','u'):
        finite=(mesa.evaluate(rho,T*np.exp(h))[key]-mesa.evaluate(rho,T*np.exp(-h))[key])/(2*h)
        complex_derivative=mesa.evaluate(rho,T*np.exp(1e-25j))[key].imag/1e-25
        np.testing.assert_allclose(complex_derivative,finite,rtol=2e-5)


def test_mesa_conservative_thermal_evolution(mesa,tmp_path):
    c=resolve({'name':'mesa_regression','eos':'mesa','mesa':{'library_path':os.environ['QSTARS_MESA_LIBRARY'],'mesa_dir':os.environ['MESA_DIR']},
       'opacity':{'model':'mesa'},'numerics':{'zones':32,'mesh':'log_radius','initial_dt_years':1e-4,'max_dt_years':1e-4},
       'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},
       'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4},
       'feeding':{'rate_msun_per_year':.01,'thermodynamics':'surface'},'wind':{'model':'prescribed','rate_msun_per_year':.008},
       'evolution':{'end_time_years':.002}})
    p,r=run(c,tmp_path)
    assert r['status']=='completed_requested_interval'
    assert abs(r['conservation']['energy_residual_relative'])<1e-9
    assert (p/'microphysics_build.json').exists()
