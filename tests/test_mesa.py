"""Opt-in compiled MESA integration tests; no mocked microphysics.

The memo unit test needs no library and runs in the default suite.
"""
import os
import types
from collections import OrderedDict
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


def test_lru_memo_hits_order_eviction_and_failures():
    # Only the memo of MesaMicrophysics._cached; no microphysics is involved.
    calls=[]
    def native(key):
        def call():
            calls.append(key);return [key]
        return call
    def failing():
        calls.append('failed');raise ValueError('native failure')
    host=types.SimpleNamespace(cache_entries=8)
    cache=OrderedDict();memo=lambda key,f=None:MesaMicrophysics._cached(host,cache,key,f or native(key))
    stored=[memo(k) for k in range(9)]
    assert calls==list(range(9)) and list(cache)==list(range(1,9))
    # The ninth distinct key evicts the oldest; a hit returns the stored object.
    assert memo(8) is stored[8] and memo(1) is stored[1] and calls==list(range(9))
    assert list(cache)==[2,3,4,5,6,7,8,1]
    memo(9)
    assert list(cache)==[3,4,5,6,7,8,1,9] and calls[-1]==9
    assert memo(0) is not stored[0] and calls[-1]==0 and 3 not in cache
    before=list(cache);n=len(calls)
    with pytest.raises(ValueError,match='native failure'):memo('bad',failing)
    assert list(cache)==before and 'bad' not in cache and calls[n:]==['failed']
    assert memo('bad')==['bad'] and calls[-1]=='bad' and list(cache)==before[1:]+['bad']


def test_raw_lru_memo_is_bitwise_neutral(mesa):
    # Any bridge ABI: raw() is shared by every MesaMicrophysics user.
    rng=np.random.default_rng(7)
    sets=[(10**rng.uniform(-12,-2,40),np.maximum(10**rng.uniform(3.7,6.5,40),7000.)) for _ in range(11)]
    def call(k,complex_step):
        rho,T=sets[k]
        return mesa.raw(rho,T*np.exp(1e-30j) if complex_step else T)
    reference={}
    for k in range(len(sets)):
        for c in (False,True):
            mesa._cache.clear();reference[(k,c)]=call(k,c).copy()
    mesa._cache.clear()
    order=[0,1,0,2,0,1,3,4,5,6,7,8,9,10,0,1,10,9,0]
    for k in order:
        for c in (False,True):
            got=call(k,c);want=reference[(k,c)]
            assert got.dtype==want.dtype and got.tobytes()==want.tobytes()
        assert len(mesa._cache)<=mesa.cache_entries==8
    # Least recently used entry is the one evicted; a hit returns the stored arrays.
    keys=list(mesa._cache)
    mesa.raw(*sets[order[-1]])
    assert list(mesa._cache)[-1]==keys[-1]
    first=mesa.raw(*sets[2]);oldest=list(mesa._cache)[0]
    assert np.shares_memory(first,mesa.raw(*sets[2]))
    assert oldest in mesa._cache and len(mesa._cache)==8
    mesa.raw(*sets[3]);mesa.raw(*sets[4])
    assert oldest not in mesa._cache


def test_mesa_ionization_and_thermodynamic_consistency(mesa):
    rho=np.array([1e-9,1e-9,1e-7]);T=np.array([5000.,1e4,1e6])
    q=mesa.evaluate(rho,T)
    assert q['free_e'][0]<.001 and .1<q['free_e'][1]<.8 and q['free_e'][2]>.8
    assert max(mesa.consistency(rho,T))<.001
    for key in ('u','p'):
        np.testing.assert_allclose(mesa.temperature(rho,q[key],key),T,rtol=1e-9)


def test_mesa_blend_counterexample_preserved(mesa):
    assert mesa.consistency(1e-4,6811.544788207988)>.1


def test_dilute_gas_pressure_counterexample_not_hidden_by_radiation(mesa):
    rho,T=1e-9,1e6
    assert mesa.consistency(rho,T)<1e-3
    screen=mesa.gas_pressure_limit(rho,T)
    assert screen['applicable'] and screen['coupling']<.001
    assert screen['screened_relative_error']>.07
    assert mesa.gas_pressure_limit(1e-7,T)['screened_relative_error']<.001


def test_external_eos_opacity_same_state_same_answer(mesa):
    if mesa.abi_version<5:pytest.skip('ABI 5 external EOS opacity required')
    rho=np.array([1e-14,1e-9,1e-7]);T=np.array([8000.,6000.,1e6])
    q=mesa.evaluate(rho,T)
    np.testing.assert_allclose(mesa.opacity_external_eos(rho,T,q['free_e'],q['eta']),mesa.opacity(rho,T),rtol=2e-14)


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


def test_native_atmosphere_against_analytic_grey_pressure(mesa):
    if mesa.abi_version<4:pytest.skip('ABI 4 atmosphere reference required')
    teff,g,kappa,tau=8000.,3.,.34,2/3
    a=mesa.native_atmosphere(teff,g,tau,mode='fixed',kappa=kappa)
    exact=tau*g/kappa+mesa.radiation_constant*teff**4/6
    assert np.isclose(a['pressure'],exact,rtol=1e-12)
    assert np.isclose(a['temperature'],teff,rtol=1e-12)
    b=mesa.native_atmosphere(teff,g,tau,mode='varying_constant_opacity',kappa=kappa)
    # Retain this counterexample: a converged native routine is not sufficient
    # to certify the vacuum boundary appropriate to a low-gravity atmosphere.
    assert b['gas_pressure']/a['gas_pressure']>70


def test_mesa_native_pressure_inversion(mesa):
    if mesa.abi_version<4:pytest.skip('ABI 4 inversion required')
    for rho,T in ((1e-12,8000.),(1e-9,5000.),(1e-7,1e6)):
        q=mesa.evaluate(rho,T)
        assert np.isclose(mesa.density_from_gas_pressure(q['pg'],T),rho,rtol=2e-7)
    with pytest.raises(ValueError):mesa.density_from_gas_pressure(-1.,8000.)


@pytest.mark.parametrize('model',['ML1','Henyey'])
def test_native_full_mlt_bounds_and_derivatives(mesa,model):
    if mesa.abi_version<5:pytest.skip('ABI 5 MLT derivatives required')
    rho,T=1e-7,1e6;q=mesa.evaluate(rho,T)
    inputs=np.array([q['chiT'],q['chiRho'],q['cp'],1.,rho,q['p'],T,.34,.3,q['grad_ad']])
    result=mesa.mlt(inputs,model)
    assert q['grad_ad']<result[0]<inputs[8] and 0<result[1]<q['cs']
    for k in (2,3,4,5,6,7,8):
        h=1e-4;plus=inputs.copy();minus=inputs.copy()
        plus[k]*=np.exp(h);minus[k]*=np.exp(-h)
        finite=(mesa.mlt(plus,model)-mesa.mlt(minus,model))/(2*h)
        perturbed=inputs.astype(complex);perturbed[k]*=np.exp(1e-25j)
        derivative=mesa.mlt(perturbed,model).imag/1e-25
        np.testing.assert_allclose(derivative,finite,rtol=3e-4,atol=1e-8)
    inputs[8]=q['grad_ad']/2
    stable=mesa.mlt(inputs,model)
    assert stable[0]==inputs[8] and stable[1]==0
