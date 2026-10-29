"""Opt-in compiled MESA ABI 7 tests: the Ball 2011 STARS opacity switch and its LRU memo.

Set QSTARS_MESA_LIBRARY to an ABI 7 build (backends/mesa/build_ball_abi7) and MESA_DIR.
QSTARS_MESA_ABI6_LIBRARY optionally names an ABI 6 build (backends/mesa/build_outer_2)
for the backward-compatibility check, which runs in a subprocess because one process
holds one MESA installation. No mocked microphysics. The raw() memo, shared by every
ABI 3-7 build, is tested in tests/test_mesa.py.
"""
import json
import os
import subprocess
import sys
from pathlib import Path
import numpy as np
import pytest
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def mesa():
    library=os.environ.get('QSTARS_MESA_LIBRARY');directory=os.environ.get('MESA_DIR')
    if not library or not directory:pytest.skip('Set QSTARS_MESA_LIBRARY and MESA_DIR for compiled MESA tests')
    m=MesaMicrophysics(Composition(.7,.3,0.),library,directory)
    if m.abi_version<7:pytest.skip('ABI 7 STARS opacity required')
    return m

def states(rng,logT_low,logT_high,n=400):
    logT=rng.uniform(logT_low,logT_high,n);logR=rng.uniform(-6.999,.999,n)
    return 10**(logR+3*logT-18),10**logT

def logR(rho,T):return np.log10(rho)-3*np.log10(T)+18


def test_stars_switch_bitwise_equals_paper_handle_outside_switch_interval(mesa):
    rng=np.random.default_rng(20111)
    for low,high,ends in ((np.log10(7000),7.,[7000.,1e7]),(np.log10(3000),3.74,[3000.,10**3.74])):
        rho,T=states(rng,low,high)
        for end in ends:
            rho=np.r_[rho,10**(np.array([-6.9,-3.,.9])+3*np.log10(end)-18)];T=np.r_[T,[end]*3]
        stars=mesa.opacity_stars(rho,T)
        np.testing.assert_array_equal(stars,mesa.opacity_paper_reference(rho,T,'gs98'))
        clamped,flags=mesa.opacity_stars(rho,T,clamp=True,flags=True)
        np.testing.assert_array_equal(clamped,stars)
        assert not np.any(flags)


def test_stars_switch_reproduces_opal_and_af94_table_nodes(mesa):
    # Between log T 3.76 and 7000 K STARS reads OPAL, the Nakauchi handle AF94.
    for name,low,high in (('gs98',3.76,np.log10(7000)),('lowT_af94_gn93',np.log10(3000),3.75)):
        path=Path(os.environ['MESA_DIR'])/'data/kap_data'/f'{name}_z0.0_x0.7.data'
        lines=path.read_text().splitlines();grid_logR=np.array([float(t) for t in lines[5].split()])
        grid=np.loadtxt(path,skiprows=7)
        rows=grid[(grid[:,0]>=low)&(grid[:,0]<=high)]
        columns=np.flatnonzero((grid_logR>-6.99)&(grid_logR<.99))
        assert len(rows) and len(columns)
        for row in rows:
            temp=10**row[0];density=10**(grid_logR[columns]+3*row[0]-18)
            actual=mesa.opacity_stars(density,np.full(len(columns),temp))
            assert np.max(abs(np.log10(actual)-row[columns+1]))<1e-9
    logT=np.linspace(3.77,3.84,8);T=10**logT;rho=10**(-4+3*logT-18)
    assert np.all(abs(mesa.opacity_stars(rho,T)/mesa.opacity_paper_reference(rho,T,'gs98')-1)>1e-3)


@pytest.mark.parametrize('clamp',[False,True])
def test_stars_complex_step_matches_finite_differences(mesa,clamp):
    T=np.array([1e6,1e5,2e4,1e4,8000.,10**3.755,5000.,4000.])
    rho=10**(np.array([-1.,-4.,-6.,-3.,-5.,-2.,-6.5,0.])+3*np.log10(T)-18)
    if clamp:
        # Off both edges at fixed T: the derivative of the clamped map.
        T=np.r_[T,4e4,5000.,1e5];rho=np.r_[rho,10**(-7.5+3*np.log10([4e4,5000.])-18),10**(1.5+15-18)]
        assert np.all(mesa.opacity_stars(rho,T,clamp=True,flags=True)[1][-3:]==[1,1,2])
    kappa=mesa.opacity_stars(rho,T,clamp=clamp)
    h,small=1e-6,1e-30
    for axis in ('rho','T'):
        up=mesa.opacity_stars(rho*np.exp(h) if axis=='rho' else rho,T*np.exp(h) if axis=='T' else T,clamp=clamp)
        down=mesa.opacity_stars(rho*np.exp(-h) if axis=='rho' else rho,T*np.exp(-h) if axis=='T' else T,clamp=clamp)
        finite=(np.log(up)-np.log(down))/(2*h)
        perturbed=mesa.opacity_stars(rho*np.exp(1j*small) if axis=='rho' else rho,T*np.exp(1j*small) if axis=='T' else T,clamp=clamp)
        np.testing.assert_array_equal(perturbed.real,kappa)
        derivative=perturbed.imag/(small*kappa)
        assert np.all(abs(derivative-finite)<=1e-6*np.maximum(1,abs(finite))),(axis,derivative,finite)
        if clamp and axis=='rho':assert np.all(derivative[-3:]==0)


def test_stars_off_table_errors_and_clamp_flags(mesa):
    T=np.array([4e4,5000.,1e5]);edges=np.array([-7.,-7.,1.])
    rho=10**(np.array([-7.5,-12.,4.])+3*np.log10(T)-18)
    with pytest.raises(ValueError,match='domain'):mesa.opacity_stars(rho,T)
    for i in range(3):
        with pytest.raises(ValueError,match='domain'):mesa.opacity_stars(rho[i],T[i])
    kappa,flags=mesa.opacity_stars(rho,T,clamp=True,flags=True)
    np.testing.assert_array_equal(flags,[1,1,2])
    # Nearest table value: the same as the state on the edge at the same T.
    inside=10**(edges+3*np.log10(T)-18)*np.array([1+1e-12,1+1e-12,1-1e-12])
    assert np.all((logR(inside,T)>=-7)&(logR(inside,T)<=1))
    np.testing.assert_allclose(kappa,mesa.opacity_stars(inside,T),rtol=1e-10)
    for bad in (2999.,1.01e7):
        with pytest.raises(ValueError,match='domain'):mesa.opacity_stars(1e-9,bad,clamp=True)
    for bad in (0.,-1e-9,np.nan,np.inf):
        with pytest.raises(ValueError):mesa.opacity_stars(bad,1e4,clamp=True)
    other=MesaMicrophysics(Composition(.75,.25,0.),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
    with pytest.raises(ValueError,match='X=.7'):other.opacity_stars(1e-9,1e4)


def test_stars_lru_memo_is_bitwise_neutral(mesa):
    rng=np.random.default_rng(7)
    sets=[(10**rng.uniform(-12,-2,40),10**rng.uniform(3.7,6.5,40)) for _ in range(11)]
    sets=[(rho,np.maximum(T,7000.)) for rho,T in sets]
    def call(k,complex_step):
        rho,T=sets[k]
        return mesa.opacity_stars(rho,T*np.exp(1e-30j) if complex_step else T,clamp=True)
    reference={}
    for k in range(len(sets)):
        for c in (False,True):
            mesa._stars_cache.clear();reference[(k,c)]=call(k,c).copy()
    mesa._stars_cache.clear()
    for k in [0,1,0,2,0,1,3,4,5,6,7,8,9,10,0,1,10,9,0]:
        for c in (False,True):
            # raw() calls in between use a separate memo and must not disturb this one.
            mesa.raw(*sets[k])
            got=call(k,c);want=reference[(k,c)]
            assert got.dtype==want.dtype and got.tobytes()==want.tobytes()
        assert len(mesa._stars_cache)<=mesa.cache_entries==8


def test_abi6_library_still_loads_and_agrees_bitwise(mesa):
    library=os.environ.get('QSTARS_MESA_ABI6_LIBRARY')
    if not library:pytest.skip('Set QSTARS_MESA_ABI6_LIBRARY to an ABI 6 build for the compatibility check')
    rho=np.array([1e-12,1e-9,1e-7,1e-4]);T=np.array([8000.,5000.,1e6,2e5])
    code=f'''
import json,sys,numpy as np
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
m=MesaMicrophysics(Composition(.7,.3,0.),{library!r},{os.environ['MESA_DIR']!r})
rho=np.array({rho.tolist()!r});T=np.array({T.tolist()!r})
try:
    m.opacity_stars(rho,T);stars='no error'
except ValueError as e:stars=str(e)
print(json.dumps(dict(abi=m.abi_version,raw=m.raw(rho,T).tobytes().hex(),
    paper=m.opacity_paper_reference(rho,T,'gs98').tobytes().hex(),stars=stars)))
'''
    env=dict(os.environ,PYTHONPATH=os.pathsep.join([str(ROOT/'src'),os.environ.get('PYTHONPATH','')]))
    r=subprocess.run([sys.executable,'-c',code],env=env,capture_output=True,text=True,timeout=300)
    assert r.returncode==0,r.stderr
    old=json.loads(r.stdout.strip().splitlines()[-1])
    assert old['abi']==6 and 'requires MESA bridge ABI 7' in old['stars']
    assert old['raw']==mesa.raw(rho,T).tobytes().hex()
    assert old['paper']==mesa.opacity_paper_reference(rho,T,'gs98').tobytes().hex()
