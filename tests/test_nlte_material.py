from dataclasses import replace
import numpy as np
import pytest
from qstars.nlte_material import hhe_resolved_material,HC,KB


def material():
    # Several cool states with independent ionization/excitation fractions.
    rho=np.array([1e-12,1e-10,1e-8]);T=np.array([4000.,6500.,10000.])
    hn=4e23*rho;hen=.45e23*rho
    p=dict(Hminus=np.array([hn*.01]),HI=np.array([hn*.39,hn*.1]),HII=hn*.5,
           HeI=np.array([hen*.7]),HeII=np.array([hen*.2]),HeIII=hen*.1)
    levels={ion:dict(wavenumber=np.array(w),ionization_wavenumber=threshold) for ion,w,threshold in [
        ('Hminus',[0.],6083.06),('HI',[0.,82259.],109678.77),
        ('HeI',[0.],198310.7),('HeII',[0.],438908.9)]}
    return hhe_resolved_material(dict(density_clump=rho,temperature=T),p,levels)


def test_extended_first_law_derivatives_independent_directions():
    m=material();dy=np.zeros_like(m.number_per_gram)
    # Free H ionization, bound-bound excitation and He ionization independently.
    dy[0]=-.002e23;dy[1]=-.098e23;dy[2]=.03e23;dy[3]=.07e23
    dy[4]=-.04e23;dy[5]=.03e23;dy[6]=.01e23
    dr=np.array([.2,-.1,.15]);dt=np.array([-.3,.2,.1]);eps=1e-5
    plus=replace(m,density=m.density*np.exp(eps*dr),temperature=m.temperature*np.exp(eps*dt),
                 number_per_gram=m.number_per_gram+eps*dy).values()
    minus=replace(m,density=m.density*np.exp(-eps*dr),temperature=m.temperature*np.exp(-eps*dt),
                  number_per_gram=m.number_per_gram-eps*dy).values()
    d=m.differential(dr,dt,dy);q=m.values()
    for key in ['pressure','internal_energy']:
        np.testing.assert_allclose((plus[key]-minus[key])/(2*eps),d[key],rtol=1e-8)
    volume_derivative=(1/(m.density*np.exp(eps*dr))-1/(m.density*np.exp(-eps*dr)))/(2*eps)
    np.testing.assert_allclose(d['first_law_heat'],d['internal_energy']+q['pressure']*volume_derivative,rtol=1e-8)
    assert q['derivative_constraint']=='fixed numbers per gram'
    assert not q['equilibrium_derivatives_available']


def test_ionization_and_negative_ion_binding_not_double_counted():
    m=material();dy=np.zeros_like(m.number_per_gram);rate=1e22
    dy[0]=-rate;dy[1]=rate # detach H- to HI, produces one electron
    d=m.differential(0.,0.,dy)
    np.testing.assert_allclose(d['internal_energy'],rate*(HC*6083.06+1.5*KB*m.temperature))
    np.testing.assert_allclose(d['pressure'],rate*m.density*KB*m.temperature)
    dy[:]=0;dy[1]=-rate;dy[3]=rate # ionize ground H
    np.testing.assert_allclose(m.differential(0.,0.,dy)['internal_energy'],rate*(HC*109678.77+1.5*KB*m.temperature))


def test_element_energy_gauge_changes_cancel_in_first_law():
    m=material();dy=np.zeros_like(m.number_per_gram);dy[1]=-1e22;dy[3]=1e22
    shifted=replace(m,atomic_energy=m.atomic_energy+np.where(m.element=='H',2e-11,7e-11))
    np.testing.assert_allclose(m.differential(.2,.1,dy)['first_law_heat'],shifted.differential(.2,.1,dy)['first_law_heat'])


def test_nuclear_creation_and_negative_electrons_rejected():
    m=material();dy=np.ones_like(m.number_per_gram)
    with pytest.raises(ValueError,match='conserve element'):m.differential(0.,0.,dy)
    y=np.zeros_like(m.number_per_gram);y[0]=1e23
    with pytest.raises(ValueError,match='negative free-electron'):replace(m,number_per_gram=y)


def test_lte_recombination_requires_population_chain_rule():
    # Independent ground-state pure-H Saha control, with no nonideal physics.
    # It tests the thermodynamic chain rule, not CMFGEN atomic rates.
    rho=1e-11;T=7000.;mass=1.6735575e-24;chi=HC*109678.77
    def state(temp):
        S=(2*np.pi*9.1093837139e-28*KB*temp/(6.62607015e-27)**2)**1.5*np.exp(-chi/(KB*temp))
        b=S/(rho/mass);x=2/(1+np.sqrt(1+4/b))
        from qstars.nlte_material import ResolvedMaterial
        return ResolvedMaterial(np.array([rho]),np.array([temp]),np.array([[(1-x)/mass],[x/mass]]),
            np.array([0,1]),np.array([0.,chi]),np.array(['H','H']),('HI:0','HII:0'))
    m=state(T);eps=1e-5;plus=state(T*np.exp(eps));minus=state(T*np.exp(-eps))
    dy=(plus.number_per_gram-minus.number_per_gram)/(2*eps)
    du=m.differential(0.,1.,dy)['internal_energy']
    numerical=(plus.values()['internal_energy']-minus.values()['internal_energy'])/(2*eps)
    np.testing.assert_allclose(du,numerical,rtol=2e-8)
    assert du[0]>3*m.values()['energy_dlnT'][0]
