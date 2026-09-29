"""Ideal resolved-atom material energy and its extended-state differential.

Population coordinates are numbers per gram, not LTE functions of rho and T.
This supplies exact *frozen-population* partials and population energy terms.
It does not supply a rate response, nonideal free energy, or a sonic closure.
"""
from dataclasses import dataclass
import numpy as np
from .constants import C

KB=1.380649e-16
HC=6.62607015e-27*C


@dataclass(frozen=True)
class ResolvedMaterial:
    density: np.ndarray
    temperature: np.ndarray
    number_per_gram: np.ndarray
    charge: np.ndarray
    atomic_energy: np.ndarray
    element: np.ndarray
    labels: tuple

    def __post_init__(self):
        rho,T,y,Z,E,el=map(np.asarray,(self.density,self.temperature,self.number_per_gram,
                                     self.charge,self.atomic_energy,self.element))
        if (rho.ndim!=1 or T.shape!=rho.shape or y.shape!=(len(Z),len(rho))
                or E.shape!=Z.shape or el.shape!=Z.shape or len(self.labels)!=len(Z)):
            raise ValueError('Inconsistent resolved material dimensions')
        if (not all(np.all(np.isfinite(a)) for a in [rho,T,y,Z,E])
                or np.any(rho<=0) or np.any(T<=0) or np.any(y<0)
                or np.any(Z@y<0) or np.any((1+Z)@y<=0)):
            raise ValueError('Invalid resolved material state or negative free-electron density')

    def values(self):
        a=(1+self.charge)@self.number_per_gram
        thermal=1.5*KB*self.temperature*a
        p=self.density*KB*self.temperature*a
        atomic=self.atomic_energy@self.number_per_gram
        return dict(pressure=p,internal_energy=thermal+atomic,atomic_energy=atomic,
                    pressure_dlnrho=p,pressure_dlnT=p,
                    energy_dlnrho=np.zeros_like(p),energy_dlnT=thermal,
                    pressure_dpopulation=self.density[None,:]*KB*self.temperature[None,:]*(1+self.charge[:,None]),
                    energy_dpopulation=1.5*KB*self.temperature[None,:]*(1+self.charge[:,None])+self.atomic_energy[:,None],
                    energy_gauge='neutral_H_He_ground',derivative_constraint='fixed numbers per gram',
                    equilibrium_derivatives_available=False,production_certified=False)

    def differential(self,dlnrho,dlnT,dpopulation):
        """Directional differential on the fixed-composition population manifold.

        dpopulation may be a perturbation or a derivative (e.g. dY/dln r).
        Element nuclei must be conserved. Electrons follow charge neutrality.
        """
        dy=np.asarray(dpopulation,float)
        if dy.shape!=self.number_per_gram.shape or not np.all(np.isfinite(dy)):
            raise ValueError('Invalid population direction')
        for elem in set(self.element):
            a=dy[self.element==elem]
            if np.any(abs(a.sum(axis=0))>1e-12*np.maximum(abs(a).sum(axis=0),1.)):
                raise ValueError('Population direction does not conserve element nuclei')
        dr,dt=np.broadcast_arrays(dlnrho,dlnT)
        if dr.shape not in [(),self.density.shape] or not np.all(np.isfinite(dr+dt)):
            raise ValueError('Invalid thermodynamic direction')
        q=self.values()
        dp=q['pressure_dlnrho']*dr+q['pressure_dlnT']*dt+np.sum(q['pressure_dpopulation']*dy,axis=0)
        du=q['energy_dlnT']*dt+np.sum(q['energy_dpopulation']*dy,axis=0)
        # du + P d(1/rho); mechanical radiation work is NOT part of this heat.
        return dict(pressure=dp,internal_energy=du,
                    first_law_heat=du-q['pressure']/self.density*dr,
                    population_energy=np.sum(q['energy_dpopulation']*dy,axis=0))


def hhe_resolved_material(state,populations,levels):
    """Pack native H-/H/He levels in a shared neutral-ground energy gauge.

    H- binding energy is negative in this gauge. Native integer nuclear masses
    remain encoded in the input density; no FreeEOS mass convention is grafted.
    Missing H- populations cannot be inferred from HI, T, or an LTE EOS.
    """
    groups=[('HI',0,'H',0.),('HII',1,'H',levels['HI']['ionization_wavenumber']),
            ('HeI',0,'He',0.),('HeII',1,'He',levels['HeI']['ionization_wavenumber']),
            ('HeIII',2,'He',levels['HeI']['ionization_wavenumber']+levels['HeII']['ionization_wavenumber'])]
    if 'Hminus' in populations:
        groups.insert(0,('Hminus',-1,'H',-levels['Hminus']['ionization_wavenumber']))
    density=np.atleast_1d(np.asarray(state['density_clump'],float))
    temperature=np.atleast_1d(np.asarray(state['temperature'],float))
    yy=[];zz=[];ee=[];elements=[];labels=[]
    for ion,z,el,offset in groups:
        pop=np.asarray(populations[ion],float)
        if pop.ndim==1:pop=pop[None,:]
        if pop.ndim!=2 or pop.shape[1]!=len(density):raise ValueError('Invalid population array')
        w=np.asarray(levels[ion]['wavenumber']) if ion in levels else np.zeros(1)
        if len(w)!=len(pop):raise ValueError('Atomic level/population dimensions differ')
        yy.extend(pop/density);zz.extend([z]*len(pop));ee.extend(HC*(offset+w))
        elements.extend([el]*len(pop));labels.extend(f'{ion}:{i}' for i in range(len(pop)))
    return ResolvedMaterial(density,temperature,np.asarray(yy),np.asarray(zz),np.asarray(ee),
                            np.asarray(elements),tuple(labels))
