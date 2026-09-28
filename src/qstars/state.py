from dataclasses import dataclass, field
import copy
import numpy as np
from .constants import G, PI

@dataclass
class State:
    time: float
    bh_mass: float
    dm: np.ndarray
    r: np.ndarray
    v: np.ndarray
    T: np.ndarray
    ledger: dict = field(default_factory=dict)
    steps: int = 0
    next_dt: float = 0.

    def copy(self):
        return copy.deepcopy(self)

    @property
    def volume(self):
        return 4*PI/3*np.diff(self.r**3)

    @property
    def rho(self):
        return self.dm/self.volume

    @property
    def weights(self):
        return .5*np.r_[self.dm[0],self.dm[:-1]+self.dm[1:],self.dm[-1]]

    @property
    def enclosed(self):
        return self.bh_mass+np.r_[0.,np.cumsum(self.dm)]

    @property
    def rc(self):
        return (.5*(self.r[:-1]**3+self.r[1:]**3))**(1/3)

    def energies(self, eos):
        internal = float(np.dot(self.dm,eos.evaluate(self.rho,self.T)["u"]))
        gravity = float(-G*np.sum(self.weights*self.enclosed/self.r))
        kinetic = float(.5*np.sum(self.weights*self.v**2))
        return dict(internal=internal,gravity=gravity,kinetic=kinetic,
                    total=internal+gravity+kinetic)
