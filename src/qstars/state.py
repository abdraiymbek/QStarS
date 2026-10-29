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
    cell_volume: np.ndarray | None = None
    face_luminosity: np.ndarray | None = None
    reference_mass_fractions: np.ndarray | None = None
    radius_compensation: np.ndarray | None = None
    # Unresolved central gas between the black hole and the resolved base
    # (Ball cavity). Zero by default: bh_mass+0.0 is exact, so mass sums,
    # gravity and energies of existing models are bitwise unchanged.
    cavity_mass: float = 0.
    photosphere_density: float | None = None
    # Quasi-static mechanics carries no radial kinetic energy as a state variable.
    kinetic_tracked: bool = True
    # JSON-serialisable driver state (controller, predictor), checkpointed with the model.
    aux: dict = field(default_factory=dict)

    def copy(self):
        return copy.deepcopy(self)

    @property
    def volume(self):
        if self.cell_volume is not None: return self.cell_volume
        return 4*PI/3*np.diff(self.r)*(self.r[1:]**2+self.r[1:]*self.r[:-1]+self.r[:-1]**2)

    @property
    def rho(self):
        return self.dm/self.volume

    @property
    def weights(self):
        return .5*np.r_[self.dm[0],self.dm[:-1]+self.dm[1:],self.dm[-1]]

    @property
    def enclosed(self):
        return self.bh_mass+self.cavity_mass+np.r_[0.,np.cumsum(self.dm)]

    @property
    def rc(self):
        return (.5*(self.r[:-1]**3+self.r[1:]**3))**(1/3)

    def energies(self, eos):
        internal = float(np.dot(self.dm,eos.evaluate(self.rho,self.T)["u"]))
        gravity = float(-G*np.sum(self.weights*self.enclosed/self.r))
        if not self.kinetic_tracked:
            return dict(internal=internal,gravity=gravity,kinetic=0.,total=internal+gravity)
        kinetic = float(.5*np.sum(self.weights*self.v**2))
        return dict(internal=internal,gravity=gravity,kinetic=kinetic,
                    total=internal+gravity+kinetic)
