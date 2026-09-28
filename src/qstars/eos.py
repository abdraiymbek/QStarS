"""Fully ionized, classical monatomic ideal gas plus LTE radiation.

The additive entropy constant is fixed by rho_ref=1 g/cm^3, T_ref=1 K.
No recombination, degeneracy, pairs, Coulomb terms, or composition evolution.
"""
from dataclasses import dataclass
import numpy as np
from .constants import A, KB, MU

@dataclass(frozen=True)
class Composition:
    X: float = 0.75
    Y: float = 0.25
    Z: float = 0.0

    def __post_init__(self):
        if not np.all(np.isfinite([self.X,self.Y,self.Z])) or min(self.X, self.Y, self.Z) < 0 or abs(self.X+self.Y+self.Z-1) > 1e-10:
            raise ValueError("X, Y, Z must be nonnegative and sum to one")

    @property
    def mu(self):
        # Metals approximation: A ~ 2Z; trace metals only.
        return 1 / (2*self.X + .75*self.Y + .5*self.Z)

class IdealGasRadiation:
    def __init__(self, composition=Composition()):
        self.composition = composition
        self.R = KB / (composition.mu * MU)

    def thermal_contrast(self,rho,T,q):
        """Stable adjacent dlnT - mean(nabla_ad)*dlnP near radiation dominance."""
        lt=np.log1p(np.diff(T)/T[:-1])
        lr=np.log1p(np.diff(rho)/rho[:-1])
        beta=q['pg']/q['p']
        one_minus_four_grad=-3*beta**2/(32-24*beta-3*beta**2)
        grad=.5*(q['grad_ad'][:-1]+q['grad_ad'][1:])
        gas_correction=np.log1p(beta[:-1]*np.expm1(lr-3*lt))
        return .5*(one_minus_four_grad[:-1]+one_minus_four_grad[1:])*lt-grad*gas_correction

    def evaluate(self, rho, T):
        rho, T = np.asarray(rho), np.asarray(T)
        if np.any(rho <= 0) or np.any(T <= 0):
            raise ValueError("EOS needs positive density and temperature")
        pg = self.R*rho*T
        pr = A*T**4/3
        p = pg+pr
        u = 1.5*self.R*T + A*T**4/rho
        s = self.R*(1.5*np.log(T)-np.log(rho)) + 4*A*T**3/(3*rho)
        cv = 1.5*self.R+4*A*T**3/rho
        beta = pg/p
        chiT, chiR = 4-3*beta, beta
        gamma3m1 = p*chiT/(rho*T*cv)
        gamma1 = chiR+chiT*gamma3m1
        grad_ad = gamma3m1/gamma1
        delta = chiT/chiR
        cp = cv+p*chiT**2/(rho*T*chiR)
        return dict(p=p, pg=pg, pr=pr, u=u, s=s, cv=cv, cp=cp,
                    gamma1=gamma1, grad_ad=grad_ad, delta=delta,
                    cs=np.sqrt(gamma1*p/rho), h=u+p/rho)

    def temperature(self, rho, value, quantity="u"):
        """Monotonic log-Newton EOS inversion, relative convergence required."""
        rho, value = np.broadcast_arrays(rho, value)
        if np.any(value <= 0) or quantity not in ("u", "p"):
            raise ValueError("EOS inversion requires positive u or p")
        if quantity == "u":
            T = np.minimum(value/(1.5*self.R), (value*rho/A)**.25)
        else:
            T = np.minimum(value/(rho*self.R), (3*value/A)**.25)
        for _ in range(40):
            f = self.evaluate(rho,T)
            deriv = T*f["cv"] if quantity == "u" else f["pg"]+4*f["pr"]
            step = np.clip((value-f[quantity])/deriv, -.5, .5)
            T *= np.exp(step)
            if np.max(abs(step)) < 2e-13:
                return T
        raise RuntimeError("EOS inversion failed")
