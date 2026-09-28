"""Restricted opacity prescriptions. No hidden cool-temperature extrapolation."""
import numpy as np

class ElectronScattering:
    source = "thomson"
    def __init__(self, composition):
        self.kappa = .2*(1+composition.X)
    def __call__(self, rho, T):
        return np.full(np.broadcast(rho,T).shape, self.kappa)

class OpacityTable:
    """Strict rectangular log10(rho), log10(T), log10(kappa) table in CGS.

    User supplies provenance. Out-of-domain queries raise, never extrapolate.
    """
    def __init__(self, path):
        from scipy.interpolate import RegularGridInterpolator
        with np.load(path, allow_pickle=False) as data:
            if "source" not in data or not str(data["source"]):
                raise ValueError("Opacity table must contain source provenance")
            self.source = str(data["source"])
            self.interp = RegularGridInterpolator((data["logrho"],data["logT"]),
                          data["logkappa"],bounds_error=True)
    def __call__(self, rho, T):
        rho,T = np.broadcast_arrays(rho,T)
        values = self.interp(np.stack([np.log10(rho),np.log10(T)],axis=-1))
        return (10**values).reshape(rho.shape)
