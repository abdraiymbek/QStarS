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
        with np.load(path, allow_pickle=False) as data:
            if "source" not in data or not str(data["source"]):
                raise ValueError("Opacity table must contain source provenance")
            self.source = str(data["source"])
            self.x,self.y,self.z=data['logrho'].copy(),data['logT'].copy(),data['logkappa'].copy()
        if self.z.shape!=(len(self.x),len(self.y)) or min(len(self.x),len(self.y))<2 or np.any(np.diff(self.x)<=0) or np.any(np.diff(self.y)<=0):
            raise ValueError('Opacity table needs strictly increasing rectangular axes')
        if not all(np.all(np.isfinite(a)) for a in (self.x,self.y,self.z)):
            raise ValueError('Opacity table contains nonfinite data')
    def __call__(self, rho, T):
        rho,T = np.broadcast_arrays(rho,T)
        x,y=np.log10(rho),np.log10(T)
        if np.any(x.real<self.x[0]) or np.any(x.real>self.x[-1]) or np.any(y.real<self.y[0]) or np.any(y.real>self.y[-1]):
            raise ValueError('Opacity query outside table domain')
        i=np.clip(np.searchsorted(self.x,x.real)-1,0,len(self.x)-2)
        j=np.clip(np.searchsorted(self.y,y.real)-1,0,len(self.y)-2)
        u=(x-self.x[i])/(self.x[i+1]-self.x[i]);v=(y-self.y[j])/(self.y[j+1]-self.y[j])
        values=(1-u)*(1-v)*self.z[i,j]+u*(1-v)*self.z[i+1,j]+(1-u)*v*self.z[i,j+1]+u*v*self.z[i+1,j+1]
        return (10**values).reshape(rho.shape)
