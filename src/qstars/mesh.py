"""Conservative overlap remap of volume and TOTAL energy in mass coordinates."""
import numpy as np
from .constants import PI,G

def remap(state,new_dm,eos):
    new_dm = np.asarray(new_dm,dtype=float)
    if np.any(new_dm <= 0) or not np.isclose(new_dm.sum(),state.dm.sum(),rtol=1e-13):
        raise ValueError("Remap must preserve total mass with positive cells")
    a,b = np.r_[0.,np.cumsum(state.dm)],np.r_[0.,np.cumsum(new_dm)]
    b[-1] = a[-1]  # Identical physical domain; suppress cumulative-roundoff sliver.
    overlap = np.maximum(0.,np.minimum(b[1:,None],a[None,1:])-np.maximum(b[:-1,None],a[None,:-1]))
    gravity = -G*.5*(state.enclosed[:-1]/state.r[:-1]+state.enclosed[1:]/state.r[1:])
    kinetic = .25*(state.v[:-1]**2+state.v[1:]**2)
    specific = eos.evaluate(state.rho,state.T)["u"]+gravity+kinetic
    E = overlap@specific
    V = overlap@(state.volume/state.dm)
    out = state.copy()
    out.dm = new_dm.copy()
    out.face_luminosity=None
    out.reference_mass_fractions=new_dm/new_dm.sum()
    out.cell_volume = V.copy()
    out.r = np.r_[state.r[0],np.cbrt(state.r[0]**3+3*np.cumsum(V)/(4*PI))]
    out.r[-1] = state.r[-1]
    out.v = np.interp(b,a,state.v)
    new_gravity = -G*.5*(out.enclosed[:-1]/out.r[:-1]+out.enclosed[1:]/out.r[1:])
    new_kinetic = .25*(out.v[:-1]**2+out.v[1:]**2)
    internal = E/new_dm-new_gravity-new_kinetic
    if np.any(internal <= 0): raise ValueError("Conservative remap would require negative internal energy")
    out.T = eos.temperature(out.rho,internal)
    ds = float(np.dot(out.dm,eos.evaluate(out.rho,out.T)["s"])-np.dot(state.dm,eos.evaluate(state.rho,state.T)["s"]))
    if out.ledger: out.ledger["remap_entropy_change"] += ds
    return out
