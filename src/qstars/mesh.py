"""Conservative overlap remap of volume and TOTAL energy in mass coordinates."""
import numpy as np
from .constants import PI,G

def overlap_integrals(old_dm,new_dm,cell_integrals):
    """Conservative overlap without subtracting large binary mass coordinates.

    Very small boundary cells lose relative precision in cumulative CGS mass
    coordinates. Decimal endpoints are used only to calculate overlap weights;
    thermodynamics and evolution remain double precision. The sweep is O(N),
    including when old and new zone counts differ.
    """
    import decimal,math
    with decimal.localcontext() as context:
        context.prec=60
        D=decimal.Decimal.from_float
        a=[decimal.Decimal(0)];b=[decimal.Decimal(0)]
        for mass in old_dm:a.append(a[-1]+D(float(mass)))
        for mass in new_dm:b.append(b[-1]+D(float(mass)))
        b[-1]=a[-1]
        out=np.zeros((len(new_dm),cell_integrals.shape[1]))
        i=0
        for j in range(len(new_dm)):
            contributions=[]
            while i<len(old_dm) and a[i]<b[j+1]:
                overlap=min(a[i+1],b[j+1])-max(a[i],b[j])
                if overlap>0:
                    contributions.append(cell_integrals[i]*float(overlap/(a[i+1]-a[i])))
                if a[i+1]<=b[j+1]:i+=1
                else:break
            for k in range(out.shape[1]):out[j,k]=math.fsum(x[k] for x in contributions)
    return out

def remap(state,new_dm,eos):
    new_dm = np.asarray(new_dm,dtype=float)
    if np.any(new_dm <= 0) or not np.isclose(new_dm.sum(),state.dm.sum(),rtol=1e-13):
        raise ValueError("Remap must preserve total mass with positive cells")
    a,b = np.r_[0.,np.cumsum(state.dm)],np.r_[0.,np.cumsum(new_dm)]
    b[-1] = a[-1]  # Identical physical domain; suppress cumulative-roundoff sliver.
    gravity = -G*.5*(state.enclosed[:-1]/state.r[:-1]+state.enclosed[1:]/state.r[1:])
    kinetic = .25*(state.v[:-1]**2+state.v[1:]**2)
    specific = eos.evaluate(state.rho,state.T)["u"]+gravity+kinetic
    integrals=overlap_integrals(state.dm,new_dm,np.column_stack([state.dm*specific,state.volume]))
    E,V=integrals.T
    out = state.copy()
    out.dm = new_dm.copy()
    out.face_luminosity=None
    out.reference_mass_fractions=new_dm/new_dm.sum()
    out.cell_volume = V.copy()
    out.r = np.r_[state.r[0],np.cbrt(state.r[0]**3+3*np.cumsum(V)/(4*PI))]
    out.r[-1] = state.r[-1]
    out.radius_compensation=np.zeros_like(out.r)
    out.v = np.interp(b,a,state.v)
    new_gravity = -G*.5*(out.enclosed[:-1]/out.r[:-1]+out.enclosed[1:]/out.r[1:])
    new_kinetic = .25*(out.v[:-1]**2+out.v[1:]**2)
    internal = E/new_dm-new_gravity-new_kinetic
    if np.any(internal <= 0): raise ValueError("Conservative remap would require negative internal energy")
    out.T = eos.temperature(out.rho,internal)
    ds = float(np.dot(out.dm,eos.evaluate(out.rho,out.T)["s"])-np.dot(state.dm,eos.evaluate(state.rho,state.T)["s"]))
    if out.ledger: out.ledger["remap_entropy_change"] += ds
    return out
