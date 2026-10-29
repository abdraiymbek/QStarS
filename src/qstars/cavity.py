"""Ball et al. 2011 moving inner boundary: Bondi radius, cavity gas, base slicing.

The inner face r0 sits at r_B/b with r_B = 2 G M_BH/c_b^2 and c_b^2 =
Gamma1 p/rho of the gas on that face (thesis 3.2-3.4; b = 1). The cavity
between the black hole and r0 holds M_cav = (8 pi/3) rho_b r0^3 for
rho ~ r^-3/2 (thesis 3.7; paper eq 46 gives 8 pi/5 for r^-1/2). STARS
evaluates the base state on the innermost meshpoint (thesis 2.2.2), so the
face mode extrapolates cells 0 and 1 half a cell hydrostatically to r0
(second order); the cell mode (variant C6) uses cell 0 itself.

Swept gas leaves the envelope with exactly the energy difference of the
discrete envelope functional (State.energies) at fixed geometry: the base
slice through mesh.remap, the cell-0 sink correction delta through sink().
No energy balance is claimed inside the cavity or the black hole.
Every function on the residual path accepts complex-step states.
"""
import math
from fractions import Fraction
import numpy as np
from .constants import G,PI
from .mesh import remap

PROFILES={'r-3/2':8*PI/3,'r-1/2':8*PI/5,'none':0.}
# Slices below this fraction of the base cell are left to the implicit sink
# delta: a sliver cut would only re-face cell 0 at rounding level.
MINIMUM_SLICE=1e-9

def cavity_mass(rho_b,r0,profile='r-3/2'):
    """Gas between the black hole and r0 (thesis 3.7; paper eq 46 for r^-1/2)."""
    if profile not in PROFILES:raise ValueError(f'Unknown cavity profile {profile}')
    return PROFILES[profile]*rho_b*r0**3

def bondi_radius(bh_mass,cs,b=1.):
    """Inner radius r0 = r_B/b, r_B = 2 G M_BH/c_b^2 (thesis 3.2, 3.4)."""
    return 2*G*bh_mass/(b*cs**2)

def boundary_state(rho,T,p,r0,m0,dm0,eos,mode='face'):
    """Gas state on the inner face from cells 0 and 1.

    face: p_b = p_0 + G m_0 dm_0/(8 pi r0^4) (hydrostatic over the lower
    half of cell 0, with the face mass and radius), and ln T, ln rho
    follow the cell 0-1 slopes in ln p: ln T_b = ln T_0 + nabla_01 ln(p_b/p_0).
    cell: the cell-0 state. Gamma1 and c_b come from the EOS at (rho_b, T_b).
    """
    if mode=='cell':
        rho_b,T_b,lnpb=rho[0],T[0],0.*rho[0]
    elif mode=='face':
        lnp01=np.log1p((p[1]-p[0])/p[0])
        if not lnp01.real<0:raise ValueError('Face extrapolation needs outward-decreasing pressure in cells 0-1')
        lnpb=np.log1p(G*m0*dm0/(8*PI*r0**4*p[0]))
        T_b=T[0]*np.exp(np.log1p((T[1]-T[0])/T[0])/lnp01*lnpb)
        rho_b=rho[0]*np.exp(np.log1p((rho[1]-rho[0])/rho[0])/lnp01*lnpb)
    else:raise ValueError(f'Unknown boundary-state mode {mode}')
    q=eos.evaluate(np.atleast_1d(rho_b),np.atleast_1d(T_b))
    return dict(rho=rho_b,T=T_b,p=q['p'][0],p_extrapolated=p[0]*np.exp(lnpb),gamma1=q['gamma1'][0],
                cs=q['cs'][0],h=q['h'][0],s=q['s'][0])

def state_boundary(state,eos,mode='face'):
    """boundary_state of an assembled State (cells 0 and 1, face 0)."""
    q=eos.evaluate(state.rho[:2],state.T[:2])
    return boundary_state(state.rho[:2],state.T[:2],q['p'],state.r[0],state.enclosed[0],state.dm[0],eos,mode)

def slice_base(state,mass,eos):
    """Move `mass` of the base gas into the transit central mass (cavity_mass).

    mass > 0 slices the bottom of the envelope with exact (Fraction) cell
    edges: whole cells leave, the cut cell keeps its (rho, u) and its
    remainder becomes cell 0; untouched cells are copied bitwise. mass < 0
    returns gas from the cavity with the cell-0 state, merged into cell 0
    (conserving E and V). |mass| below MINIMUM_SLICE of the base cell is
    left to the implicit sink. Returns (state, info) with mesh.remap's slab
    info: mass, volume, entropy, and the functional energy difference
    'energy' = E(before) - E(after).
    """
    if not np.isfinite(mass):raise ValueError('Swept mass must be finite')
    if abs(mass)<=MINIMUM_SLICE*state.dm[0]:
        return state.copy(),dict(mass=0.,volume=0.,internal=0.,entropy=0.,overlap_energy=0.,energy=0.,
                                 reconstruction_defect=0.,copied_cells=len(state.dm),fixed_state_cells=0,
                                 reconstructed_cells=0)
    if mass>0:
        edges=[Fraction(0)]
        for x in state.dm:edges.append(edges[-1]+Fraction(float(x)))
        cut=Fraction(float(mass))
        if not cut<edges[-1]:raise ValueError('Swept mass exceeds the envelope')
        k=max(i for i in range(len(state.dm)) if edges[i]<=cut)
        remainder=float(edges[k+1]-cut)
        if not remainder>0:raise ValueError('Base slice leaves no remainder')
        new_dm=np.r_[remainder,state.dm[k+1:]]
        return remap(state,new_dm,eos,offset=float(mass),copy_unchanged=True,return_slab=True)
    returned=-float(mass)
    new_dm=np.r_[state.dm[0]+returned,state.dm[1:]]
    return remap(state,new_dm,eos,offset=-returned,copy_unchanged=True,return_slab=True)

def sink(state,delta,eos):
    """Remove delta from cell 0 at fixed geometry with the cell-0 enthalpy.

    dm_0 -> dm_0 - delta, U_0 -> dm_0 u_0 - delta h_0 (negative delta
    returns gas with the cell-0 state), the transit central mass gains
    delta, so the enclosed masses of faces 1..N do not change. Returns the
    new state; its energy difference from `state` is what delta carries.
    """
    q=eos.evaluate(state.rho[:1],state.T[:1])
    out=state.copy()
    if delta==0:return out
    mass=state.dm[0]-delta
    if not mass>0:raise ValueError('Sink correction exceeds the base cell')
    internal=(state.dm[0]*q['u'][0]-delta*q['h'][0])/mass
    if not internal>0:raise ValueError('Sink correction exhausts the base-cell internal energy')
    out.dm=state.dm.copy();out.dm[0]=mass
    out.cell_volume=state.volume.copy()
    out.T=state.T.copy();out.T[0]=eos.temperature(np.atleast_1d(mass/state.volume[0]),np.atleast_1d(internal))[0]
    out.cavity_mass=state.cavity_mass+delta
    out.reference_mass_fractions=out.dm/out.dm.sum()
    return out

def sink_energy(state,delta,eos):
    """Closed form of E(state) - E(sink(state, delta)): delta h_0 plus the re-weighted potential."""
    q=eos.evaluate(state.rho[:1],state.T[:1])
    m,r,dm=state.enclosed,state.r,state.dm
    return float(delta*q['h'][0]-.5*G*delta*((m[0]-dm[0]+delta)/r[0]+m[1]/r[1]))

def defect_work(state,mass):
    """Energy gained when `mass` leaves the central mass at fixed radii: G mass sum(w/r).

    Exact, because the discrete potential energy is linear in M0.
    """
    return float(G*mass*math.fsum(state.weights/state.r))

def piston_work(p0,r0,r0bar,w0,m0,dr=None):
    """Inner-face work of the secant identity: p0 (4 pi/3)(r0^3 - r0bar^3) + G w0 m0 (1/r0bar - 1/r0).

    dr is the solver's displacement (default r0 - r0bar); the step's cell
    volumes use the same secant area times dr, so the identity is exact.
    """
    dr=r0-r0bar if dr is None else dr
    return float(p0*4*PI/3*(r0*r0+r0*r0bar+r0bar*r0bar)*dr+G*w0*m0*dr/(r0*r0bar))

def bookings(energy,entropy,swept,accreted,epsilon):
    """Ledger increments for one step's swept gas (spec 4, booking 1).

    The swept gas (slice plus sink) leaves with `energy` and `entropy`.
    The split by mass fraction is a declared labelling choice: the rest mass
    swallowed this step (accreted) is advected into the black hole, the rest
    of the swept mass (the cavity change) is cavity_swept. eps*accreted
    leaves the gravitating mass (mass_defect).
    """
    share=accreted/swept if swept else (1. if accreted else 0.)
    return dict(advected_energy=energy*share,cavity_swept_energy=energy-energy*share,
                advected_entropy=entropy*share,cavity_swept_entropy=entropy-entropy*share,
                cavity_swept_mass=swept,accreted_rest_mass=accreted,mass_defect=epsilon*accreted)
