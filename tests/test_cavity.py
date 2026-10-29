"""Gate 6 Stage 4, U6 and U8: the Ball cavity boundary (cavity.py).

U6: the face extrapolation to r0 is second order on an analytic loaded
polytrope (the cell mode is first order). U8: every booking equals the
State.energies difference of its sub-operator to 1e-13 of the internal
energy, and the swept volume agrees with outer_interface.moving_flux.
Analytic microphysics only; no relaxation.
"""
import math
from fractions import Fraction
import numpy as np
import pytest
from qstars import cavity,ledger,quasistatic
from qstars.config import resolve
from qstars.constants import A,G,PI,MSUN,RSUN
from qstars.eos import IdealGasRadiation,Composition
from qstars.opacity import ElectronScattering
from qstars.outer_interface import InterfaceState,moving_flux

EOS=IdealGasRadiation(Composition(.7,.3,0.))
OPACITY=ElectronScattering(EOS.composition)

def toy_config(zones):
    c=resolve(dict(composition=dict(X=.7,Y=.3,Z=0.),numerics=dict(zones=zones)))
    c['engine'].update(model='ball2012_eq3_15',coupling=1.5e-4)
    c['transport']['convection']='mlt_ml1';c['boundary']['thermal']='stars_photosphere'
    return c

@pytest.fixture(scope='module')
def seed():
    """Discretely hydrostatic loaded-polytrope seed (100 Msun, M_BH = 1 Msun, 200 Rsun)."""
    c=toy_config(64)
    s,info=quasistatic.seed(c,EOS,OPACITY,64,1*MSUN,100*MSUN,200*RSUN)
    ledger.initialize(s,EOS,cavity=True)
    return s,info

def scale(s):
    return s.energies(EOS)['internal']

# ---------------------------------------------------------------- formulas

def test_profiles_bondi_radius_and_bookings_split():
    assert cavity.cavity_mass(2.,3.,'r-3/2')==pytest.approx(8*PI/3*2*27,rel=1e-15)
    assert cavity.cavity_mass(2.,3.,'r-1/2')==pytest.approx(8*PI/5*2*27,rel=1e-15)
    assert cavity.cavity_mass(2.,3.,'none')==0.
    with pytest.raises(ValueError):cavity.cavity_mass(2.,3.,'r-2')
    assert cavity.bondi_radius(MSUN,1e8,b=2.)==pytest.approx(G*MSUN/1e16,rel=1e-15)
    b=cavity.bookings(10.,4.,5.,2.,.1)
    assert b['advected_energy']==4. and b['cavity_swept_energy']==6. and b['advected_entropy']==pytest.approx(1.6)
    assert b['cavity_swept_mass']==5. and b['accreted_rest_mass']==2. and b['mass_defect']==pytest.approx(.2)
    # Relaxation (no rest mass swallowed) books everything as cavity sweep; an exact
    # zero swept mass with accretion (the cavity returned what was swallowed) advects it all.
    assert cavity.bookings(3.,1.,-2.,0.,.1)['cavity_swept_energy']==3.
    assert cavity.bookings(3.,1.,0.,2.,.1)['advected_energy']==3.

def test_boundary_state_modes_and_complex_step(seed):
    s,_=seed
    face=cavity.state_boundary(s,EOS,'face');cell=cavity.state_boundary(s,EOS,'cell')
    q=EOS.evaluate(s.rho[:2],s.T[:2])
    assert cell['rho']==s.rho[0] and cell['T']==s.T[0] and cell['p']==pytest.approx(q['p'][0],rel=1e-15)
    # Face: the half-cell hydrostatic pressure, then the cell 0-1 slopes in ln p.
    pb=q['p'][0]+G*s.enclosed[0]*s.dm[0]/(8*PI*s.r[0]**4)
    assert face['p_extrapolated']==pytest.approx(pb,rel=1e-14)
    slope=np.log(s.T[1]/s.T[0])/np.log(q['p'][1]/q['p'][0])
    assert np.log(face['T']/s.T[0])==pytest.approx(slope*np.log(pb/q['p'][0]),rel=1e-12)
    assert face['cs']==pytest.approx(np.sqrt(face['gamma1']*face['p']/face['rho']),rel=1e-15)
    with pytest.raises(ValueError,match='outward-decreasing'):  # pressure inversion across cells 0-1
        cavity.boundary_state(s.rho[:2],s.T[:2],np.array([q['p'][0],2*q['p'][0]]),s.r[0],s.enclosed[0],s.dm[0],EOS)
    # Complex steps through the extrapolation equal central differences.
    rho,T,p=s.rho[:2].astype(complex),s.T[:2].astype(complex),q['p'][:2]
    for k,x in enumerate((rho,T)):
        h=1e-30*x[0].real;y=x.copy();y[0]+=1j*h
        args=(y,T,p) if k==0 else (rho,y,p)
        exact=cavity.boundary_state(*args,s.r[0],s.enclosed[0],s.dm[0],EOS)['cs'].imag/h
        up,down=[a.real.copy() for a in (rho,T)],[a.real.copy() for a in (rho,T)]
        up[k][0]*=1+1e-6;down[k][0]*=1-1e-6
        f=lambda r_,T_:cavity.boundary_state(r_,T_,p,s.r[0],s.enclosed[0],s.dm[0],EOS)['cs']
        assert exact==pytest.approx((f(*up)-f(*down))/(2e-6*x[0].real),rel=1e-6)

# ---------------------------------------------------------------- U6

def loaded_cells(sol,alpha,rho_c,K,h,cells=2):
    """Cells of width h (in xi) from the inner face of an analytic loaded n = 3 polytrope."""
    xi=sol['xi0']+h*np.arange(cells+1)
    phi=sol['solution'].sol(xi)[1]
    dm=4*PI*alpha**3*rho_c*np.diff(phi);r=alpha*xi
    V=4*PI/3*np.diff(r**3);rho=dm/V;p=K*rho**(4/3)
    return rho,EOS.temperature(rho,p,'p'),p,r[0],4*PI*alpha**3*rho_c*phi[0],dm[0]

def test_face_extrapolation_is_second_order_on_loaded_polytrope():
    sol=quasistatic.loaded_polytrope(5e-3)
    M,R=1e3*MSUN,1e3*RSUN
    alpha=R/sol['xi_s'];rho_c=M/(4*PI*alpha**3*sol['phi_s']);K=PI*G*alpha**2*rho_c**(2/3)
    p_c=K*rho_c**(4/3);T_c=EOS.temperature(np.array([rho_c]),np.array([p_c]),'p')[0]
    exact=EOS.evaluate(np.array([rho_c]),np.array([T_c]))
    bh=4*PI*alpha**3*rho_c*sol['xi0']/6
    errors={'face':[],'cell':[]}
    for k in range(4):
        h=.2*sol['xi0']/2**k
        rho,T,p,r0,m0,dm0=loaded_cells(sol,alpha,rho_c,K,h)
        for mode in errors:
            B=cavity.boundary_state(rho,T,p,r0,m0,dm0,EOS,mode)
            errors[mode].append([abs(B['rho']/rho_c-1),abs(B['T']/T_c-1),abs(B['p']/p_c-1),
                                 abs(cavity.bondi_radius(bh,B['cs'])/cavity.bondi_radius(bh,exact['cs'][0])-1)])
    face,cell=np.array(errors['face']),np.array(errors['cell'])
    face_order,cell_order=np.log2(face[:-1]/face[1:]),np.log2(cell[:-1]/cell[1:])
    assert np.all(face_order>=1.8),face_order
    assert np.all(abs(cell_order-1)<=.2),cell_order
    assert np.all(face[-1]<.05*cell[-1])
    # r0 = r_B by construction of the loaded polytrope (Gamma = 4/3): the ideal-gas Gamma1 sets the rest.
    gamma=exact['gamma1'][0]
    assert cavity.bondi_radius(bh,exact['cs'][0])==pytest.approx(alpha*sol['xi0']*(4/3)/gamma,rel=1e-9)

# ---------------------------------------------------------------- U8 bookings

@pytest.mark.parametrize('fraction',[1e-3,-1e-3,2e-2,-2e-2])
def test_sink_booking_equals_functional_difference(seed,fraction):
    s,_=seed;delta=fraction*s.dm[0]
    out=cavity.sink(s,delta,EOS)
    q=EOS.evaluate(s.rho[:1],s.T[:1])
    assert out.dm[0]==s.dm[0]-delta and out.cavity_mass==s.cavity_mass+delta and out.volume[0]==s.volume[0]
    assert math.fsum(out.dm)+out.cavity_mass==pytest.approx(math.fsum(s.dm)+s.cavity_mass,rel=1e-15)
    assert out.dm[0]*EOS.evaluate(out.rho[:1],out.T[:1])['u'][0]==pytest.approx(s.dm[0]*q['u'][0]-delta*q['h'][0],rel=1e-13)
    assert out.r.tobytes()==s.r.tobytes() and out.T[1:].tobytes()==s.T[1:].tobytes()
    np.testing.assert_array_equal(out.enclosed[1:],s.enclosed[1:]) if delta==0 else \
        np.testing.assert_allclose(out.enclosed[1:],s.enclosed[1:],rtol=1e-15)
    functional=s.energies(EOS)['total']-out.energies(EOS)['total']
    assert abs(functional-cavity.sink_energy(s,delta,EOS))<=1e-13*scale(s)

def test_defect_work_equals_functional_difference(seed):
    s,_=seed
    for mass in (1e-6*MSUN,1e-3*MSUN):
        t=s.copy();t.cavity_mass-=mass
        assert abs(t.energies(EOS)['total']-s.energies(EOS)['total']-cavity.defect_work(s,mass))<=1e-13*scale(s)
        assert cavity.defect_work(s,mass)==pytest.approx(G*mass*np.sum(s.weights/s.r),rel=1e-14)

@pytest.mark.parametrize('epsilon',[1e-4,-1e-4,1e-1,-1e-1])
def test_piston_and_outer_boundary_work_equal_functional_differences(seed,epsilon):
    """Moving face 0 (or face N) by epsilon of the adjacent cell's thickness at a chosen
    face pressure p*: that cell's internal energy takes the p* dV work, the potential
    re-faces; E(after) - E(before) is the secant booking."""
    s,_=seed;E0=s.energies(EOS)['total']
    for face in (0,-1):
        t=s.copy();cell=0 if face==0 else -1
        t.r=s.r.copy();t.r[face]=s.r[face]+epsilon*(s.r[1]-s.r[0] if face==0 else s.r[-1]-s.r[-2])
        dV=4*PI/3*(t.r[face]**3-s.r[face]**3)
        sign=-1 if face==0 else 1  # face 0 moving out shrinks cell 0
        t.cell_volume=s.volume.copy();t.cell_volume[cell]=s.volume[cell]+sign*dV
        pstar=EOS.evaluate(s.rho,s.T)['p'][cell]*1.1
        u=EOS.evaluate(s.rho[[cell]],s.T[[cell]])['u'][0]-sign*pstar*dV/s.dm[cell]
        t.T=s.T.copy();t.T[cell]=EOS.temperature(np.array([s.dm[cell]/t.cell_volume[cell]]),np.array([u]))[0]
        change=t.energies(EOS)['total']-E0
        if face==0:
            booked=cavity.piston_work(pstar,t.r[0],s.r[0],s.weights[0],s.enclosed[0])
            assert booked==pytest.approx(cavity.piston_work(pstar,t.r[0],s.r[0],s.weights[0],s.enclosed[0],dr=t.r[0]-s.r[0]),rel=1e-15)
        else:
            booked=-pstar*dV-G*s.weights[-1]*s.enclosed[-1]*(1/t.r[-1]-1/s.r[-1])
        assert abs(change-booked)<=1e-13*scale(s)

def test_slice_booking_and_swept_volume_against_moving_flux(seed):
    s,_=seed;E0=s.energies(EOS)['total']
    mass=1e-3*s.dm[0]
    cut,info=cavity.slice_base(s,mass,EOS)
    # The slab is the exact remainder of the rounded cut cell: exact to an ulp of dm_0.
    assert abs(info['mass']-mass)<=2*np.spacing(s.dm[0]) and cut.cavity_mass==s.cavity_mass+info['mass']
    assert math.fsum(cut.dm)+cut.cavity_mass==pytest.approx(math.fsum(s.dm)+s.cavity_mass,rel=1e-15)
    assert abs(cut.energies(EOS)['total']+info['energy']-E0)<=1e-13*scale(s)
    assert cut.T[0]==s.T[0] and cut.T[1:].tobytes()==s.T[1:].tobytes() and cut.r[1:].tobytes()==s.r[1:].tobytes()
    # The base face moved through gas at rest: the Eulerian swept volume and mass
    # through outer_interface.moving_flux at the mid radius (midpoint quadrature in time).
    dt=1e6;u=(cut.r[0]-s.r[0])/dt;rmid=.5*(cut.r[0]+s.r[0])
    q=EOS.evaluate(s.rho[:1],s.T[:1])
    face=InterfaceState(radius=rmid,enclosed_mass=s.enclosed[0],density=s.rho[0],velocity=0.,
                        gas_pressure=q['pg'][0],gas_internal_energy=q['u'][0]-A*s.T[0]**4/s.rho[0],
                        radiation_energy=A*s.T[0]**4,radiation_flux=0.,radiation_pressure_rr=A*s.T[0]**4/3,
                        potential=-G*s.enclosed[0]/rmid,energy_gauge='newtonian_zero_at_infinity')
    flux=moving_flux(face,face_velocity=u)
    assert (cut.r[0]-s.r[0])/s.r[0]<2e-3
    assert abs(-flux['mass']*dt/info['mass']-1)<=2e-7
    assert abs(-flux['mass']*dt/s.rho[0]/info['volume']-1)<=2e-7
    # Returned gas: merged into cell 0 with the cell-0 state, out of the cavity.
    back,binfo=cavity.slice_base(s,-mass,EOS)
    assert back.cavity_mass==s.cavity_mass-mass and back.dm[0]==pytest.approx(s.dm[0]+mass,rel=1e-15)
    assert abs(back.energies(EOS)['total']+binfo['energy']-E0)<=1e-13*scale(s) and binfo['energy']<0
    with pytest.raises(ValueError,match='exceeds the cavity'):cavity.slice_base(s,-2*s.cavity_mass,EOS)
    # A sliver below MINIMUM_SLICE of the base cell is left to the implicit sink.
    same,zero=cavity.slice_base(s,.5*cavity.MINIMUM_SLICE*s.dm[0],EOS)
    assert zero['mass']==0. and same.r.tobytes()==s.r.tobytes() and same.cavity_mass==s.cavity_mass
    # Whole cells leave with an exact Fraction remainder.
    edges=[Fraction(0)]
    for x in s.dm:edges.append(edges[-1]+Fraction(float(x)))
    big=float(edges[3])+.25*s.dm[3]
    cut,info=cavity.slice_base(s,big,EOS)
    assert len(cut.dm)==len(s.dm)-3 and cut.dm[0]==pytest.approx(.75*s.dm[3],rel=1e-12)
    assert abs(cut.energies(EOS)['total']+info['energy']-E0)<=1e-13*scale(s)
