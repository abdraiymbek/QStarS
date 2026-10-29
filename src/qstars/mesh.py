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

def remap(state,new_dm,eos,kinetic=None,offset=0.,copy_unchanged=False,return_slab=False):
    """Conservative remap onto new_dm (mass-ordered from the base).

    Defaults reproduce the original remap bitwise. Options (see _remap):
    kinetic (None: state.kinetic_tracked) includes the face kinetic energy;
    offset > 0 slices the base into a transit central mass (cavity_mass;
    the slab is the exact Decimal remainder, offset is checked against it
    to 1e-13), offset < 0 returns gas with the cell-0 state;
    copy_unchanged copies identical boundary-aligned cells bitwise and keeps
    the (rho, T) of the cut remainder (new cell 0 of a slice or of returned
    gas, when it lies inside one old cell); return_slab also
    returns the slab integrals and functional energy differences.
    """
    kinetic = state.kinetic_tracked if kinetic is None else bool(kinetic)
    if not (kinetic and state.kinetic_tracked and offset == 0 and not copy_unchanged and not return_slab):
        return _remap(state,new_dm,eos,kinetic,float(offset),copy_unchanged,return_slab)
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

def _sweep(a,b,cell_integrals):
    """Overlap sweep between Decimal boundary lists with a[0]=b[0], a[-1]=b[-1].

    Same arithmetic as overlap_integrals; also returns, for each new cell, the
    old cell containing it entirely (-1 when it spans an old face). The
    caller holds the 60-digit Decimal context.
    """
    import math
    out=np.zeros((len(b)-1,cell_integrals.shape[1]));inside=np.full(len(b)-1,-1)
    i=0
    for j in range(len(b)-1):
        contributions=[];owners=[]
        while i<len(a)-1 and a[i]<b[j+1]:
            overlap=min(a[i+1],b[j+1])-max(a[i],b[j])
            if overlap>0:
                contributions.append(cell_integrals[i]*float(overlap/(a[i+1]-a[i])));owners.append(i)
            if a[i+1]<=b[j+1]:i+=1
            else:break
        if len(owners)==1:inside[j]=owners[0]
        for k in range(out.shape[1]):out[j,k]=math.fsum(x[k] for x in contributions)
    return out,inside

def _total_energy(state,eos,kinetic):
    import copy as _copy
    s=_copy.copy(state);s.kinetic_tracked=kinetic
    return s.energies(eos)["total"]

def _cell_parts(state,eos,kinetic):
    """Per-cell split of State.energies: dm*(u + face-average potential + face kinetic)."""
    q=eos.evaluate(state.rho,state.T)
    m,r=state.enclosed,state.r
    gravity=-G*.5*(m[:-1]/r[:-1]+m[1:]/r[1:])
    kin=.25*(state.v[:-1]**2+state.v[1:]**2) if kinetic else np.zeros(len(state.dm))
    return q,gravity,kin

def _returned_gas(state,mass):
    """Prepend gas of the given mass below the base with the cell-0 (rho, T).

    The gas comes out of the central mass (cavity_mass, as transit), so the
    enclosed masses of all existing faces are unchanged. The cavity cannot
    return more than it holds.
    """
    if mass>state.cavity_mass*(1+1e-13):raise ValueError("Returned gas exceeds the cavity mass")
    ext=state.copy()
    volume=mass*state.volume[0]/state.dm[0]
    inner=state.r[0]**3-3*volume/(4*PI)
    if not inner>0:raise ValueError("Returned gas does not fit inside the inner radius")
    ext.dm=np.r_[mass,state.dm];ext.cell_volume=np.r_[volume,state.volume]
    ext.r=np.r_[np.cbrt(inner),state.r];ext.v=np.r_[state.v[0],state.v];ext.T=np.r_[state.T[0],state.T]
    ext.cavity_mass=state.cavity_mass-mass
    if state.face_luminosity is not None:ext.face_luminosity=np.r_[state.face_luminosity[0],state.face_luminosity]
    if state.radius_compensation is not None:ext.radius_compensation=np.r_[0.,state.radius_compensation]
    return ext

def _remap(state,new_dm,eos,kinetic,offset,copy_unchanged,return_slab):
    """Generalised remap. Mass coordinates are exact 60-digit Decimal sums.

    Reconstructed cells conserve volume and their overlap total energy
    u + face-average potential (+ face kinetic), with T by EOS inversion.
    With copy_unchanged, cells that are bitwise identical and aligned with
    the surface (or, without a slice, with the base) are copied bitwise, and
    new cell 0 of a slice (or of returned gas) keeps its old cell's (rho, T)
    when it lies inside that cell: the fixed-(rho, u) cut of a base slice.
    All other cells conserve E and V, including cells that only split an old
    cell, so an ordinary rezone books nothing. A sliced slab (offset > 0) joins
    cavity_mass, so the enclosed masses of the remaining gas do not change;
    gas returned with offset < 0 comes out of cavity_mass with the cell-0
    state. The slab energy is the functional difference E(before)-E(after)
    of State.energies. remap_entropy_change books everything except the slab
    entropy; energy bookings are left to the caller.
    """
    import bisect,decimal
    new_dm=np.asarray(new_dm,dtype=float)
    slicing=offset>0
    old=_returned_gas(state,-offset) if offset<0 else state
    n,m=len(old.dm),len(new_dm)
    if (np.any(new_dm<=0) or not (offset if slicing else 0.)<old.dm.sum() or
            not np.isclose(new_dm.sum()+(offset if slicing else 0.),old.dm.sum(),rtol=1e-13)):
        raise ValueError("Remap must preserve total mass with positive cells")
    suffix=prefix=0
    if copy_unchanged:
        while suffix<min(n,m) and new_dm[m-1-suffix]==old.dm[n-1-suffix]:suffix+=1
        if not slicing:
            while prefix<min(n,m)-suffix and new_dm[prefix]==old.dm[prefix]:prefix+=1
    top,ntop=n-suffix,m-suffix
    if offset==0 and prefix==top==ntop:
        out=state.copy()
        info=dict(mass=0.,volume=0.,internal=0.,entropy=0.,overlap_energy=0.,energy=0.,reconstruction_defect=0.,
                  copied_cells=m,fixed_state_cells=0,reconstructed_cells=0)
        return (out,info) if return_slab else out
    if not slicing and prefix==top==ntop:
        # Returned gas placed as its own cell below an otherwise identical mesh.
        out=old.copy();out.reference_mass_fractions=new_dm/new_dm.sum()
        q,gravity,kin=_cell_parts(old,eos,kinetic)
        slab=-np.array([old.dm[0]*(q["u"][0]+gravity[0]+kin[0]),old.volume[0],old.dm[0]*q["s"][0],old.dm[0]*q["u"][0]])
        return _remap_result(state,out,eos,kinetic,return_slab,slab,offset,prefix+suffix,0,0)
    if not slicing and not (prefix<top and prefix<ntop):
        raise ValueError("Remap must preserve total mass with positive cells")
    q,gravity,kin=_cell_parts(old,eos,kinetic)
    integrands=np.column_stack([old.dm*(q["u"]+gravity+kin),old.volume,old.dm*q["s"],old.dm*q["u"]])
    with decimal.localcontext() as context:
        context.prec=60
        D=decimal.Decimal.from_float
        a=[decimal.Decimal(0)]
        for x in old.dm:a.append(a[-1]+D(float(x)))
        if slicing:
            base=a[top]
            for x in new_dm[:ntop]:base-=D(float(x))
            if not base>0:raise ValueError("Base slice must remove a positive mass")
            b=[decimal.Decimal(0),base]
        else:b=[a[prefix]]
        for x in new_dm[prefix:ntop]:b.append(b[-1]+D(float(x)))
        b[-1]=a[top]
        integrals,inside=_sweep(a[prefix:top+1],b,integrands[prefix:top])
        inside[inside>=0]+=prefix
        if slicing:
            slab,integrals,inside=integrals[0],integrals[1:],inside[1:]
            owner=min(bisect.bisect_right(a,base)-1,n-1)
            if base==a[owner]:r_base=old.r[owner]
            else:
                fraction=float((base-a[owner])/(a[owner+1]-a[owner]))
                r_base=np.cbrt(old.r[owner]**3+3*fraction*old.volume[owner]/(4*PI))
            slab_mass=float(base)
        else:
            slab=-integrands[0] if offset<0 else np.zeros(4)
            r_base=old.r[prefix];slab_mass=offset
        faces=np.array([float(x) for x in b[len(b)-(ntop-prefix)-1:-1]])
        old_faces=np.array([float(x) for x in a])
    E,V=integrals[:,0],integrals[:,1]
    out=old.copy()
    out.dm=new_dm.copy()
    if slicing:out.cavity_mass=old.cavity_mass+slab_mass
    out.reference_mass_fractions=new_dm/new_dm.sum()
    volume=np.empty(m);volume[:prefix]=old.volume[:prefix];volume[ntop:]=old.volume[top:];volume[prefix:ntop]=V
    out.cell_volume=volume
    r=np.empty(m+1);r[:prefix]=old.r[:prefix];r[ntop:]=old.r[top:]
    r[prefix]=r_base
    r[prefix+1:ntop]=np.cbrt(r_base**3+3*np.cumsum(V)[:-1]/(4*PI))
    out.r=r
    # Faces kept bitwise keep their velocity, luminosity and sub-ULP compensation.
    kept=slice(0,prefix+1 if not slicing else prefix)
    v=np.empty(m+1);v[kept]=old.v[kept];v[ntop:]=old.v[top:]
    moved=slice(len(range(m+1)[kept]),ntop)
    v[moved]=np.interp(faces[moved.start-prefix:],old_faces,old.v)
    out.v=v
    if copy_unchanged and old.face_luminosity is not None:
        L=np.empty(m+1);L[kept]=old.face_luminosity[kept];L[ntop:]=old.face_luminosity[top:]
        L[moved]=np.interp(faces[moved.start-prefix:],old_faces,old.face_luminosity)
        out.face_luminosity=L
    else:out.face_luminosity=None
    compensation=np.zeros(m+1)
    if old.radius_compensation is not None:
        compensation[kept]=old.radius_compensation[kept];compensation[ntop:]=old.radius_compensation[top:]
    out.radius_compensation=compensation
    T=np.empty(m);T[:prefix]=old.T[:prefix];T[ntop:]=old.T[top:]
    # Only the cut remainder keeps its (rho, u): new cell 0 of a slice, or of
    # returned gas, lying inside one old cell. Every other changed cell is rebuilt.
    fixed=np.zeros(ntop-prefix,bool)
    if copy_unchanged and (slicing or offset<0) and prefix==0 and fixed.size and inside[0]>=0:fixed[0]=True
    T[prefix:ntop][fixed]=old.T[inside[fixed]]
    if not np.all(fixed):
        mm,rr=out.enclosed,out.r
        new_gravity=-G*.5*(mm[:-1]/rr[:-1]+mm[1:]/rr[1:])
        new_kinetic=.25*(out.v[:-1]**2+out.v[1:]**2) if kinetic else np.zeros(m)
        cells=np.arange(prefix,ntop)[~fixed]
        internal=E[~fixed]/new_dm[cells]-new_gravity[cells]-new_kinetic[cells]
        if np.any(internal<=0):raise ValueError("Conservative remap would require negative internal energy")
        T[cells]=eos.temperature(new_dm[cells]/V[~fixed],internal)
    out.T=T
    return _remap_result(state,out,eos,kinetic,return_slab,slab,slab_mass,prefix+suffix,int(np.sum(fixed)),int(np.sum(~fixed)))

def _remap_result(state,out,eos,kinetic,return_slab,slab,slab_mass,copied,fixed,reconstructed):
    before,after=_total_energy(state,eos,kinetic),_total_energy(out,eos,kinetic)
    s_before=float(np.dot(state.dm,eos.evaluate(state.rho,state.T)["s"]))
    s_after=float(np.dot(out.dm,eos.evaluate(out.rho,out.T)["s"]))
    if out.ledger:out.ledger["remap_entropy_change"]+=s_after+float(slab[2])-s_before
    if not return_slab:return out
    info=dict(mass=slab_mass,volume=float(slab[1]),internal=float(slab[3]),entropy=float(slab[2]),
              overlap_energy=float(slab[0]),energy=before-after,reconstruction_defect=after+float(slab[0])-before,
              copied_cells=int(copied),fixed_state_cells=int(fixed),reconstructed_cells=int(reconstructed))
    return out,info

# Thesis eq 2.19 with the Table 2.1 defaults. ct1 = 10**(10 c1) with c1 = 9.99
# exceeds every stellar pressure, so p + ct1 is a constant and is dropped (thesis
# 2.3.1); only the denominator goes, the c5 ln(p + c9) term stays (Fig. 2.1).
STARS_Q = dict(c2=0.,c3=.05,c4=.5,c5=.15,c6=.02,c7=.45,c8=1e18,c9=1e15,c10=3e19,ct10=2e4)

def stars_q(p,T,m,r,total_mass,constants=STARS_Q):
    """STARS mesh-spacing function Q (cgs). It decreases outward; meshpoints are equally spaced in Q."""
    c=constants
    return (c["c4"]*np.log(p)+c["c5"]*np.log(p+c["c9"])+c["c2"]*np.log(p+c["c10"])
            +c["c7"]*np.log(T/(T+c["ct10"]))-np.log((m/total_mass)**(2/3)/c["c6"]+1)
            -c["c3"]*np.log(r**2/c["c8"]+1))

def q_targets(state,eos,n,constants=STARS_Q,mass_above_base=None):
    """Target cell masses for n cells equally spaced in Q on this state.

    Q is evaluated on faces: m and r are exact there; ln p and ln T are
    interpolated linearly in mass between cell centres and extrapolated
    from the two outermost cells at the boundary faces. Cell k has target
    dm_k * dQ / |Q_{k+1} - Q_k| with dQ = |Q_N - Q_0| / n. With
    mass_above_base (positions, g, measured from the base of this state) the
    targets are interpolated in ln dm* at those positions. Evaluate it on a
    converged state: a just-split, piecewise-constant profile has a
    staircase Q and would ask to undo the split.
    """
    q=eos.evaluate(state.rho,state.T)
    dm=state.dm
    def faces(x):
        inner=(dm[1:]*x[:-1]+dm[:-1]*x[1:])/(dm[:-1]+dm[1:])
        first=x[0]-(x[1]-x[0])*dm[0]/(dm[0]+dm[1])
        last=x[-1]+(x[-1]-x[-2])*dm[-1]/(dm[-2]+dm[-1])
        return np.r_[first,inner,last]
    m=state.enclosed
    Q=stars_q(np.exp(faces(np.log(q["p"]))),np.exp(faces(np.log(state.T))),m,state.r,m[-1],constants)
    dQ=-np.diff(Q)
    if not np.all(dQ>0):raise ValueError("STARS mesh function is not monotone on this state")
    target=dm*((Q[0]-Q[-1])/n)/dQ
    if mass_above_base is None:return target
    return np.exp(np.interp(mass_above_base,np.cumsum(dm)-.5*dm,np.log(target)))

def split_merge(state,eos,target,split=2.,merge=.5):
    """Split cells above split*target and merge cells below merge*target.

    One pass from the base groups the cells, measuring a group by the number
    of targets it spans, x = sum(dm/target). A cell above split is split
    into 2**k equal parts (k >= 1, dm/2**k nearest its target). A cell below
    merge joins the neighbour with the smaller x (on the left, the group
    already formed there; ties go outward), and a merged group still below
    merge absorbs the next cell, so a slice remainder or returned gas cannot
    survive as a sliver. A merged group above split is cut into 2**k
    equal-mass parts. info counts groups cut into several parts as splits
    and groups of several cells as merges.

    Every changed group conserves its volume and its total energy
    u + face-average potential (+ face kinetic). Children of one cell keep
    its density and share one specific internal energy; children of a merged
    group take overlap shares of its volume and internal energy. The change
    of the re-faced potential (+ kinetic) energy, info["potential_reallocation"],
    is taken out of u per unit mass of the group, and T follows by EOS
    inversion. Kept faces keep r, v, face luminosity and compensation
    bitwise; new faces interpolate v and luminosity linearly in mass. The
    functional energy difference (EOS inversion only) is booked to
    remap_energy_defect, the entropy change to remap_entropy_change;
    info["max_event_defect"] is the largest |EOS-inversion defect| of one
    changed group (an event; plan reconstruction_defect_per_event_relative).
    Returns (state, info).

    The reallocation is a quadrature error of the face-average potential.
    It changes the entropy of a group by -dPhi/T to first order: m/r is
    convex in mass inside a uniform-density cell, so splits cool, and merges
    across a density drop cool too. info["merge_mixing_entropy"] isolates the
    thermodynamic part of merged groups (overlap shares mixed at their summed
    internal energy and volume), which is never negative.
    """
    import math
    dm,V,r=state.dm,state.volume,state.r
    n=len(dm);target=np.broadcast_to(np.asarray(target,dtype=float),(n,))
    kinetic=state.kinetic_tracked
    x=dm/target
    groups=[];k=0  # [first cell, number of cells]
    while k<n:
        if n>1 and x[k]<merge:
            left=math.fsum(x[groups[-1][0]:k]) if groups else np.inf
            if k+1<n and x[k+1]<=left:
                groups.append([k,2]);k+=2
            else:
                groups[-1][1]+=1;k+=1
        else:
            groups.append([k,1]);k+=1
        g=groups[-1]
        while g[1]>1 and k<n and math.fsum(x[g[0]:k])<merge:
            g[1]+=1;k+=1
    plan=[]
    for first,count in groups:
        mass,cells=(dm[first],x[first]) if count==1 else (math.fsum(dm[first:first+count]),math.fsum(x[first:first+count]))
        plan.append((first,count,2**max(1,int(np.floor(np.log2(cells)+.5))) if cells>split else 1,mass))
    splits=sum(1 for p in plan if p[2]>1);merges=sum(1 for p in plan if p[1]>1)
    if not splits and not merges:
        return state.copy(),dict(splits=0,merges=0,energy_defect=0.,entropy_change=0.,merge_mixing_entropy=0.,
                                 potential_reallocation=0.,max_event_defect=0.)
    q,gravity,kin=_cell_parts(state,eos,kinetic)
    U=dm*q["u"];cell_energy=dm*(q["u"]+gravity+kin)
    new_dm=[];new_V=[];new_T=[];r_new=[r[0]];v_new=[state.v[0]]
    L=state.face_luminosity;L_new=[L[0]] if L is not None else None
    rc=state.radius_compensation;rc_new=[rc[0]] if rc is not None else None
    pending=[]  # (first new cell, parts, first old cell, cells, mass, internal-energy shares of merged groups)
    for first,count,parts,mass in plan:
        top=first+count
        if count==1 and parts==1:
            new_dm.append(dm[first]);new_V.append(V[first]);new_T.append(state.T[first])
        elif count==1:
            child,volume=dm[first]/parts,V[first]/parts
            pending.append((len(new_dm),parts,first,1,mass,None))
            for i in range(1,parts):
                fraction=i/parts
                r_new.append(np.cbrt(r[first]**3+3*i*volume/(4*PI)))
                v_new.append(state.v[first]+(state.v[top]-state.v[first])*fraction)
                if L is not None:L_new.append(L[first]+(L[top]-L[first])*fraction)
                if rc is not None:rc_new.append(0.)
                new_dm.append(child);new_V.append(volume);new_T.append(np.nan)
            new_dm.append(child);new_V.append(volume);new_T.append(np.nan)
        else:
            # Overlap shares of the merged cells in parts equal-mass children.
            a=np.r_[0.,np.cumsum(dm[first:top])];b=np.r_[np.arange(parts)*(mass/parts),a[-1]]
            w=np.clip(np.minimum(a[None,1:],b[1:,None])-np.maximum(a[None,:-1],b[:-1,None]),0.,None)/dm[None,first:top]
            volumes=w@V[first:top];shares=w@U[first:top]
            volumes[-1]=math.fsum(V[first:top])-math.fsum(volumes[:-1])
            shares[-1]=math.fsum(U[first:top])-math.fsum(shares[:-1])
            pending.append((len(new_dm),parts,first,count,mass,shares))
            faces=b[1:-1];cumulative=np.cumsum(volumes)[:-1]
            for i in range(parts-1):
                r_new.append(np.cbrt(r[first]**3+3*cumulative[i]/(4*PI)))
                v_new.append(np.interp(faces[i],a,state.v[first:top+1]))
                if L is not None:L_new.append(np.interp(faces[i],a,L[first:top+1]))
                if rc is not None:rc_new.append(0.)
            for i in range(parts):
                new_dm.append(mass/parts);new_V.append(volumes[i]);new_T.append(np.nan)
        r_new.append(r[top]);v_new.append(state.v[top])
        if L is not None:L_new.append(L[top])
        if rc is not None:rc_new.append(rc[top])
    out=state.copy()
    out.dm=np.array(new_dm);out.cell_volume=np.array(new_V);out.r=np.array(r_new);out.v=np.array(v_new)
    out.T=np.array(new_T)
    if L is not None:out.face_luminosity=np.array(L_new)
    if rc is not None:out.radius_compensation=np.array(rc_new)
    out.reference_mass_fractions=out.dm/out.dm.sum()
    mm,rr=out.enclosed,out.r
    new_gravity=-G*.5*(mm[:-1]/rr[:-1]+mm[1:]/rr[1:])
    new_kinetic=.25*(out.v[:-1]**2+out.v[1:]**2) if kinetic else np.zeros(len(out.dm))
    rho=out.dm/out.cell_volume;u=np.full(len(out.dm),np.nan)
    reallocation=0.;mixed_rows=[];mixed_u=[];mixed_old=0.
    for start,parts,first,count,mass,shares in pending:
        c=slice(start,start+parts);cells=slice(first,first+count)
        potential=np.sum(out.dm[c]*(new_gravity[c]+new_kinetic[c]))
        energy=cell_energy[first] if count==1 else math.fsum(cell_energy[cells])
        u[c]=(energy-potential)/mass
        reallocation+=float(potential-np.sum(dm[cells]*(gravity[cells]+kin[cells])))
        if shares is not None:
            if parts>1:u[c]+=shares/out.dm[c]-math.fsum(shares)/mass
            mixed_rows.extend(range(start,start+parts));mixed_u.extend(shares/out.dm[c])
            mixed_old+=math.fsum(dm[cells]*q["s"][cells])
    changed=np.isnan(out.T)
    if np.any(u[changed]<=0):raise ValueError("Conservative split/merge would require negative internal energy")
    out.T[changed]=eos.temperature(rho[changed],u[changed])
    energy_defect=out.energies(eos)["total"]-state.energies(eos)["total"]
    q_out=eos.evaluate(out.rho,out.T)
    # Per event (one changed group): the EOS-inversion defect sum dm (u(rho, T) - u_target).
    event_defect=max(abs(math.fsum(out.dm[c]*(q_out["u"][c]-u[c]))) for c in (slice(p[0],p[0]+p[1]) for p in pending))
    entropy_change=float(np.dot(out.dm,q_out["s"])-np.dot(dm,q["s"]))
    mixing=0.
    if mixed_rows:
        rows=np.array(mixed_rows);mixed=eos.temperature(rho[rows],np.array(mixed_u))
        mixing=float(math.fsum(out.dm[rows]*eos.evaluate(rho[rows],mixed)["s"])-mixed_old)
    if out.ledger:
        out.ledger["remap_entropy_change"]+=entropy_change
        out.ledger["remap_energy_defect"]=out.ledger.get("remap_energy_defect",0.)+energy_defect
    return out,dict(splits=splits,merges=merges,energy_defect=energy_defect,entropy_change=entropy_change,
                    merge_mixing_entropy=mixing,potential_reallocation=reallocation,max_event_defect=float(event_defect))
