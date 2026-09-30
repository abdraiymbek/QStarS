"""Auditable CMFGEN component interface, not a certified wind boundary.

No routine in this module selects mdot or authorizes production coupling.
Frequency data are comoving. Extinction and emissivity in the native snapshot
refer to material inside microclumps; returned exchange/force use volume means.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib, re, struct
import numpy as np
from .constants import C


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def _numbers(line):
    try:return [float(re.sub(r'(?<=[0-9.])([+-]\d{3})$', r'E\1', x.replace('D','E').replace('d','e'))) for x in line.split()]
    except ValueError:return None


def nonlinear_corrections(text):
    """Retain raw largest corrections: upstream stopping can discard ten levels.

    Percent decrease is printed relative to the old value, whereas MAXCH uses
    the new value. Missing raw records fail closed; neither a native zero nor
    a converged majority of levels certifies the complete atomic state.
    """
    rows=[];current=None;increase=None;decrease=None
    for line in text.splitlines():
        m=re.search(r'Current great iteration count is\s+(\d+)',line)
        if m:current=int(m[1]);increase=decrease=None
        m=re.search(r'Maximum % (increase|decrease) at depth\s+\d+ is\s+(\S+)',line)
        if m:
            value=_numbers(m[2])
            if m[1]=='increase':increase=value[0] if value else None
            else:decrease=value[0] if value else None
        m=re.search(r'Maximm changes as returned by SOLVEBA_V13 is\s+(\S+)',line)
        if m:
            value=_numbers(m[1]);reported=value[0] if value else None
            raw=None
            if increase is not None and decrease is not None and np.isfinite([increase,decrease]).all():
                down=100*decrease/(100-decrease) if decrease<99.999 else 1e7
                raw=max(abs(increase),down)
            complete=max(raw,abs(reported)) if raw is not None and reported is not None and np.isfinite(reported) else None
            rows.append(dict(iteration=current,reported_percent=reported,raw_percent=raw,all_level_percent=complete))
    return rows


def vector_section(lines, label, count):
    """Strict labeled vector; reject missing/truncated/nonfinite output."""
    matches=[i for i,s in enumerate(lines) if s.strip()==label]
    if len(matches)!=1:raise ValueError(f'Expected exactly one {label!r} section')
    values=[]
    for line in lines[matches[0]+1:]:
        if not line.strip():continue
        row=_numbers(line)
        if row is None:break
        values.extend(row)
        if len(values)>=count:break
    a=np.array(values)
    if a.size!=count or not np.all(np.isfinite(a)):raise ValueError(f'Malformed {label!r} vector')
    return a


def read_rvtj(path):
    lines=Path(path).read_text().splitlines()
    nd=[int(s.split(':',1)[1]) for s in lines if s.strip().startswith('ND:')]
    if len(nd)!=1 or nd[0]<2:raise ValueError('Invalid RVTJ ND')
    n=nd[0]
    labels={'radius':('Radius (10^10 cm)',1e10),'velocity':('Velocity (km/s)',1e5),
            'temperature':('Temperature (10^4K)',1e4),'electron_density':('Electron density',1),
            'density_clump':('Mass Density (gm/cm^3)',1),'filling_factor':('Clumping Factor',1),
            'hydrogen_density':('Hydrogen Density',1),'helium_density':('Helium Density',1),
            'chi_rosseland':('Rosseland Mean Opacity',1e-10),'chi_flux':('Flux Mean Opacity',1e-10)}
    q={key:vector_section(lines,label,n)[::-1]*scale for key,(label,scale) in labels.items()}
    for name in ['J','H','K']:
        label=name+' moment of radiation field'
        if any(s.strip()==label for s in lines):q[name]=vector_section(lines,label,n)[::-1]
    if not np.all(np.diff(q['radius'])>0):raise ValueError('RVTJ radius must decrease in native ordering')
    if np.any(q['density_clump']<=0) or np.any(q['temperature']<=0):raise ValueError('Nonpositive material state')
    f=q['filling_factor']
    if np.any((f<=0)|(f>1)):raise ValueError('Invalid microclump filling factor')
    q['density']=q['density_clump']*f
    q['mass_flux']=4*np.pi*q['radius']**2*q['density']*q['velocity']
    q['kappa_flux']=q['chi_flux']/q['density_clump']
    q['kappa_rosseland']=q['chi_rosseland']/q['density_clump']
    return q


def read_hhe_populations(path):
    """Read full H/He level populations and top ions from the RVTJ format."""
    lines=Path(path).read_text().splitlines()
    n=int(next(s.split(':',1)[1] for s in lines if s.strip().startswith('ND:')))
    q={}
    ions=[('HI','Number of Hydrogen levels:'),('HeI','Number of Helium I levels:'),('HeII','Number of Helium II levels:')]
    if any(s.strip().startswith('Number of Hminus levels:') for s in lines):
        ions.insert(0,('Hminus','Number of Hminus levels:'))
    for ion,label in ions:
        ii=[i for i,s in enumerate(lines) if s.strip().startswith(label)]
        if len(ii)!=1:raise ValueError(f'Missing {ion} populations')
        i=ii[0];nl=int(lines[i].split(':',1)[1]);vals=[]
        expected={'HI':'HI','HeI':'HeI','HeII':'He2','Hminus':'H0'}[ion]+' populations'
        if i<1 or lines[i-1].strip()!=expected:
            raise ValueError('RVTJ ion label mismatch; unpatched native H- output uses incorrect fixed ion indices')
        if nl<1:raise ValueError(f'Empty {ion} atom')
        for line in lines[i+2:]:
            row=_numbers(line)
            if row is None:break
            vals.extend(row)
            if len(vals)>=nl*n:break
        a=np.asarray(vals)
        if a.size!=nl*n or not np.all(np.isfinite(a)) or np.any(a<0):raise ValueError(f'Invalid {ion} populations')
        q[ion]=a.reshape(n,nl)[::-1].T
    q['HII']=vector_section(lines,'DHYD population',n)[::-1]
    q['HeIII']=vector_section(lines,'DHe2 population',n)[::-1]
    if np.any(q['HII']<0) or np.any(q['HeIII']<0):raise ValueError('Negative top ion density')
    return q


def read_cooling(path, n):
    """Native GENCOOL net = cooling minus heating, including enabled advection.

    This electron budget is a diagnostic, not automatically the total gas source.
    The native percentage is 200*net/sum(abs(individual channels)).
    """
    lines=Path(path).read_text().splitlines();result={}
    for key,label in [('net_cooling','Net Cooling Rate [ergs/cm**3/sec]'),('net_percent','Net Cooling as percentage of total.')]:
        vals=[]
        for i,s in enumerate(lines):
            if s.strip()==label:
                row=_numbers(lines[i+1])
                if row is None:raise ValueError('Malformed cooling block')
                vals.extend(row)
        a=np.array(vals)
        if a.size!=n or not np.all(np.isfinite(a)):raise ValueError('Missing/truncated cooling budget')
        result[key]=a[::-1]
    return result


@dataclass
class SpectralSnapshot:
    radius: np.ndarray
    velocity: np.ndarray
    temperature: np.ndarray
    density_clump: np.ndarray
    filling_factor: np.ndarray
    frequency: np.ndarray
    weight_hz: np.ndarray
    extinction: np.ndarray
    scattering: np.ndarray
    emissivity: np.ndarray
    J: np.ndarray
    H: np.ndarray
    K: np.ndarray

    def moments(self, *, coherent_scattering):
        if not coherent_scattering:
            raise NotImplementedError("Inelastic scattering energy redistribution is not exported by ABI 1")
        w=self.weight_hz
        integrate=lambda a: np.einsum('f,fr->r',w,a)
        J,H,K=map(integrate,[self.J,self.H,self.K])
        f=self.filling_factor;rho=self.density_clump*f
        force_density=4*np.pi/C*f*integrate(self.extinction*self.H)
        heating=4*np.pi*f*integrate((self.extinction-self.scattering)*self.J-self.emissivity)
        E,F,P=4*np.pi*J/C,4*np.pi*H,4*np.pi*K/C
        beta=self.velocity/C
        if np.any(abs(beta)>=.01):raise ValueError('Outside audited nonrelativistic interface: |v|/c >= 0.01')
        g2=1/(1-beta**2)
        Elab=g2*(E+2*beta*F/C+beta**2*P)
        Flab=g2*((1+beta**2)*F+beta*C*(E+P))
        Plab=g2*(P+2*beta*F/C+beta**2*E)
        return dict(radius=self.radius,velocity=self.velocity,temperature=self.temperature,density=rho,
          J=J,H=H,K=K,E_comoving=E,F_comoving=F,P_comoving=P,
          E_lab=Elab,F_lab=Flab,P_lab=Plab,luminosity_comoving=4*np.pi*self.radius**2*F,
          luminosity_lab=4*np.pi*self.radius**2*Flab,heating_comoving=heating,
          radiative_force_density=force_density,radiative_acceleration=force_density/rho,
          # Lorentz transformation of the radiation four-force, not a second work source.
          energy_to_matter_lab=np.sqrt(g2)*(heating+self.velocity*force_density),
          kappa_flux=np.divide(force_density*C,rho*F,out=np.full_like(F,np.nan),where=F!=0),
          production_certified=False,wind_eigenvalue_solved=False)


def read_spectral(path):
    """Read opt-in QSCMF001 stream, little-endian int32/float64, with exact length."""
    path=Path(path)
    with path.open('rb') as f:
        if f.read(8)!=b'QSCMF001':raise ValueError('Unknown CMFGEN spectral ABI')
        dims=f.read(8)
        if len(dims)!=8:raise ValueError('Truncated spectral header')
        n,nf=struct.unpack('<ii',dims)
    if n<2 or nf<2 or n>100000 or nf>10000000:raise ValueError('Invalid spectral dimensions')
    expected=16+8*(5*n+nf*(2+6*n))
    if path.stat().st_size!=expected:raise ValueError('Incomplete spectral snapshot or incompatible numeric kind')
    a=np.memmap(path,dtype='<f8',mode='r',offset=16)
    r,v,T,rho,fill=a[:5*n].reshape(5,n)[:,::-1]
    data=a[5*n:].reshape(nf,2+6*n);nu=data[:,0];weights=data[:,1]
    arrays=data[:,2:].reshape(nf,6,n)[:,:,::-1]
    if not np.all(np.isfinite(a)):raise ValueError('Nonfinite spectral data')
    if np.any(r<=0) or np.any(nu<=0) or np.any(np.diff(r)<=0) or np.any(np.diff(nu)>=0) or np.any(weights<=0):raise ValueError('Invalid spectral grid/quadrature')
    if np.any(T<=0) or np.any(rho<=0) or np.any((fill<=0)|(fill>1)):raise ValueError('Invalid spectral material state')
    # Negative extinction from inversions is retained; it must not be clipped.
    return SpectralSnapshot(r,v,T,rho,fill,nu,weights,*[arrays[:,i,:] for i in range(6)])


def production_outer(*args, **kwargs):
    raise RuntimeError('CMFGEN component integration is not a validated material-energy/sonic/envelope closure; production wind remains gated')


def read_atomic_levels(path, count):
    """CMFGEN oscillator header: g and excitation wavenumber, cm^-1.

    The selected model takes the first count full levels in the atomic file.
    No photoionization fit or collision rate is reimplemented here.
    """
    lines=Path(path).read_text().splitlines()
    ni=next(i for i,s in enumerate(lines) if '!Number of energy levels' in s)
    available=int(lines[ni].split('!')[0]);ion=float(next(s.split('!')[0] for s in lines if '!Ionization energy' in s))
    start=next(i for i,s in enumerate(lines) if '!Number of transitions' in s)+1
    vals=[]
    for line in lines[start:]:
        tokens=line.split()
        if len(tokens)<4:continue
        try:g=float(tokens[1]);energy=float(tokens[2])
        except ValueError:continue
        vals.append((g,energy))
        if len(vals)==count:break
    a=np.array(vals)
    if count>available or a.shape!=(count,2) or np.any(a[:,0]<=0) or not np.all(np.isfinite(a)):
        raise ValueError('Invalid atomic level table')
    return dict(statistical_weight=a[:,0],wavenumber=a[:,1]-a[0,1],ionization_wavenumber=ion-a[0,1],source_sha256=sha256(path))


def hhe_material(state,populations,levels):
    """Ideal translational + resolved excitation/ionization energy diagnostic.

    Neutral ground H and He have zero atomic energy. No molecules, degeneracy,
    Coulomb terms, radiation energy or turbulence are included. This is not an
    LTE FreeEOS call and supplies no equilibrium derivatives for a sonic solve.
    """
    from .nlte_material import hhe_resolved_material
    p=populations;s=state
    hm=p['Hminus'].sum(axis=0) if 'Hminus' in p else 0.
    nh=p['HI'].sum(axis=0)+p['HII']+hm;nhe=p['HeI'].sum(axis=0)+p['HeII'].sum(axis=0)+p['HeIII']
    ne=p['HII']+p['HeII'].sum(axis=0)+2*p['HeIII']-hm
    resolved=hhe_resolved_material(s,p,levels).values()
    pressure=resolved['pressure'];thermal=1.5*pressure
    atomic=resolved['atomic_energy']*s['density_clump']
    residuals={'hydrogen':np.max(abs(nh/s['hydrogen_density']-1)),
               'helium':np.max(abs(nhe/s['helium_density']-1)),
               'charge':np.max(abs(ne/s['electron_density']-1))}
    return dict(gas_pressure=pressure,gas_internal_energy=(thermal+atomic)/s['density_clump'],
        translational_energy_density_clump=thermal,atomic_energy_density_clump=atomic,
        energy_gauge='neutral_H_He_ground',conservation_relative={k:float(v) for k,v in residuals.items()},
        equilibrium_derivatives_available=False,frozen_population_derivatives_available=True,production_certified=False)
