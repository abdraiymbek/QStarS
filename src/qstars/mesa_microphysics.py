"""Direct MESA EOS/opacity calls via a versioned, separately built C interface.

Initial scope: fixed H/He composition (Z=0). No interpolated surrogate EOS,
no new ionization physics. The opacity and EOS routines receive identical
composition, density and temperature. One instance per process; not threaded.
"""
import ctypes as ct
import hashlib,json,os
from pathlib import Path
import numpy as np
from .constants import A,KB,MU,C,PI
from .eos import IdealGasRadiation

class MesaMicrophysics(IdealGasRadiation):
    source='mesa_24_08_1'
    _loaded_library=None
    def __init__(self,composition,library_path,mesa_dir):
        super().__init__(composition)
        if composition.Z != 0: raise ValueError('MESA interface currently requires Z=0; isotope abundances are explicit H1 and He4')
        path=Path(library_path).resolve()
        manifest=json.loads((path.parent/'build.json').read_text())
        if manifest['mesa_tag'] != '24.08.1': raise ValueError('MESA interface requires audited release 24.08.1')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['library_sha256']:
            raise ValueError('MESA shared library does not match its build manifest')
        identity=(str(path),str(Path(mesa_dir).resolve()))
        if self._loaded_library is not None and identity != self._loaded_library:
            raise ValueError('Different MESA installations require separate processes')
        os.environ['MESA_DIR']=identity[1]
        os.environ.setdefault('OMP_NUM_THREADS','1')
        self.lib=ct.CDLL(str(path))
        self.lib.qstars_api_version.restype=ct.c_int
        self.abi_version=self.lib.qstars_api_version()
        if self.abi_version not in (3,4,5,6):raise ValueError('Incompatible QStarS/MESA bridge ABI')
        self.lib.qstars_constants.argtypes=[ct.POINTER(ct.c_double)]
        arad=ct.c_double();self.lib.qstars_constants(ct.byref(arad));self.radiation_constant=arad.value
        self.lib.qstars_init.argtypes=[ct.POINTER(ct.c_int)]
        array=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        self.lib.qstars_eval.argtypes=[ct.c_int,array,array,ct.c_double,array,array,array,ct.POINTER(ct.c_int)]
        err=ct.c_int();self.lib.qstars_init(ct.byref(err))
        if err.value:raise RuntimeError(f'MESA initialization error {err.value}')
        type(self)._loaded_library=identity
        self.manifest=manifest
        self._cache=None
        self._mlt_cache=None

    def raw(self,rho,T):
        rho,T=np.broadcast_arrays(rho,T)
        rr=np.ascontiguousarray(rho.real.ravel(),dtype=np.float64)
        tt=np.ascontiguousarray(T.real.ravel(),dtype=np.float64)
        if np.any(~np.isfinite(rr)) or np.any(~np.isfinite(tt)) or np.any(rr<=0) or np.any(tt<=0):
            raise ValueError('MESA requires finite, positive rho and T')
        # Deliberately bounded interface; these limits do not certify every point.
        if min(rr)<1e-14 or max(rr)>1e3 or min(tt)<3e3 or max(tt)>1e8:
            raise ValueError('MESA interface outside audited query rectangle: rho=1e-14..1e3, T=3000..1e8')
        key=(rr.tobytes(),tt.tobytes())
        if self._cache is None or key!=self._cache[0]:
            out=np.empty((21,len(rr)),order='F');dr=np.empty_like(out);dt=np.empty_like(out)
            ierr=ct.c_int()
            self.lib.qstars_eval(len(rr),rr,tt,self.composition.X,out,dr,dt,ct.byref(ierr))
            if ierr.value:raise ValueError(f'MESA EOS/opacity evaluation error {ierr.value}')
            if not np.all(np.isfinite(out)):raise ValueError('MESA returned nonfinite microphysics')
            self._cache=(key,out,dr,dt)
        _,out,dr,dt=self._cache
        if np.iscomplexobj(rho) or np.iscomplexobj(T):
            out=out+1j*(dr*(rho.imag.ravel()/rr)+dt*(T.imag.ravel()/tt))
        return out.reshape((21,)+rho.shape)

    def evaluate(self,rho,T):
        rho,T=np.broadcast_arrays(rho,T)
        x=self.raw(rho,T)
        if np.any(x[18].real>1e-6):raise ValueError('MESA neutral ideal backstop is outside supported ionization physics')
        pg=np.exp(x[0]);pr=self.radiation_constant*T**4/3;p=pg+pr
        u=np.exp(x[1]);s=np.exp(x[2]);cv,cp,gamma1,grad,chir,chit=x[3:9]
        return dict(p=p,pg=pg,pr=pr,u=u,s=s,cv=cv,cp=cp,gamma1=gamma1,
                    grad_ad=grad,delta=chit/chir,cs=np.sqrt(gamma1*p/rho),h=u+p/rho,
                    chiT=chit,chiRho=chir,free_e=np.exp(x[9]),eta=x[10],mu=x[11],
                    conductive_fraction=1-x[20]/x[19])

    def temperature(self,rho,value,quantity='u'):
        rho,value=np.broadcast_arrays(rho,value)
        if quantity not in ('u','p') or np.any(value<=0):raise ValueError('EOS inversion needs positive u or p')
        # Only an initial guess comes from the analytic EOS. Every iterate is MESA.
        T=IdealGasRadiation(self.composition).temperature(rho,value,quantity)
        T=np.clip(T,3001.,9.99e7)
        for _ in range(70):
            q=self.evaluate(rho,T)
            deriv=T*q['cv'] if quantity=='u' else q['p']*q['chiT']
            step=np.clip((value-q[quantity])/deriv,-.3,.3)
            T*=np.exp(step)
            if np.max(abs(step))<1e-10:return T
        raise RuntimeError('MESA temperature inversion failed')

    def thermal_contrast(self,rho,T,q):
        grad=.5*(q['grad_ad'][:-1]+q['grad_ad'][1:])
        return np.log1p(np.diff(T)/T[:-1])-grad*np.log1p(np.diff(q['p'])/q['p'][:-1])

    def opacity(self,rho,T):return self.raw(rho,T)[19]

    def gas_pressure_limit(self,rho,T,q=None):
        """Independent dilute classical-gas screen, unaffected by Prad.

        H/He only: use charge 2 and rho/m_u as conservative upper estimates
        in the Coulomb parameter. Gamma<0.01 and eta<-10 make Coulomb and
        degeneracy corrections far smaller than the 1% acceptance threshold.
        The T cut excludes molecular/neutral-fluid interaction regimes.
        This is a validity screen, never a correction to the library output.
        """
        rho,T=np.broadcast_arrays(rho,T)
        if q is None:q=self.evaluate(rho,T)
        radius=(3*MU/(4*PI*rho))**(1/3)
        charge=1.602176634e-19*.1*C
        coupling=4*charge**2/(radius*KB*T)
        dilute=(T>=4e4)&(coupling<.01)&(q['eta']<-10)
        expected=rho*KB*T/(MU*q['mu'])
        discrepancy=abs(q['pg']/expected-1)
        return dict(applicable=dilute,coupling=coupling,relative_error=discrepancy,
                    screened_relative_error=np.where(dilute,discrepancy,0.))

    def _require_stellar_bridge(self):
        if self.abi_version<4:
            raise ValueError('Atmosphere and full MLT require a freshly built ABI 4 MESA bridge')

    def opacity_paper_reference(self,rho,T,mixture='gs98'):
        """OPAL/AF94 with the Nakauchi 7000 K switch; diagnostic only.

        Uses MESA's processed tables and interpolation, not the original
        authors' files. At Z=0 the two metal-mixture labels should largely
        agree; retaining both measures table realization differences.
        """
        if self.abi_version<6:raise RuntimeError('Paper opacity reference requires MESA bridge ABI 6')
        if mixture not in ('gs98','gn93') or self.composition.X!=.7:
            raise ValueError('Paper opacity supports X=.7, Z=0 and gs98/gn93 only')
        rho,T=np.broadcast_arrays(rho,T)
        if np.iscomplexobj(rho) or np.iscomplexobj(T):raise ValueError('Paper opacity requires real states')
        if np.any(~np.isfinite(rho+T)) or np.any(rho<=0) or np.any(T<=0):
            raise ValueError('Paper opacity requires finite positive states')
        rr=np.ascontiguousarray(rho.ravel(),dtype=np.float64)
        tt=np.ascontiguousarray(T.ravel(),dtype=np.float64)
        kap=np.empty(len(rr));fracs=np.empty((3,len(rr)),order='F');err=ct.c_int()
        array=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        f=self.lib.qstars_opacity_paper
        f.argtypes=[ct.c_int,array,array,ct.c_double,ct.c_int,array,array,ct.POINTER(ct.c_int)]
        f(len(rr),rr,tt,self.composition.X,1 if mixture=='gs98' else 2,kap,fracs,ct.byref(err))
        if err.value:raise ValueError(f'Paper opacity table/domain error {err.value}')
        target=np.stack([tt<7000,tt>=7000,np.zeros(len(tt),dtype=bool)])
        if not np.all(np.isfinite(kap)) or np.any(kap<=0) or not np.allclose(fracs,target,atol=1e-12,rtol=0):
            raise ValueError('Paper opacity has unexpected table blend or invalid value')
        return kap.reshape(rho.shape)

    def opacity_external_eos(self,rho,T,free_e,eta):
        """MESA opacity reference with caller-supplied EOS state; no derivatives.

        This does not expand the validated density/temperature domain of any
        underlying opacity table. Inspect table boundary use separately.
        """
        if self.abi_version<5:raise ValueError('External-EOS opacity requires bridge ABI 5')
        shape=np.broadcast_shapes(np.shape(rho),np.shape(T),np.shape(free_e),np.shape(eta))
        args=[np.ascontiguousarray(np.broadcast_to(v,shape).ravel(),dtype=np.float64) for v in (rho,T,free_e,eta)]
        if any(np.any(~np.isfinite(v)) for v in args) or any(np.any(v<=0) for v in args[:3]):
            raise ValueError('Invalid opacity state')
        out=np.empty_like(args[0]);array=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        f=self.lib.qstars_opacity_external;f.argtypes=[ct.c_int,array,array,ct.c_double,array,array,array,ct.POINTER(ct.c_int)]
        err=ct.c_int();f(out.size,args[0],args[1],self.composition.X,args[2],args[3],out,ct.byref(err))
        if err.value or np.any(~np.isfinite(out)) or np.any(out<=0):raise ValueError(f'MESA opacity failed: {err.value}')
        return out.reshape(shape)

    def density_from_gas_pressure(self,pgas,T):
        """Use MESA's native Pgas,T inversion without clipping negative pressure."""
        self._require_stellar_bridge()
        if not np.isfinite(pgas+T) or pgas<=0 or T<=0:
            raise ValueError('MESA density inversion needs finite positive gas pressure and T')
        f=self.lib.qstars_density_pt
        f.argtypes=[ct.c_double]*3+[ct.POINTER(ct.c_double),ct.POINTER(ct.c_int)]
        rho=ct.c_double();err=ct.c_int()
        f(float(pgas),float(T),self.composition.X,ct.byref(rho),ct.byref(err))
        if err.value or not np.isfinite(rho.value) or rho.value<=0:
            raise ValueError(f'MESA Pgas,T inversion failed: {err.value}')
        return rho.value

    def native_atmosphere(self,teff,gravity,tau=2/3,mode='varying',kappa=.34,tolerance=1e-9):
        """Auditable reference to MESA's unmodified Eddington atmosphere.

        The varying modes retain MESA's built-in nonzero starting gas pressure.
        This method is not by itself a validated QStarS evolution boundary.
        """
        self._require_stellar_bridge()
        modes={'fixed':0,'varying_constant_opacity':1,'varying':2}
        if mode not in modes:raise ValueError('Unsupported atmosphere reference mode')
        if not all(np.isfinite(v) and v>0 for v in (teff,gravity,tau,kappa,tolerance)):
            raise ValueError('Atmosphere inputs must be finite and positive')
        array=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        f=self.lib.qstars_atmosphere
        f.argtypes=[ct.c_int]+[ct.c_double]*6+[array,ct.POINTER(ct.c_int)]
        out=np.zeros(2);err=ct.c_int()
        f(modes[mode],tau,teff,gravity,self.composition.X,kappa,tolerance,out,ct.byref(err))
        if err.value or not np.all(np.isfinite(out)):
            raise ValueError(f'MESA atmosphere integration failed: {err.value}')
        T,p=np.exp(out)
        return dict(temperature=float(T),pressure=float(p),gas_pressure=float(p-self.radiation_constant*T**4/3))

    def mlt(self,inputs,model='Henyey',alpha=2.):
        """Full native MLT including blob radiative losses and AD derivatives.

        Input leading axis: chiT, chiRho, cp, g, rho, P, T, kappa, grad_rad,
        grad_ad. Fixed-composition Schwarzschild criterion equals Ledoux here.
        Models and Henyey constants are explicit MESA prescriptions.
        """
        # ABI 4 existed in exploratory builds with two different MLT layouts.
        # Refuse both rather than risk an incorrect native argument list.
        if self.abi_version<5:raise ValueError('Native MLT derivatives require bridge ABI 5')
        modes={'ML1':0,'Henyey':1}
        x=np.asarray(inputs)
        if model not in modes or not np.isfinite(alpha) or alpha<=0:raise ValueError('Invalid MLT model or alpha')
        if x.shape[0]!=10 or not np.all(np.isfinite(x)) or np.any(x[:8].real<=0):
            raise ValueError('Invalid MLT thermodynamic inputs')
        real=np.asfortranarray(x.real.reshape(10,-1),dtype=np.float64)
        key=(model,alpha,real.tobytes())
        if self._mlt_cache is None or self._mlt_cache[0]!=key:
            out=np.zeros((4,real.shape[1]),order='F');jac=np.zeros((4,10,real.shape[1]),order='F')
            array=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
            f=self.lib.qstars_mlt;f.argtypes=[ct.c_int,ct.c_int,ct.c_double,array,array,array,ct.POINTER(ct.c_int)]
            err=ct.c_int();f(modes[model],real.shape[1],alpha,real,out,jac,ct.byref(err))
            if err.value or not np.all(np.isfinite(out)) or not np.all(np.isfinite(jac)):
                raise ValueError(f'MESA MLT evaluation failed: {err.value}')
            if np.any(out[0]<np.minimum(real[8],real[9])-1e-10) or np.any(out[0]>real[8]+1e-10):
                raise ValueError('MESA MLT returned a gradient outside radiative/adiabatic limits')
            self._mlt_cache=(key,out,jac)
        _,out,jac=self._mlt_cache
        if np.iscomplexobj(x):out=out+1j*np.einsum('oin,in->on',jac,x.imag.reshape(10,-1))
        return out.reshape((4,)+x.shape[1:])

    def consistency(self,rho,T):
        """Finite-difference first-law audit; fractions identify EOS blends.

        Normalize by |dU|+P/rho, so cancellation in the first law does not
        spuriously amplify the error. This measures the supplied EOS, not
        global finite-volume energy conservation.
        """
        rho,T=np.broadcast_arrays(rho,T);q=self.evaluate(rho,T);h=1e-4
        defects=[]
        for axis in ('rho','T'):
            qp=self.evaluate(rho*np.exp(h) if axis=='rho' else rho,T*np.exp(h) if axis=='T' else T)
            qm=self.evaluate(rho*np.exp(-h) if axis=='rho' else rho,T*np.exp(-h) if axis=='T' else T)
            du=(qp['u']-qm['u'])/(2*h);ds=(qp['s']-qm['s'])/(2*h)
            rhs=du-q['p']/rho if axis=='rho' else du
            defects.append(abs(T*ds-rhs)/(abs(du)+q['p']/rho))
        return np.maximum(*defects)
