"""Direct MESA EOS/opacity calls via a versioned, separately built C interface.

Initial scope: fixed H/He composition (Z=0). No interpolated surrogate EOS,
no new ionization physics. The opacity and EOS routines receive identical
composition, density and temperature. One instance per process; not threaded.
"""
import ctypes as ct
import hashlib,json,os
from pathlib import Path
import numpy as np
from .constants import A
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
        if self.lib.qstars_api_version()!=3:raise ValueError('Incompatible QStarS/MESA bridge ABI')
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
