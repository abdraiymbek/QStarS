"""Auditable direct FreeEOS reference, not an enabled evolution EOS.

Gas EOS1 is supplied by FreeEOS; the blackbody radiation contribution uses
the selected radiation constant explicitly. The zero of internal energy is
the upstream FreeEOS convention. No stitching to the MESA EOS is performed.
"""
import ctypes as ct
import hashlib,json
from pathlib import Path
import numpy as np

class FreeEOSReference:
    def __init__(self,library_path,hydrogen=.7,radiation_constant=None):
        path=Path(library_path).resolve();self.manifest=json.loads((path.parent/'build.json').read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=self.manifest['library_sha256']:
            raise ValueError('FreeEOS library differs from recorded build')
        if not 0<=hydrogen<=1:raise ValueError('Invalid hydrogen mass fraction')
        self.hydrogen=hydrogen;self.lib=ct.CDLL(str(path))
        self.lib.qstars_freeeos_version.restype=ct.c_int
        if self.lib.qstars_freeeos_version()!=1:raise ValueError('Unknown FreeEOS bridge ABI')
        array=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
        self.lib.qstars_freeeos_eval.argtypes=[ct.c_int,ct.c_int,array,array,ct.c_double,array,ct.POINTER(ct.c_int)]
        self.lib.qstars_freeeos_constants.argtypes=[ct.POINTER(ct.c_double)]
        arad=ct.c_double();self.lib.qstars_freeeos_constants(ct.byref(arad))
        self.native_radiation_constant=arad.value
        self.radiation_constant=arad.value if radiation_constant is None else float(radiation_constant)
        if not np.isfinite(self.radiation_constant) or self.radiation_constant<=0:raise ValueError('Invalid radiation constant')

    def raw(self,values,T,mode=2):
        values,T=np.broadcast_arrays(values,T)
        if np.iscomplexobj(values) or np.iscomplexobj(T):raise ValueError('Reference EOS accepts real states only')
        v=np.ascontiguousarray(values.ravel(),dtype=np.float64);t=np.ascontiguousarray(T.ravel(),dtype=np.float64)
        if mode not in (1,2) or np.any(~np.isfinite(v+t)) or np.any(v<=0) or np.any(t<=0):
            raise ValueError('Reference EOS needs finite positive states and mode 1 or 2')
        out=np.zeros((24,len(v)),order='F');err=ct.c_int()
        self.lib.qstars_freeeos_eval(mode,len(v),v,t,self.hydrogen,out,ct.byref(err))
        if err.value or not np.all(np.isfinite(out)):raise ValueError(f'FreeEOS evaluation failed: {err.value}')
        return out.reshape((24,)+values.shape)

    def density_from_gas_pressure(self,pgas,T):
        return np.exp(self.raw(pgas,T,mode=1)[9])

    def evaluate(self,rho,T):
        rho,T=np.broadcast_arrays(rho,T);x=self.raw(rho,T)
        pg=np.exp(x[0]);pr=self.radiation_constant*T**4/3;p=pg+pr
        u=x[3]+3*pr/rho;s=x[6]+4*pr/(rho*T)
        cv=x[5]/T+12*pr/(rho*T)
        chir=pg*x[1]/p;chit=(pg*x[2]+4*pr)/p
        cp=cv+p*chit**2/(rho*T*chir)
        gamma1=chir+p*chit**2/(rho*T*cv)
        grad=p*chit/(rho*T*cv*gamma1)
        return dict(pg=pg,pr=pr,p=p,u=u,s=s,cv=cv,cp=cp,chiRho=chir,chiT=chit,
                    grad_ad=grad,gamma1=gamma1,free_e=x[15],eta=x[16],mu=x[17],
                    cs=np.sqrt(gamma1*p/rho),h=u+p/rho,gas_u=x[3],gas_s=x[6],
                    h_ionized=x[18],he_once=x[19],he_twice=x[20])
