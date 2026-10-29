"""DRYRUN microphysics probe (no model built): old vs new photosphere.solve on a grid of (L, R) at 1e4 Msun,
against a brute-force scan of the residual. The grid is a fixed lattice, not chosen from any table."""
from rp_common import *
import importlib.util
from qstars.mass_transport import DomainError
spec=importlib.util.spec_from_file_location('old_photosphere',HERE+'/snap/qstars/photosphere.py')
old=importlib.util.module_from_spec(spec);old.__package__='qstars';spec.loader.exec_module(old)
config=vbt.build_config(PLAN,'C1',LIB,MESA_DIR);eos,opacity=vbt.physics(config)
M=1e4*MSUN
def call(f,L,Rr,guess=None):
    try:return ('root',f(L,Rr,M,eos,opacity,guess=guess))
    except DomainError as e:return ('gamma' if 'Gamma_ph >= 1' in str(e) else 'nobracket',str(e)[:60])
    except Exception as e:return ('error',type(e).__name__+': '+str(e)[:60])
def brute(L,Rr):
    xs=np.linspace(np.log(1.2e-14),np.log(1e-3),1500);vals=[];gam=[]
    for x in xs:
        try:
            o=photosphere._row(np.exp(x),L,Rr,M,eos,opacity);vals.append(float(np.real(o['residual'])));gam.append(float(np.real(o['eddington_factor'])))
        except ValueError:vals.append(np.nan);gam.append(np.nan)
    vals,gam=np.array(vals),np.array(gam)
    roots=[(float(np.exp(xs[k])),gam[k]) for k in range(len(xs)-1) if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]<=0<vals[k+1]]
    down=[float(np.exp(xs[k])) for k in range(len(xs)-1) if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]>0>=vals[k+1]]
    ok=np.isfinite(gam)
    return roots,down,(float(np.nanmin(gam)) if ok.any() else None),int(ok.sum())
tally={};bad=[]
for lgL in np.arange(7.0,10.01,0.25):
    for lgR in np.arange(2.7,4.11,0.1):
        L,Rr=10**lgL*LSUN,10**lgR*RSUN
        teff=float(photosphere.effective_temperature(L,Rr))
        if not 3000<teff<1e7:continue
        o,n=call(old.solve,L,Rr),call(photosphere.solve,L,Rr)
        roots,down,gmin,nok=brute(L,Rr)
        key=(o[0],n[0],'brute_roots=%d'%len(roots));tally[key]=tally.get(key,0)+1
        if o[0]=='root' and (n[0]!='root' or n[1]!=o[1]):bad.append(('old root not reproduced bitwise',lgL,lgR,o,n))
        if n[0]=='gamma' and (roots or (gmin is not None and gmin<1)):bad.append(('new says Gamma>=1 everywhere but brute force finds Gamma<1 or a root',lgL,lgR,roots,gmin))
        if n[0]=='root':
            if not any(abs(np.log(n[1]/r[0]))<0.02 for r in roots):bad.append(('new root not among brute-force roots',lgL,lgR,n,roots))
            if len(roots)>1:bad.append(('several roots; new picked %.4g'%n[1],lgL,lgR,roots))
        if n[0] in ('nobracket','error'):bad.append(('new: '+n[0],lgL,lgR,n,roots,gmin))
for k in sorted(tally):print('DRYRUN grid old=%s new=%s %s : %d'%(k+(tally[k],)))
print('DRYRUN anomalies:',len(bad))
for b in bad[:25]:print('  ',b)
