"""DRYRUN experiments on seed surface variants (scratch; nothing written under benchmarks/results)."""
from dry import *
from qstars.transport import luminosities
from qstars.numerics import StepRejected
import copy,sys
def grey_T(L,R,rho,dm,opacity,T0):
    lnT=np.log(T0)
    for _ in range(60):
        kap=float(np.atleast_1d(opacity(np.atleast_1d(rho),np.atleast_1d(np.exp(lnT))))[0])
        tau=.5*kap*dm/(4*PI*R**2)+2/3
        new=.25*np.log(3*L*(tau+2/3)/(4*PI*R**2*A*C))
        if abs(new-lnT)<1e-13:break
        lnT=.5*(lnT+new)
    return float(np.exp(lnT))
def variant(s,info,eos,opacity,kind):
    s=s.copy();L=info['law_luminosity'];R=s.r[-1]
    if kind=='A':
        s.T[-1]=grey_T(L,R,s.rho[-1],s.dm[-1],opacity,6000.)
    elif kind=='B':
        # layer: T_k = min(T_eos, T_grey(tau_k)) marching inward with tau from the cells above
        tau=2/3;teff4=L/(4*PI*R**2*SIGMA)
        for k in range(len(s.dm)-1,-1,-1):
            T=min(s.T[k],8000.)
            for _ in range(60):
                kap=float(np.atleast_1d(opacity(np.atleast_1d(s.rho[k]),np.atleast_1d(T)))[0])
                half=.5*kap*s.dm[k]/(4*PI*s.r[k+1]**2)
                Tn=(.75*teff4*(tau+half+2/3))**.25
                if abs(Tn/T-1)<1e-12:break
                T=np.sqrt(T*Tn)
            if T>=s.T[k]:break
            s.T[k]=T;tau+=2*half
        print('DRYRUN variant B replaced',len(s.dm)-1-k,'cells')
    return s
if __name__=='__main__':
    kind=sys.argv[1];R=float(sys.argv[2]);N=int(sys.argv[3]);nsteps=int(sys.argv[4])
    config,eos,opacity=setup({200:'C1',400:'C0',800:'C2'}[N])
    s,info=quasistatic.seed(config,eos,opacity,N,5*MSUN,1e4*MSUN,R*RSUN)
    s=variant(s,info,eos,opacity,kind)
    bc=copy.deepcopy(config);bc['transport']['convection']='off'
    L,_=luminosities(s,eos,opacity,bc,0.)
    print('DRYRUN',kind,R,N,'T_last',s.T[-4:],'L_N/L_law',L[-1]/info['law_luminosity'],flush=True)
    ledger.initialize(s,eos,cavity=True)
    dt=.01*YEAR;acc=rej=0;t0=time.perf_counter()
    for i in range(nsteps):
        t=time.perf_counter()
        try:new,d=quasistatic.bondi_step(s,dt,eos,opacity,config,mode='relax')
        except StepRejected as e:
            rej+=1;print('DRYRUN rej dt=%.3g yr %.1fs: %s'%(dt/YEAR,time.perf_counter()-t,str(e)[:160]),flush=True);dt*=.5
            if rej>40:break
            continue
        acc+=1
        change=max(np.max(abs(np.log(new.r/s.r))),np.max(abs(np.log(new.T/s.T)))) if len(new.dm)==len(s.dm) else np.inf
        print('DRYRUN acc %d dt=%.3g yr %.1fs it=%d pass=%d N=%d change=%.3g L*/LBH-1=%.3g L*=%.4g R=%.5g r0=%.5g Gam=%.3f spl=%d mrg=%d'%(acc,dt/YEAR,time.perf_counter()-t,d['newton_iterations'],d['passes'],len(new.dm),change,d['luminosity']/d['bh_luminosity']-1,d['luminosity']/LSUN,d['radius']/RSUN,d['inner_radius']/RSUN,d.get('photosphere_eddington_factor',-1),d['splits'],d['merges']),flush=True)
        s=new;dt=min(2*dt,1e8*YEAR)
    print('DRYRUN total acc',acc,'rej',rej,'wall %.1f'%(time.perf_counter()-t0))
