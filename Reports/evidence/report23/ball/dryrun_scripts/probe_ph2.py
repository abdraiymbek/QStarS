from dry import *
config,eos,opacity=setup()
M=1e4*MSUN
for L,R,g in ((2.195e8,3030,None),(3.48e8,3500,None),(3.48e8,3030,None),(6.2e8,2000,None),(3.5e9,1000,None),(4.0e7,6000,None),(2.195e8,3030,1e-9),(3.48e8,3500,1e-8),(4.0e7,6000,1e-13)):
    t=time.perf_counter()
    try:
        rho=photosphere.solve(L*LSUN,R*RSUN,M,eos,opacity,guess=g)
        out=photosphere.evaluate(rho,L*LSUN,R*RSUN,M,eos,opacity)
        print(L,R,g,'rho_ph',rho,'Gamma',float(out['eddington_factor'].real),'res',float(out['residual'].real),'%.3fs'%(time.perf_counter()-t))
    except Exception as e:
        print(L,R,g,type(e).__name__,e,'%.3fs'%(time.perf_counter()-t))
