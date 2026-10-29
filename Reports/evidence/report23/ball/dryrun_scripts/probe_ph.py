from dry import *
config,eos,opacity=setup()
M=1e4*MSUN
for L,R in ((2.195e8,3030),(3.48e8,3500),(3.48e8,3030),(6.2e8,2000),(3.5e9,1000),(4.0e7,6000)):
    try:
        rho=photosphere.solve(L*LSUN,R*RSUN,M,eos,opacity)
        out=photosphere.evaluate(rho,L*LSUN,R*RSUN,M,eos,opacity)
        print(L,R,'rho_ph',rho,'Gamma',float(out['eddington_factor'].real))
    except Exception as e:
        print(L,R,type(e).__name__,e)
    # scan
    teff=float(photosphere.effective_temperature(L*LSUN,R*RSUN))
    res=[]
    for lr in np.linspace(-14,-6,33):
        try:
            o=photosphere.evaluate(10**lr,L*LSUN,R*RSUN,M,eos,opacity)
            res.append((lr,float(o['residual'].real),float(o['eddington_factor'].real)))
        except Exception as e:
            res.append((lr,type(e).__name__,None))
    print('  Teff',teff,[(round(a,2),(round(b,3) if not isinstance(b,str) else b),(round(c,3) if c is not None else None)) for a,b,c in res])
