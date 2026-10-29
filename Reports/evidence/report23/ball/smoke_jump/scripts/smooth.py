from common import *
state,config,eos,opacity=load()
for rho in (2.3e-11,3.66e-11,5.0e-11):
    lnT=np.linspace(np.log(6400.),np.log(13500.),4001);T=np.exp(lnT);r=np.full_like(T,rho)
    k=np.real(opacity(r,T));q=eos.evaluate(r,T)
    for name,f in (('ln kappa',np.log(k)),('ln p',np.log(np.real(q['p']))),('u',np.real(q['u'])),('grad_ad',np.real(q['grad_ad'])),('cp',np.real(q['cp']))):
        d1=np.diff(f)/np.diff(lnT);d2=np.diff(d1)
        # a jump shows as an isolated spike of |d2| far above its neighbours
        spike=np.max(abs(d2))/np.median(abs(d2)+1e-300)
        print(f'rho={rho:.3g} {name:8s}: max|d/dlnT|={np.max(abs(d1)):.4g}  max|second difference|/median={spike:.3g}  at T={T[1+int(np.argmax(abs(d2)))]:.0f}')
    fl=opacity.clamp_flags(r,T);print('   clamp flags any:',bool(np.any(fl)))
