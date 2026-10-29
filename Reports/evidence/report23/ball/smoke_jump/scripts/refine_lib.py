from common import *
from qstars.constants import PI
def smooth_refine(post,factor,first=140,last=10**9):
    """Split cells first.. into `factor` equal-mass children with ln T and ln rho log-linear in mass between
    parent centres (parent volume kept, so parent faces stay put). Diagnostic only: energy is not conserved."""
    dm,V,T,r=post.dm,post.volume,post.T,post.r
    mu=np.cumsum(dm)-.5*dm;lnT=np.log(T);lnrho=np.log(dm/V)
    new_dm=[];new_T=[];new_V=[];new_r=[r[0]];L=post.face_luminosity;new_L=[L[0]]
    for k in range(len(dm)):
        if k<first or k>last:
            new_dm.append(dm[k]);new_T.append(T[k]);new_V.append(V[k]);new_r.append(r[k+1]);new_L.append(L[k+1]);continue
        lo=mu[k]-.5*dm[k];child=dm[k]/factor
        c=lo+(np.arange(factor)+.5)*child
        t=np.exp(np.interp(c,mu,lnT));rh=np.exp(np.interp(c,mu,lnrho))
        v=child/rh;v*=V[k]/v.sum()
        cum=np.cumsum(v)
        for i in range(factor):
            new_dm.append(child);new_T.append(t[i]);new_V.append(v[i])
            if i<factor-1:
                new_r.append(np.cbrt(r[k]**3+3*cum[i]/(4*PI)));new_L.append(L[k]+(L[k+1]-L[k])*(i+1)/factor)
        new_r.append(r[k+1]);new_L.append(L[k+1])
    out=post.copy()
    out.dm=np.array(new_dm);out.T=np.array(new_T);out.cell_volume=np.array(new_V);out.r=np.array(new_r)
    out.v=np.zeros(len(out.r));out.face_luminosity=np.array(new_L)
    if post.radius_compensation is not None:out.radius_compensation=np.zeros(len(out.r))
    out.reference_mass_fractions=out.dm/out.dm.sum()
    return out
