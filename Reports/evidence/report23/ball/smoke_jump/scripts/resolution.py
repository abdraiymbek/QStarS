"""Truncation-error scale of the outer opacity peak: the same state brought to thermal equilibrium (a ladder of
steps up to 1 yr, no further split/merge) on meshes whose outer cells are merged (the event), unchanged, split in
2, 4 and 8; T compared at fixed mass below the surface."""
from common import *
from locate_lib import internals
from qstars.mesh import split_merge
MSUN=1.98847e33
state,config,eos,opacity=load()
from refine_lib import smooth_refine
def refine(factor,first=140):
    return lambda post:smooth_refine(post,factor,first)
def _old_refine(factor,first=140):
    def hook(post):
        target=post.dm.copy();target[first:]=post.dm[first:]/(factor*1.0001)
        out,info=split_merge(post,eos,target,split=1.5,merge=.01)
        return out
    return hook
def ladder(rm,factor):
    st=state;hook=refine(factor) if factor>1 else None;first=True
    for dt in [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.]:
        I=internals(st,dt,eos,opacity,config,remesh=(rm and first),post_hook=hook if first else None)
        st,_=qs._assemble(I['P'],I['x'],I['info'],I['cut_info'],dict(I['remesh'],regrid=I['regrid'],effect=None),I['passes'],I['pred']['swept'])
        first=False
    return st
res={}
for label,rm,f in (('merge',True,1),('pre-event',False,1),('outer x2',False,2),('outer x4',False,4),('outer x8',False,8)):
    try:st=ladder(rm,f)
    except Exception as e:print(label,'FAILED',e);continue
    dm=st.dm;M=dm.sum();below=M-(np.cumsum(dm)-.5*dm)
    res[label]=(below,st.T,st.rho)
    print(f'{label:12s} n={len(dm)} L_N/Lsun={st.face_luminosity[-1]/3.828e33:.6g} R/Rsun={st.r[-1]/6.957e10:.6g}')
probe=np.array([.25,.15,.10,.08,.0684,.0601,.05,.04,.03,.02,.01])*MSUN
ref=list(res)[-1]
print(f"{'m_below/Msun':>12} "+' '.join(f'{k:>17}' for k in res))
for m in probe:
    vals=[np.exp(np.interp(m,b[::-1],np.log(T[::-1]))) for b,T,_ in res.values()]
    print(f'{m/MSUN:12.4g} '+' '.join(f'{v:8.1f}({np.log(v/vals[-1]):+.4f})' for v in vals))
print(f'(parenthesis: ln(T/T[{ref}]) at the same mass below the surface; log-linear interpolation between cell centres)')
np.savez(SCR/'resolution_profiles.npz',**{f'{k}_{q}':v for k,(b,T,rho) in res.items() for q,v in (('below',b),('T',T),('rho',rho))})
