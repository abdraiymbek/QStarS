from common import *
from locate_lib import internals
from refine_lib import smooth_refine
MSUN=1.98847e33
state,config,eos,opacity=load()
for f,first,last in ((2,140,173),(2,140,170),(4,140,170),(4,140,168),(8,140,168)):
  for ladder in ([1.],[1e-4,1e-3,1e-2,1e-1,1.]):
    st=state;ok=True
    for i,dt in enumerate(ladder):
        try:
            I=internals(st,dt,eos,opacity,config,remesh=False,post_hook=(lambda p:smooth_refine(p,f,first,last)) if i==0 else None)
            st,_=qs._assemble(I['P'],I['x'],I['info'],I['cut_info'],dict(I['remesh'],regrid=I['regrid'],effect=None),I['passes'],I['pred']['swept'])
        except Exception as e:
            print(f'x{f} cells {first}-{last} ladder {ladder}: failed at dt={dt}: {e}');ok=False;break
    if ok:
        dm=st.dm;M=dm.sum();below=M-(np.cumsum(dm)-.5*dm)
        Ts=[np.exp(np.interp(m*MSUN,below[::-1],np.log(st.T[::-1]))) for m in (.1,.08,.0684,.0601,.05,.04)]
        print(f'x{f} cells {first}-{last} ladder {ladder}: ok n={len(dm)} T at m_below 0.1,0.08,0.0684,0.0601,0.05,0.04 = '+', '.join(f'{t:.1f}' for t in Ts))
        np.savez(SCR/f'refined_x{f}_{first}_{last}.npz',below=below,T=st.T,rho=st.rho)
        break
