from common import *
from qstars import cavity,photosphere
from qstars.constants import PI,C,G
from qstars.transport import _mixing_length_faces,luminosities
from qstars.mlt import closure,radiative_gradient
MSUN=1.98847e33
def internals(old,dt_years,eos,opacity,config,remesh=True,post_hook=None):
    s=qs.settings(config);q_=s['quasistatic'];zones=config['numerics']['zones']
    dt=dt_years*YEAR
    pred=qs.predict(old,dt,eos,config,'evolve',None,s)
    swept=pred['swept']
    for passes in range(1,4):
        cut,ci=cavity.slice_base(old,swept,eos)
        premesh,regrid=qs.regrid_base(old,cut,eos,zones,q_['base_mass_fraction'])
        if remesh:post,rm=qs.remesh_base(old,premesh,eos,zones,q_['split'],q_['merge'],protect=regrid['region_cells'])
        else:post,rm=premesh,dict(splits=0,merges=0,energy_defect=0.,entropy_change=0.,merge_mixing_entropy=0.,potential_reallocation=0.,passes=0,max_event_defect=0.)
        if post_hook:post=post_hook(post)
        P=qs.Problem(old,post,dt,eos,opacity,config,'evolve',None,pred,s)
        x,info=P.solve()
        if not info['success']:raise RuntimeError(info['message'])
        delta=P.delta(x)
        if abs(delta)<=q_['pass_fraction']*post.dm[0]:break
        swept+=delta
    return dict(pred=pred,cut=cut,cut_info=ci,premesh=premesh,regrid=regrid,post=post,remesh=rm,P=P,x=x,info=info,passes=passes)

def table(st,eos,opacity,config,lo=158,out=None):
    q=eos.evaluate(st.rho,st.T);kap=np.real(opacity(st.rho,st.T))
    n=len(st.dm);L=st.face_luminosity
    faces,mass=_mixing_length_faces(st,q,kap)
    grad_rad=np.real(radiative_gradient(L[1:-1],faces['kappa'],faces['p'],faces['T'],mass))
    cl=closure(config['transport']['convection'],grad_rad,**faces,alpha=config['transport']['mixing_length_alpha'])
    tau=np.cumsum((kap*st.dm/(4*PI*(.5*(st.r[1:]+st.r[:-1]))**2))[::-1])[::-1]  # optical depth to the outer face (cell-integrated)
    Medd=kap*L[1:]/(4*PI*C*G*st.enclosed[1:])
    try:flags=opacity.clamp_flags(st.rho,st.T)
    except Exception:flags=np.zeros(n)
    nabla=np.real(np.log(q['p'][1:]/q['p'][:-1]));nabla=np.real(np.log(st.T[1:]/st.T[:-1]))/nabla
    M=st.dm.sum()
    print(f"{'k':>4} {'mbelow/Msun':>11} {'dm/Msun':>9} {'T':>9} {'rho':>9} {'kappa':>9} {'tau_out':>9} {'Gam_loc':>7} {'clamp':>5} | face k+1: {'L/Lsun':>9} {'nab_rad':>8} {'nab_ad':>7} {'nabla':>8} {'conv':>4} {'effic':>8} {'A':>8}")
    mu=np.cumsum(st.dm)
    for k in range(lo,n):
        s=f"{k:4d} {(M-mu[k]+st.dm[k]*.5)/MSUN:11.4g} {st.dm[k]/MSUN:9.4g} {st.T[k]:9.5g} {st.rho[k]:9.4g} {kap[k]:9.4g} {tau[k]:9.4g} {Medd[k]:7.3f} {int(np.asarray(flags)[k]):5d} |"
        if k<n-1:
            s+=f" {L[k+1]/3.828e33:9.4g} {grad_rad[k]:8.4g} {np.real(faces['grad_ad'][k]):7.4f} {nabla[k]:8.4g} {int(cl['convective'][k]):4d} {np.real(cl.get('efficiency',np.zeros(n-1))[k]):8.3g} {np.real(cl.get('A',np.zeros(n-1))[k]):8.3g}"
        else:s+=f" L_N={L[-1]/3.828e33:9.4g} rho_ph={st.photosphere_density}"
        print(s)
