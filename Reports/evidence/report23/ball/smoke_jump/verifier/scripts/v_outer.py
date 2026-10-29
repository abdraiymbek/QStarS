from v_common import *
from qstars.mesh import q_targets
from qstars.constants import G,C,PI
state,config,eos,opacity=load()
n=len(state.dm);zones=config['numerics']['zones']
t=q_targets(state,eos,zones)
q=eos.evaluate(state.rho,state.T)
kap=np.atleast_1d(opacity(state.rho,state.T))
L=state.face_luminosity
m=state.enclosed
mbs=state.dm.sum()-np.cumsum(state.dm)+.5*state.dm
print('cell  dm/target  mbs[Msun]   T        rho        kappa   dlnT(k->k+1) dlnkap  Gamma_loc(face k+1)')
for k in range(150,n):
    Gam=kap[k]*L[k+1]/(4*PI*C*G*m[k+1])
    dT=np.log(state.T[k+1]/state.T[k]) if k+1<n else np.nan
    dk=np.log(kap[k+1]/kap[k]) if k+1<n else np.nan
    print(f'{k:4d} {state.dm[k]/t[k]:8.3f} {mbs[k]/MSUN:9.4f} {state.T[k]:9.1f} {state.rho[k]:10.3e} {kap[k]:7.3f} {dT:+8.4f} {dk:+8.4f} {Gam:7.3f}')
# targets as the remesh sees them in the step (targets(old,cut)) -> approximate with q_targets of the state itself
