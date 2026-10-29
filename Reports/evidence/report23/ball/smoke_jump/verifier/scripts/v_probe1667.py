"""The no-merge DRYRUN's final state (step 1667, M_BH 55.0): outer structure and single steps vs dt (unpatched code)."""
from v_common import *
ck=OUT/'DRYRUN_registered_controller_merge0_600/checkpoint.npz'
state,config=qs.load_checkpoint(ck)
eos,opacity=vbt.physics(config)
n=len(state.dm)
kap=np.atleast_1d(opacity(state.rho,state.T))
print('step',state.steps,'Mbh',state.bh_mass/MSUN,'zones',n,'merge setting',config['quasistatic']['merge'])
L=state.face_luminosity;m=state.enclosed
from qstars.constants import G,C,PI
for k in range(n-16,n):
    Gam=kap[k]*L[k+1]/(4*PI*C*G*m[k+1])
    print(f'{k:4d} T={state.T[k]:9.1f} rho={state.rho[k]:.3e} kap={kap[k]:.3f} dlnT={np.log(state.T[min(k+1,n-1)]/state.T[k]):+.4f} dlnkap={np.log(kap[min(k+1,n-1)]/kap[k]):+.3f} Gamma={Gam:.3f}')
for dt in (1e-3,1e-2,.1,.3,.5,.7,1.,2.):
    try:
        new,d,m_=attempt(state,dt,eos,opacity,config)
        dl=np.log(new.T/state.T) if len(new.T)==n else None
        print(f'dt={dt:6.3g} ok it={d["newton_iterations"]} m/s={d["merges"]}/{d["splits"]} delta={m_["delta_actual"]:.4g} roundoff_exit={d["roundoff_exit"]} cell={int(np.argmax(abs(dl))) if dl is not None else None}',flush=True)
    except Exception as e:
        print(f'dt={dt:6.3g} FAIL {str(e)[:150]}',flush=True)
