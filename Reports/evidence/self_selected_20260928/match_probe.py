import os,json,time
import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp
from qstars.constants import MSUN
from qstars.eos import Composition
from qstars.freeeos_reference import FreeEOSReference
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.spherical_wind import FreeEOSOpacity,DiffusionWind,match_base
mesa=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
micro=FreeEOSOpacity(FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,mesa.radiation_constant),mesa)
w=DiffusionWind(micro,1e4*MSUN,1.152e15,3e-11,16000.)
eps=1e-7;x=-.03
sol=solve_ivp(w.rhs,[-eps,x],-eps*np.array([w.slope,w.sonic_temperature_slope]),method='DOP853',rtol=1e-10,atol=1e-12,max_step=.001)
r,v,T,rho,q,L,*_=w.state(x,sol.y[:,-1])
base=dict(mass=w.mass,radius=r,density=rho,temperature=T,luminosity=L+w.mdot*(v*v/2+q.h))
rows=[]
for factors in [[1.001,1.01,1.001],[1.01,.7,.98],[1.001,1.2,1.01],[1.01,1.1,1.03],[1.003,.8,1.005],[1.01,.9,1.01]]:
 row=dict(factors=factors);tick=time.perf_counter()
 try:
  p,s=match_base(micro,**base,sonic_guess=np.array([w.rs,w.rhos,w.ts])*factors,rtol=1e-9,match_tolerance=1e-8,max_nfev=250)
  row.update(status='matched',mdot_error=p.mdot/w.mdot-1,report=s)
 except Exception as exc:row.update(status='failed',error=str(exc),trials=getattr(exc,'trials',[]))
 row['seconds']=time.perf_counter()-tick;rows.append(row)
 print(factors,row['status'],row.get('error',''),row.get('mdot_error',''),len(row.get('trials',row.get('report',{}).get('trials',[]))),flush=True)
 Path('Reports/evidence/self_selected_20260928/match_probe.json').write_text(json.dumps(rows,indent=2)+'\n')
