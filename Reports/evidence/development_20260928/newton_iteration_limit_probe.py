import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.optimize import OptimizeResult
from qstars import numerics,ledger
from qstars.config import resolve
from qstars.evolution import physics,attempt_step
from qstars.structure import initialize
from qstars.constants import YEAR
c=resolve({'numerics':{'zones':2048,'mesh':'log_radius_cosine','solver_tolerance':1e-10},'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4},'feeding':{'rate_msun_per_year':.01,'thermodynamics':'surface'},'wind':{'model':'prescribed','rate_msun_per_year':.008}})
e,op=physics(c);s,_=initialize(c,e);ledger.initialize(s,e)
def debug(fun,x,tolerance,max_iterations=40,jacobian=None):
 n=2048
 for it in range(220):
  f=fun(x);j=jacobian(fun,x);dx=spsolve(j,-f)
  print('iteration',it,'residuals',[(float(max(abs(a))),int(np.argmax(abs(a)))) for a in (f[:n],f[n:2*n],f[2*n:])],'dx',max(abs(dx)),flush=True)
  norm=np.linalg.norm(f)
  for k in range(25):
   try:g=fun(x+dx*2.**(-k))
   except ValueError:continue
   if np.linalg.norm(g)<norm*(1-1e-4*2.**(-k)):
    print('step',k,flush=True);x+=dx*2.**(-k);break
  else:break
 return OptimizeResult(x=x,message='debug',nfev=0)
numerics.newton=debug
try:attempt_step(s,1e-5*YEAR,e,op,c)
except Exception as exc:print(type(exc),exc)
