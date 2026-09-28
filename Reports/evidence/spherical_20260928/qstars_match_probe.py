exec(open('/tmp/qstars_wind_probe.py').read().split('rows=[]')[0].replace('spherical_wind_probe','spherical_matching_probe'))
from scipy.integrate import solve_ivp
micro=FreeEOSOpacity(eos,mesa);w=DiffusionWind(micro,1e4*MSUN,1.152e15,3e-11,16000.)
x=-.01;e=1e-7
sol=solve_ivp(w.rhs,[-e,x],-e*np.array([w.slope,w.sonic_temperature_slope]),method='DOP853',rtol=1e-10,atol=1e-12,max_step=.001)
r,v,T,rho,q,L,*_=w.state(x,sol.y[:,-1]);flux=w.boundary_flux(x,sol.y[:,-1]);io.json_write(out/'base.json',flux)
try:
 z,s=match_base(micro,mass=w.mass,radius=r,density=rho,temperature=T,luminosity=L+w.mdot*(v*v/2+q.h),sonic_guess=[w.rs*1.0001,w.rhos*1.001,w.ts*1.0001],rtol=1e-9)
 s.update(mass_flux_error=z.mdot/w.mdot-1,sonic_radius_error=z.rs/w.rs-1);io.json_write(out/'match.json',s);print({k:v for k,v in s.items() if k!='trials'},flush=True)
except Exception as e:
 io.json_write(out/'failure.json',{'reason':str(e),'traceback':traceback.format_exc()});print(traceback.format_exc(),flush=True)
io.json_write(out/'provenance_check.json',{'source_snapshot_verified':io.finish_provenance(out)})
