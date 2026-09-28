#!/usr/bin/env python3
"""Compare internal profiles and a separate continuum hydrostatic solution.

This campaign concerns the bounded hot-envelope solver; it cannot certify
full quasi-star evolution or production photospheres.
"""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from qstars import io
from qstars.config import resolve
from qstars.evolution import run,physics
from qstars.structure import initialize,continuum_polytrope
from qstars.constants import PI
from qstars.validation_metrics import interpolant_norms

def profile_difference(a,b,eos):
    """Whole-mass overlap L2 plus smooth interpolation on the common center domain."""
    qa=np.r_[0,np.cumsum(a.dm)/a.dm.sum()];qb=np.r_[0,np.cumsum(b.dm)/b.dm.sum()]
    qa[-1]=qb[-1]=1.
    overlap=np.maximum(0,np.minimum(qa[1:,None],qb[None,1:])-np.maximum(qa[:-1,None],qb[None,:-1]))
    ca=(qa[:-1]+qa[1:])/2;cb=(qb[:-1]+qb[1:])/2
    lo=max(ca[0],cb[0]);hi=min(ca[-1],cb[-1])
    sample=np.unique(np.r_[np.linspace(lo,hi,1001),ca[(ca>=lo)&(ca<=hi)],cb[(cb>=lo)&(cb<=hi)]])
    ea=eos.evaluate(a.rho,a.T);eb=eos.evaluate(b.rho,b.T)
    fields={'density':(a.rho,b.rho),'temperature':(a.T,b.T),'pressure':(ea['p'],eb['p']),
            'specific_internal_energy':(ea['u'],eb['u']),'radius':(a.rc,b.rc)}
    out={'interpolation_excluded_mass_fraction':float(lo+1-hi),'fields':{}}
    for name,(x,y) in fields.items():
        lx,ly=np.log(x),np.log(y)
        linf,l2=interpolant_norms(ca,lx,cb,ly)
        out['fields'][name]={'mass_overlap_log_L2':float(np.sqrt(np.sum(overlap*(lx[:,None]-ly[None,:])**2))),
                            'interpolated_log_Linf':linf,
                            'interpolated_mass_log_L2':l2}
    # Entropy has an arbitrary zero; normalize differences by local heat capacity.
    sa=PchipInterpolator(ca,ea['s'])(sample);sb=PchipInterpolator(cb,eb['s'])(sample)
    cv=PchipInterpolator(cb,eb['cv'])(sample)
    out['entropy_difference_over_cv_Linf']=float(max(abs(sa-sb)/cv))
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--end-time',type=float,default=.1)
    ap.add_argument('--zones',type=int,nargs='+',default=[64,128,256,512]);ap.add_argument('--parent',default='results')
    ap.add_argument('--max-dt',type=float,default=.001)
    ap.add_argument('--mechanics',choices=['midpoint','backward_euler'],default='midpoint')
    ap.add_argument('--mesh',choices=['log_radius','log_radius_cosine','log_radius_surface'],default='log_radius_cosine')
    ap.add_argument('--max-wall-seconds',type=float,default=600.)
    a=ap.parse_args();out=io.create_directory(a.parent,'spatial_validation')
    base=resolve({'name':'spatial','numerics':{'mesh':a.mesh,'mechanics':a.mechanics,'nonlinear_solver':'newton_flux','solver_tolerance':1e-10,'initial_dt_years':a.max_dt,'max_dt_years':a.max_dt},
          'engine':{'model':'prescribed','rest_mass_rate_msun_per_year':2e-4},
          'feeding':{'rate_msun_per_year':.01,'thermodynamics':'surface'},'wind':{'model':'prescribed','rate_msun_per_year':.008},
          'boundary':{'pressure_cgs':1e8,'thermal':'luminosity','luminosity_erg_s':1e42},
          'evolution':{'end_time_years':a.end_time,'max_wall_seconds':a.max_wall_seconds},'output':{'profile_every_steps':100}})
    rows=[];states=[]
    eos,_=physics(base);continuum=continuum_polytrope(base,eos)
    z=np.linspace(np.log(base['initial']['inner_radius_cm']/base['initial']['radius_cm']),0,30001)
    cq=continuum.sol(z)[0];cq[0]=0.;cq[-1]=1.
    inv=PchipInterpolator(cq,np.exp(z)*base['initial']['radius_cm'])
    np.savez_compressed(out/'continuum.npz',r=np.exp(z)*base['initial']['radius_cm'],q=cq,
                        h_dimensionless=continuum.sol(z)[1],K=continuum.K)
    for n in a.zones:
        c=copy.deepcopy(base);c['name']=f'{a.mesh}_{n}';c['numerics']['zones']=n
        s,err=initialize(c,eos)
        edges=np.r_[0,np.cumsum(s.dm)/s.dm.sum()];edges[-1]=1.
        rr=inv(edges);vv=4*PI/3*np.diff(rr)*(rr[1:]**2+rr[1:]*rr[:-1]+rr[:-1]**2)
        density_error=np.log(s.rho/(s.dm/vv))
        p,r=run(c,out);final,_=io.load_checkpoint(p/'checkpoint.npz');states.append(final)
        row={'zones':n,'path':str(p),'initial_density_log_Linf':float(max(abs(density_error))),
             'initial_density_mass_log_L2':float(np.sqrt(np.dot(s.dm,density_error**2)/s.dm.sum())),
             'numerical_kinetic_heating_erg':final.ledger.get('numerical_kinetic_heating',0.),**r}
        rows.append(row);io.json_write(out/'campaign.json',{'rows':rows,'complete':False})
        print(n,r['status'],r['runtime_seconds'],flush=True)
    for i,row in enumerate(rows[:-1]):
        if row['status']==rows[-1]['status']=='completed_requested_interval':
            row['final_comparison_to_finest']=profile_difference(states[i],states[-1],eos)
    summary={'scope':'Hot pressure-confined n=3 initialized envelope, finite thermal evolution and two streams',
        'end_time_years':a.end_time,'max_dt_years':a.max_dt,'rows':rows,'complete':True,'v1_validated':False}
    io.json_write(out/'campaign.json',summary)
    print(out)
if __name__=='__main__':main()
