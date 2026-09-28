#!/usr/bin/env python3
"""Report 10: attribution experiments, never fitted wind inputs or stellar runs."""
import argparse,csv,hashlib,itertools,json,os,resource,time,traceback
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from qstars import io
from qstars.constants import G,C,MSUN,RSUN,YEAR
from qstars.eos import Composition
from qstars.mesa_microphysics import MesaMicrophysics
from qstars.freeeos_reference import FreeEOSReference
from qstars.spherical_wind import FreeEOSOpacity,DiffusionWind

class MassNumberAbundances:
    """An explicit abundance-convention experiment, NOT the authors' EOS.

    Evaluate unmodified FreeEOS at rho'=f*rho and X'=X*A_H/f. Then nuclei per
    physical gram are X/1 and Y/4. Specific extensive quantities scale by f.
    Quantum partition masses and all FreeEOS corrections remain unchanged.
    This is not a fitted correction or a selectable production EOS.
    """
    def __init__(self,library,arad):
        self.hydrogen=.7;self.radiation_constant=arad
        self.factor=.7*1.007825031898+.3*4.002603254130/4
        self.eos=FreeEOSReference(library,.7*1.007825031898/self.factor,arad)
    def raw(self,rho,T):
        x=self.eos.raw(rho*self.factor,T).copy()
        x[3:9]*=self.factor;x[9]-=np.log(self.factor)
        x[12]*=self.factor;x[15]*=self.factor;x[17]/=self.factor
        return x


def euler_outward(w,epsilon,outer_ratio=100,sonic_offset=1e-7):
    """Paper section 2.2.3: Euler in physical r,v,T; algebraic rho and L.

    No log-variable Euler shortcut. Linear event interpolation locates the
    endpoint within the last Euler step. The sonic start uses the verified
    critical slope; exact authors' start handling is unspecified.
    """
    r=w.rs*(1+sonic_offset)
    v=w.vs*(1+w.slope*sonic_offset)
    T=w.ts*(1+w.sonic_temperature_slope*sonic_offset)
    knots=[];termination=None
    def state(r,v,T):return w.state(np.log(r/w.rs),[np.log(v/w.vs),np.log(T/w.ts)])
    def events(a):
        r,v,T,rho,q,L,*_=a
        return np.array([3*L/(4*np.pi*r*r*C*w.micro.arad*T**4)-.75,v/np.sqrt(q.ct2)-1.5])
    a=state(r,v,T)
    for i in range(1000000):
        knots.append((r,v,T));old=events(a)
        vlog=a[7]/a[8];g=a[6]
        dr=min(epsilon*r/max(abs(vlog),abs(g)),w.rs*outer_ratio-r)
        vnew=v+dr*v/r*vlog;Tnew=T+dr*T/r*g;rnew=r+dr
        anew=state(rnew,vnew,Tnew);new=events(anew)
        event=0 if old[0]<0<=new[0] else (1 if old[1]>0>=new[1] else None)
        if event is not None:
            def func(t):return events(state(r+t*dr,v+t*(vnew-v),T+t*(Tnew-T)))[event]
            frac=brentq(func,0,1,xtol=1e-12)
            r+=frac*dr;v+=frac*(vnew-v);T+=frac*(Tnew-T)
            knots.append((r,v,T));termination=['paper_photosphere','paper_stalled_mach_1p5'][event];break
        r,v,T,a=rnew,vnew,Tnew,anew
        if r>=w.rs*outer_ratio:termination='radial_bound';break
    if termination is None:raise RuntimeError('Euler step budget exceeded; no no-root inference')
    raw=np.array(knots).T
    # Preserve every step in raw file; evaluate thermodynamics on a bounded
    # subset for comparison, not for claiming independent quadrature accuracy.
    take=np.unique(np.r_[np.arange(0,len(knots),max(1,len(knots)//5000)),len(knots)-1])
    rr,vv,tt=raw[:,take]
    f=w.profile(np.log(rr/w.rs),np.array([np.log(vv/w.vs),np.log(tt/w.ts)]))
    return f,raw,dict(steps=i+1,termination=termination,epsilon=epsilon)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='results');a=p.parse_args()
    start=time.perf_counter();out=io.create_directory(a.output,'outer_discrepancy')
    design=dict(no_input_fitting=True,baseline=[100*MSUN,50*RSUN,9.71e-9,2.5e5],
        numerical_profile_log_target=.001,endpoint_relative_target=.001,
        density_rounding_half_width=.005e-9,
        conditional_printed_rounding_half_widths=[.5*MSUN,.5*RSUN,.005e-9,.05e5],
        rounding_interpretation='Deterministic sensitivity design only: not probability bounds or proof authors rounded nominal model inputs',
        euler_epsilons=[1e-3,3e-4,1e-4,3e-5],
        stage_b_allowed=False)
    io.provenance(out,design);(out/'driver.py').write_bytes(Path(__file__).read_bytes());print(out,flush=True)
    mesa=MesaMicrophysics(Composition(.7,.3,0),os.environ['QSTARS_MESA_LIBRARY'],os.environ['MESA_DIR'])
    eos=FreeEOSReference(os.environ['QSTARS_FREEEOS_LIBRARY'],.7,mesa.radiation_constant)
    io.json_write(out/'microphysics.json',dict(mesa=mesa.manifest,freeeos=eos.manifest,
        table_files=[dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in [
            Path(os.environ['MESA_DIR'])/'data/kap_data'/name for name in [
                'gs98_z0.0_x0.7.data','gn93_z0.0_x0.7.data','lowT_af94_gn93_z0.0_x0.7.data',
                'lowT_fa05_gs98_z0.0_x0.7.data']]],
        original_authors_microphysics_recovered=False,
        note='MESA AF94 file header says Alexander & Ferguson (1995); retained as distributed, not certified identical to 1994 paper files.'))
    micro={key:FreeEOSOpacity(eos,mesa,key) for key in ['default','paper_gs98','paper_gn93']}
    micro['mass_number_diagnostic']=FreeEOSOpacity(MassNumberAbundances(os.environ['QSTARS_FREEEOS_LIBRARY'],mesa.radiation_constant),mesa,'paper_gn93')
    rows=[];profiles={};winds={};comparisons=[]
    baseline=np.array(design['baseline'])
    opts=dict(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,samples=1001,stalled_mach=1.5)
    def run(label,args=baseline,family='default',changes=None,derivative_step=1e-4):
        t=time.perf_counter();options=dict(opts,**(changes or {}))
        try:
            w=DiffusionWind(micro[family],*args,derivative_step=derivative_step)
            f,s=w.integrate(**options);np.savez_compressed(out/(label+'.npz'),**f)
            profiles[label]=f;winds[label]=w
            s.update(status='completed_component' if s['solver_success'] else 'numerical_failure',
                photo_rsun=s['endpoint_radius']/RSUN,mdot_msun_year=w.mdot*YEAR/MSUN)
        except Exception as exc:
            s=dict(status='numerical_or_domain_failure',reason=str(exc),traceback=traceback.format_exc(),physical_impossibility_claimed=False)
        s.update(label=label,inputs=list(args),family=family,options=options,seconds=time.perf_counter()-t)
        rows.append(s);io.json_write(out/'cases.json',rows)
        print(label,s['status'],s.get('photo_rsun',s.get('reason')),flush=True)
    # No stopping early when a comparison improves. Every predeclared case runs.
    for family in micro:run(family,family=family)
    for rho in [9.65e-9,9.85e-9,10.5e-9]:
        run('paper_'+str(rho),np.array([100*MSUN,50*RSUN,rho,2.5e5]),family='paper_gn93')
    for label,changes in [('tight',{'rtol':1e-10}),('radial_step',{'max_log_step':.005}),
                          ('radau',{'method':'Radau'}),('sonic_offset',{'sonic_offset':1e-8}),
                          ('strict_diffusion',{'endpoint':'diffusion_limit'}),
                          ('report9_stall',{'stalled_mach':np.sqrt(1.01)})]:
        run(label,family='paper_gn93',changes=changes)
    run('critical_derivative',family='paper_gn93',derivative_step=1e-5)
    # Narrow density precision and separate one-variable effects, then all 16
    # corners of a declared, explicitly conditional rounding interpretation.
    half=np.array(design['conditional_printed_rounding_half_widths'])
    for k in range(4):
        for sign in [-1,1]:
            delta=np.zeros(4);delta[k]=sign*half[k]
            run(f'round_axis{k}_{sign:+}',baseline+delta,family='paper_gn93')
    for signs in itertools.product([-1,1],repeat=4):
        run('round_corner_'+''.join('p' if s>0 else 'm' for s in signs),baseline+half*signs,family='paper_gn93')
    for eps in design['euler_epsilons']:
        t=time.perf_counter();label=f'euler_{eps:g}'
        try:
            f,raw,s=euler_outward(winds['paper_gn93'],eps)
            np.savez_compressed(out/(label+'.npz'),**f)
            np.savez_compressed(out/(label+'_all_steps.npz'),radius=raw[0],velocity=raw[1],temperature=raw[2])
            profiles[label]=f
            s.update(status='completed_component',photo_rsun=f['radius'][-1]/RSUN,
                endpoint_velocity_over_escape=f['velocity'][-1]/f['escape_speed'][-1],label=label)
        except Exception as exc:
            s=dict(status='numerical_or_domain_failure',label=label,reason=str(exc),traceback=traceback.format_exc(),physical_impossibility_claimed=False)
        s['seconds']=time.perf_counter()-t;rows.append(s);io.json_write(out/'cases.json',rows)
        print(label,s['status'],s.get('photo_rsun',s.get('reason')),flush=True)
    base=profiles['paper_gn93']
    for label,f in profiles.items():
        x=np.linspace(max(np.log(f['radius'][0]),np.log(base['radius'][0])),min(np.log(f['radius'][-1]),np.log(base['radius'][-1])),1001)
        errors={key:float(max(abs(PchipInterpolator(np.log(f['radius']),np.log(f[key]))(x)-PchipInterpolator(np.log(base['radius']),np.log(base[key]))(x)))) for key in ['density','temperature','velocity','opacity','gas_pressure','radiative_luminosity','hydrogen_ionized']}
        comparisons.append(dict(label=label,max_profile_log_errors=errors,endpoint_relative_change=float(f['radius'][-1]/base['radius'][-1]-1)))
    io.json_write(out/'comparisons.json',comparisons)
    # Opacity substitutions at IDENTICAL states distinguish microphysics from
    # nonlinear changes to the wind solution.
    take=np.arange(0,len(base['radius']),max(1,len(base['radius'])//301))
    states={key:base[key][take] for key in ['radius','density','temperature','opacity']}
    for family,u in micro.items():
        states[family+'_opacity']=np.array([u(r,t).kappa for r,t in zip(states['density'],states['temperature'])])
    np.savez_compressed(out/'same_state_microphysics.npz',**states)
    # The abundance conversion is checked thermodynamically, not just compared
    # by its endpoint. Constants use FreeEOS's own isotope mass source.
    u=micro['mass_number_diagnostic'];h=1e-4;defects=[]
    for rho,T in zip(states['density'],states['temperature']):
        q=u(rho,T)
        for axis in [0,1]:
            qp=u(rho*np.exp(h) if axis==0 else rho,T*np.exp(h) if axis==1 else T)
            qm=u(rho*np.exp(-h) if axis==0 else rho,T*np.exp(-h) if axis==1 else T)
            du=((qp.gas_u+3*qp.pr/(rho*np.exp(h) if axis==0 else rho))-(qm.gas_u+3*qm.pr/(rho*np.exp(-h) if axis==0 else rho)))/(2*h)
            rhs=du-(q.pg+q.pr)/rho if axis==0 else du
            defects.append(abs(T*(qp.s-qm.s)/(2*h)-rhs)/(abs(du)+(q.pg+q.pr)/rho))
    # Coverage and plot precision are kept separate; no uncovered sample passes.
    data=list(csv.DictReader(open('benchmarks/nakauchi2017_figure2.csv')))
    coverage=[]
    for rho,label in [(9.65e-9,'paper_9.65e-09'),(9.71e-9,'paper_gn93'),(9.85e-9,'paper_9.85e-09'),(10.5e-9,'paper_1.05e-08')]:
        if label not in profiles:continue
        f=profiles[label]
        subset=[d for d in data if float(d['sonic_density'])==rho]
        outside=[d for d in subset if not f['radius'][0]/RSUN<=float(d['radius_rsun'])<=f['radius'][-1]/RSUN]
        coverage.append(dict(sonic_density=rho,total_samples=len(subset),uncovered_samples=len(outside),radial_range_rsun=[float(f['radius'][0]/RSUN),float(f['radius'][-1]/RSUN)]))
    io.json_write(out/'coverage.json',coverage)
    summary=dict(outcome='OUTER FAIL',v1_validated=False,stage_b_run=False,
        original_benchmark_reproduced=False,discrepancy_explained=False,
        abundance_first_law_max_defect=float(max(defects)),
        cases=len(rows),completed=sum(r['status']=='completed_component' for r in rows),
        elapsed_seconds=time.perf_counter()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
        source_snapshot_verified=io.finish_provenance(out))
    io.json_write(out/'campaign.json',summary);print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
