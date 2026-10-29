#!/usr/bin/env python3
"""Gate 6 domain sweep: MESA EOS and the ABI-7 STARS opacity on Ball-like states.

Registered plan benchmarks/ball2011_track_plan.json, execution step 1; the
acceptance of spec Stage 1. No envelope is solved, relaxed or evolved.

State families:
  base         printed (rho0, T0) of Table 1 and thesis Table 3.1;
  polytrope    n = 3 loaded polytropes (the reference.ball_loaded equations,
               small-xi0 branch) with T from the MESA EOS at (rho, p): the
               registered seeds (M_BH/M* = 5e-4, R_seed 1000-12000 Rsun) and
               printed rows 0-4 of Table 1 (M_BH/M*, R*), from r0 out to the
               point where the ideal-gas-plus-radiation T equals Teff;
  photosphere  thesis eq 2.16 at every printed row and seed: the roots rho_ph of
               P_g,EOS(rho, Teff) = (2 g / 3 kappa)(1 - kappa L / (4 pi c G M*)).

Acceptance (spec Stage 1): every state inside the EOS query rectangle, MESA
neutral-ideal fraction exactly 0 and first-law defect <= 1%. The minimum
opacity log R, states outside the registered opacity rectangle (where the
replication-mode STARS clamp would act), the other always-enforced MESA screens
and photosphere existence (risk R2) are reported and do not decide acceptance.
"""
import os
for _name in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_name]='1'
import argparse,csv,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from qstars import io,reference
from qstars.constants import G,C,PI,MSUN,RSUN,LSUN,SIGMA
from qstars.eos import Composition,IdealGasRadiation
from qstars.mesa_microphysics import MesaMicrophysics

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/'benchmarks/ball2011_track_plan.json'
REGISTRATION=ROOT/'benchmarks/ball2011_track_plan.registration.json'
TABLES={'paper':ROOT/'benchmarks/ball2011_table1.csv','thesis':ROOT/'benchmarks/ball2012_thesis_table3_1.csv'}
EOS_RECTANGLE=dict(rho=[1e-14,1e3],T=[3e3,1e8])   # MesaMicrophysics.raw query rectangle
OPACITY_RECTANGLE=dict(logR=[-7.,1.],T=[3e3,1e7])  # plan model.opacity.domain
DEFECT_LIMIT=.01                                   # spec Stage 1 acceptance
SCREEN_LIMIT=.01                                   # plan: dilute-gas and conduction screens
COLUMNS=['frac_ideal','first_law_defect','dilute_gas_defect','conduction_fraction','gamma1','kappa_stars','kappa_default']

def sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def finite(x):return float(x) if x is not None and np.isfinite(x) else None

def verify_registration():
    """Refuse to run unless the plan and both targets match their registered sha256."""
    record=json.loads(REGISTRATION.read_text());plan=sha256(PLAN)
    if plan!=record['plan_sha256']:
        raise SystemExit(f'Refusing to run: {PLAN} sha256 {plan} differs from the registered {record["plan_sha256"]}')
    for name,digest in record['target_files_sha256'].items():
        if sha256(ROOT/name)!=digest:raise SystemExit(f'Refusing to run: {name} differs from its registered sha256')
    return dict(plan=str(PLAN.relative_to(ROOT)),plan_sha256=plan,registration_sha256=sha256(REGISTRATION),
                registered_utc=record['registered_utc'],target_files_sha256=record['target_files_sha256'])

def printed_rows(plan):
    """Printed rows in cgs; M* is reduced by the booked mass defect eps/(1-eps) dM_BH."""
    model=plan['model'];eps=model['accretion_law']['epsilon']
    total,initial=model['total_mass_msun'],model['initial_bh_mass_msun']
    rows=[]
    for source,path in TABLES.items():
        with open(path,newline='') as f:
            for i,r in enumerate(csv.DictReader(f)):
                mbh=float(r['bh_mass_msun'])
                L=float(r['luminosity_1e8_lsun'])*1e8*LSUN;R=float(r['r_outer_1e4_rsun'])*1e4*RSUN
                rows.append(dict(label=f'{source}_row{i}',source=source,role='printed row',bh_mass_msun=mbh,
                    total_mass_msun=total-eps/(1-eps)*(mbh-initial),luminosity=L,radius=R,
                    teff=(L/(4*PI*R*R*SIGMA))**.25,rho0=float(r['rho_inner_cgs']),T0=float(r['T_inner_1e5_K'])*1e5))
    return rows

def loaded_profile(x0,n=3,rtol=1e-10):
    """reference.ball_loaded (paper eqs 24-26, 30, 32, 33) with dense output."""
    def rhs(x,y):return [-y[1]/x**2,x**2*max(y[0],0)**n]
    def surface(x,y):return y[0]
    surface.terminal=True;surface.direction=-1
    sol=solve_ivp(rhs,[x0,1e4],[1.,x0/(2*n)+2*x0**3/3],events=surface,
                  rtol=rtol,atol=rtol*.01,max_step=.2,dense_output=True)
    if len(sol.t_events[0])!=1:raise RuntimeError('Loaded polytrope has no surface')
    return sol

def polytrope(q,total_mass,radius,teff,ideal,limit,samples):
    """Small-xi0 loaded polytrope with M_BH/M* = q scaled to (M*, R*): rho, p from r0 to T_ideal = Teff."""
    if q>=limit['bh_fraction']:
        return None,dict(status='no_loaded_polytrope',bh_fraction=q,
                         reason=f'M_BH/M* = {q:.5g} exceeds the n = 3 loaded-polytrope maximum {limit["bh_fraction"]:.5g}')
    x0=brentq(lambda x:reference.ball_loaded(x)['bh_fraction']-q,1e-6,limit['xi0'],xtol=1e-15,rtol=1e-13)
    sol=loaded_profile(x0);xi1=sol.t[-1];ys=sol.y[1,-1]
    if abs(x0/6/ys/reference.ball_loaded(x0)['bh_fraction']-1)>1e-10:
        raise RuntimeError('Dense loaded polytrope differs from reference.ball_loaded')
    # theta(xi) at rho = rho_c theta^3, p = K rho^(4/3); alpha^2 = K rho_c^(-2/3) / (pi G).
    alpha=radius/xi1;rho_c=total_mass/(4*PI*alpha**3*ys);K=PI*G*alpha**2*rho_c**(2/3)
    xi=np.unique(np.r_[np.geomspace(x0,xi1,20001)[:-1],np.linspace(x0,xi1,20001)[:-1]])
    theta=sol.sol(xi)[0];keep=theta>0;xi,theta=xi[keep],theta[keep]
    T_ideal=ideal.temperature(rho_c*theta**3,K*rho_c**(4/3)*theta**4,'p')
    inside=np.flatnonzero(T_ideal>=teff)
    if inside[-1]!=len(inside)-1:raise RuntimeError('Ideal temperature is not monotone in the polytrope')
    lnp=np.log(K*rho_c**(4/3)*theta[inside]**4)
    # Equal steps in ln p from the base to the truncation point.
    xs=np.interp(-np.linspace(lnp[0],lnp[-1],samples),-lnp,xi[inside])
    theta=sol.sol(xs)[0];rho=rho_c*theta**3;p=K*rho_c**(4/3)*theta**4
    return (rho,p),dict(status='built',bh_fraction=q,xi0=x0,xi1=xi1,r0_over_R=x0/xi1,alpha_cm=alpha,
                        rho_base=rho[0],p_base=p[0],rho_outer=rho[-1],p_outer=p[-1],truncation_teff=teff)

def temperatures(mesa,rho,p):
    """MESA T at (rho, p); a per-state fallback isolates the failing states."""
    try:return mesa.temperature(rho,p,'p'),[None]*len(rho)
    except (ValueError,RuntimeError):
        T=np.full(len(rho),np.nan);errors=[]
        for i in range(len(rho)):
            try:T[i]=mesa.temperature(rho[i:i+1],p[i:i+1],'p')[0];errors.append(None)
            except (ValueError,RuntimeError) as e:errors.append(f'temperature inversion: {e}')
        return T,errors

def audit(mesa,rho,T):
    """Per-state EOS and opacity audit; a state outside a domain is recorded, never raised."""
    n=len(rho);out={k:np.full(n,np.nan) for k in COLUMNS}
    out['clamp_flag']=np.full(n,-1);out['error']=[None]*n
    out['inside_eos_rectangle']=((rho>=EOS_RECTANGLE['rho'][0])&(rho<=EOS_RECTANGLE['rho'][1])&
                                 (T>=EOS_RECTANGLE['T'][0])&(T<=EOS_RECTANGLE['T'][1]))
    for i in range(n):
        if not out['inside_eos_rectangle'][i]:
            out['error'][i]='outside the EOS query rectangle';continue
        r,t=rho[i:i+1],T[i:i+1]
        try:
            x=mesa.raw(r,t);out['frac_ideal'][i]=x[18,0];out['kappa_default'][i]=x[19,0]
            q=mesa.evaluate(r,t)
            out['gamma1'][i]=q['gamma1'][0];out['conduction_fraction'][i]=q['conductive_fraction'][0]
            out['first_law_defect'][i]=mesa.consistency(r,t)[0]
            out['dilute_gas_defect'][i]=mesa.gas_pressure_limit(r,t,q)['screened_relative_error'][0]
        except (ValueError,RuntimeError) as e:out['error'][i]=f'EOS: {e}'
    for i in np.flatnonzero((T>=OPACITY_RECTANGLE['T'][0])&(T<=OPACITY_RECTANGLE['T'][1])):
        kappa,flag=mesa.opacity_stars(rho[i],T[i],clamp=True,flags=True)
        out['kappa_stars'][i]=kappa;out['clamp_flag'][i]=flag
    return out

def photosphere(mesa,case,grid=301):
    """Every root in ln rho of thesis eq 2.16 at T = Teff inside the EOS query rectangle."""
    T,L,M,R=case['teff'],case['luminosity'],case['total_mass_msun']*MSUN,case['radius']
    g=G*M/R**2;ledd=4*PI*C*G*M
    def parts(rho):
        rho=np.atleast_1d(rho);temp=np.full(len(rho),T)
        try:pg=np.exp(mesa.raw(rho,temp)[0])
        except ValueError:
            pg=np.full(len(rho),np.nan)
            for i in range(len(rho)):
                try:pg[i]=np.exp(mesa.raw(rho[i:i+1],temp[i:i+1])[0,0])
                except ValueError:pass
        kappa,flag=mesa.opacity_stars(rho,temp,clamp=True,flags=True)
        gamma=kappa*L/ledd
        return pg,kappa,flag,gamma,2*g/(3*kappa)*(1-gamma)
    def residual(lnrho):
        pg,_,_,gamma,target=parts(np.exp(lnrho))
        return np.where(gamma<1,np.log(pg)-np.log(np.where(gamma<1,target,1.)),np.inf)
    lnrho=np.linspace(np.log(EOS_RECTANGLE['rho'][0]),np.log(1e-4),grid);values=residual(lnrho)
    roots=[]
    for i in np.flatnonzero(np.isfinite(values[:-1])&np.isfinite(values[1:])&(np.sign(values[:-1])!=np.sign(values[1:]))):
        x=brentq(lambda s:residual(np.array([s]))[0],lnrho[i],lnrho[i+1],xtol=1e-13,rtol=1e-14)
        pg,kappa,flag,gamma,target=(v[0] for v in parts(np.array([np.exp(x)])))
        roots.append(dict(rho=float(np.exp(x)),logR=float(x/np.log(10)-3*np.log10(T)+18),kappa=float(kappa),
                          clamp_flag=int(flag),gamma_ph=float(gamma),gas_pressure=float(pg),
                          relative_residual=float(pg/target-1)))
    gamma=parts(np.exp(lnrho))[3]
    return dict(label=case['label'],role=case['role'],teff=T,gravity=g,luminosity=L,total_mass_msun=case['total_mass_msun'],
                radius=R,roots=roots,exists=bool(roots),gamma_on_grid=[finite(np.nanmin(gamma)),finite(np.nanmax(gamma))],
                note=None if roots else 'no root of eq 2.16 inside the EOS query rectangle (Gamma_ph >= 1 or rho_ph < 1e-14)')

def summarize(s,mask):
    def ext(key,fn):
        v=s[key][mask];v=v[np.isfinite(v)]
        return float(fn(v)) if len(v) else None
    ratio=s['kappa_stars'][mask]/s['kappa_default'][mask];ratio=ratio[np.isfinite(ratio)]
    off_T=mask&np.isfinite(s['T'])&((s['T']<OPACITY_RECTANGLE['T'][0])|(s['T']>OPACITY_RECTANGLE['T'][1]))
    return dict(states=int(mask.sum()),errors=int(sum(e is not None for e,m in zip(s['error'],mask) if m)),
        rho_range=[ext('rho',np.min),ext('rho',np.max)],T_range=[ext('T',np.min),ext('T',np.max)],
        logR_range=[ext('logR',np.min),ext('logR',np.max)],
        max_frac_ideal=ext('frac_ideal',np.max),max_first_law_defect=ext('first_law_defect',np.max),
        max_dilute_gas_defect=ext('dilute_gas_defect',np.max),max_conduction_fraction=ext('conduction_fraction',np.max),
        opacity_T_outside_domain=int(off_T.sum()),opacity_clamped_low=int((mask&(s['clamp_flag']==1)).sum()),
        opacity_clamped_high=int((mask&(s['clamp_flag']==2)).sum()),
        kappa_stars_over_default=[float(ratio.min()),float(ratio.max())] if len(ratio) else None)

def main():
    p=argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--library',default=os.environ.get('QSTARS_MESA_LIBRARY',str(ROOT/'backends/mesa/build_ball_abi7/libqstars_mesa.dylib')))
    p.add_argument('--mesa-dir',default=os.environ.get('MESA_DIR'),help='default: mesa_dir of the build manifest')
    p.add_argument('--output',default=None,help='parent directory; default benchmarks/results/<UTC date>-ball-track')
    p.add_argument('--samples',type=int,default=400,help='states per loaded polytrope')
    a=p.parse_args()
    registration=verify_registration()
    plan=json.loads(PLAN.read_text());model=plan['model'];initial=plan['initial_model']
    library=Path(a.library).resolve()
    mesa_dir=a.mesa_dir or json.loads((library.parent/'build.json').read_text())['mesa_dir']
    composition=Composition(**model['composition']);ideal=IdealGasRadiation(composition)
    mesa=MesaMicrophysics(composition,library,mesa_dir)
    if mesa.abi_version<7:raise SystemExit('The registered STARS opacity switch requires a MESA bridge ABI 7 build')
    parent=Path(a.output) if a.output else ROOT/'benchmarks/results'/(datetime.now(timezone.utc).strftime('%Y-%m-%d')+'-ball-track')
    start=time.perf_counter();out=io.create_directory(parent,'domain_sweep')
    io.provenance(out,dict(kind='ball2011_domain_sweep',plan=registration,library=str(library),mesa_dir=str(mesa_dir),
        samples_per_polytrope=a.samples,eos_rectangle=EOS_RECTANGLE,opacity_rectangle=OPACITY_RECTANGLE,
        first_law_defect_limit=DEFECT_LIMIT,screen_limit=SCREEN_LIMIT,
        threads={k:os.environ[k] for k in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS')}))
    (out/'driver.py').write_bytes(Path(__file__).read_bytes())
    io.json_write(out/'plan_registration.json',registration)
    io.json_write(out/'microphysics.json',dict(mesa=mesa.manifest,abi_version=mesa.abi_version,library_sha256=sha256(library)))
    print(out,flush=True)

    rows=printed_rows(plan)
    seed_L=next(r['luminosity'] for r in rows if r['label']=='paper_row0')
    seeds=[dict(label=f'seed_{int(R)}rsun',role='registered seed' if R in initial['seed_radii_rsun'] else 'registered fallback seed',
                bh_mass_msun=model['initial_bh_mass_msun'],total_mass_msun=model['total_mass_msun'],
                luminosity=seed_L,radius=R*RSUN,teff=(seed_L/(4*PI*(R*RSUN)**2*SIGMA))**.25)
           for R in sorted(initial['seed_radii_rsun']+initial['fallback_seed_radii_rsun'])]
    table={k:[] for k in ('family','case','rho','T','p','error')}
    def add(family,case,rho,T,p,errors):
        for k,v in zip(table,([family]*len(rho),[case]*len(rho),rho,T,p,errors)):table[k].extend(v)

    for r in rows:add('base',r['label'],[r['rho0']],[r['T0']],[np.nan],[None])
    limit=reference.ball_limit();polytropes=[]
    for c in seeds+[r for r in rows if r['source']=='paper']:
        built,info=polytrope(c['bh_mass_msun']/c['total_mass_msun'],c['total_mass_msun']*MSUN,c['radius'],c['teff'],ideal,limit,a.samples)
        info.update(label=c['label'],role=c['role'],radius_rsun=c['radius']/RSUN,total_mass_msun=c['total_mass_msun'])
        if built is not None:
            rho,pres=built;T,errors=temperatures(mesa,rho,pres)
            info.update(T_base=finite(T[0]),T_outer=finite(T[-1]),temperature_failures=sum(e is not None for e in errors))
            if 'rho0' in c:info.update(printed_rho0=c['rho0'],printed_T0=c['T0'])
            add('polytrope',c['label'],rho,T,pres,errors)
        polytropes.append(info);print('polytrope',c['label'],info['status'],flush=True)
    atmospheres=[]
    for c in rows+seeds:
        ph=photosphere(mesa,c);atmospheres.append(ph)
        for root in ph['roots']:add('photosphere',c['label'],[root['rho']],[c['teff']],[np.nan],[None])
        print('photosphere',c['label'],[(f"{x['rho']:.3e}",f"logR={x['logR']:.3f}",f"Gamma={x['gamma_ph']:.4f}") for x in ph['roots']],flush=True)

    rho=np.array(table['rho'],dtype=float);T=np.array(table['T'],dtype=float);n=len(rho)
    family=np.array(table['family']);case=np.array(table['case']);ok=np.isfinite(T)
    s=dict(rho=rho,T=T,logR=np.where(ok,np.log10(rho)-3*np.log10(np.where(ok,T,1.))+18,np.nan),
           clamp_flag=np.full(n,-1),inside_eos_rectangle=np.zeros(n,dtype=bool),error=list(table['error']))
    for k in COLUMNS:s[k]=np.full(n,np.nan)
    checked=audit(mesa,rho[ok],T[ok])
    for k in COLUMNS+['clamp_flag','inside_eos_rectangle']:s[k][ok]=checked[k]
    for j,i in enumerate(np.flatnonzero(ok)):s['error'][i]=s['error'][i] or checked['error'][j]

    with open(out/'states.csv','w',newline='') as f:
        w=csv.writer(f);keys=['rho','T','logR']+COLUMNS
        w.writerow(['family','case']+keys+['clamp_flag','p_total','error'])
        for i in range(n):
            w.writerow([family[i],case[i]]+[repr(float(s[k][i])) for k in keys]+[int(s['clamp_flag'][i]),repr(float(table['p'][i])),s['error'][i] or ''])

    fi=s['frac_ideal'];de=s['first_law_defect'];errors=sum(e is not None for e in s['error'])
    accepted=bool(errors==0 and s['inside_eos_rectangle'].all() and np.all(fi==0) and np.all(de<=DEFECT_LIMIT))
    families={name:summarize(s,family==name) for name in ('base','polytrope','photosphere')}
    i_min=int(np.nanargmin(s['logR']))
    screens={k:finite(np.nanmax(s[k])) for k in ('dilute_gas_defect','conduction_fraction')}
    summary=dict(
        acceptance=dict(criterion='every state inside the EOS query rectangle (rho 1e-14..1e3 g/cm3, T 3000..1e8 K), '
                        'MESA neutral-ideal fraction exactly 0 and first-law defect <= 1% (spec Stage 1)',
                        accepted=accepted,states=n,errors=errors,max_frac_ideal=finite(np.nanmax(fi)),max_first_law_defect=finite(np.nanmax(de))),
        minimum_logR=dict(value=float(s['logR'][i_min]),family=str(family[i_min]),case=str(case[i_min]),
                          rho=float(rho[i_min]),T=float(T[i_min]),design_estimate=-6.2),
        opacity_domain=dict(rectangle=OPACITY_RECTANGLE,states_needing_clamp=int((s['clamp_flag']>0).sum()),
                            states_with_T_outside=int((ok&((T<OPACITY_RECTANGLE['T'][0])|(T>OPACITY_RECTANGLE['T'][1]))).sum()),
                            cases_with_T_outside=sorted({str(c) for c,t in zip(case,T) if np.isfinite(t) and not OPACITY_RECTANGLE['T'][0]<=t<=OPACITY_RECTANGLE['T'][1]}),
                            note='reported only: the replication-mode STARS clamp acts in log R at fixed T; off-table T is always an error'),
        always_enforced_screens=dict(limit=SCREEN_LIMIT,maxima=screens,exceeded=[k for k,v in screens.items() if v is not None and v>SCREEN_LIMIT]),
        photospheres_missing=[ph['label'] for ph in atmospheres if not ph['exists']],
        base_gamma1_mesa_over_ideal={r['label']:finite(s['gamma1'][i]/ideal.evaluate(np.array([r['rho0']]),np.array([r['T0']]))['gamma1'][0])
                                     for i,r in zip(np.flatnonzero(family=='base'),rows)},
        families=families,cases={name:summarize(s,case==name) for name in dict.fromkeys(case)},
        polytropes=polytropes,photospheres=atmospheres,runtime_seconds=time.perf_counter()-start)
    io.json_write(out/'domain_sweep.json',summary)
    print(json.dumps({k:summary[k] for k in ('acceptance','minimum_logR','opacity_domain','always_enforced_screens','photospheres_missing')},indent=1))
    return 0 if accepted else 1

if __name__=='__main__':raise SystemExit(main())
