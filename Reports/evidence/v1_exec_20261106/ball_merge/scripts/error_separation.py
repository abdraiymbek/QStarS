"""Goal 2: separate and bound the two error sources at the step-1255 merge event, on common physical coordinates
(mass below the surface; ln T linear in it between cell centres, the convention of the controller's mapped metric).

A  temporal: merge suppressed (merge threshold 0.45 in memory; no merge fires). Controller change vs dt, the change at
   the front probe, and a step-doubling estimate of the backward-Euler local error (one dt step vs two dt/2 steps).
B  structural: registered merge (0.5) vs merge suppressed at the same dt: J(dt) = lnT_merged - lnT_nomerge.
C  spatial: the outer 14 cells (160..173) refined x f (f = 1, 2, 4, 8) in a one-off remesh, step dt = 1 yr
   (>> thermal time of the outer cells, << evolution time: the result is the quasi-static structure on that mesh).
   Then from each refined state S_f a further step with and without merging the pair of cells that straddles the
   old 168/169 face (at f = 1 this is the registered event) gives the merge jump at resolution f, both as the
   registered controller metric and on common coordinates.
All in memory; remesh_base and Problem.solve monkeypatched in this process only. Writes data/error_separation.json.
Environment: BM_TOL (Newton tolerance for C; default the registered 1e-12).
"""
from bm_common import *
from qstars.mesh import split_merge
state,config,eos,opacity=load()
S=sc.settings(config)
N0=len(state.dm)
PROBE_CELLS=list(range(N0-10,N0-2))           # 164..171
probe=mbs(state)[PROBE_CELLS]                 # mass below the surface, g
P167=PROBE_CELLS.index(167)
lnT0=lnT_at(state,probe)
region=np.linspace(0.,mbs(state)[N0-20],400)   # outer 0.33 Msun, common grid for max |jump|
out=dict(probe_cells=PROBE_CELLS,probe_mass_below_surface_msun=(probe/MSUN).tolist(),
         delta_target=S['delta_target'],reject_threshold=S['delta_target']*S['reject_ratio'])

def front_change(new,ref=None):
    """ln T change at the probes and the max over the outer region, relative to ref (default the accepted state)."""
    a=lnT_at(new,probe);b=lnT0 if ref is None else lnT_at(ref,probe)
    ra=lnT_at(new,region);rb=lnT_at(state if ref is None else ref,region)
    k=int(np.argmax(abs(ra-rb)))
    return (a-b),float((ra-rb)[k]),float(region[k]/MSUN)

# ------------------------------------------------------------------ A temporal
nomerge=with_merge(config,.45)
A=[]
for dt in [1e-3,1e-2,.1,1.,3.,10.,30.,100.,FIRST_REJECTION_DT,2*FIRST_REJECTION_DT]:
    new,d,m=attempt(state,dt,eos,opacity,nomerge)
    assert d['merges']==0 and d['splits']==0
    dp,mx,where=front_change(new)
    A.append(dict(dt_years=dt,delta_actual=m['delta_actual'],component=max((v,kk) for kk,v in m.items() if kk!='delta_actual')[1],
                  profile_dlnT=m['profile_dlnT'],dlnT_probe167=float(dp[P167]),max_outer_dlnT=mx,at_msun=where,
                  newton_iterations=d['newton_iterations']))
    print('A',A[-1],flush=True)
# step doubling: one step dt vs two dt/2 steps (both merge suppressed); BE local error of the full step ~ 2 x difference
for row in A:
    dt=row['dt_years']
    if dt<1:continue
    one,_,_=attempt(state,dt,eos,opacity,nomerge)
    half,dh,_=attempt(state,dt/2,eos,opacity,nomerge)
    two,dh2,_=attempt(half,dt/2,eos,opacity,nomerge)
    assert dh['merges']==0 and dh2['merges']==0
    dp,mx,where=front_change(one,two)
    row.update(stepdoubling_diff_probe167=float(dp[P167]),stepdoubling_diff_max_outer=mx,
               be_local_error_probe167=2*float(dp[P167]),be_local_error_max_outer=2*mx)
    print('A2',dt,row['be_local_error_probe167'],row['be_local_error_max_outer'],flush=True)
big=[r for r in A if r['dt_years']>=10]
c_ctrl=np.polyfit([r['dt_years'] for r in big],[r['profile_dlnT'] for r in big],1)
c_be=[r['be_local_error_max_outer']/r['dt_years']**2 for r in A if 'be_local_error_max_outer' in r]
out['A_temporal']=dict(rows=A,controller_linear_fit_dt_ge_10yr=dict(slope_per_year=c_ctrl[0],intercept=c_ctrl[1]),
                       be_local_error_over_dt2=c_be)

# ------------------------------------------------------------------ B structural
B=[]
for dt in [1e-8,1e-6,1e-5,1e-4,1e-3,1e-2,.1,1.,10.,100.,FIRST_REJECTION_DT]:
    nm,dn,mn=attempt(state,dt,eos,opacity,nomerge)
    mg,dg,mm=attempt(state,dt,eos,opacity,config)
    assert dg['merges']==1 and dn['merges']==0
    jp,jmax,where=front_change(mg,nm)
    B.append(dict(dt_years=dt,delta_actual_merged=mm['delta_actual'],delta_actual_nomerge=mn['delta_actual'],
                  jump_probe167=float(jp[P167]),jump_probes=jp.tolist(),jump_max_outer=jmax,at_msun=where))
    print('B',{k:v for k,v in B[-1].items() if k!='jump_probes'},flush=True)
out['B_structural']=B

# ------------------------------------------------------------------ C spatial (refinement) and merge jump vs resolution
TOL=float(os.environ.get('BM_TOL','1e-12'))
orig_remesh=qs.remesh_base;orig_solve=qs.Problem.solve
mode=dict(spec=None);stash={}
def custom(old,cut,eos_,zones,split=2.,merge=.5,max_passes=40,protect=0):
    spec=mode['spec']
    st,info=orig_remesh(old,cut,eos_,zones,split,spec['registered_merge'],max_passes,protect)
    info=dict(info)
    n=len(st.dm)
    if spec.get('refine'):
        x=np.ones(n);x[n-14:]={1:1.,2:2.5,4:4.5,8:9.}[spec['refine']]
        if spec['refine']>1:
            st,i2=split_merge(st,eos_,st.dm/x,2.,0.);info['splits']+=i2['splits']
    if spec.get('merge_pair'):
        n=len(st.dm);x=np.ones(n);a,b=spec['merge_pair'];x[n+a]=x[n+b]=.3
        st,i3=split_merge(st,eos_,st.dm/x,1e9,.5);info['merges']+=i3['merges']
        assert i3['merges']==1,i3
    return st,info
orig_effect=qs._remesh_effect
def solve(self):
    x,info=orig_solve(self)
    if not stash.get('in_effect'):stash['last']=(self,x,info)   # the accepted pass, not the counterfactual solve
    return x,info
def effect(*a,**k):
    stash['in_effect']=True
    try:return orig_effect(*a,**k)
    finally:stash['in_effect']=False
def floor(p,x,trials=8):
    """Residual max|f| at x and its spread under relative perturbations of x at 4 ulp (roundoff floor test)."""
    rng=np.random.default_rng(1);base=float(np.max(abs(p.evaluate(x))))
    vals=[float(np.max(abs(p.evaluate(x*(1+4e-16*rng.standard_normal(x.size)))))) for _ in range(trials)]
    return dict(at_x=base,perturbed_min=min(vals),perturbed_max=max(vals))
qs.remesh_base=custom;qs.Problem.solve=solve;qs._remesh_effect=effect
def tolcfg(c):
    c=copy.deepcopy(c);c['quasistatic']['tolerance']=TOL;return c
def run_spec(st,dt,spec,cfg):
    mode['spec']=spec;stash.clear()
    try:
        new,d,m=attempt(st,dt,eos,opacity,cfg)
        p,x,info=stash['last']
        return new,dict(ok=True,zones=d['zones'],it=d['newton_iterations'],residual=info['residual'],last_step=info['last_step'],floor=floor(p,x),
                        roundoff_exit=info['roundoff_exit'],delta_actual=m['delta_actual'],merges=d['merges'],splits=d['splits'])
    except qs.StepRejected as e:
        p,x,info=stash['last'];n=p.n
        f=p.evaluate(x);k=int(np.argmax(abs(f)))
        kind=('momentum',k) if k<n else ('energy',k-n) if k<2*n else ('closure',k-2*n) if k<3*n-1 else ('boundary',k-3*n+1)
        T=p.post.T*np.exp(x[n+1:2*n+1]);ghost=p.post.copy();ghost.T=T
        return ghost,dict(ok=False,zones=n,residual=info['residual'],last_step=info['last_step'],it=info['iterations'],floor=floor(p,x),
                          worst_row=kind,message=str(e)[:160])
C=[];refined={}
for f in (1,2,4,8):
    Sf,info=run_spec(state,1.,dict(registered_merge=0.,refine=f),tolcfg(config))
    dp=lnT_at(Sf,probe)-lnT0
    row=dict(factor=f,build=info,dlnT_probes_vs_174=dp.tolist(),dlnT_probe167_vs_174=float(dp[P167]))
    refined[f]=(Sf,info['ok'])
    if info['ok']:
        for dt in (1.,FIRST_REJECTION_DT):
            a,ia=run_spec(Sf,dt,dict(registered_merge=0.),tolcfg(config))
            b,ib=run_spec(Sf,dt,dict(registered_merge=0.,merge_pair=(-5*f-1,-5*f)),tolcfg(config))
            jp,jmax,where=front_change(b,a) if ia['ok'] and ib['ok'] else (np.full(len(probe),np.nan),np.nan,np.nan)
            row[f'dt{dt:.4g}']=dict(nomerge=ia,merge=ib,jump_probe167=float(jp[P167]),jump_max_outer=jmax,at_msun=where)
    C.append(row);print('C',json.dumps(row,default=str)[:900],flush=True)
qs.remesh_base=orig_remesh;qs.Problem.solve=orig_solve;qs._remesh_effect=orig_effect
# Richardson on the probe-167 sequence (refinement of the outer cells by 1, 2, 4, 8)
v=[r['dlnT_probe167_vs_174'] for r in C]
d1,d2,d3=v[1]-v[0],v[2]-v[1],v[3]-v[2]
p12=math.log(d1/d2,2) if d1*d2>0 else None;p23=math.log(d2/d3,2) if d2*d3>0 else None
pe=p23 if p23 else 2.
ext=v[3]+d3/(2**pe-1)
out['C_spatial']=dict(newton_tolerance=TOL,rows=C,observed_order=[p12,p23],extrapolated_dlnT_probe167_vs_174=ext,
                      spatial_error_174_at_probe167=-ext,spatial_error_merged_at_probe167=-ext+B[[r['dt_years'] for r in B].index(1.)]['jump_probe167'])
print('Richardson',out['C_spatial']['observed_order'],ext)
jdump(out,DATA/('error_separation.json' if TOL==1e-12 else f'error_separation_tol{TOL:g}.json'))
