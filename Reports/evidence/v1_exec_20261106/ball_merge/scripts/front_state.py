"""Front resolution at the final checkpoints of the goal-3 DRYRUNs, and for a stopped front-refined run the Newton end
state of its last rejected attempt repeated with the same in-memory front refinement (residual floor test).
usage: front_state.py [RUN C W FMAX DT]   (writes data/front_state.json, appending per call)"""
from bm_common import *
from qstars.constants import G,C as CLIGHT,PI
def outer(state,opacity,cells=30):
    kap=opacity(state.rho,state.T);n=len(state.dm)
    gf=np.diff(np.log(kap));L=state.face_luminosity;m=state.enclosed
    gam=kap*0.5*(L[:-1]+L[1:])/(4*PI*CLIGHT*G*0.5*(m[:-1]+m[1:]))
    sl=slice(n-cells,n)
    return dict(zones=n,bh_mass_msun=state.bh_mass/MSUN,steps=state.steps,max_abs_dlnkappa_face=float(np.max(abs(gf[n-cells:]))),
                faces_dlnkappa_over_0_3=int(np.sum(abs(gf[n-cells:])>.3)),max_eddington=float(np.max(gam[sl])),
                cells_gamma_over_1=int(np.sum(gam[sl]>1)),max_dlnT_cell=float(np.max(abs(np.diff(np.log(state.T[sl]))))),
                min_abs_dlnP_face=float(np.min(abs(np.diff(np.log(state.eos_p))))) if hasattr(state,'eos_p') else None)
res_path=DATA/'front_state.json'
res=json.loads(res_path.read_text()) if res_path.exists() else {}
_,_,eos,opacity=load()
res['registered_step1255']=outer(qs.load_checkpoint(RUN/'checkpoint.npz')[0],opacity)
for d in sorted(BIG.glob('E*')):
    if (d/'checkpoint.npz').exists() and (d/'DRYRUN_report.json').exists():
        st,_=qs.load_checkpoint(d/'checkpoint.npz')
        q=eos.evaluate(st.rho,st.T);o=outer(st,opacity);n=len(st.dm)
        o['min_abs_dlnP_outer_faces']=float(np.min(abs(np.diff(np.log(q['p'][n-30:])))))
        o['faces_closure_scale_floored_outer30']=int(np.sum(.25*abs(np.diff(np.log(q['p'][n-30:])))<1e-3))
        res[d.name]=dict(res.get(d.name,{}),final=o)
if len(sys.argv)>5:
    name,CF,W,FMAX,DT=sys.argv[1],float(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])
    st,config=qs.load_checkpoint(BIG/name/'checkpoint.npz')
    def front_factor(old):
        kap=np.log(opacity(old.rho,old.T));gf=abs(np.diff(kap))
        g=np.maximum(np.r_[0.,gf],np.r_[gf,0.])
        pad=np.r_[np.zeros(W),g,np.zeros(W)]
        g=np.max(np.lib.stride_tricks.sliding_window_view(pad,2*W+1),axis=1)
        return np.clip(CF/np.maximum(g,1e-300),1./FMAX,1.)
    def targets_front(old,cut,eos_,zones):
        lnt=np.log(qs.q_targets(old,eos_,zones))+np.log(front_factor(old))
        centre_old=np.cumsum(old.dm)-.5*old.dm;centre=np.cumsum(cut.dm)-.5*cut.dm
        total_old,total=old.dm.sum(),cut.dm.sum()
        inner=np.interp(centre,centre_old,lnt)
        outer_=np.interp((total-centre)[::-1],(total_old-centre_old)[::-1],lnt[::-1])[::-1]
        return np.exp(np.where(centre<=.5*total,inner,outer_))
    qs.targets=targets_front
    stash={};orig=qs.Problem.solve
    def solve(self):
        x,info=orig(self);stash.setdefault('main',(self,x,info));return x,info
    qs.Problem.solve=solve
    rows=[]
    for dt in (DT,DT/2,DT/4,DT/8):
        stash.clear()
        try:attempt(st,dt,eos,opacity,config);ok=True
        except qs.StepRejected:ok=False
        p,x,info=stash['main'];n=p.n;f=p.evaluate(x);k=int(np.argmax(abs(f)))
        rng=np.random.default_rng(1)
        pert=[float(np.max(abs(p.evaluate(x*(1+4e-16*rng.standard_normal(x.size)))))) for _ in range(8)]
        row=dict(dt_years=dt,ok=ok,iterations=info['iterations'],residual=info['residual'],last_step=info['last_step'],
                 perturbed_residual_min=min(pert),perturbed_residual_max=max(pert),
                 worst_row=('closure face %d, %d cells from surface, scale %.2e'%(k-2*n+1,n-(k-2*n+1),p.contrast_scale[k-2*n])) if 2*n<=k<3*n-1 else int(k))
        rows.append(row);print(row,flush=True)
    res[name]=dict(res.get(name,{}),last_attempt_probe=rows)
jdump(res,res_path)
print(json.dumps({k:v.get('final',v) for k,v in res.items()},indent=0,default=str)[:3000])
