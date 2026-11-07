"""Repeat the last rejected attempt of a front-refined DRYRUN from its final checkpoint, with the same in-memory
front targets, with and without merges, and report the event and where the change sits.
usage: probe_event_stop.py RUN C W FMAX DT   (appends to data/event_stop_probe.json)"""
from bm_common import *
name,CF,W,FMAX,DT=sys.argv[1],float(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])
st,config=qs.load_checkpoint(BIG/name/'checkpoint.npz');_,_,eos,opacity=load()
def front_factor(old):
    kap=np.log(opacity(old.rho,old.T));gf=abs(np.diff(kap))
    g=np.maximum(np.r_[0.,gf],np.r_[gf,0.])
    pad=np.r_[np.zeros(W),g,np.zeros(W)]
    g=np.max(np.lib.stride_tricks.sliding_window_view(pad,2*W+1),axis=1)
    return np.clip(CF/np.maximum(g,1e-300),1./FMAX,1.)
def targets_front(old,cut,eos_,zones):
    lnt=np.log(qs.q_targets(old,eos_,zones))+np.log(front_factor(old))
    co=np.cumsum(old.dm)-.5*old.dm;c=np.cumsum(cut.dm)-.5*cut.dm;To,T=old.dm.sum(),cut.dm.sum()
    inner=np.interp(c,co,lnt);outer=np.interp((T-c)[::-1],(To-co)[::-1],lnt[::-1])[::-1]
    return np.exp(np.where(c<=.5*T,inner,outer))
qs.targets=targets_front
n=len(st.dm);kap=opacity(st.rho,st.T);fac=front_factor(st)
rows=[]
for label,cfg in (('as run',config),('merges off',with_merge(config,0.)),('splits off',(lambda c:(c['quasistatic'].__setitem__('split',1e9),c)[1])(copy.deepcopy(config)))):
    for dt in (DT,DT/8):
        try:
            new,d,m=attempt(st,dt,eos,opacity,cfg)
            k=int(np.argmax(abs(lnT_at(new,mbs(st))-np.log(st.T))))
            r=dict(case=label,dt=dt,ok=True,delta_actual=m['delta_actual'],merges=d['merges'],splits=d['splits'],zones=d['zones'],
                   it=d['newton_iterations'],max_change_cell_from_surface=n-k,mass_below_surface_msun=float(mbs(st)[k]/MSUN))
        except qs.StepRejected as e:r=dict(case=label,dt=dt,ok=False,message=str(e)[:120])
        rows.append(r);print(r,flush=True)
outer=[dict(from_surface=n-k,dlnkappa=float(np.log(kap[k+1]/kap[k])) if k+1<n else None,factor=float(fac[k]),
            x=None) for k in range(n-16,n)]
p=DATA/'event_stop_probe.json';res=json.loads(p.read_text()) if p.exists() else {}
res[name]=dict(state=dict(steps=st.steps,bh_mass_msun=st.bh_mass/MSUN,zones=n),attempts=rows,outer=outer)
jdump(res,p)
print([(o['from_surface'],round(o['dlnkappa'] or 0,3),round(o['factor'],3)) for o in outer])
