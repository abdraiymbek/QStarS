#!/usr/bin/env python3
"""Exact frozen recurrence attribution and rejected repair controls, never a wind."""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.linalg import solve, solve_banded
from qstars.cmfgen_backend import _numbers, sha256


def read_seed(path):
    lines=path.read_text().splitlines(); head=lines[0].split();n=int(head[1])
    x=np.array([_numbers(s) for s in lines[2:2+n]])
    if x.shape!=(n,41) or not np.isfinite(x).all():raise ValueError('Invalid complete recurrence')
    return head,x


def frozen_operator(x):
    n=len(x);a,b,c,rhs=x[:,30:34].T
    A=np.diag(b)+np.diag(a[1:],-1)+np.diag(c[:-1],1)
    BJ=np.zeros((n,n));BH=np.zeros((n,n-1));F=np.zeros((n-1,n));G=np.zeros((n-1,n))
    BJ[0,0]=x[0,22]
    for i in range(1,n-1):
        BJ[i,i]=x[i,22]-x[i-1,27]+x[i,27]
        BJ[i,i-1]=-x[i-1,27];BJ[i,i+1]=x[i,27]
        BH[i,i-1]=-x[i-1,25];BH[i,i]=x[i,25]
    for i in range(n-1):
        F[i,i]=-x[i,24]-x[i,26];F[i,i+1]=x[i,23]-x[i,26]
        G[i,i:i+2]=x[i,27]
    s=rhs-BJ@x[:,28]-BH@x[:-1,29]
    return A,BJ,BH,F,G,s


def node_candidate(x,A,rhs):
    """Identical two-node midpoint reconstruction to the experimental native patch."""
    n=len(x);beta=.5*(x[:-1,1]+x[1:,1])*3.33564e-6;sig=.5*(x[:-1,2]+x[1:,2])
    cn=(sig+1)/(1-beta**2)-1;ck=beta*(sig+1)/(1-beta**2)
    fac=.5*x[:-1,36]/(1+x[:-1,36])
    EL=fac*(cn*x[:-1,9]+ck*x[:-1,7]);ER=fac*(cn*x[1:,9]+ck*x[1:,7])
    PL=fac*(cn*x[:-1,40]+ck*x[:-1,38]);PR=fac*(cn*x[1:,40]+ck*x[1:,38])
    out=A.copy();source=rhs.copy();p=x[:,28]
    for i in range(1,n-1):
        out[i,i-1]+=x[i-1,26]-EL[i-1];out[i,i+1]+=ER[i]-x[i,26]
        out[i,i]+=x[i-1,26]-ER[i-1]+EL[i]-x[i,26]
        old=-x[i-1,27]*(p[i-1]+p[i])+x[i,27]*(p[i]+p[i+1])
        new=-PL[i-1]*p[i-1]-PR[i-1]*p[i]+PL[i]*p[i]+PR[i]*p[i+1]
        source[i]+=new-old
    return solve(out,source)


def analyze(directory,output,candidate=None):
    output.mkdir(parents=True,exist_ok=False);head,x=read_seed(directory/'QSTARS_MOM_SEED')
    A,BJ,BH,F,G,s=frozen_operator(x);rhs=x[:,33];raw=x[:,34];n=len(x)
    band=np.zeros((3,n));band[0,1:]=x[:-1,32];band[1]=x[:,31];band[2,:-1]=x[1:,30]
    independent=solve_banded((1,1),band,rhs)
    componentwise=max(abs(A@raw-rhs)/(abs(A)@abs(raw)+abs(rhs)+1e-300))
    i=int(np.argmin(raw));j=x[:,28];h=x[:-1,29]
    terms=[x[i,19]*x[i,4]/x[i,15]*x[i,0]**2,
           -x[i-1,25]*h[i-1]+x[i,25]*h[i],x[i,22]*j[i],
           -x[i-1,27]*(j[i-1]+j[i]),x[i,27]*(j[i]+j[i+1])]
    # This counterfactual is attribution only. It is NOT an accepted H overwrite.
    angular_h=.5*(x[:-1,39]*j[:-1]+x[1:,39]*j[1:])
    replay=solve(A,s+BJ@j+BH@angular_h)
    # A shared, outward-upwind beta*x flux for the zeroth equation, including
    # removal of its expanded beta' reaction. Native boundaries remain untouched;
    # this rejected diagnostic is not a complete conservative replacement.
    adv=A.copy();beta=x[:,1]*3.33564e-6
    for k in range(1,n-1):
        adv[k,k]-=x[k,19]*(x[k,14]-1+(x[k,15]-x[k,3]/x[k,6])/x[k,15])
        adv[k,k]+=.5*(beta[k-1]+beta[k]);adv[k,k+1]-=.5*(beta[k]+beta[k+1])
    implicit=solve(adv,rhs);node=node_candidate(x,A,rhs)
    result=dict(frequency_1e15Hz=float(head[2]),row_1based=i+1,
      native_raw_J=float(raw[i]),independent_relative_difference=float(max(abs(raw-independent))/max(abs(independent))),
      componentwise_backward_residual=float(componentwise),seed_row_P_J=float(x[i,14]),
      seed_rhs=float(rhs[i]),rhs_terms=dict(zip(['local','previous_H_divergence','previous_J','minus_previous_EPS','plus_previous_EPS'],map(float,terms))),
      previous_face_H=float(h[i]),neighbour_previous_J=j[i:i+2].tolist(),
      attribution_only_positive_angular_H_replay_min=float(min(replay)),
      rejected_implicit_J_advection_min=float(min(implicit)),rejected_node_reconstruction_min=float(min(node)),
      source_sha256=sha256(directory/'QSTARS_MOM_SEED'),driver_sha256=sha256(__file__),
      scope='Frozen failed checkpoint. Exact seed in signed previous-frequency H/EPS forcing; positive J-advection diagonal. Counterfactual H overwrite is not a physical repair. No nonlinear, thermal, energy or full relaxed response claim.',production_certified=False)
    if candidate:
        p=candidate/'QSTARS_NEGATIVE_J';lines=p.read_text().splitlines();hd=lines[0].split();m=int(hd[1])
        z=np.array([_numbers(r) for r in lines[2:2+m]]);mat=np.diag(z[:,1])+np.diag(z[1:,0],-1)+np.diag(z[:-1,2],1)
        b=np.zeros((3,m));b[0,1:]=z[:-1,2];b[1]=z[:,1];b[2,:-1]=z[1:,0]
        sol=solve_banded((1,1),b,z[:,3]);local=abs(mat@sol-z[:,3])/(abs(mat)@abs(sol)+abs(z[:,3])+1e-300)
        result['native_node_candidate_failure']=dict(frequency_1e15Hz=float(hd[2]),row_1based=int(np.argmin(sol))+1,
          min_raw_J=float(min(z[:,4])),min_independent_J=float(min(sol)),componentwise_backward_residual=float(max(local)),
          seed_row_rhs=float(z[np.argmin(sol),3]),sha256=sha256(p),
          scope='Strict unclipped full-cycle request rejects first negative state before reaching the He I edge. Tiny absolute magnitude is retained, not floored or certified away.')
        np.savez_compressed(output/'node_failure.npz',state=z,matrix=mat,independent=sol)
    np.savez_compressed(output/'seed.npz',state=x,matrix=A,independent=independent)
    (output/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('output',type=Path);p.add_argument('--candidate',type=Path)
    analyze(**vars(p.parse_args()))
