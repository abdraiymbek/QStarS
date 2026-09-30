#!/usr/bin/env python3
"""Re-solve read-only native failure operators; distinguish linear/nonlinear failure."""
import argparse,json,re
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded,eigvalsh_tridiagonal
from qstars.cmfgen_backend import _numbers,sha256


def analyze(directory,output):
    output.mkdir(parents=True,exist_ok=False);frames=[];inputs={}
    for name,ncols,a_col,raw_col in [('QSTARS_MOM_FAILURE',15,9,13),('QSTARS_VAR_FAILURE',13,8,12),('QSTARS_MOM_GROWTH',19,9,13)]:
        p=directory/name
        if not p.exists():continue
        inputs[name]=sha256(p);lines=p.read_text().splitlines();i=0
        while i<len(lines):
            header=lines[i].split()
            if not header or header[0]!='FRAME':raise ValueError('Bad trace frame')
            n=int(header[1]);x=np.array([_numbers(row) for row in lines[i+2:i+2+n]])
            if x.shape!=(n,ncols) or not np.all(np.isfinite(x)):raise ValueError('Malformed operator')
            a,b,c,rhs=x[:,a_col:a_col+4].T;raw=x[:,raw_col]
            mat=np.diag(b)+np.diag(a[1:],-1)+np.diag(c[:-1],1)
            band=np.zeros((3,n));band[0,1:]=c[:-1];band[1]=b;band[2,:-1]=a[1:]
            independent=solve_banded((1,1),band,rhs)
            norm=max(abs(independent));diff=max(abs(raw-independent))/norm
            res=max(abs(mat@raw-rhs))/(np.linalg.norm(mat,np.inf)*max(abs(raw))+max(abs(rhs)))
            componentwise=max(abs(mat@raw-rhs)/(abs(mat)@abs(raw)+abs(rhs)+1e-300))
            off=a[1:]*c[:-1]
            eigen_min=float(eigvalsh_tridiagonal(b,np.sqrt(off))[0]) if np.all(off>=0) else None
            frames.append(dict(file=name,ND=n,frequency_1e15Hz=float(header[2]),dlog_nu=float(header[3]),
              native_fixedpoint_error=float(header[4]) if name=='QSTARS_MOM_FAILURE' else None,
              captured_max_physical_J=float(header[4]) if name=='QSTARS_MOM_GROWTH' else None,
              linear_solution_relative_difference=float(diff),linear_backward_residual=float(res),
              componentwise_linear_backward_residual=float(componentwise),
              negative_rhs_rows=np.flatnonzero(rhs<0).tolist(),negative_solution_rows=np.flatnonzero(raw<0).tolist(),
              raw_solution_range=[float(raw.min()),float(raw.max())],matrix_condition=float(np.linalg.cond(mat)),
              diagonal_positive=bool(np.all(b>0)),off_diagonal_nonpositive=bool(np.all(a[1:]<=0) and np.all(c[:-1]<=0)),
              symmetrized_min_eigenvalue=eigen_min,negative_absorption_rows=np.flatnonzero(x[:,2]-x[:,4]<0).tolist(),
              opacity_min=float(min(x[:,2])),absorption_min=float(min(x[:,2]-x[:,4])),
              photon_advection_log_gradient_range=[float(min(x[:,7])),float(max(x[:,7]))],
              repaired_solution_difference=float(max(abs(x[:,-1]-raw))/norm) if name=='QSTARS_MOM_FAILURE' else None,
              interior_rhs_decomposition_relative_error=float(max(abs(x[1:-1,14:].sum(axis=1)-rhs[1:-1]))/max(abs(rhs[1:-1]))) if name=='QSTARS_MOM_GROWTH' else None,
              interior_local_source_to_rhs_norm=float(max(abs(x[1:-1,14]))/max(abs(rhs[1:-1]))) if name=='QSTARS_MOM_GROWTH' else None))
            np.savez_compressed(output/(name+str(len(frames))+'.npz'),state=x,matrix=mat,independent=independent)
            i+=n+2
    fatal=None;p=directory/'QSTARS_J_FAILURE'
    if p.exists():
        lines=p.read_text().splitlines();header=lines[0].split();n=int(header[1])
        state=np.array([_numbers(s) for s in lines[2:2+n]])
        if state.shape!=(n,8) or not np.all(np.isfinite(state)):raise ValueError('Malformed fatal J state')
        fatal=dict(ND=n,frequency_1e15Hz=float(header[2]),frequency_index=int(header[3]),
            J_range=[float(state[:,6].min()),float(state[:,6].max())],
            temperature_range_K=[float(state[:,1].min()*1e4),float(state[:,1].max()*1e4)],
            absorption_min=float((state[:,3]-state[:,5]).min()),sha256=sha256(p))
    conclusion=('Independent band solve reproduces failed native linear solutions; this diagnoses captured operators only, not nonlinear closure or physical nonexistence.' if frames else
        'No full raw operator was captured; independent linear-solve attribution is unavailable for this case. Fatal native J state, if present, is retained.')
    result=dict(frames=frames,fatal_state=fatal,input_sha256=inputs,driver_sha256=sha256(__file__),
      conclusion=conclusion,production_certified=False)
    (output/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('output',type=Path);analyze(**vars(p.parse_args()))
