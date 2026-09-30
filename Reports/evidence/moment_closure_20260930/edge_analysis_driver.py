#!/usr/bin/env python3
"""Audit projected shared-face moments; do not confuse them with formal rays."""
import argparse,json
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import _numbers,sha256


def analyze(path,output):
    lines=path.read_text().splitlines();rows=[];first=None;t=0
    while t<len(lines):
        head=lines[t].split();n=int(head[1]);x=np.array([_numbers(s) for s in lines[t+2:t+n+2]])
        if x.shape!=(n,41) or not np.isfinite(x).all():raise ValueError('Invalid frame')
        j=x[:,34]
        h=x[:-1,23]*j[1:]-x[:-1,24]*j[:-1]+x[:-1,25]*x[:-1,29]+x[:-1,27]*(x[:-1,28]+x[1:,28])-x[:-1,26]*(j[:-1]+j[1:])
        k=x[:-1,12]*(j[:-1]+j[1:]);nn=x[:-1,11]*(j[:-1]+j[1:])+x[:-1,10]*h
        # For I>=0 on [-1,1], both moment matrices of (1+mu)I and
        # (1-mu)I are PSD. Their Schur bounds give the minimum feasible j.
        with np.errstate(divide='ignore',invalid='ignore'):
            required=np.maximum((h+k)**2/(k+nn)-h,(h-k)**2/(k-nn)+h)
        maximum=np.maximum(j[:-1],j[1:]);bad=np.where((required>maximum*(1+1e-8))|(k<abs(nn)))[0]
        if first is None and len(bad):
            i=int(bad[0]);first=dict(frequency_1e15Hz=float(head[2]),face_1based=i+1,
              neighbour_gJ=j[i:i+2].tolist(),gH=float(h[i]),gK=float(k[i]),gN=float(nn[i]),
              minimum_feasible_gJ=float(required[i]),required_to_maximum_node_ratio=float(required[i]/maximum[i]))
        rows.append(dict(frequency_1e15Hz=float(head[2]),min_raw_gJ=float(min(j)),failed_faces=(bad+1).tolist()))
        t+=n+2
    result=dict(frames=len(rows),first_in_selected_window=first,sequence=rows,
       input_sha256=sha256(path),driver_sha256=sha256(__file__),
       validity_boundary='Projected staggered face K/N/H cannot be moments of nonnegative angular intensity for any monotone gJ interpolation at the reported face. This does not audit the underlying formal-ray I against its own face interpolation, nor locate the earliest failure outside the selected window.',
       production_certified=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sequence'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('path',type=Path);p.add_argument('output',type=Path);analyze(**vars(p.parse_args()))
