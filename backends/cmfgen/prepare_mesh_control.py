#!/usr/bin/env python3
"""Fresh-grid transfer controls from completed H/He populations (not a wind solve)."""
import argparse,json,re,shutil
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from prepare_hhe import make_hhe
from prepare_reference import keyword
from qstars.cmfgen_backend import sha256

def prepare(source,template,destination,atomic,factor=1):
    if factor not in (1,2,4,8):raise ValueError('Control supports dyadic factors 1, 2, 4, 8')
    make_hhe(template,destination,atomic,1)
    for name in ['VADAT','MODEL_SPEC','RVSIG_COL']:
        shutil.copy2(source/name,destination/name)
    for name in ['HI','HeI','He2']:
        shutil.copy2(source/(name+'OUT'),destination/(name+'_IN'))
    text=(destination/'VADAT').read_text()
    (destination/'VADAT').write_text(''.join(s for s in text.splitlines(True) if '[STOP_ON_WARN]' not in s))
    for k in ['FIX_T','COMP_F','USE_J_REL','INCL_REL','INCL_ADV_TRANS']:keyword(destination/'VADAT',k,'T')
    for k,v in [('DO_LAM_IT','T'),('DO_LAM_AUTO','F'),('DO_T_AUTO','F')]:keyword(destination/'IN_ITS',k,v)
    lines=(source/'RVSIG_COL').read_text().splitlines()
    idx=next(i for i,s in enumerate(lines) if '!Number of depth points' in s);n=int(lines[idx].split()[0])
    rows=[]
    for s in lines[idx+1:]:
        if not s.strip() or s.lstrip().startswith('!'):continue
        rows.append([float(x.replace('D','E')) for x in s.split()[:3]])
        if len(rows)==n:break
    rv=np.array(rows)
    if rv.shape!=(n,3):raise ValueError('Malformed reference radius/velocity grid')
    if factor>1:
        x=np.log(rv[::-1,0]);spline=CubicHermiteSpline(x,np.log(rv[::-1,1]),rv[::-1,2]+1)
        xx=np.r_[np.concatenate([np.linspace(a,b,factor,endpoint=False) for a,b in zip(x[:-1],x[1:])]),x[-1]]
        rv=np.column_stack((np.exp(xx),np.exp(spline(xx)),spline(xx,1)-1))[::-1]
    if not np.all(np.diff(rv[:,1])<0):raise ValueError('Interpolated flow is not monotonic')
    n=len(rv);nc=15
    keyword(destination/'MODEL_SPEC','ND',str(n));keyword(destination/'MODEL_SPEC','NP',str(n+nc))
    (destination/'RVSIG_COL').write_text(f'{n} !Number of depth points\n'+''.join(f'{r:.14e} {v:.14e} {sigma:.14e}\n' for r,v,sigma in rv))
    (destination/'QSTARS_EXPORT').touch()
    (destination/'changes.json').write_text(json.dumps(dict(scope='frozen incoming H/He state transfer; no thermal/spatial equilibrium certification',
       source=str(source),factor=factor,radii=n,core_rays=nc,material_interpolation='native departure-coefficient initialization from completed *_OUT',
       flow_interpolation='cubic Hermite log(v)(log(r)), prescribed dlog(v)/dlog(r)=sigma+1; geometric interval bisection',
       source_hashes={n:sha256(source/n) for n in ['HIOUT','HeIOUT','He2OUT','RVSIG_COL']},
       production=False),indent=2)+'\n')
    return destination
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('template',type=Path);p.add_argument('destination',type=Path)
    p.add_argument('--atomic',type=Path,required=True);p.add_argument('--factor',type=int,default=1);print(prepare(**vars(p.parse_args())))
