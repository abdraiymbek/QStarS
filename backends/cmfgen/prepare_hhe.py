#!/usr/bin/env python3
"""Declared H/He prescribed-flow component derived from supplied Zeta Pup grid."""
import argparse,re,json
from pathlib import Path
from prepare_reference import prepare,keyword

def make_hhe(source,destination,atomic,iterations=40):
    prepare(source,destination,atomic,iterations)
    p=destination/'MODEL_SPEC';lines=[]
    for line in p.read_text().splitlines():
        m=re.search(r'\[([^]]+)_ISF\]',line)
        if m and m[1] not in {'HI','HeI','He2'}:continue
        lines.append(line)
    p.write_text('\n'.join(lines)+'\n')
    p=destination/'VADAT'
    for key in re.findall(r'\[([^]]+/X)\]',p.read_text()):
        keyword(p,key,{'HYD/X':'9.333333333333333','HE/X':'1.0'}.get(key,'0.0'))
    for key,value in {'DO_HYDRO':'F','DO_CL':'F','INC_AD':'F','INC_ADV':'F','FIX_T':'F','FIX_T_AUTO':'F','INC_XRAYS':'F'}.items():keyword(p,key,value)
    keyword(destination/'IN_ITS','DO_LAM_IT','T')
    keyword(destination/'IN_ITS','DO_LAM_AUTO','T')
    (destination/'component_definition.json').write_text(json.dumps(dict(
        purpose='NLTE H/He transfer/population component, not a quasi-star or wind eigenvalue',
        input_flow='supplied Zeta Pup radius/velocity and mdot; hydro updates disabled',
        mass_fractions={'H':.7,'He':.3},mass_convention='CMFGEN integer H=1, He=4',
        microclumping=False,adiabatic_cooling=False,population_advection=False,
        atomic_levels={'HI':30,'HeI':69,'HeII':30},superlevels_equal_full_levels=True,
        mdot_is_prescribed=True,production_certified=False),indent=2)+'\n')
    return destination
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('destination',type=Path)
    p.add_argument('--atomic',type=Path,required=True);p.add_argument('--iterations',type=int,default=40)
    print(make_hhe(**vars(p.parse_args())))
