#!/usr/bin/env python3
"""Prepare a supplied CMFGEN reference without executing its shell scripts."""
import argparse, hashlib, json, re, shutil
from pathlib import Path

def keyword(path, key, value):
    text=path.read_text(); pattern=r'(?m)^.*?\['+re.escape(key)+r'\].*$'
    line=f'{value:<22} [{key}]'
    text,n=re.subn(pattern,line,text)
    if n>1:raise ValueError(f'Duplicate keyword {key}')
    path.write_text(text if n else text+'\n'+line+'\n')

def prepare(source, destination, atomic, iterations=1):
    source=source.resolve();atomic=atomic.resolve();destination.mkdir(parents=True,exist_ok=False)
    for name in ['VADAT','IN_ITS','MODEL_SPEC','RVSIG_COL','HYDRO_DEFAULTS','ROSSELAND_LTE_TAB','RDINR','IT_SPECIFIER','ADJUST_R_DEFAULTS','RVTJ','GREY_SCL_FACIN','GREY_SCL_FAC_IN']:
        if (source/name).is_file():shutil.copy2(source/name,destination/name)
    for p in source.glob('*_IN'):shutil.copy2(p,destination/p.name)
    # Official cpmod semantics: use converged departure coefficients as restart guesses.
    for p in source.glob('*OUT'):
        if re.fullmatch(r'[A-Za-z0-9]+OUT',p.name):shutil.copy2(p,destination/(p.name[:-3]+'_IN'))
    if (source/'GAMMAS').exists():shutil.copy2(source/'GAMMAS',destination/'GAMMAS_IN')
    links={}
    for batch in ['batch.sh','batch_ins.sh']:
        if not (source/batch).exists():continue
        for line in (source/batch).read_text().splitlines():
            m=re.match(r'\s*ln\s+-sf\s+\$(?:ATOMIC|atomic)/([^\s]+)\s+([\w.-]+)\s*(?:#.*)?$',line)
            if m:
                target=(atomic/m[1]).resolve();target.relative_to(atomic)
                links[m[2]]=target
    active=set()
    for line in (source/'MODEL_SPEC').read_text().splitlines():
        m=re.match(r'\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\[([^]]+)_ISF\]',line)
        if m and int(m[3])>0:active.add(m[4])
    def required(name):
        if name in {'TWO_PHOT_DATA','XRAY_PHOT_FITS','HYD_L_DATA','GBF_N_DATA'}:return True
        return any(name.startswith(ion+'_') or name.startswith('PHOT'+ion+'_') or name=='DIE'+ion for ion in active)
    missing=[str(p) for n,p in links.items() if required(n) and not p.is_file()]
    if missing:raise FileNotFoundError('Missing active atomic data: '+str(missing))
    skipped=[dict(link=n,path=str(p),reason='inactive model ion') for n,p in links.items() if not p.is_file()]
    (destination/'inactive_missing_links.json').write_text(json.dumps(skipped,indent=2)+'\n')
    manifest=[]
    for name,path in links.items():
        if not path.is_file():continue
        (destination/name).symlink_to(path)
        manifest.append(dict(link=name,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    # Supplied batch files use a temperature/departure file as T_IN.
    if (destination/'He2_IN').is_file():(destination/'T_IN').symlink_to('He2_IN')
    elif (destination/'HI_IN').is_file():(destination/'T_IN').symlink_to('HI_IN')
    keyword(destination/'IN_ITS','NUM_ITS',str(iterations))
    (destination/'atomic_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (destination/'reference_origin.json').write_text(json.dumps(dict(source=str(source),iterations=iterations,
      mdot_role='prescribed_reference_only',production_wind=False),indent=2)+'\n')
    return destination

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('source',type=Path);a.add_argument('destination',type=Path)
    a.add_argument('--atomic',type=Path,required=True);a.add_argument('--iterations',type=int,default=1)
    opts=a.parse_args();print(prepare(**vars(opts)))
