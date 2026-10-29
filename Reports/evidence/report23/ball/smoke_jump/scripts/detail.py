from common import *
from locate_lib import internals,table
state,config,eos,opacity=load()
print('=== outer cells: accepted state (step 1255) ===')
table(state,eos,opacity,config,lo=156)
for dt in [1.]:
    for rm in (False,True):
        I=internals(state,dt,eos,opacity,config,remesh=rm)
        new,_=qs._assemble(I['P'],I['x'],I['info'],I['cut_info'],dict(I['remesh'],regrid=I['regrid'],effect=None),I['passes'],I['pred']['swept'])
        print(f'=== outer cells: solution dt={dt} yr, merge={"on" if rm else "off"} ===')
        table(new,eos,opacity,config,lo=156)
