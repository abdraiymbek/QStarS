from common import *
from locate_lib import internals
state,config,eos,opacity=load()
def flags(st):
    t=qs.face_transport(st,eos,opacity,config);return np.asarray(t['convective_face']),st
acc,_=flags(state)
out={}
for rm in (False,True):
    I=internals(state,1.,eos,opacity,config,remesh=rm)
    new,_=qs._assemble(I['P'],I['x'],I['info'],I['cut_info'],dict(I['remesh'],regrid=I['regrid'],effect=None),I['passes'],I['pred']['swept'])
    f,_=flags(new);out[rm]=f
# faces are interior faces 1..N-1; compare the pre-event solution with the accepted state face by face,
# and the event solution with the pre-event one on the faces both meshes share (all but the removed face 169)
print('accepted: convective faces',int(acc.sum()),'of',len(acc))
print('pre-event solution vs accepted: faces that changed regime',np.flatnonzero(out[False]!=acc).tolist())
pre=out[False];ev=out[True]
shared=np.r_[pre[:168],pre[169:]]  # pre-event interior face index 168 (between cells 168 and 169) disappears in the merge
print('event solution vs pre-event (shared faces): faces that changed regime',np.flatnonzero(ev!=shared).tolist(),'| convective faces',int(ev.sum()),'of',len(ev))
