#!/usr/bin/env python3
"""Experimental conservative common-face K/N reconstruction and matched VAR.

Use gamma*r^2 moments at the two nodes with the same midpoint weights.
Do not collapse their distinct coefficients into a single N_mid/J_sum ratio.
No atomic data or physical terms change. This is not production certification.
"""
import argparse, difflib, hashlib, json
from pathlib import Path


def apply(root, output):
    changes = {}
    name = 'plane/mom_jrel_v9.f'; old = (root/name).read_text(); s = old
    if 'QSTARS_NODE_MOMENT' in old: raise ValueError('Already patched')
    if 'QSTARS_RAW' not in old or 'QSTARS_OPERATOR' not in old: raise ValueError('Apply read-only moment tracing first')
    s = s.replace('ALLOCATABLE :: EPS(:),EPS_PREV(:)',
                  'ALLOCATABLE :: EPS(:),EPS_PREV(:),EPS_R(:),EPS_PREV_R(:)', 1)
    s = s.replace('ALLOCATE (EPS_PREV(ND),STAT=IOS)',
                  'ALLOCATE (EPS_PREV(ND),STAT=IOS)\n\t  IF(IOS .EQ. 0)ALLOCATE (EPS_R(ND),EPS_PREV_R(ND),STAT=IOS)', 1)
    s = s.replace('EPS_PREV(:)=0.0_LDP',
                  'EPS_PREV(:)=0.0_LDP\n\t  EPS_R(:)=0.0_LDP;EPS_PREV_R(:)=0.0_LDP', 1)
    start = s.index('\t    EPS(I)=DELTAH(I)'); end = s.index('\n\t  END DO', start)
    s = s[:start]+'''! QSTARS_NODE_MOMENT: one shared face, distinct node contributions.
\t    IF(N_TYPE .NE. 'N_ON_J')THEN
\t      WRITE(6,*)'QSTARS_NODE_MOMENT requires N_ON_J'
\t      STOP
\t    END IF
\t    EPS(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J(I)+
\t1       CON_dKdNUH(I)*K_ON_J(I))/(P_H(I)+W(I))
\t    EPS_R(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J(I+1)+
\t1       CON_dKdNUH(I)*K_ON_J(I+1))/(P_H(I)+W(I))
\t    EPS_PREV(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J_PREV(I)+
\t1       CON_dKdNUH(I)*K_ON_J_PREV(I))/(P_H(I)+W(I))
\t    EPS_PREV_R(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J_PREV(I+1)+
\t1       CON_dKdNUH(I)*K_ON_J_PREV(I+1))/(P_H(I)+W(I))'''+s[end:]
    s = s.replace('TC(I)=-HU(I)+EPS(I)', 'TC(I)=-HU(I)+EPS_R(I)', 1)
    s = s.replace('HU(I-1)-EPS(I-1)+EPS(I)', 'HU(I-1)-EPS_R(I-1)+EPS(I)', 1)
    # Replace every flux contribution, including the read-only trace terms.
    s = s.replace('EPS_PREV(I-1)*(GAM_RSQJNU_PREV(I-1)+GAM_RSQJNU_PREV(I))',
                  '(EPS_PREV(I-1)*GAM_RSQJNU_PREV(I-1)+EPS_PREV_R(I-1)*GAM_RSQJNU_PREV(I))')
    s = s.replace('EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1))',
                  '(EPS_PREV(I)*GAM_RSQJNU_PREV(I)+EPS_PREV_R(I)*GAM_RSQJNU_PREV(I+1))')
    s = s.replace('EPS(I)*(XM(I)+XM(I+1))', '(EPS(I)*XM(I)+EPS_R(I)*XM(I+1))')
    # Capture and fail before native J replacement. No radiation clipping.
    marker='! Check that no negative mean intensities have been computed.\n'
    s = s.replace(marker, marker+'''\tIF(MINVAL(XM) .LT. 0.0_LDP)THEN
\t  WRITE(6,*)'QSTARS_NODE_MOMENT negative J',FREQ,MINLOC(XM),MINVAL(XM)
\t  OPEN(NEWUNIT=QSTARS_LU,FILE='QSTARS_NEGATIVE_J',STATUS='UNKNOWN',POSITION='APPEND')
\t  WRITE(QSTARS_LU,*)'FRAME',ND,FREQ,dLOG_NU,ICOUNT
\t  WRITE(QSTARS_LU,*)'a b c rhs rawJ'
\t  DO QSTARS_I=1,ND
\t    WRITE(QSTARS_LU,'(5ES25.16)')QSTARS_OPERATOR(QSTARS_I,:),QSTARS_RAW(QSTARS_I)
\t  END DO
\t  CLOSE(QSTARS_LU)
\t  STOP
\tEND IF
''', 1)
    changes[name] = (old,s)
    name='plane/var_jrel_v5.f'; old=(root/name).read_text();s=old
    start=s.index('\t    EPS_A(I)=DELTAH(I)');end=s.index('\n\t  END DO',start)
    s=s[:start]+'''! QSTARS_NODE_MOMENT: identical shared face and separate node weights.
\t    EPS_A(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J(I)+
\t1       CON_dKdNUH(I)*K_ON_J(I))/(P_H(I)+W(I))*GAM_RSQ(I)
\t    EPS_B(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J(I+1)+
\t1       CON_dKdNUH(I)*K_ON_J(I+1))/(P_H(I)+W(I))*GAM_RSQ(I+1)
\t    EPS_PREV_A(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J_PREV(I)+
\t1       CON_dKdNUH(I)*K_ON_J_PREV(I))/(P_H(I)+W(I))*GAM_RSQ(I)
\t    EPS_PREV_B(I)=0.5_LDP*DELTAH(I)*(CON_dNdNUH(I)*N_ON_J_PREV(I+1)+
\t1       CON_dKdNUH(I)*K_ON_J_PREV(I+1))/(P_H(I)+W(I))*GAM_RSQ(I+1)'''+s[end:]
    changes[name]=(old,s);output.mkdir(parents=True,exist_ok=True);patch=[];manifest={}
    for name,(old,new) in changes.items():
        if old==new: raise ValueError('Pinned anchor absent: '+name)
        (root/name).write_text(new)
        manifest[name]={'original':hashlib.sha256(old.encode()).hexdigest(),'patched':hashlib.sha256(new.encode()).hexdigest()}
        patch.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    (output/'node_moment.patch').write_text(''.join(patch))
    (output/'node_moment_manifest.json').write_text(json.dumps(dict(sources=manifest,scope='Experimental shared midpoint flux reconstruction; fixed-factor VAR propagation updated. Full relaxed formal-factor response and energy/convergence certification still required.'),indent=2)+'\n')
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path)
    print(json.dumps(apply(**vars(p.parse_args())),indent=2))
