#!/usr/bin/env python3
"""Read-only native high-J radiation recurrence capture, once per iteration."""
import argparse,difflib,hashlib,json
from pathlib import Path


def apply(root,output):
    output.mkdir(parents=True,exist_ok=True);records={};patch=[]
    name='plane/mom_jrel_v9.f';p=root/name;old=p.read_text();s=old
    s=s.replace('\tINTEGER QSTARS_LU,QSTARS_I','\tLOGICAL, SAVE :: QSTARS_BIG_WRITTEN=.FALSE.\n\tINTEGER QSTARS_LU,QSTARS_I',1)
    s=s.replace('\tIF(INIT)THEN','\tIF(INIT)QSTARS_BIG_WRITTEN=.FALSE.\n\tIF(INIT)THEN',1)
    anchor='\t  QSTARS_RAW=XM\n'
    body='''! QSTARS_GROWTH_TRACE: unmodified recurrence and raw linear solution.
\t  IF(.NOT. QSTARS_BIG_WRITTEN .AND. MAXVAL(ABS(XM/GAM_REL/R**2)) .GT. 1.0E20_LDP)THEN
\t    OPEN(NEWUNIT=QSTARS_LU,FILE='QSTARS_MOM_GROWTH',STATUS='UNKNOWN',POSITION='APPEND')
\t    WRITE(QSTARS_LU,*)'FRAME',ND,FREQ,dLOG_NU,MAXVAL(ABS(XM/GAM_REL/R**2))
\t    WRITE(QSTARS_LU,*)'R V chi eta scattering K/J H/J dlnGJ P_J a b c rhs rawJ local prevH prevJ minusEPS plusEPS'
\t    DO QSTARS_I=1,ND
\t      I=QSTARS_I
\t      IF(I .GT. 1 .AND. I .LT. ND)THEN
\t        WRITE(QSTARS_LU,'(19ES25.16)')R(I),V(I),CHI(I),ETA(I),ESEC(I),K_ON_J(I),H_ON_J(I),
\t1          dlnGRSQJdlnR(I),P_J(I),QSTARS_OPERATOR(I,:),QSTARS_RAW(I),
\t1          DTAUONQ(I)*SOURCE(I)*R(I)*R(I),VB(I)*GAM_RSQHNU_PREV(I-1)+VC(I)*GAM_RSQHNU_PREV(I),
\t1          PSIPREV(I)*GAM_RSQJNU_PREV(I),
\t1          -EPS_PREV(I-1)*(GAM_RSQJNU_PREV(I-1)+GAM_RSQJNU_PREV(I)),
\t1          EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1))
\t      ELSE
\t        WRITE(QSTARS_LU,'(19ES25.16)')R(I),V(I),CHI(I),ETA(I),ESEC(I),K_ON_J(I),H_ON_J(I),
\t1          dlnGRSQJdlnR(I),P_J(I),QSTARS_OPERATOR(I,:),QSTARS_RAW(I),0.0_LDP,0.0_LDP,0.0_LDP,0.0_LDP,0.0_LDP
\t      END IF
\t    END DO
\t    CLOSE(QSTARS_LU)
\t    QSTARS_BIG_WRITTEN=.TRUE.
\t  END IF
'''
    s=s.replace(anchor,anchor+body,1);records[name]=(old,s)
    name='new_main/mod_subs/comp_j_blank.f';p=root/name;old=p.read_text();s=old
    s=s.replace('\tINTEGER I\n','\tINTEGER QSTARS_LU\n\tINTEGER I\n',1)
    anchor="\tIF(ABS(RJ(1)) .GT. 1.0E+30_LDP)THEN\n"
    body='''! QSTARS_GROWTH_TRACE: native fatal state, unchanged stopping threshold.
\t  OPEN(NEWUNIT=QSTARS_LU,FILE='QSTARS_J_FAILURE',STATUS='UNKNOWN',POSITION='APPEND')
\t  WRITE(QSTARS_LU,*)'FRAME',ND,FL,FREQ_INDX,SECTION
\t  WRITE(QSTARS_LU,*)'R T V chi eta scattering J K'
\t  DO I=1,ND
\t    WRITE(QSTARS_LU,'(8ES25.16)')R(I),T(I),V(I),CHI_CLUMP(I),ETA_CLUMP(I),CHI_SCAT_CLUMP(I),RJ(I),K_MOM(I)
\t  END DO
\t  CLOSE(QSTARS_LU)
'''
    s=s.replace(anchor,anchor+body,1);records[name]=(old,s)
    manifest={}
    for name,(old,s) in records.items():
        if s==old or 'QSTARS_GROWTH_TRACE' in old:raise ValueError('Pinned anchor absent or already traced')
        (root/name).write_text(s);manifest[name]={k:hashlib.sha256(t.encode()).hexdigest() for k,t in [('original',old),('traced',s)]}
        patch.extend(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    (output/'growth_trace.patch').write_text(''.join(patch));(output/'growth_trace_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);print(json.dumps(apply(**vars(p.parse_args())),indent=2))
