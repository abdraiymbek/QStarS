#!/usr/bin/env python3
"""Read-only failure operator snapshots; original native failure behavior retained."""
import argparse,hashlib,json
from pathlib import Path

def apply(root):
    records={}
    for name in ['plane/mom_jrel_v9.f','plane/var_jrel_v5.f']:
        p=root/name;original=p.read_text();s=original
        if 'QSTARS_FAILURE_TRACE' in s:raise ValueError('Trace already installed')
        if name.endswith('mom_jrel_v9.f'):
            # Only the MOM routine has ICOUNT, not its module/allocation routine.
            marker='\tINTEGER ICOUNT\n'
            declarations='\tINTEGER QSTARS_LU,QSTARS_I\n\tREAL(KIND=LDP), ALLOCATABLE :: QSTARS_OPERATOR(:,:),QSTARS_RAW(:)\n'
            s=s.replace(marker,declarations+marker,1)
            marker='\tDO WHILE(.NOT. ACCURATE)\n'
            s=s.replace(marker,'\tALLOCATE(QSTARS_OPERATOR(ND,4),QSTARS_RAW(ND))\n'+marker,1)
            marker='\t  CALL THOMAS(TA,TB,TC,XM,ND,1)\n'
            save='\t  QSTARS_OPERATOR(:,1)=TA; QSTARS_OPERATOR(:,2)=TB\n\t  QSTARS_OPERATOR(:,3)=TC; QSTARS_OPERATOR(:,4)=XM\n'
            s=s.replace(marker,save+marker+'\t  QSTARS_RAW=XM\n',1)
            marker="\t      WRITE(LUER,*)'Error in MOM_J_REL_V9: excessive iteration count.'\n"
            body="""! QSTARS_FAILURE_TRACE: record only, no change to operator/solution.
\t      OPEN(NEWUNIT=QSTARS_LU,FILE='QSTARS_MOM_FAILURE',STATUS='UNKNOWN',POSITION='APPEND')
\t      WRITE(QSTARS_LU,*)'FRAME',ND,FREQ,dLOG_NU,MAX_ER
\t      WRITE(QSTARS_LU,*)'R V chi eta scattering K/J H/J dlnGJ P_J a b c rhs rawJ repairedJ'
\t      DO QSTARS_I=1,ND
\t        WRITE(QSTARS_LU,'(16ES25.16)')R(QSTARS_I),V(QSTARS_I),CHI(QSTARS_I),ETA(QSTARS_I),
\t1          ESEC(QSTARS_I),K_ON_J(QSTARS_I),H_ON_J(QSTARS_I),dlnGRSQJdlnR(QSTARS_I),P_J(QSTARS_I),
\t1          QSTARS_OPERATOR(QSTARS_I,:),QSTARS_RAW(QSTARS_I),XM(QSTARS_I)
\t      END DO
\t      CLOSE(QSTARS_LU)
"""
            s=s.replace(marker,body+marker,1)
        else:
            marker='\tINTEGER ICNT\n'
            s=s.replace(marker,'\tINTEGER QSTARS_LU,QSTARS_I\n'+marker,1)
            marker="\t    WRITE(LUER,*)'Error in VAR_JREL_V2: RJ blowing up'\n"
            body="""! QSTARS_FAILURE_TRACE: native fatal state and unmodified operator.
\t    OPEN(NEWUNIT=QSTARS_LU,FILE='QSTARS_VAR_FAILURE',STATUS='UNKNOWN',POSITION='APPEND')
\t    WRITE(QSTARS_LU,*)'FRAME',ND,FREQ,dLOG_NU
\t    WRITE(QSTARS_LU,*)'R V chi eta scattering K/J H/J dlnGJ a b c rhs rawJ'
\t    DO QSTARS_I=1,ND
\t      WRITE(QSTARS_LU,'(13ES25.16)')R(QSTARS_I),V(QSTARS_I),CHI(QSTARS_I),ETA(QSTARS_I),
\t1        ESEC(QSTARS_I),K_ON_J(QSTARS_I),H_ON_J(QSTARS_I),dlnGRSQJdlnR(QSTARS_I),
\t1        TA_SAV(QSTARS_I),TB_SAV(QSTARS_I),TC_SAV(QSTARS_I),XM_SAV(QSTARS_I),XM(QSTARS_I)
\t    END DO
\t    CLOSE(QSTARS_LU)
"""
            s=s.replace(marker,body+marker,1)
        if s==original or 'QSTARS_FAILURE_TRACE' not in s:raise ValueError('Pinned native anchors absent')
        records[name]={'original':hashlib.sha256(original.encode()).hexdigest(),'traced':hashlib.sha256(s.encode()).hexdigest()}
        p.write_text(s)
    (root/'qstars_moment_trace.json').write_text(json.dumps({'scope':'Read-only failure snapshots. Original native operators, clipping, stopping and coefficients retained.','sources':records},indent=2)+'\n')
    return records

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('root',type=Path);print(json.dumps(apply(**vars(a.parse_args())),indent=2))
