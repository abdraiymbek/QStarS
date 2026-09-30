#!/usr/bin/env python3
"""Pinned CMFGEN H-minus FF initialization, last-depth and Jacobian repair.

No table, coefficient, physical process or rate is changed. Normal continuum,
electron-balance and cooling paths use the same neutral-H ground population.
The strict stopping change retains every level in the convergence criterion.
"""
import argparse, difflib, hashlib, json
from pathlib import Path


def apply(root, output):
    changes = {}
    name = 'newsubs/do_h0_ff.f'
    p = root/name; old = p.read_text(); s = old
    start = s.index('\tIF(FIRST_TIME)THEN')
    end = s.index('\n\tLOG_NU=LOG(CONT_FREQ)', start)
    load = s[start:end]
    # Keep native parsing/interpolation/units. Close the file after loading.
    load = load.replace('\t  FIRST_TIME=.FALSE.', '\t  CLOSE(LUIN)\n\t  FIRST_TIME=.FALSE.')
    s = s[:start]+'\tCALL INIT_H0_FF_DATA(LUIN)\n'+s[end:]
    anchor = '\tLOG_NU=LOG(CONT_FREQ)'
    pos = s.index(anchor, s.index('SUBROUTINE DO_H0_FF_COOL'))
    s = s[:pos]+'\tCALL INIT_H0_FF_DATA(7)\n'+s[pos:]
    s += '''
! QSTARS_HMINUS_FF_REPAIR: lazy loading for every independent entry point.
\tSUBROUTINE INIT_H0_FF_DATA(LUIN)
\tUSE H0_FF_DATA
\tIMPLICIT NONE
\tINTEGER LUIN,I,J
\tREAL(KIND=LDP) T1,BOLTZMANN_CONSTANT
\tEXTERNAL BOLTZMANN_CONSTANT
\tCHARACTER(LEN=200) STRING
\tLOGICAL, SAVE :: FIRST_TIME=.TRUE.
'''+load+'''
\tEND
!
! Exact derivatives of native bilinear log-table interpolation and Planck
! emissivity. Derivatives are at fixed frequency and independent populations.
! Ion argument is the neutral-H ground population, not H-minus population.
\tSUBROUTINE DO_H0_FF_VAR(VETA,VCHI,ION_DEN,ED,TEMP,NU,EQION,ND,NVAR,LAST_ONLY)
\tUSE H0_FF_DATA
\tIMPLICIT NONE
\tINTEGER ND,NVAR,EQION,I,J,K,KSTART
\tLOGICAL LAST_ONLY
\tREAL(KIND=LDP) VETA(NVAR,ND),VCHI(NVAR,ND),ION_DEN(ND),ED(ND),TEMP(ND),NU
\tREAL(KIND=LDP) W,U,C0,C1,SLOPE,BASE,OP,EM,B,X,E,DB
\tREAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
\tCOMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
\tCALL INIT_H0_FF_DATA(7)
\tJ=1
\tDO WHILE(J .LT. NNU-1)
\t  IF(LOG(NU) .LE. LOG_NU_TAB(J+1))EXIT
\t  J=J+1
\tEND DO
\tW=(LOG(NU)-LOG_NU_TAB(J))/(LOG_NU_TAB(J+1)-LOG_NU_TAB(J))
\tKSTART=1
\tIF(LAST_ONLY)KSTART=ND
\tDO K=KSTART,ND
\t  I=1
\t  DO WHILE(I .LT. NT-1)
\t    IF(LOG(TEMP(K)) .LT. LOG_T_TAB(I+1))EXIT
\t    I=I+1
\t  END DO
\t  C0=(1.0_LDP-W)*CROSS(I,J)+W*CROSS(I,J+1)
\t  C1=(1.0_LDP-W)*CROSS(I+1,J)+W*CROSS(I+1,J+1)
\t  SLOPE=(C1-C0)/(LOG_T_TAB(I+1)-LOG_T_TAB(I))
\t  U=(LOG(TEMP(K))-LOG_T_TAB(I))/(LOG_T_TAB(I+1)-LOG_T_TAB(I))
\t  BASE=TEMP(K)*EXP((1.0_LDP-U)*C0+U*C1)
\t  OP=BASE*ED(K)*ION_DEN(K)
\t  X=HDKT*NU/TEMP(K); E=EXP(-X)
\t  B=TWOHCSQ*NU**3*E/(1.0_LDP-E)
\t  DB=B*X/TEMP(K)/(1.0_LDP-E)
\t  EM=OP*B
\t  VCHI(EQION,K)=VCHI(EQION,K)+BASE*ED(K)
\t  VETA(EQION,K)=VETA(EQION,K)+BASE*ED(K)*B
\t  VCHI(NVAR-1,K)=VCHI(NVAR-1,K)+BASE*ION_DEN(K)
\t  VETA(NVAR-1,K)=VETA(NVAR-1,K)+BASE*ION_DEN(K)*B
\t  VCHI(NVAR,K)=VCHI(NVAR,K)+OP*(1.0_LDP+SLOPE)/TEMP(K)
\t  VETA(NVAR,K)=VETA(NVAR,K)+EM*(1.0_LDP+SLOPE)/TEMP(K)+OP*DB
\tEND DO
\tEND
'''
    changes[name] = (old,s)
    name='new_main/subs/var_op_v10.f';p=root/name;old=p.read_text()
    s=old.replace('\tIF( IONFF .AND. PHOT_ID .EQ. 1)THEN', '''! QSTARS_HMINUS_FF_REPAIR: same ground population as GENOPAETA_V10.
\tIF(Z .EQ. 0.0_LDP)THEN
\t  POP_SUM=DI_S(1,:)
\t  CALL DO_H0_FF_VAR(PETA,PCHI,POP_SUM,ED,T,NU,GS_ION_EQ,ND,NT,LST_DEPTH_ONLY)
\tELSE IF( IONFF .AND. PHOT_ID .EQ. 1)THEN''',1)
    changes[name]=(old,s)
    name='newsubs/genopaeta_v10.f';p=root/name;old=p.read_text()
    s=old.replace('\t    K=ND_LOC\n\t    COR_FAC(K)=DI(1,K)', '\t    K=ND\n\t    COR_FAC(K)=DI(1,K)',1)
    s=s.replace('EMHNUKT(K),NU,I,K)','EMHNUKT(K),NU,I,1)',1)
    changes[name]=(old,s)
    name='newsubs/comp_free_free_v2.f';p=root/name;old=p.read_text()
    s=old.replace('CALL DO_H0_FF(ETA_FF,CHI_FF,POP_SUM(1,ID),ED,T,EMHNUKT,FREQ,I,ND)',
       'GFF_VAL=ATM(ID+1)%XzV(1,:)\n\t    CALL DO_H0_FF(ETA_FF,CHI_FF,GFF_VAL,ED,T,EMHNUKT,FREQ,I,ND)',1)
    anchor='\t    ELSE\n\t      GFF_VAL=GFF_STORE(:,NINT(ATM(ID)%ZXzV))'
    s=s.replace(anchor,'''\t    ELSE IF(ION_ID(ID) .EQ. 'H0')THEN
\t      GFF_VAL=ATM(ID+1)%XzV(1,:)
\t      EQION=ATM(ID+1)%EQXzV
\t      CALL DO_H0_FF_VAR(VETA_FF,VCHI_FF,GFF_VAL,ED,T,FREQ,EQION,ND,NT,.FALSE.)
\t    ELSE IF(NINT(ATM(ID)%ZXzV) .GT. 0)THEN
\t      GFF_VAL=GFF_STORE(:,NINT(ATM(ID)%ZXzV))''',1)
    changes[name]=(old,s)
    name='new_main/subs/solveba_v13.f';p=root/name;old=p.read_text()
    a=old.index('\tIF(DO_LEVEL_CHK)THEN\n\t  MAXCH=');b=old.index('\tEND IF',a)+len('\tEND IF')
    s=old[:a]+'''! QSTARS_STRICT_LEVEL_STOP: do not discard ten largest corrections.
! Keep DO_LEVEL_CHK for native scaling, but MAXCH remains the full maximum.'''+old[b:]
    changes[name]=(old,s)
    output.mkdir(parents=True,exist_ok=True)
    manifest={};patch=[]
    for name,(old,new) in changes.items():
        if old==new:raise ValueError('Pinned anchor absent: '+name)
        if 'QSTARS_HMINUS_FF_REPAIR' in old or 'QSTARS_STRICT_LEVEL_STOP' in old:raise ValueError('Already repaired')
        (root/name).write_text(new)
        manifest[name]={'original':hashlib.sha256(old.encode()).hexdigest(),'patched':hashlib.sha256(new.encode()).hexdigest()}
        patch.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    (output/'hminus_ff.patch').write_text(''.join(patch))
    (output/'hminus_ff_patch_manifest.json').write_text(json.dumps({'source_hashes':manifest,'scope':'native FF table and coefficients unchanged; ground-H consistent opacity/cooling/EHB derivatives; full-level nonlinear stopping'},indent=2)+'\n')
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);print(json.dumps(apply(**vars(p.parse_args())),indent=2))
