#!/usr/bin/env python3
"""Opt-in read-only CMFGEN frequency diagnostic; never changes physical arrays."""
import argparse,hashlib,difflib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
f=a.root/'new_main/cmfgen_sub.f';old=f.read_text()
assert 'QSTAR_EXPORT' not in old
oldhash=hashlib.sha256(f.read_bytes()).hexdigest()
if oldhash != '092f8526661b1f9a5eaeb7a875f07f55623f427ad07c383f2659f2fb67143374':
    raise ValueError('Source does not match the inspected 18-Jun-2025 release')
new=old.replace('\tREAL(KIND=LDP) LUM_SCL_FAC','\tINTEGER QSTAR_LU\n\tLOGICAL QSTAR_EXPORT\n\tREAL(KIND=LDP) LUM_SCL_FAC',1)
marker='! Evaluate (CMF) observd X-ray luminosities.'
hook='''! QStarS diagnostic ABI 1: comoving spectral data; cgs, outer to inner.
! An empty QSTARS_EXPORT file enables this read-only snapshot.
\tIF(LST_ITERATION)THEN
\t  IF(ML .EQ. 1)THEN
\t    INQUIRE(FILE='QSTARS_EXPORT',EXIST=QSTAR_EXPORT)
\t    IF(QSTAR_EXPORT)THEN
\t      OPEN(NEWUNIT=QSTAR_LU,FILE='QSTARS_SPECTRAL.bin',
\t1      STATUS='REPLACE',ACCESS='STREAM',FORM='UNFORMATTED')
\t      WRITE(QSTAR_LU)'QSCMF001',ND,NCF
\t      WRITE(QSTAR_LU)R*1.0E10_LDP,V*1.0E5_LDP,
\t1       T*1.0E4_LDP,DENSITY,CLUMP_FAC
\t    END IF
\t  END IF
\t  IF(QSTAR_EXPORT)THEN
\t    WRITE(QSTAR_LU)NU(ML)*1.0E15_LDP,FQW(ML),
\t1      CHI*1.0E-10_LDP,CHI_SCAT*1.0E-10_LDP,
\t1      ETA*1.0E-10_LDP,RJ,H_MOM,K_MOM
\t    IF(ML .EQ. NCF)CLOSE(QSTAR_LU)
\t  END IF
\tEND IF
!
'''
assert old.count(marker)==1
new=new.replace(marker,hook+marker,1);f.write_text(new)
out=Path(__file__).parent
(out/'export.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='new_main/cmfgen_sub.f',tofile='new_main/cmfgen_sub.f')))
(out/'export_patch_manifest.json').write_text(json.dumps(dict(original_sha256=oldhash,patched_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),abi=1,physical_arrays_modified=False),indent=2)+'\n')
