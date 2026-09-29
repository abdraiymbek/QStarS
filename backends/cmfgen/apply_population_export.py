#!/usr/bin/env python3
"""Repair read-only RVTJ ion indexing when H- precedes HI in the ion list."""
import argparse,difflib,json
from pathlib import Path
from qstars.cmfgen_backend import sha256

p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
f=a.root/'new_main/cmfgen_sub.f';old=f.read_text();oldhash=sha256(f)
expected=json.loads((Path(__file__).parent/'export_patch_manifest.json').read_text())['patched_sha256']
if oldhash!=expected:raise ValueError('Unrecognized CMFGEN source')
start=old.index('\t    IF(ATM(1)%XzV_PRES)THEN',old.index("' Hydrogen Density'"))
end=old.index('\t  CLOSE(UNIT=LU_POP)',start)
block='''! QStarS: output by ion name, not fixed ion index. H0 shifts HI and He.
	    DO ID=1,NUM_IONS
	      IF(.NOT. ATM(ID)%XzV_PRES)CYCLE
	      SELECT CASE(TRIM(ION_ID(ID)))
	      CASE('H0')
	        WRITE(LU_POP,'(A)')' H0 populations'
	        WRITE(LU_POP,'(A,I4)')' Number of Hminus levels:',ATM(ID)%NXzV_F
	      CASE('HI')
	        WRITE(LU_POP,'(A)')' HI populations'
	        WRITE(LU_POP,'(A,I4)')' Number of Hydrogen levels:',ATM(ID)%NXzV_F
	      CASE('HeI')
	        WRITE(LU_POP,'(A)')' HeI populations'
	        WRITE(LU_POP,'(A,I4)')' Number of Helium I levels:',ATM(ID)%NXzV_F
	      CASE('He2')
	        WRITE(LU_POP,'(A)')' He2 populations'
	        WRITE(LU_POP,'(A,I4)')' Number of Helium II levels:',ATM(ID)%NXzV_F
	      CASE DEFAULT
	        CYCLE
	      END SELECT
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(ID))//
	1             ' oscillator date:',ATM(ID)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(ID)%XzV_F
	      SELECT CASE(TRIM(ION_ID(ID)))
	      CASE('HI')
	        WRITE(LU_POP,'(A)')' DHYD population'
	      CASE('HeI')
	        WRITE(LU_POP,'(A)')' DHeI population'
	      CASE('He2')
	        WRITE(LU_POP,'(A)')' DHe2 population'
	      CASE DEFAULT
	        CYCLE
	      END SELECT
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(ID)%DXzV_F
	    END DO
!
'''
new=old[:start]+block+old[end:];f.write_text(new)
out=Path(__file__).parent
(out/'population_export.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),
    fromfile='new_main/cmfgen_sub.f',tofile='new_main/cmfgen_sub.f')))
(out/'population_export_patch_manifest.json').write_text(json.dumps(dict(original_sha256=oldhash,
    patched_sha256=sha256(f),physical_arrays_modified=False,scope='RVTJ output ion names including H-'),indent=2)+'\n')
