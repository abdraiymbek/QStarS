#!/usr/bin/env python3
"""Compile actual native Rayleigh call sites; independent opacity/Jacobian check.

Use the same gfortran environment as the backend. This is a component test,
not a test of Rayleigh's published fit accuracy or a wind validation.
"""
import argparse,json,subprocess,os,platform
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import sha256

def check(root,output):
    root=root.resolve();output.mkdir(parents=True,exist_ok=False)
    op=root/'new_main/mod_subs/comp_opac.f';va=root/'new_main/mod_subs/do_var_cont.f'
    s=op.read_text();opacity=s[s.index('\t  CHI_RAY(1:ND)=0.0_LDP'):s.index('\t  CHI_SCAT(1:ND)=')]
    s=va.read_text();start=s.index('\t      IF(SPECIES_PRES(1))THEN',s.index("CALL TUNE(1,'VC_VCHI')"))
    derivative=s[start:s.index('! Now must do line terms.',start)].rstrip('!\n')
    # Exercise H0 present, nonuniform depths and 30 HI levels (wrong stride cannot hide).
    program='''program probe
 use set_kind_module
 implicit none
 integer, parameter :: ND=3
 type ion
   real(kind=LDP),allocatable :: XzV_F(:,:), AXzV_F(:,:), EDGEXZV_F(:)
   integer :: NXzV_F, EQXzV
 end type
 type(ion) :: ATM(2)
 character(3) :: ION_ID(2)
 logical :: SPECIES_PRES(1),INCL_RAY_SCAT
 real(kind=LDP) :: CHI_RAY(ND), CONT_FREQ, VJ(40,1,ND), dj_LOC(5,1,ND)
 integer :: ID,I,J,K,L
 ION_ID=[character(3)::'H0','HI']; SPECIES_PRES=.true.; INCL_RAY_SCAT=.true.
 do ID=1,2
   ATM(ID)%NXzV_F=1
   if(ID==2)ATM(ID)%NXzV_F=30
   I=ATM(ID)%NXzV_F
   allocate(ATM(ID)%XzV_F(I,ND),ATM(ID)%AXzV_F(I,I),ATM(ID)%EDGEXZV_F(I))
   ATM(ID)%XzV_F=0.7_LDP; ATM(ID)%AXzV_F=0; ATM(ID)%EDGEXZV_F=0.182365_LDP
   ATM(ID)%EQXzV=ID
 end do
 ATM(2)%XzV_F(1,:)=[2.0e10_LDP,5.0e11_LDP,9.0e12_LDP]
 ATM(2)%EDGEXZV_F(1)=3.288087_LDP
 CONT_FREQ=0.6_LDP
'''+opacity+'''
 write(*,'(3ES25.16)') CHI_RAY
 VJ=0; dj_LOC=1; J=1
 do K=1,ND
   L=K
'''+derivative+'''
 end do
 write(*,'(3ES25.16)') VJ(2,1,:)
 write(*,'(3ES25.16)') VJ(1,1,:)
end program
'''
    src=output/'probe.f90';src.write_text(program)
    # The native snippets are fixed-form; remove continuation column markers in
    # this small call and retain ordinary indentation in a free-form harness.
    text=src.read_text().replace('ATM(ID)%EDGEXZV_F,\n\t1','ATM(ID)%EDGEXZV_F, &\n     &')
    src.write_text(text)
    exe=output/'probe';ray=root/'newsubs/rayleigh_scat.f'
    command=['gfortran','-ffixed-line-length-none','-I'+str(root/'mod'),str(src),str(ray),'-o',str(exe)]
    env=dict(os.environ)
    if platform.system()=='Darwin':env['MACOSX_DEPLOYMENT_TARGET']='14.0'
    r=subprocess.run(command,capture_output=True,text=True,env=env)
    (output/'build.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Native probe compilation failed; see build.log')
    r=subprocess.run([str(exe)],capture_output=True,text=True,check=True)
    (output/'native.txt').write_text(r.stdout+r.stderr);data=np.loadtxt(output/'native.txt')
    x=(.6/3.288087)**2;cross=6.65e-15*x*x*(1.26+5.068*x+708.3*x**5)
    target=cross*np.array([2e10,5e11,9e12])
    errors=dict(opacity=float(max(abs(data[0]/target-1))),jacobian=float(max(abs(data[1]/cross-1))),
                spurious_hminus_derivative=float(max(abs(data[2]))/cross))
    result=dict(status='PASS' if max(errors.values())<1e-12 else 'FAIL',errors=errors,limit=1e-12,
        source_hashes={str(p.relative_to(root)):sha256(p) for p in [op,va,ray]},
        scope='actual upstream call-site extraction with H- present; low-frequency HI analytic formula and linear derivative',
        test_source_sha256=sha256(src),executable_sha256=sha256(exe),production_certified=False)
    (output/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);check(**vars(p.parse_args()))
