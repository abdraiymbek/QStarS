#!/usr/bin/env python3
"""Native cool H-minus FF regression: independent entry, energy and Jacobian.

Compile actual opacity/variation/EHB sources with bound-free contributions
isolated to zero. This isolates FF implementation, never certifies an atom or
relaxed non-LTE response. No scientific run changes bound-free processes.
"""
import argparse, json, os, platform, subprocess
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import sha256

MOCK = '''module mod_lev_dis_blk
 use set_kind_module
 logical :: MOD_DO_LEV_DIS=.false.
 real(kind=LDP) :: PHOT_DIS_PARAMETER=0.0001_LDP
 real(kind=LDP) :: X_LEV_DIS(3),A_LEV_DIS(3),B_LEV_DIS(3)
end module
module mod_cmfgen
 use set_kind_module
 type ion
  logical :: XzV_PRES
  real(kind=LDP),allocatable :: XzV(:,:),GIONXzV_F(:)
  integer :: EQXzV,NXzV
  real(kind=LDP) :: ZXzV
 end type
 integer,parameter :: NUM_IONS=2
 type(ion) :: ATM(3)
 character(3) :: ION_ID(2)
 real(kind=LDP) :: ED(3),T(3)
end module
'''
PROGRAM = '''program probe
 use set_kind_module
 use mod_cmfgen
 implicit none
 integer,parameter :: ND=3,NT=40
 real(kind=LDP) :: chi(ND),eta(ND),cool(ND),bp(ND),jp(ND),em(ND),nu
 real(kind=LDP) :: vc(NT,ND),ve(NT,ND),vcn(NT,ND),ven(NT,ND)
 real(kind=LDP) :: hs(1,ND),di(30,ND),dl(30,ND),edge(1),g(1),x
 real(kind=LDP) :: chiE(ND),etaE(ND),vcE(NT,ND),veE(NT,ND)
 real(kind=LDP) :: CHIBF,CHIFF,HDKT,TWOHCSQ
 common/constants/CHIBF,CHIFF,HDKT,TWOHCSQ
 integer :: map(1),i,k,mode,eq
 logical :: imp(NT)
 read(*,*)mode,nu,T,ED,di(1,:)
 CHIBF=1;CHIFF=1;HDKT=4.79924307336622_LDP;TWOHCSQ=0.014745003464_LDP
 chi=0;eta=0;cool=0;hs=1;edge=3.288;g=2;dl=0;map=1;imp=.true.;eq=2
 do k=1,ND
  x=HDKT*nu/T(k);em(k)=exp(-x)
  bp(k)=TWOHCSQ*nu**3*em(k)/(1-em(k));jp(k)=0.37_LDP*bp(k)
 end do
 if(mode==1)call DO_H0_FF_COOL(cool,di(1,:),ED,T,bp,jp,nu,ND)
 call DO_H0_FF(eta,chi,di(1,:),ED,T,em,nu,7,ND)
 if(mode/=1)call DO_H0_FF_COOL(cool,di(1,:),ED,T,bp,jp,nu,ND)
 do k=1,ND
  write(*,'(A,3ES25.16)')'DATA ',chi(k),eta(k),cool(k)
 end do
 if(mode/=0)return
 di(2:,:)=0
 vc=0;ve=0;vcn=0;ven=0
 call VAR_OP_V10(vc,ve,vcn,ven,hs,hs,hs,hs,1,hs,edge,1,map, &
  di,dl,dl,30,1,1,ED,T,em,imp,nu,0.0_LDP,1,.true.,1,eq,NT,ND,.false.)
 do k=1,ND
  write(*,'(A,6ES25.16)')'DERIV ',vc(eq,k),ve(eq,k),vc(NT-1,k),ve(NT-1,k),vc(NT,k),ve(NT,k)
 end do
 ION_ID=[character(3)::'H0','HI']
 do i=1,3
  ATM(i)%XzV_PRES=.false.;ATM(i)%NXzV=30;ATM(i)%EQXzV=eq;ATM(i)%ZXzV=1
  allocate(ATM(i)%XzV(30,ND),ATM(i)%GIONXzV_F(30))
  ATM(i)%XzV=0;ATM(i)%GIONXzV_F=2
 end do
 ATM(1)%XzV_PRES=.true.;ATM(1)%ZXzV=0
 ATM(2)%XzV(1,:)=di(1,:);ATM(2)%XzV(2,:)=0.2_LDP*di(1,:)
 call COMP_FREE_FREE_V2(chiE,etaE,vcE,veE,nu,nu,.true.,.true.,.true.,ND,NT)
 do k=1,ND
  write(*,'(A,8ES25.16)')'EHB ',chiE(k),etaE(k),vcE(eq,k),veE(eq,k), &
    vcE(NT-1,k),veE(NT-1,k),vcE(NT,k),veE(NT,k)
 end do
end program
subroutine SUB_PHOT_GEN(ID,a,nu,edge,n,phot,dis)
 use set_kind_module
 integer :: ID,n,phot
 real(kind=LDP) :: a(n),nu,edge(n)
 logical :: dis
 a=0
end subroutine
subroutine GFF_VEC(a,nu,t,z,n)
 use set_kind_module
 integer :: n
 real(kind=LDP) :: a(n),nu,t(n),z
 a=1
end subroutine
subroutine FF_RES_GAUNT(a,nu,t,id,g,z,n)
 use set_kind_module
 integer :: n,id
 real(kind=LDP) :: a(n),nu,t(n),g(*),z
end subroutine
'''


def check(root, output, table):
    root=root.resolve();output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    (output/'H0_FF').symlink_to(table.resolve())
    (output/'mock.f90').write_text(MOCK);(output/'probe.f90').write_text(PROGRAM)
    sources=['newsubs/do_h0_ff.f','new_main/subs/var_op_v10.f','newsubs/comp_free_free_v2.f']
    cmd=['gfortran','-ffixed-line-length-none','-fcheck=all','-fbacktrace','-I'+str(root/'mod'),
         str(output/'mock.f90'),*[str(root/n) for n in sources],str(output/'probe.f90'),
         '-L'+str(root/'lib'),'-lsubs','-lmod']
    if platform.system()=='Darwin':cmd+=['-framework','Accelerate']
    cmd+=['-o',str(output/'probe')]
    env=dict(os.environ,MACOSX_DEPLOYMENT_TARGET='14.0')
    r=subprocess.run(cmd,cwd=output,env=env,capture_output=True,text=True)
    (output/'build.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Probe build failed; retained build.log')
    cases=[]; maxerr={}; failed=[]
    baseT=np.array([.33,.61,1.13]);baseED=np.array([2e9,3e11,5e12]);baseH=np.array([4e10,7e12,9e13])
    def run(mode,nu,T,ED,H,label):
        inp=' '.join(str(x) for x in [mode,nu,*T,*ED,*H])+'\n'
        r=subprocess.run([str(output/'probe')],input=inp,cwd=output,capture_output=True,text=True)
        (output/(label+'.txt')).write_text('INPUT '+inp+r.stdout+r.stderr)
        data={}
        for tag in ['DATA','DERIV','EHB']:
            rows=[line.split()[1:] for line in r.stdout.splitlines() if line.startswith(tag+' ')]
            if rows:data[tag]=np.array(rows,float)
        cases.append(dict(label=label,exit_code=r.returncode,input=inp.strip()))
        if r.returncode:failed.append(label)
        return data
    for j,nu in enumerate([.03,.2,.7,2.0]):
        first=run(1,nu,baseT,baseED,baseH,f'cool_first_{j}')
        ordinary=run(2,nu,baseT,baseED,baseH,f'op_first_{j}')
        if 'DATA' in first:
            maxerr['entry_order']=max(maxerr.get('entry_order',0),float(np.max(abs(first['DATA']/ordinary['DATA']-1))))
        d=run(0,nu,baseT,baseED,baseH,f'derivative_{j}')
        a=ordinary['DATA'];x=4.79924307336622*nu/baseT
        B=.014745003464*nu**3*np.exp(-x)/(1-np.exp(-x))
        for key,err in [('kirchhoff',np.max(abs(a[:,1]/a[:,0]/B-1))),('energy',np.max(abs(a[:,2]/(1.256637061e-9*a[:,0]*.63*B)-1)))]:
            maxerr[key]=max(maxerr.get(key,0),float(err))
        if 'DERIV' not in d:continue
        for i,param in enumerate(['H','ED','T']):
            args=[baseT.copy(),baseED.copy(),baseH.copy()];idx=2-i;h=1e-5
            args[idx]*=(1+h);plus=run(2,nu,*args,f'plus_{param}_{j}')['DATA'][:,:2]
            args[idx]=[baseT,baseED,baseH][idx]*(1-h);minus=run(2,nu,*args,f'minus_{param}_{j}')['DATA'][:,:2]
            fd=(plus-minus)/(2*h*[baseH,baseED,baseT][i][:,None])
            err=float(np.max(abs(d['DERIV'][:,2*i:2*i+2]-fd)/np.maximum(abs(fd),1e-300)))
            maxerr[param+'_derivative']=max(maxerr.get(param+'_derivative',0),err)
        if 'EHB' in d:
            target=np.column_stack([a[:,:2],d['DERIV']])
            maxerr['EHB_consistency']=max(maxerr.get('EHB_consistency',0),float(np.max(abs(d['EHB']/target-1))))
    result=dict(status='PASS' if not failed and max(maxerr.values())<1e-6 else 'FAIL',errors=maxerr,failed_native_cases=failed,
      threshold=1e-6,cases=cases,source_hashes={n:sha256(root/n) for n in sources},table_sha256=sha256(table),
      executable_sha256=sha256(output/'probe'),driver_sha256=sha256(__file__),
      scope='actual native FF sources, T=3300/6100/11300 K, four frequencies; bound-free isolated in regression only; no relaxed material/population or wind certificate',production_certified=False)
    (output/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','errors','failed_native_cases']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('table',type=Path);check(**vars(p.parse_args()))
