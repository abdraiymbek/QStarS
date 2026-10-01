! Nonlocal response check for QS_FRM_ALI (build_13 native library).
! Same synthetic wind as formal_ali_probe.f90. The radiation field is the
! formal-solution ALI fixed point; the response is VAR_JREL_V5 (moment
! linearization) evaluated about that formal state through QS_FORMAL_MOM_STATE,
! exactly as DO_VAR_CONT uses it. Central finite differences of the complete
! formal-ALI sweep, for an additive eta or chi perturbation at depth D0 applied
! at every frequency (the convention TX accumulates), give the exact discrete
! response it is compared with. Scattering coefficient is held in chi changes.
program formal_ali_jacobian_probe
use set_kind_module
use mod_ray_mom_store
implicit none
#ifndef NDV
#define NDV 40
#endif
integer,parameter::nd=NDV,nc=10,np=nd+nc,nm=2
integer::nf,d0,mode,lu
real(kind=ldp),allocatable::nu(:)
real(kind=ldp)::r(nd),v(nd),sigma(nd),p(np),dop(nd),t(nd),ne(nd)
real(kind=ldp)::eta(nd),chi(nd),esec(nd),theta(nd),ta(nd),tc(nd),fedd(nd),rj(nd),rsqh(nd),iplus(np)
real(kind=ldp)::tx(nd,nd,nm),tvx(nd-1,nd,nm),ki(nd,nd,nm),work(nd,nd),rhsc(nd-1,nd)
real(kind=ldp)::txt(nd),txdt(nd),tvxt(nd),tvxdt(nd)
real(kind=ldp)::jbase(nd),jp(nd),jm(nd),resp_var(nd,2),resp_fd(nd,2),delta(2)
real(kind=ldp)::bnue,dbb,dlnu,hib,hob,vinf,v0,rmax,kap_hi,kap_lo,nu_edge,front_r,front_fac,t1,t2,res,eps
real(kind=ldp)::ng(nd,4),q1,q2,q3,a11,a12,a22,cc1,cc2,na,nb,det,tmp(nd)
integer::i,k,ifr,it,nhist,iout
logical::first,newf,conv,do_tx(nm)
character(len=6)::method='LOGMON'
real(kind=ldp),parameter::tol=1.0e-13_ldp
read(*,*)kap_hi,kap_lo,front_r,front_fac,vinf,nf,d0,eps
allocate(nu(nf))
rmax=30.0_ldp;v0=10.0_ldp
do i=1,nd
  r(i)=rmax**(real(nd-i,ldp)/(nd-1))
enddo
r(nd)=1.0_ldp
do i=1,nd
  t1=1.0_ldp-0.999_ldp/r(i)
  v(i)=v0+(vinf-v0)*t1
  sigma(i)=(vinf-v0)*0.999_ldp/r(i)/v(i)-1.0_ldp
  t(i)=2.0_ldp*(1.0_ldp/r(i))**0.5_ldp+0.3_ldp
  ne(i)=1.0_ldp/r(i)**2/(v(i)/v0)
enddo
dop=1.0e4_ldp
call impar(p,r,r(nd),nc,nd,np)
do k=1,nf
  nu(k)=3.0_ldp*(0.5_ldp)**(real(k-1,ldp)/(nf-1))
enddo
nu_edge=2.2_ldp
iout=nf
do k=1,nf
  if(nu(k)<nu_edge)then
    iout=k;exit
  endif
enddo
qs_frm_ali=.true.
! Perturbation sizes: eps times the value at D0 at the output frequency.
call model(nu(iout),0,0.0_ldp)
delta(1)=eps*eta(d0);delta(2)=eps*(chi(d0)-esec(d0))
call sweep(0,0.0_ldp,jbase,.true.)
do mode=1,2
  call sweep(mode,delta(mode),jp,.false.)
  call sweep(mode,-delta(mode),jm,.false.)
  resp_fd(:,mode)=(jp-jm)/(2*delta(mode))
enddo
open(newunit=lu,file='jacobian.txt',status='replace')
write(lu,'(A,ES24.15,A,I4)')'# nu_out ',nu(iout),' d0 ',d0
write(lu,'(A)')'# r J dJ/deta(FD) dJ/deta(VAR) dJ/dchi(FD) dJ/dchi(VAR)'
do i=1,nd
  write(lu,'(6ES24.15)')r(i),jbase(i),resp_fd(i,1),resp_var(i,1),resp_fd(i,2),resp_var(i,2)
enddo
close(lu)
do mode=1,2
  t1=maxval(abs(resp_var(:,mode)-resp_fd(:,mode)))/maxval(abs(resp_fd(:,mode)))
  t2=abs(resp_var(d0,mode)-resp_fd(d0,mode))/abs(resp_fd(d0,mode))
  write(*,'(A,I2,A,ES12.4,A,ES12.4)')'RESULT mode ',mode,' max|VAR-FD|/max|FD| ',t1,' local rel ',t2
enddo
contains
subroutine model(f,pm,dv)
  real(kind=ldp),intent(in)::f,dv
  integer,intent(in)::pm
  real(kind=ldp)::kap,x,b1
  integer::j
  kap=kap_lo;if(f>nu_edge)kap=kap_hi
  do j=1,nd
    x=4.8_ldp*f/t(j)
    chi(j)=kap*ne(j)**1.5_ldp+2.0_ldp*ne(j)
    if(r(j)<front_r)chi(j)=chi(j)+front_fac*kap*ne(j)**1.5_ldp
    esec(j)=2.0_ldp*ne(j)
    eta(j)=(chi(j)-esec(j))*f**3/(exp(x)-1)
  enddo
  if(pm==1)eta(d0)=eta(d0)+dv
  if(pm==2)chi(d0)=chi(d0)+dv
  x=4.8_ldp*f/t(nd);bnue=f**3/(exp(x)-1)
  b1=4.8_ldp*f/t(nd-1);b1=f**3/(exp(b1)-1)
  dbb=(bnue-b1)/(r(nd-1)-r(nd))
end subroutine
subroutine sweep(pm,dv,jout,with_var)
  integer,intent(in)::pm
  real(kind=ldp),intent(in)::dv
  real(kind=ldp),intent(out)::jout(nd)
  logical,intent(in)::with_var
  rj=0
  do ifr=1,nf
    first=(ifr==1)
    dlnu=0;if(.not.first)dlnu=log(nu(ifr-1)/nu(ifr))
    call model(nu(ifr),pm,dv)
    newf=.true.;nhist=0;conv=.false.;it=0
    do while(.not.conv)
      ta=eta+esec*rj
      call cmf_formal_rel_v4(ta,chi,esec,v,sigma,r,p,tc,fedd,hib,hob,iplus,nu(ifr),dlnu,bnue,dbb,&
           'DIFFUSION',.false.,dop,0.5_ldp,0.0_ldp,method,first,newf,nc,np,nd)
      newf=.false.;it=it+1
      res=0
      do i=1,nd
        t2=max(tc(i),rj(i));if(t2>0)res=max(res,abs(tc(i)-rj(i))/t2)
      enddo
      if(res<tol .or. it>=5000)then
        rj=tc;conv=.true.
      else
        do i=1,nd
          t1=esec(i)*qs_lstar(i);t2=tc(i)
          if(t1>0 .and. t1<1)then
            t2=(tc(i)-t1*rj(i))/(1-t1);if(t2<0)t2=tc(i)
          endif
          rj(i)=t2
        enddo
        nhist=nhist+1;ng(:,1:3)=ng(:,2:4);ng(:,4)=rj
        if(nhist>=4 .and. it>=5)then
          nhist=0;a11=0;a12=0;a22=0;cc1=0;cc2=0
          do i=1,nd
            t1=1/ng(i,4)**2
            q1=ng(i,4)-2*ng(i,3)+ng(i,2);q2=ng(i,4)-ng(i,3)-ng(i,2)+ng(i,1);q3=ng(i,4)-ng(i,3)
            a11=a11+t1*q1*q1;a12=a12+t1*q1*q2;a22=a22+t1*q2*q2;cc1=cc1+t1*q1*q3;cc2=cc2+t1*q2*q3
          enddo
          det=a11*a22-a12*a12
          if(det>0)then
            na=(cc1*a22-cc2*a12)/det;nb=(cc2*a11-cc1*a12)/det
            tmp=(1-na-nb)*ng(:,4)+na*ng(:,3)+nb*ng(:,2)
            if(all(tmp>0))rj=tmp
          endif
        endif
      endif
    enddo
    call qs_formal_mom_state(r,v,nu(ifr),'N_ON_J',dbb,chi(nd),'DIFFUSION',first,nd)
    if(with_var)then
      theta=esec/chi
      do_tx=.true.
      call var_jrel_v5(eta,chi,esec,theta,v,sigma,r,tx,tvx,txt,txdt,tvxt,tvxdt,ki,work,rhsc,&
           first,nu(ifr),dlnu,'NONE','NONE','NONE','DIFFUSION','HONJ',0.0_ldp,1.0_ldp,dbb,1.0_ldp,1.0_ldp,&
           .true.,.true.,do_tx,method,nd,nm,nm)
      if(ifr==iout)then
        resp_var(:,1)=tx(:,d0,2);resp_var(:,2)=tx(:,d0,1)
      endif
    endif
    if(ifr==iout)jout=rj
  enddo
end subroutine
end program
