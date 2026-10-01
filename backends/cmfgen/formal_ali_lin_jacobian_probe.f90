! Local ALI-Jacobian check for QS_ALI_LIN (build_15 native library).
! Same synthetic wind as formal_ali_jacobian_probe.f90. The radiation field is
! the formal-solution ALI fixed point; the response is the local formal-solution
! ALI operator installed in VAR_JREL_V5 under QS_ALI_LIN (approximate Jacobian;
! drives Newton only). Central finite differences of the converged formal-ALI
! map at the output frequency, for a fractional eta or chi perturbation applied
! at depth i (at the output frequency only), give the exact local response:
!     Psi_i = L*_i/(1-esec_i*L*_i)          vs  dJ_i/deta_i (FD)
!     -Psi_i*S_i                            vs  dJ_i/dchi_i (FD)
! Scattering coefficient is held in chi changes. The perturbation is applied at
! the output frequency only so the FD measures the strictly local map.
program formal_ali_lin_jacobian_probe
use set_kind_module
use mod_ray_mom_store
implicit none
#ifndef NDV
#define NDV 40
#endif
integer,parameter::nd=NDV,nc=10,np=nd+nc,nm=2
integer::nf,d0,mode,lu,ipd
real(kind=ldp),allocatable::nu(:)
real(kind=ldp)::r(nd),v(nd),sigma(nd),p(np),dop(nd),t(nd),ne(nd)
real(kind=ldp)::eta(nd),chi(nd),esec(nd),theta(nd),ta(nd),tc(nd),fedd(nd),rj(nd),rsqh(nd),iplus(np)
real(kind=ldp)::tx(nd,nd,nm),tvx(nd-1,nd,nm),ki(nd,nd,nm),work(nd,nd),rhsc(nd-1,nd)
real(kind=ldp)::txt(nd),txdt(nd),tvxt(nd),tvxdt(nd)
real(kind=ldp)::jbase(nd),jp(nd),jm(nd),psi(nd),mps(nd),fd_eta(nd),fd_chi(nd)
real(kind=ldp)::eta0(nd),chi0(nd),esec0(nd)
real(kind=ldp)::bnue,dbb,dlnu,hib,hob,vinf,v0,rmax,kap_hi,kap_lo,nu_edge,front_r,front_fac,t1,t2,res,eps
real(kind=ldp)::ng(nd,4),q1,q2,q3,a11,a12,a22,cc1,cc2,na,nb,det,tmp(nd)
real(kind=ldp)::rmin_e,rmax_e,rmin_c,rmax_c,cmax,odmax
integer::i,k,ifr,it,nhist,iout,nsum,nsgn
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
qs_ali_lin=.true.
! Base model coefficients at the output frequency set the perturbation scale.
call model(nu(iout),0,1,0.0_ldp)
eta0=eta;chi0=chi;esec0=esec
call sweep(0,1,0.0_ldp,jbase,.true.)
! psi/mps hold the ALI-branch TX diagonal at the output frequency.
! Finite-difference diagonal of the converged formal-ALI map at iout.
do ipd=1,nd
  call sweep(1,ipd,eps,jp,.false.)
  call sweep(1,ipd,-eps,jm,.false.)
  fd_eta(ipd)=(jp(ipd)-jm(ipd))/(2.0_ldp*eps*eta0(ipd))
  call sweep(2,ipd,eps,jp,.false.)
  call sweep(2,ipd,-eps,jm,.false.)
  fd_chi(ipd)=(jp(ipd)-jm(ipd))/(2.0_ldp*eps*(chi0(ipd)-esec0(ipd)))
enddo
open(newunit=lu,file='ali_lin_jacobian.txt',status='replace')
write(lu,'(A,ES24.15)')'# nu_out ',nu(iout)
write(lu,'(A)')'# r J Psi=dJ/deta(TX) dJ/deta(FD) -Psi*S=dJ/dchi(TX) dJ/dchi(FD) Psi/FDeta chiFD/chiTX'
cmax=0.0_ldp
do i=1,nd
  cmax=max(cmax,abs(fd_eta(i)),abs(fd_chi(i)))
enddo
do i=1,nd
  t1=0.0_ldp;t2=0.0_ldp
  if(abs(fd_eta(i))>1.0e-6_ldp*cmax)t1=psi(i)/fd_eta(i)
  if(abs(fd_chi(i))>1.0e-6_ldp*cmax)t2=fd_chi(i)/mps(i)
  write(lu,'(8ES24.15)')r(i),jbase(i),psi(i),fd_eta(i),mps(i),fd_chi(i),t1,t2
enddo
close(lu)
! Off-diagonal leakage check: the ALI response must be diagonal in TX(:,:,1:2).
odmax=0.0_ldp
do i=1,nd
  do k=1,nd
    if(k/=i)odmax=max(odmax,abs(tx(i,k,1)),abs(tx(i,k,2)))
  enddo
enddo
rmin_e=1.0e30_ldp;rmax_e=0.0_ldp;nsum=0
rmin_c=1.0e30_ldp;rmax_c=0.0_ldp;nsgn=0
do i=1,nd
  if(abs(fd_eta(i))>1.0e-6_ldp*cmax)then
    t1=psi(i)/fd_eta(i)
    rmin_e=min(rmin_e,t1);rmax_e=max(rmax_e,t1);nsum=nsum+1
  endif
  if(abs(fd_chi(i))>1.0e-6_ldp*cmax)then
    t1=fd_chi(i)/mps(i)
    rmin_c=min(rmin_c,t1);rmax_c=max(rmax_c,t1)
    if(sign(1.0_ldp,fd_chi(i))==sign(1.0_ldp,mps(i)))nsgn=nsgn+1
  endif
enddo
write(*,'(A,ES12.4)')'RESULT offdiag_max ',odmax
write(*,'(A,I3,A,ES12.4,A,ES12.4)')'RESULT Psi/FD_deta n',nsum,' min ',rmin_e,' max ',rmax_e
write(*,'(A,I3,A,ES12.4,A,ES12.4)')'RESULT FDdchi/(-Psi*S) n',nsgn,' min ',rmin_c,' max ',rmax_c
contains
subroutine model(f,pm,pm_depth,frac)
  real(kind=ldp),intent(in)::f,frac
  integer,intent(in)::pm,pm_depth
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
  if(pm==1)eta(pm_depth)=eta(pm_depth)*(1.0_ldp+frac)
  if(pm==2)chi(pm_depth)=chi(pm_depth)+frac*(chi(pm_depth)-esec(pm_depth))
  x=4.8_ldp*f/t(nd);bnue=f**3/(exp(x)-1)
  b1=4.8_ldp*f/t(nd-1);b1=f**3/(exp(b1)-1)
  dbb=(bnue-b1)/(r(nd-1)-r(nd))
end subroutine
subroutine sweep(pm,pm_depth,frac,jout,with_var)
  integer,intent(in)::pm,pm_depth
  real(kind=ldp),intent(in)::frac
  real(kind=ldp),intent(out)::jout(nd)
  logical,intent(in)::with_var
  integer::pm_eff
  rj=0
  do ifr=1,nf
    first=(ifr==1)
    dlnu=0;if(.not.first)dlnu=log(nu(ifr-1)/nu(ifr))
    pm_eff=0;if(ifr==iout)pm_eff=pm
    call model(nu(ifr),pm_eff,pm_depth,frac)
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
        do i=1,nd
          psi(i)=tx(i,i,2);mps(i)=tx(i,i,1)
        enddo
      endif
    endif
    if(ifr==iout)jout=rj
  enddo
end subroutine
end program
