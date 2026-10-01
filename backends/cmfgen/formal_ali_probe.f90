! Native-library fixture for the QS_FRM_ALI transfer operator (build_13).
! Synthetic spherical beta-law wind with a sharp continuum edge and strong
! coherent scattering. Two independent blue-to-red sweeps on identical input:
!   NATIVE: formal factors + MOM_JREL_V9 alternation, no J/H/XM checks.
!   ALI:    formal solution + local ALI + Ng, as in COMP_J_BLANK (QS_FRM_ALI).
! For ALI, the full Lambda operator is assembled column by column at every
! frequency, and the fixed point is compared with a dense direct solve.
program formal_ali_probe
use set_kind_module
use mod_ray_mom_store
use mod_space_grid_v2
implicit none
#ifndef NDV
#define NDV 40
#endif
integer,parameter::nd=NDV,nc=10,np=nd+nc
integer::iprof,lu_prof
real(kind=ldp)::jnative(nd)
integer::nf
real(kind=ldp)::r(nd),v(nd),sigma(nd),p(np),dop(nd),t(nd),ne(nd)
real(kind=ldp)::eta(nd),chi(nd),esec(nd),ta(nd),tc(nd),fedd(nd),rj(nd),rsqh(nd),iplus(np)
real(kind=ldp),allocatable::nu(:)
real(kind=ldp)::front_r,front_fac,bnue,dbb,dlnu,hib,hob,beta_law,vinf,v0,rmax,kap_hi,kap_lo,nu_edge
real(kind=ldp)::a(nd,nd),c0(nd),mat(nd,nd),rhs(nd),jd(nd),jali(nd),t1,t2,res
real(kind=ldp)::ng(nd,4),q1,q2,q3,a11,a12,a22,cc1,cc2,na,nb,det,tmp(nd)
real(kind=ldp)::min_native,min_formal,min_ray,max_fp_err,max_lstar_excess,min_lambda,min_c0,min_lstar_ratio
real(kind=ldp)::max_lstar_ratio,worst_nu_native
integer::i,k,ifr,it,nhist,ipiv(nd),info,max_it,n_native_neg,ip,iz,sweep,n_ng
logical::first,newf,conv
character(len=6)::method='LOGMON'
real(kind=ldp),parameter::tol=1.0e-12_ldp
! Inputs: edge opacities, radial opacity front (contrast inside front_r),
! terminal speed (km/s) and frequency count.
read(*,*)kap_hi,kap_lo,front_r,front_fac,vinf,nf
allocate(nu(nf))
! Geometry and flow: R in 1e10 cm, V in km/s, sigma = dlnV/dlnR - 1.
rmax=30.0_ldp;v0=10.0_ldp;beta_law=1.0_ldp
do i=1,nd
  r(i)=rmax**(real(nd-i,ldp)/(nd-1))
enddo
r(nd)=1.0_ldp
do i=1,nd
  t1=1.0_ldp-0.999_ldp/r(i)
  v(i)=v0+(vinf-v0)*t1**beta_law
  sigma(i)=(vinf-v0)*beta_law*t1**(beta_law-1)*0.999_ldp/r(i)/v(i)-1.0_ldp
  t(i)=2.0_ldp*(1.0_ldp/r(i))**0.5_ldp+0.3_ldp          ! 10^4 K
  ne(i)=1.0_ldp/r(i)**2/(v(i)/v0)
enddo
dop=1.0e4_ldp
call impar(p,r,r(nd),nc,nd,np)
do k=1,nf
  nu(k)=3.0_ldp*(0.5_ldp)**(real(k-1,ldp)/(nf-1))           ! 3 -> 1.5 (10^15 Hz)
enddo
nu_edge=2.2_ldp
iprof=nf
do k=1,nf
  if(nu(k)<nu_edge)then
    iprof=k;exit
  endif
enddo
max_fp_err=0;max_lstar_excess=-huge(1.0_ldp);min_lambda=huge(1.0_ldp);min_c0=huge(1.0_ldp)
min_native=huge(1.0_ldp);min_formal=huge(1.0_ldp);min_ray=huge(1.0_ldp);max_it=0;n_native_neg=0
min_lstar_ratio=huge(1.0_ldp);max_lstar_ratio=0;worst_nu_native=0;n_ng=0
do sweep=1,2
  qs_frm_ali=(sweep==2)
  rj=0
  do ifr=1,nf
    first=(ifr==1)
    dlnu=0;if(.not.first)dlnu=log(nu(ifr-1)/nu(ifr))
    call model(nu(ifr))
    newf=.true.
    if(sweep==1)then
      ! Native: alternate formal factors and moment solution until f settles.
      do it=1,30
        ta=eta+esec*rj
        call cmf_formal_rel_v4(ta,chi,esec,v,sigma,r,p,tc,fedd,hib,hob,iplus,nu(ifr),dlnu,bnue,dbb,&
             'DIFFUSION',.false.,dop,0.5_ldp,0.0_ldp,method,first,newf,nc,np,nd)
        newf=.false.
        tmp=fedd
        call mom_jrel_v9(eta,chi,esec,v,sigma,r,rj,rsqh,hib,hob,dop,0.5_ldp,nu(ifr),dlnu,dbb,&
             'NONE','NONE','NONE',0.0_ldp,'DIFFUSION','HONJ',method,.true.,'N_ON_J',.true.,.true.,first,nd)
        if(minval(rj)<min_native)then
          min_native=minval(rj);worst_nu_native=nu(ifr)
        endif
        if(it>1 .and. maxval(abs(tmp-fedd))<1.0e-4_ldp)exit
      enddo
      if(minval(rj)<0)n_native_neg=n_native_neg+1
      if(ifr==iprof)jnative=rj
    else
      ! ALI sweep, same algorithm as COMP_J_BLANK with QS_FRM_ALI.
      nhist=0;conv=.false.;it=0
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
            t1=esec(i)*qs_lstar(i)
            t2=tc(i)
            if(t1>0 .and. t1<1)then
              t2=(tc(i)-t1*rj(i))/(1-t1)
              if(t2<0)t2=tc(i)
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
              if(all(tmp>0))then
                rj=tmp;n_ng=n_ng+1
              endif
            endif
          endif
        endif
      enddo
      max_it=max(max_it,it);jali=rj
      if(ifr==iprof)then
        open(newunit=lu_prof,file='profile.txt',status='replace')
        write(lu_prof,'(A,ES24.15)')'# nu ',nu(ifr)
        do i=1,nd
          write(lu_prof,'(3ES24.15)')r(i),jali(i),jnative(i)
        enddo
        close(lu_prof)
      endif
      min_formal=min(min_formal,minval(jali))
      do ip=1,np
        do iz=1,ray(ip)%nz
          min_ray=min(min_ray,ray(ip)%i_p(iz),ray(ip)%i_m(iz))
        enddo
      enddo
      ! Full operator J = A*eta_total + c0, with the same stored blue-side rays.
      ta=0
      call cmf_formal_rel_v4(ta,chi,esec,v,sigma,r,p,c0,fedd,hib,hob,iplus,nu(ifr),dlnu,bnue,dbb,&
           'DIFFUSION',.false.,dop,0.5_ldp,0.0_ldp,method,first,newf,nc,np,nd)
      do k=1,nd
        ta=0;ta(k)=eta(k)
        call cmf_formal_rel_v4(ta,chi,esec,v,sigma,r,p,tc,fedd,hib,hob,iplus,nu(ifr),dlnu,bnue,dbb,&
             'DIFFUSION',.false.,dop,0.5_ldp,0.0_ldp,method,first,newf,nc,np,nd)
        a(:,k)=(tc-c0)/eta(k)
        t1=a(k,k)
        max_lstar_excess=max(max_lstar_excess,(qs_lstar(k)-t1)/t1)
        min_lstar_ratio=min(min_lstar_ratio,qs_lstar(k)/t1);max_lstar_ratio=max(max_lstar_ratio,qs_lstar(k)/t1)
      enddo
      min_lambda=min(min_lambda,minval(a)/maxval(a));min_c0=min(min_c0,minval(c0))
      do k=1,nd
        mat(:,k)=-a(:,k)*esec(k)
        mat(k,k)=mat(k,k)+1
      enddo
      rhs=matmul(a,eta)+c0
      call dgesv(nd,1,mat,nd,ipiv,rhs,nd,info)
      if(info/=0)stop 'direct solve failed'
      jd=rhs
      max_fp_err=max(max_fp_err,maxval(abs(jali-jd)/jd))
      ! Restore stored rays to the converged field for the next frequency.
      ta=eta+esec*jali
      call cmf_formal_rel_v4(ta,chi,esec,v,sigma,r,p,tc,fedd,hib,hob,iplus,nu(ifr),dlnu,bnue,dbb,&
           'DIFFUSION',.false.,dop,0.5_ldp,0.0_ldp,method,first,newf,nc,np,nd)
      rj=tc
    endif
  enddo
enddo
write(*,'(A,ES24.15)')'RESULT native_min_moment_J ',min_native
write(*,'(A,ES24.15)')'RESULT native_worst_frequency ',worst_nu_native
write(*,'(A,I8)')'RESULT native_frequencies_with_negative_J ',n_native_neg
write(*,'(A,ES24.15)')'RESULT ali_min_J ',min_formal
write(*,'(A,ES24.15)')'RESULT ali_min_ray_intensity ',min_ray
write(*,'(A,ES24.15)')'RESULT ali_max_rel_fixed_point_error_vs_direct ',max_fp_err
write(*,'(A,ES24.15)')'RESULT lambda_min_entry_over_max ',min_lambda
write(*,'(A,ES24.15)')'RESULT lambda_min_affine_part ',min_c0
write(*,'(A,2ES24.15)')'RESULT lstar_over_diag_min_max ',min_lstar_ratio,max_lstar_ratio
write(*,'(A,I8)')'RESULT ali_max_formal_solutions ',max_it
write(*,'(A,I8)')'RESULT ng_steps ',n_ng
contains
subroutine model(f)
  real(kind=ldp),intent(in)::f
  real(kind=ldp)::kap,x
  integer::j
  kap=kap_lo;if(f>nu_edge)kap=kap_hi
  do j=1,nd
    x=4.8_ldp*f/t(j)
    chi(j)=kap*ne(j)**1.5_ldp+2.0_ldp*ne(j)
    if(r(j)<front_r)chi(j)=chi(j)+front_fac*kap*ne(j)**1.5_ldp
    esec(j)=2.0_ldp*ne(j)
    eta(j)=(chi(j)-esec(j))*f**3/(exp(x)-1)
  enddo
  x=4.8_ldp*f/t(nd);bnue=f**3/(exp(x)-1)
  t1=4.8_ldp*f/t(nd-1);t1=f**3/(exp(t1)-1)
  dbb=(bnue-t1)/(r(nd-1)-r(nd))
end subroutine
end program
