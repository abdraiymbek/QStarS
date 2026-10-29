! SPDX-License-Identifier: GPL-3.0-only
! Differentiated adaptation of CMFGEN formal characteristics (DLM/DJH);
! upstream GPLv3 provenance: backend README.
! Linearized form of qs_formal_tangent%tangent_map: the primal ray quantities
! and all local partial derivatives are computed ONCE per frequency
! (lin_setup); each directional application (lin_apply) is then a linear sweep
! with no exponentials or dual arithmetic along the rays. The node-level
! extrapolation and Steffen (LOGMON) coefficient derivatives are evaluated with
! the validated dual code (O(ND) per application). Same branches, same
! truncated small-dtau polynomial, same moments as the dual reference.
module qs_formal_tangent_lin
use set_kind_module
use omp_lib, only: omp_in_parallel,omp_get_max_threads
use qs_dual
use qs_formal_tangent, only: mon_coeff
use mod_space_grid_v2
use ext_rel_grid_v4, only: r_ext,log_r_ext,nd_ext,nd_add,idmin,idmax
use mod_ray_mom_store, only: qs_lstar
implicit none
private
public :: lin_setup,lin_apply,lin_scattering,lin_ready,lin_apply_multi,lin_em_build,asig_ok,lin_scattering_multi,&
          lin_setup_omp
type ray_lin
 integer :: n=0
 integer,allocatable :: kk(:)          ! interval index, or -node when coincident
 real(kind=ldp),allocatable :: z(:),crv(:),erv(:)
 real(kind=ldp),allocatable :: ctm(:),ctp(:),srcm(:),srcp(:),imv(:),ipv(:)
 real(kind=ldp),allocatable :: eem(:),wum(:),wdm(:),deem(:),dwum(:),dwdm(:),dsm(:)
 real(kind=ldp),allocatable :: eep(:),wup(:),wdp(:),deep(:),dwup(:),dwdp(:),dsp(:)
end type
type(ray_lin),allocatable,save :: rl(:)
real(kind=ldp),save :: s_nuon=0,s_b=0,s_dbb=0,s_chind=1
real(kind=ldp),allocatable,save :: s_eta(:),s_chi(:),s_scat(:)
integer,save :: s_nc=0,s_nray=0,s_nd=0
character(len=12),save :: s_bnd=' '
logical,save :: lin_ready=.false.
! Build 30 (Report 22): explicit scattering operator of the current frequency, A(:,j) = J response of the
! linear formal map to the emissivity s_scat(j) at node j (exactly the GMRES product of lin_scattering
! applied to e_j). Built once per frequency by lin_em_build and shared by all directions; invalidated by
! lin_setup. When asig_ok, gmres_scat uses A*v instead of a ray sweep (same map; summation order only).
real(kind=ldp),allocatable,save :: asig(:,:)
logical,save :: asig_ok=.false.
! Build 31 (Report 22): lin_setup's ray loop in parallel (per-ray data are independent; the primal moments
! are still summed serially in ray order afterwards, so every result is unchanged). Off: QSTARS_JV_NOSETUP_OMP.
logical,save :: lin_setup_omp=.true.
contains

subroutine dweights(t,ee,wu,wd,dee,dwu,dwd)
real(kind=ldp),intent(in)::t
real(kind=ldp),intent(out)::ee,wu,wd,dee,dwu,dwd
type(dual)::td,e0,e1,e2,eed,wdd,wud
td=q(t,1.0_ldp)
if(t>0.01_ldp)then
 eed=exp(-td);e0=1.0_ldp-eed;e1=td-e0
else
 e2=td*td*td*(1.0_ldp-0.25_ldp*td*(1.0_ldp-0.2_ldp*td*(1.0_ldp-td/6.0_ldp*&
    (1.0_ldp-td/7.0_ldp))))/3.0_ldp
 e1=0.5_ldp*(td*td-e2);e0=td-e1;eed=1.0_ldp-e0
endif
wdd=e1/td;wud=e0-wdd
ee=eed%v;dee=eed%d;wu=wud%v;dwu=wud%d;wd=wdd%v;dwd=wdd%d
end subroutine

! Node values (dual) exactly as tangent_map: extended extinction and emissivity.
subroutine nodes(eta,chi,scat,deta,dchi,dscat,ce,et,cc,ec)
real(kind=ldp),intent(in)::eta(:),chi(:),scat(:),deta(:),dchi(:),dscat(:)
type(dual),intent(out)::ce(nd_ext),et(nd_ext),cc(nd_ext,4),ec(nd_ext,4)
type(dual)::ch(size(eta)),em(size(eta)),es(size(eta)),lc(nd_ext),le(nd_ext),alpha,pow,z
real(kind=ldp)::lr
integer::i,nd
nd=size(eta)
do i=1,nd
 ch(i)=q(chi(i),dchi(i));em(i)=q(eta(i),deta(i));es(i)=q(scat(i),dscat(i))
enddo
ce(nd_add+1:nd_ext)=ch;et(nd_add+1:nd_ext)=em
if(nd_add>0)then
 lr=log(r_ext(nd_add+idmin)/r_ext(nd_add+idmax))
 pow=log(es(idmax)/es(idmin))/lr
 if(pow%v<2.0_ldp)pow=q(2.0_ldp,0.0_ldp)
 if(chi(idmin)<=scat(idmin).or.chi(idmax)<=scat(idmax))then
  do i=1,nd_add
   ce(i)=ch(idmin)*exp(pow*log(r_ext(nd_add+idmin)/r_ext(i)))
  enddo
 else
  alpha=log((ch(idmax)-es(idmax))/(ch(idmin)-es(idmin)))/lr
  if(alpha%v<2.0_ldp)alpha=q(2.0_ldp,0.0_ldp)
  do i=1,nd_add
   z=q(log(r_ext(nd_add+idmin)/r_ext(i)),0.0_ldp)
   ce(i)=(ch(idmin)-es(idmin))*exp(alpha*z)+es(idmin)*exp(pow*z)
  enddo
 endif
 alpha=log(em(idmax)/em(idmin))/lr
 if(alpha%v<3.5_ldp)alpha=q(3.5_ldp,0.0_ldp)
 do i=1,nd_add
  et(i)=em(idmin)*exp(alpha*log(r_ext(nd_add+idmin)/r_ext(i)))
  if(et(i)%v<=1.0e-280_ldp)et(i)=q(1.0e-280_ldp,0.0_ldp)
 enddo
endif
if(ray_points_inserted)then
 do i=1,nd_ext
  lc(i)=log(ce(i));le(i)=log(et(i))
 enddo
 call mon_coeff(lc,log_r_ext,cc);call mon_coeff(le,log_r_ext,ec)
else
 cc=q(0.0_ldp,0.0_ldp);ec=q(0.0_ldp,0.0_ldp)
endif
end subroutine

! Primal pass at the converged state; eta is the TOTAL emissivity eta+sigma*J.
subroutine lin_setup(eta,chi,scat,b,dbb,dlnu,nc,boundary,prim_moms)
real(kind=ldp),intent(in)::eta(:),chi(:),scat(:),b,dbb,dlnu
integer,intent(in)::nc
character(len=*),intent(in)::boundary
real(kind=ldp),intent(out),optional::prim_moms(4,size(eta))
type(dual)::ce(nd_ext),et(nd_ext),cc(nd_ext,4),ec(nd_ext,4)
real(kind=ldp),allocatable::zero(:)
real(kind=ldp)::nuon,z,ds,t,gradv
integer::nd,np,ir,j,n,i,k
nd=size(eta);np=np_sav
if(nd<3.or.nd/=nd_sav)error stop 'lin_setup: native geometry mismatch'
if(boundary/='DIFFUSION'.and.boundary/='ZERO_FLUX')error stop 'lin_setup: unsupported inner boundary'
if(nd_ext/=nd+nd_add)error stop 'lin_setup: inconsistent native grid'
allocate(zero(nd));zero=0
asig_ok=.false.
call nodes(eta,chi,scat,zero,zero,zero,ce,et,cc,ec)
nuon=0;if(dlnu>0)nuon=1.0_ldp/dlnu
s_nuon=nuon;s_b=b;s_dbb=dbb;s_chind=chi(nd);s_nc=nc;s_bnd=boundary;s_nd=nd
if(allocated(s_eta))deallocate(s_eta,s_chi,s_scat)
allocate(s_eta(nd),s_chi(nd),s_scat(nd));s_eta=eta;s_chi=chi;s_scat=scat
j=np-1;if(nd_add>0)j=np
s_nray=j
if(.not.allocated(rl))allocate(rl(size(ray)))
if(present(prim_moms))prim_moms=0
gradv=dbb/chi(nd)
!$omp parallel do schedule(dynamic,4) default(shared) private(ir,n,i,k,z,ds,t) &
!$omp if(lin_setup_omp.and..not.omp_in_parallel())
do ir=1,j
 n=ray(ir)%nz
 if(rl(ir)%n/=n)then
  if(allocated(rl(ir)%kk))deallocate(rl(ir)%kk,rl(ir)%z,rl(ir)%crv,rl(ir)%erv,rl(ir)%ctm,rl(ir)%ctp,&
     rl(ir)%srcm,rl(ir)%srcp,rl(ir)%imv,rl(ir)%ipv,rl(ir)%eem,rl(ir)%wum,rl(ir)%wdm,rl(ir)%deem,&
     rl(ir)%dwum,rl(ir)%dwdm,rl(ir)%dsm,rl(ir)%eep,rl(ir)%wup,rl(ir)%wdp,rl(ir)%deep,rl(ir)%dwup,&
     rl(ir)%dwdp,rl(ir)%dsp)
  allocate(rl(ir)%kk(n),rl(ir)%z(n),rl(ir)%crv(n),rl(ir)%erv(n),rl(ir)%ctm(n),rl(ir)%ctp(n),&
     rl(ir)%srcm(n),rl(ir)%srcp(n),rl(ir)%imv(n),rl(ir)%ipv(n),rl(ir)%eem(n),rl(ir)%wum(n),rl(ir)%wdm(n),&
     rl(ir)%deem(n),rl(ir)%dwum(n),rl(ir)%dwdm(n),rl(ir)%dsm(n),rl(ir)%eep(n),rl(ir)%wup(n),rl(ir)%wdp(n),&
     rl(ir)%deep(n),rl(ir)%dwup(n),rl(ir)%dwdp(n),rl(ir)%dsp(n))
  rl(ir)%n=n
 endif
 associate(a=>rl(ir))
 if(ray_points_inserted)then
  k=1
  do i=1,n
   do while(k<nd_ext-1.and.ray(ir)%r_ray(i)<r_ext(k+1))
    k=k+1
   enddo
   if(ray(ir)%r_ray(i)==r_ext(k))then
    a%kk(i)=-k;a%z(i)=0;a%crv(i)=ce(k)%v;a%erv(i)=et(k)%v
   else if(ray(ir)%r_ray(i)==r_ext(k+1))then
    a%kk(i)=-(k+1);a%z(i)=0;a%crv(i)=ce(k+1)%v;a%erv(i)=et(k+1)%v
   else
    z=log(ray(ir)%r_ray(i)/r_ext(k))
    a%kk(i)=k;a%z(i)=z
    a%crv(i)=exp(((cc(k,1)%v*z+cc(k,2)%v)*z+cc(k,3)%v)*z+cc(k,4)%v)
    a%erv(i)=exp(((ec(k,1)%v*z+ec(k,2)%v)*z+ec(k,3)%v)*z+ec(k,4)%v)
   endif
  enddo
 else
  do i=1,n
   a%kk(i)=-i;a%z(i)=0;a%crv(i)=ce(i)%v;a%erv(i)=et(i)%v
  enddo
 endif
 ! Minus direction (outer boundary inward).
 do i=1,n
  a%ctm(i)=a%crv(i)+(nuon+3.0_ldp)*ray(ir)%b_m(i)
  a%srcm(i)=(a%erv(i)+nuon*ray(ir)%b_m(i)*ray(ir)%i_m_prev(i))/a%ctm(i)
 enddo
 a%imv(1)=0
 do i=1,n-1
  ds=0.5_ldp*(ray(ir)%s_m(i+1)-ray(ir)%s_m(i))
  a%dsm(i)=ds
  t=ds*(a%ctm(i)+a%ctm(i+1))
  call dweights(t,a%eem(i),a%wum(i),a%wdm(i),a%deem(i),a%dwum(i),a%dwdm(i))
  a%imv(i+1)=a%imv(i)*a%eem(i)+a%wum(i)*a%srcm(i)+a%wdm(i)*a%srcm(i+1)
 enddo
 ! Plus direction (inner boundary outward).
 do i=1,n
  a%ctp(i)=a%crv(i)+(nuon+3.0_ldp)*ray(ir)%b_p(i)
  a%srcp(i)=(a%erv(i)+nuon*ray(ir)%b_p(i)*ray(ir)%i_p_prev(i))/a%ctp(i)
 enddo
 if(ir>nc.or.boundary=='ZERO_FLUX')then
  a%ipv(n)=a%imv(n)
 else
  a%ipv(n)=b+ray(ir)%mu_p(n)*gradv
 endif
 do i=n-1,1,-1
  ds=0.5_ldp*(ray(ir)%s_p(i)-ray(ir)%s_p(i+1))
  a%dsp(i)=ds
  t=ds*(a%ctp(i)+a%ctp(i+1))
  call dweights(t,a%eep(i),a%wup(i),a%wdp(i),a%deep(i),a%dwup(i),a%dwdp(i))
  a%ipv(i)=a%ipv(i+1)*a%eep(i)+a%wup(i)*a%srcp(i+1)+a%wdp(i)*a%srcp(i)
 enddo
 end associate
enddo
!$omp end parallel do
if(present(prim_moms))then
 do ir=1,j
  associate(a=>rl(ir))
  do i=1,min(nd,np-ir+1)
   k=ray(ir)%lnk(i)
   prim_moms(1,i)=prim_moms(1,i)+a%ipv(k)*jqw_p(i,ir)+a%imv(k)*jqw_m(i,ir)
   prim_moms(2,i)=prim_moms(2,i)+a%ipv(k)*hqw_p(i,ir)+a%imv(k)*hqw_m(i,ir)
   prim_moms(3,i)=prim_moms(3,i)+a%ipv(k)*kqw_p(i,ir)+a%imv(k)*kqw_m(i,ir)
   prim_moms(4,i)=prim_moms(4,i)+a%ipv(k)*nqw_p(i,ir)+a%imv(k)*nqw_m(i,ir)
  enddo
  end associate
 enddo
endif
lin_ready=.true.
end subroutine

! Directional derivative of the formal map about the lin_setup state.
! deta is the change of TOTAL emissivity; prev_p/m are previous-frequency
! ray derivatives (zero for a pure current-frequency emissivity product).
subroutine lin_apply(deta,dchi,dscat,db,ddbb,prev_p,prev_m,dmoms,out_p,out_m,jonly)
real(kind=ldp),intent(in)::deta(:),dchi(:),dscat(:),db,ddbb,prev_p(:,:),prev_m(:,:)
real(kind=ldp),intent(out)::dmoms(4,size(deta)),out_p(:,:),out_m(:,:)
logical,intent(in),optional::jonly
type(dual)::ce(nd_ext),et(nd_ext),cc(nd_ext,4),ec(nd_ext,4)
real(kind=ldp)::dcr(maxval(rl(1:s_nray)%n)),der(size(dcr)),dct(size(dcr)),dsrc(size(dcr))
real(kind=ldp)::dim_(size(dcr)),dip(size(dcr))
real(kind=ldp)::ddt,dgrad,nuon
integer::nd,np,ir,n,i,k,kk
logical::jo
if(.not.lin_ready)error stop 'lin_apply before lin_setup'
jo=.false.;if(present(jonly))jo=jonly
nd=size(deta);np=np_sav
if(nd/=s_nd)error stop 'lin_apply: size mismatch'
call nodes(s_eta,s_chi,s_scat,deta,dchi,dscat,ce,et,cc,ec)
nuon=s_nuon
dmoms=0
dgrad=ddbb/s_chind-s_dbb*dchi(nd)/(s_chind*s_chind)
! Rays are independent; moments are reduced. Inside an enclosing parallel
! region (several directions) this runs serially (no nested parallelism).
!$omp parallel do schedule(dynamic,4) default(shared) &
!$omp private(ir,n,i,k,kk,ddt,dcr,der,dct,dsrc,dim_,dip) reduction(+:dmoms) if(.not.omp_in_parallel())
do ir=1,s_nray
 associate(a=>rl(ir))
 n=a%n
 do i=1,n
  kk=a%kk(i)
  if(kk<0)then
   dcr(i)=ce(-kk)%d;der(i)=et(-kk)%d
  else
   dcr(i)=a%crv(i)*(((cc(kk,1)%d*a%z(i)+cc(kk,2)%d)*a%z(i)+cc(kk,3)%d)*a%z(i)+cc(kk,4)%d)
   der(i)=a%erv(i)*(((ec(kk,1)%d*a%z(i)+ec(kk,2)%d)*a%z(i)+ec(kk,3)%d)*a%z(i)+ec(kk,4)%d)
  endif
 enddo
 ! Minus direction.
 do i=1,n
  dsrc(i)=(der(i)+nuon*ray(ir)%b_m(i)*prev_m(i,ir)-a%srcm(i)*dcr(i))/a%ctm(i)
 enddo
 dim_(1)=0
 do i=1,n-1
  ddt=a%dsm(i)*(dcr(i)+dcr(i+1))
  dim_(i+1)=dim_(i)*a%eem(i)+a%imv(i)*a%deem(i)*ddt+a%dwum(i)*ddt*a%srcm(i)+a%wum(i)*dsrc(i)&
           +a%dwdm(i)*ddt*a%srcm(i+1)+a%wdm(i)*dsrc(i+1)
 enddo
 ! Plus direction.
 do i=1,n
  dsrc(i)=(der(i)+nuon*ray(ir)%b_p(i)*prev_p(i,ir)-a%srcp(i)*dcr(i))/a%ctp(i)
 enddo
 if(ir>s_nc.or.s_bnd=='ZERO_FLUX')then
  dip(n)=dim_(n)
 else
  dip(n)=db+ray(ir)%mu_p(n)*dgrad
 endif
 do i=n-1,1,-1
  ddt=a%dsp(i)*(dcr(i)+dcr(i+1))
  dip(i)=dip(i+1)*a%eep(i)+a%ipv(i+1)*a%deep(i)*ddt+a%dwup(i)*ddt*a%srcp(i+1)+a%wup(i)*dsrc(i+1)&
        +a%dwdp(i)*ddt*a%srcp(i)+a%wdp(i)*dsrc(i)
 enddo
 if(.not.jo)then
  out_p(1:n,ir)=dip(1:n);out_m(1:n,ir)=dim_(1:n)
 endif
 do i=1,min(nd,np-ir+1)
  k=ray(ir)%lnk(i)
  dmoms(1,i)=dmoms(1,i)+dip(k)*jqw_p(i,ir)+dim_(k)*jqw_m(i,ir)
  if(.not.jo)then
   dmoms(2,i)=dmoms(2,i)+dip(k)*hqw_p(i,ir)+dim_(k)*hqw_m(i,ir)
   dmoms(3,i)=dmoms(3,i)+dip(k)*kqw_p(i,ir)+dim_(k)*kqw_m(i,ir)
   dmoms(4,i)=dmoms(4,i)+dip(k)*nqw_p(i,ir)+dim_(k)*nqw_m(i,ir)
  endif
 enddo
 end associate
enddo
!$omp end parallel do
end subroutine

! Coherent-scattering tangent with the forcing from ONE full linear application
! and the implicit (I - A diag(sigma)) solve by right-preconditioned GMRES whose
! products are emissivity-only linear applications. The final state is the sum
! of the two linear applications (exact by linearity). eta is THERMAL here.
subroutine lin_scattering(scat,jbase,deta,dchi,dscat,db,ddbb,prev_p,prev_m,dmoms,out_p,out_m,&
                          closure_error,products,tol)
real(kind=ldp),intent(in)::scat(:),jbase(:),deta(:),dchi(:),dscat(:),db,ddbb,prev_p(:,:),prev_m(:,:)
real(kind=ldp),intent(out)::dmoms(4,size(scat)),out_p(:,:),out_m(:,:),closure_error
integer,intent(out)::products
real(kind=ldp),intent(in),optional::tol
real(kind=ldp)::m1(4,size(scat)),m2(4,size(scat)),zero(size(scat)),rhs(size(scat)),tl
real(kind=ldp),allocatable::p0(:,:),o2p(:,:),o2m(:,:)
integer::nd
nd=size(scat);zero=0
tl=1.0e-12_ldp;if(present(tol))tl=tol
allocate(p0(size(prev_p,1),size(prev_p,2)),o2p(size(prev_p,1),size(prev_p,2)),o2m(size(prev_p,1),size(prev_p,2)))
p0=0;o2p=0;o2m=0
call lin_apply(deta+dscat*jbase,dchi,dscat,db,ddbb,prev_p,prev_m,m1,out_p,out_m)
rhs=m1(1,:)
call gmres_scat(rhs,products)
call lin_apply(scat*rhs,zero,zero,0.0_ldp,0.0_ldp,p0,p0,m2,o2p,o2m)
dmoms=m1+m2
out_p=out_p+o2p;out_m=out_m+o2m
closure_error=maxval(abs(dmoms(1,:)-rhs))/max(maxval(abs(rhs)),tiny(1.0_ldp))
contains
subroutine gmres_scat(solution,count)
real(kind=ldp),intent(inout)::solution(nd)
integer,intent(out)::count
real(kind=ldp)::basis(nd,nd+1),h(nd+1,nd),w(nd),diag(nd),g(nd+1),cs(nd),sn(nd),y(nd),mm(4,nd)
real(kind=ldp)::beta,hn,den,tmp,proj,tolerance
integer::ii,jj,pass,nused
diag=1.0_ldp-scat*qs_lstar(1:nd)
where(diag<=0.0_ldp)diag=1.0_ldp
solution=solution/diag
beta=sqrt(sum(solution*solution))
count=0
if(beta==0.0_ldp)then
 solution=0;return
endif
basis=0;h=0;g=0;g(1)=beta;basis(:,1)=solution/beta
tolerance=tl*beta
nused=nd
do jj=1,nd
 if(asig_ok)then
  w=matmul(asig,basis(:,jj))
 else
  call lin_apply(scat*basis(:,jj),zero,zero,0.0_ldp,0.0_ldp,p0,p0,mm,o2p,o2m,.true.)
  w=mm(1,:)
 endif
 count=count+1
 w=(basis(:,jj)-w)/diag
 do pass=1,2
  do ii=1,jj
   proj=dot_product(basis(:,ii),w);h(ii,jj)=h(ii,jj)+proj
   w=w-proj*basis(:,ii)
  enddo
 enddo
 hn=sqrt(sum(w*w));h(jj+1,jj)=hn
 if(hn>tiny(1.0_ldp))basis(:,jj+1)=w/hn
 do ii=1,jj-1
  tmp=cs(ii)*h(ii,jj)+sn(ii)*h(ii+1,jj)
  h(ii+1,jj)=-sn(ii)*h(ii,jj)+cs(ii)*h(ii+1,jj);h(ii,jj)=tmp
 enddo
 den=sqrt(h(jj,jj)**2+h(jj+1,jj)**2)
 if(den==0.0_ldp)error stop 'scattering GMRES singular'
 cs(jj)=h(jj,jj)/den;sn(jj)=h(jj+1,jj)/den
 h(jj,jj)=den;h(jj+1,jj)=0
 g(jj+1)=-sn(jj)*g(jj);g(jj)=cs(jj)*g(jj)
 if(abs(g(jj+1))<=tolerance.or.hn<=tiny(1.0_ldp))then
  nused=jj;exit
 endif
enddo
y=0
do ii=nused,1,-1
 y(ii)=(g(ii)-dot_product(h(ii,ii+1:nused),y(ii+1:nused)))/h(ii,ii)
enddo
solution=matmul(basis(:,1:nused),y(1:nused))
end subroutine
end subroutine

! Build 30 (Report 22): m directional derivatives in ONE ray sweep. Exactly the statements of lin_apply
! with the direction index innermost (coefficients read once per point, recurrences vectorized across
! directions); moments are accumulated over rays in the order ir = 1..s_nray, i.e. the serial (direction-
! parallel mode) order of lin_apply, so each direction's result equals a serial lin_apply bit for bit.
! noprev: previous-frequency ray derivatives are zero (prev_p/prev_m not referenced).
! jonly: only J moments; out_p/out_m not referenced.
subroutine lin_apply_multi(m,deta,dchi,dscat,db,ddbb,noprev,prev_p,prev_m,dmoms,out_p,out_m,jonly)
integer,intent(in)::m
real(kind=ldp),intent(in)::deta(:,:),dchi(:,:),dscat(:,:),db(:),ddbb(:),prev_p(:,:,:),prev_m(:,:,:)
logical,intent(in)::noprev,jonly
real(kind=ldp),intent(out)::dmoms(:,:,:),out_p(:,:,:),out_m(:,:,:)
type(dual)::ce(nd_ext),et(nd_ext),cc(nd_ext,4),ec(nd_ext,4)
real(kind=ldp),allocatable::ced(:,:),etd(:,:),ccd(:,:,:),ecd(:,:,:),dcr(:,:),der(:,:),dsrc(:,:),dim_(:,:),dip(:,:),&
  dm(:,:,:)
real(kind=ldp)::ddt(m),dgrad(m),nuon,z,bm,bp
integer::nd,np,ir,n,i,k,kk,c,nz
if(.not.lin_ready)error stop 'lin_apply_multi before lin_setup'
nd=size(deta,1);np=np_sav
if(nd/=s_nd)error stop 'lin_apply_multi: size mismatch'
nz=maxval(rl(1:s_nray)%n)
allocate(ced(m,nd_ext),etd(m,nd_ext),ccd(m,nd_ext,4),ecd(m,nd_ext,4),dcr(m,nz),der(m,nz),dsrc(m,nz),dim_(m,nz),&
  dip(m,nz),dm(m,4,nd))
do c=1,m
 call nodes(s_eta,s_chi,s_scat,deta(:,c),dchi(:,c),dscat(:,c),ce,et,cc,ec)
 do i=1,nd_ext
  ced(c,i)=ce(i)%d;etd(c,i)=et(i)%d
  do k=1,4
   ccd(c,i,k)=cc(i,k)%d;ecd(c,i,k)=ec(i,k)%d
  enddo
 enddo
 dgrad(c)=ddbb(c)/s_chind-s_dbb*dchi(nd,c)/(s_chind*s_chind)
enddo
nuon=s_nuon
dm=0
do ir=1,s_nray
 associate(a=>rl(ir))
 n=a%n
 do i=1,n
  kk=a%kk(i)
  if(kk<0)then
   dcr(:,i)=ced(:,-kk);der(:,i)=etd(:,-kk)
  else
   z=a%z(i)
   dcr(:,i)=a%crv(i)*(((ccd(:,kk,1)*z+ccd(:,kk,2))*z+ccd(:,kk,3))*z+ccd(:,kk,4))
   der(:,i)=a%erv(i)*(((ecd(:,kk,1)*z+ecd(:,kk,2))*z+ecd(:,kk,3))*z+ecd(:,kk,4))
  endif
 enddo
 ! Minus direction.
 if(noprev)then
  do i=1,n
   dsrc(:,i)=(der(:,i)+nuon*ray(ir)%b_m(i)*0.0_ldp-a%srcm(i)*dcr(:,i))/a%ctm(i)
  enddo
 else
  do i=1,n
   bm=nuon*ray(ir)%b_m(i)
   dsrc(:,i)=(der(:,i)+bm*prev_m(:,i,ir)-a%srcm(i)*dcr(:,i))/a%ctm(i)
  enddo
 endif
 dim_(:,1)=0
 do i=1,n-1
  ddt=a%dsm(i)*(dcr(:,i)+dcr(:,i+1))
  dim_(:,i+1)=dim_(:,i)*a%eem(i)+a%imv(i)*a%deem(i)*ddt+a%dwum(i)*ddt*a%srcm(i)+a%wum(i)*dsrc(:,i)&
           +a%dwdm(i)*ddt*a%srcm(i+1)+a%wdm(i)*dsrc(:,i+1)
 enddo
 ! Plus direction.
 if(noprev)then
  do i=1,n
   dsrc(:,i)=(der(:,i)+nuon*ray(ir)%b_p(i)*0.0_ldp-a%srcp(i)*dcr(:,i))/a%ctp(i)
  enddo
 else
  do i=1,n
   bp=nuon*ray(ir)%b_p(i)
   dsrc(:,i)=(der(:,i)+bp*prev_p(:,i,ir)-a%srcp(i)*dcr(:,i))/a%ctp(i)
  enddo
 endif
 if(ir>s_nc.or.s_bnd=='ZERO_FLUX')then
  dip(:,n)=dim_(:,n)
 else
  dip(:,n)=db(1:m)+ray(ir)%mu_p(n)*dgrad
 endif
 do i=n-1,1,-1
  ddt=a%dsp(i)*(dcr(:,i)+dcr(:,i+1))
  dip(:,i)=dip(:,i+1)*a%eep(i)+a%ipv(i+1)*a%deep(i)*ddt+a%dwup(i)*ddt*a%srcp(i+1)+a%wup(i)*dsrc(:,i+1)&
        +a%dwdp(i)*ddt*a%srcp(i)+a%wdp(i)*dsrc(:,i)
 enddo
 if(.not.jonly)then
  do i=1,n
   out_p(:,i,ir)=dip(:,i);out_m(:,i,ir)=dim_(:,i)
  enddo
 endif
 do i=1,min(nd,np-ir+1)
  k=ray(ir)%lnk(i)
  dm(:,1,i)=dm(:,1,i)+dip(:,k)*jqw_p(i,ir)+dim_(:,k)*jqw_m(i,ir)
  if(.not.jonly)then
   dm(:,2,i)=dm(:,2,i)+dip(:,k)*hqw_p(i,ir)+dim_(:,k)*hqw_m(i,ir)
   dm(:,3,i)=dm(:,3,i)+dip(:,k)*kqw_p(i,ir)+dim_(:,k)*kqw_m(i,ir)
   dm(:,4,i)=dm(:,4,i)+dip(:,k)*nqw_p(i,ir)+dim_(:,k)*nqw_m(i,ir)
  endif
 enddo
 end associate
enddo
do c=1,m
 do i=1,nd
  dmoms(1:4,i,c)=dm(c,1:4,i)
 enddo
enddo
end subroutine

! Build 30 (Report 22): explicit scattering operator of the current frequency (see asig). Columns in blocks
! of up to 16 with lin_apply_multi (unit emissivity s_scat(j) at node j, no previous-frequency term, J only);
! blocks run in parallel in a top-level region (never nested). Called after lin_setup, before directions.
! Build 31: block size ceil(nd/(2*threads)) (at most 16), so the blocks fill every thread twice (build 30:
! 11 blocks of 16 for nd=167, i.e. 2 rounds on 8 threads with 3 busy in the second). Columns are computed
! by the same statements whatever the block size.
subroutine lin_em_build()
integer::mb
integer::nb,jb,j0,j1,m,j,nd
real(kind=ldp),allocatable::de(:,:),z0(:,:),dm(:,:,:),zz(:),pa(:,:,:),pb(:,:,:),oa(:,:,:),ob(:,:,:)
if(.not.lin_ready)error stop 'lin_em_build before lin_setup'
nd=s_nd
if(allocated(asig))then
 if(size(asig,1)/=nd)deallocate(asig)
endif
if(.not.allocated(asig))allocate(asig(nd,nd))
mb=16
if(.not.omp_in_parallel())mb=max(1,min(16,(nd+2*omp_get_max_threads()-1)/(2*omp_get_max_threads())))
nb=(nd+mb-1)/mb
!$omp parallel do schedule(dynamic) default(shared) private(jb,j0,j1,m,j,de,z0,dm,zz,pa,pb,oa,ob) &
!$omp if(.not.omp_in_parallel())
do jb=1,nb
 j0=(jb-1)*mb+1;j1=min(nd,jb*mb);m=j1-j0+1
 allocate(de(nd,m),z0(nd,m),dm(4,nd,m),zz(m),pa(1,1,1),pb(1,1,1),oa(1,1,1),ob(1,1,1))
 de=0;z0=0;zz=0;pa=0;pb=0
 do j=j0,j1
  de(j,j-j0+1)=s_scat(j)
 enddo
 call lin_apply_multi(m,de,z0,z0,zz,zz,.true.,pa,pb,dm,oa,ob,.true.)
 do j=j0,j1
  asig(:,j)=dm(1,:,j-j0+1)
 enddo
 deallocate(de,z0,dm,zz,pa,pb,oa,ob)
enddo
!$omp end parallel do
asig_ok=.true.
end subroutine

! Build 31 (Report 22): lin_scattering for m directions at once, with the explicit operator of this frequency
! (lin_em_build first). The forcing sweep and the final sweep are ONE lin_apply_multi each (the statements of
! the serial lin_apply with the direction index innermost, rays in order); each direction's GMRES is the asig
! branch of gmres_scat verbatim (gmres_asig). Same linear map, algorithm, preconditioner and tolerance as
! lin_scattering. Ray arrays are (m,nzmax,nray) here: prev_p/prev_m in, out_p/out_m out (zeroed first).
subroutine lin_scattering_multi(m,scat,jbase,deta,dchi,dscat,db,ddbb,prev_p,prev_m,dmoms,out_p,out_m,&
                                closure_error,products)
integer,intent(in)::m
real(kind=ldp),intent(in)::scat(:),jbase(:),deta(:,:),dchi(:,:),dscat(:,:),db(:),ddbb(:),prev_p(:,:,:),prev_m(:,:,:)
real(kind=ldp),intent(out)::dmoms(:,:,:),out_p(:,:,:),out_m(:,:,:),closure_error(:)
integer,intent(out)::products(:)
real(kind=ldp),allocatable::ef(:,:),es(:,:),z2(:,:),zv(:),m2(:,:,:),o2p(:,:,:),o2m(:,:,:),rhs(:,:),pz(:,:,:)
integer::nd,c
if(.not.asig_ok)error stop 'lin_scattering_multi without the explicit scattering operator'
nd=size(scat)
allocate(ef(nd,m),es(nd,m),z2(nd,m),zv(m),m2(4,nd,m),rhs(nd,m),pz(1,1,1),&
         o2p(size(out_p,1),size(out_p,2),size(out_p,3)),o2m(size(out_m,1),size(out_m,2),size(out_m,3)))
do c=1,m
 ef(:,c)=deta(:,c)+dscat(:,c)*jbase
enddo
out_p=0;out_m=0
call lin_apply_multi(m,ef,dchi,dscat,db,ddbb,.false.,prev_p,prev_m,dmoms,out_p,out_m,.false.)
do c=1,m
 rhs(:,c)=dmoms(1,:,c)
 call gmres_asig(scat,rhs(:,c),products(c))
 es(:,c)=scat*rhs(:,c)
enddo
z2=0;zv=0;pz=0;o2p=0;o2m=0
call lin_apply_multi(m,es,z2,z2,zv,zv,.true.,pz,pz,m2,o2p,o2m,.false.)
dmoms=dmoms+m2
out_p=out_p+o2p;out_m=out_m+o2m
do c=1,m
 closure_error(c)=maxval(abs(dmoms(1,:,c)-rhs(:,c)))/max(maxval(abs(rhs(:,c))),tiny(1.0_ldp))
enddo
end subroutine

! The asig branch of lin_scattering's gmres_scat, statement for statement (module-level, thread-safe).
subroutine gmres_asig(scat,solution,count)
real(kind=ldp),intent(in)::scat(:)
real(kind=ldp),intent(inout)::solution(:)
integer,intent(out)::count
real(kind=ldp)::basis(size(scat),size(scat)+1),h(size(scat)+1,size(scat)),w(size(scat)),diag(size(scat)),&
                g(size(scat)+1),cs(size(scat)),sn(size(scat)),y(size(scat))
real(kind=ldp)::beta,hn,den,tmp,proj,tolerance,tl
integer::ii,jj,pass,nused,nd
nd=size(scat)
tl=1.0e-12_ldp
diag=1.0_ldp-scat*qs_lstar(1:nd)
where(diag<=0.0_ldp)diag=1.0_ldp
solution=solution/diag
beta=sqrt(sum(solution*solution))
count=0
if(beta==0.0_ldp)then
 solution=0;return
endif
basis=0;h=0;g=0;g(1)=beta;basis(:,1)=solution/beta
tolerance=tl*beta
nused=nd
do jj=1,nd
 w=matmul(asig,basis(:,jj))
 count=count+1
 w=(basis(:,jj)-w)/diag
 do pass=1,2
  do ii=1,jj
   proj=dot_product(basis(:,ii),w);h(ii,jj)=h(ii,jj)+proj
   w=w-proj*basis(:,ii)
  enddo
 enddo
 hn=sqrt(sum(w*w));h(jj+1,jj)=hn
 if(hn>tiny(1.0_ldp))basis(:,jj+1)=w/hn
 do ii=1,jj-1
  tmp=cs(ii)*h(ii,jj)+sn(ii)*h(ii+1,jj)
  h(ii+1,jj)=-sn(ii)*h(ii,jj)+cs(ii)*h(ii+1,jj);h(ii,jj)=tmp
 enddo
 den=sqrt(h(jj,jj)**2+h(jj+1,jj)**2)
 if(den==0.0_ldp)error stop 'scattering GMRES singular'
 cs(jj)=h(jj,jj)/den;sn(jj)=h(jj+1,jj)/den
 h(jj,jj)=den;h(jj+1,jj)=0
 g(jj+1)=-sn(jj)*g(jj);g(jj)=cs(jj)*g(jj)
 if(abs(g(jj+1))<=tolerance.or.hn<=tiny(1.0_ldp))then
  nused=jj;exit
 endif
enddo
y=0
do ii=nused,1,-1
 y(ii)=(g(ii)-dot_product(h(ii,ii+1:nused),y(ii+1:nused)))/h(ii,ii)
enddo
solution=matmul(basis(:,1:nused),y(1:nused))
end subroutine
end module
