! SPDX-License-Identifier: GPL-3.0-only
! Differentiated adaptation of CMFGEN formal characteristics and Steffen
! interpolation (DLM/DJH); upstream GPLv3 provenance: backend README.
! Analytic directional derivative of build-17 QS_FRM_ALI on its native ray grid.
! Component interface only: does not replace VAR_JREL_V5 or certify a wind.
! Fixed R,V,ray geometry and frequency grid; DIFFUSION/ZERO_FLUX boundaries.
! Carries the full blue-side ray derivative, including inserted-point LOGMON
! interpolation and thick outer extrapolation. At limiter ties the map is
! piecewise differentiable; the selected primal branch is differentiated.
module qs_dual
use set_kind_module
implicit none
private
public :: dual,q,operator(+),operator(-),operator(*),operator(/),exp,log,abs,dmin
 type dual
 real(kind=ldp)::v=0,d=0
 end type
interface operator(+)
 module procedure add_dd,add_dr,add_rd
end interface
interface operator(-)
 module procedure sub_dd,sub_dr,sub_rd
 module procedure neg
end interface
interface operator(*)
 module procedure mul_dd,mul_dr,mul_rd
end interface
interface operator(/)
 module procedure div_dd,div_dr,div_rd
end interface
interface exp
 module procedure dexp
end interface
interface log
 module procedure dlog
end interface
interface abs
 module procedure dabs
end interface
contains
pure function q(v,d) result(z)
real(kind=ldp),intent(in)::v,d
type(dual)::z
z%v=v;z%d=d
end function
pure function add_dd(a,b) result(z)
type(dual),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a%v+b%v
z%d=a%d+b%d
end function
pure function add_dr(a,b) result(z)
type(dual),intent(in)::a
real(kind=ldp),intent(in)::b
type(dual)::z
z%v=a%v+b
z%d=a%d+0.0_ldp
end function
pure function add_rd(a,b) result(z)
real(kind=ldp),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a+b%v
z%d=0.0_ldp+b%d
end function
pure function sub_dd(a,b) result(z)
type(dual),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a%v-b%v
z%d=a%d-b%d
end function
pure function sub_dr(a,b) result(z)
type(dual),intent(in)::a
real(kind=ldp),intent(in)::b
type(dual)::z
z%v=a%v-b
z%d=a%d-0.0_ldp
end function
pure function sub_rd(a,b) result(z)
real(kind=ldp),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a-b%v
z%d=0.0_ldp-b%d
end function
pure function mul_dd(a,b) result(z)
type(dual),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a%v*b%v
z%d=a%d*b%v+a%v*b%d
end function
pure function mul_dr(a,b) result(z)
type(dual),intent(in)::a
real(kind=ldp),intent(in)::b
type(dual)::z
z%v=a%v*b
z%d=a%d*b+a%v*0.0_ldp
end function
pure function mul_rd(a,b) result(z)
real(kind=ldp),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a*b%v
z%d=0.0_ldp*b%v+a*b%d
end function
pure function div_dd(a,b) result(z)
type(dual),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a%v/b%v
z%d=(a%d-z%v*b%d)/b%v
end function
pure function div_dr(a,b) result(z)
type(dual),intent(in)::a
real(kind=ldp),intent(in)::b
type(dual)::z
z%v=a%v/b
z%d=(a%d-z%v*0.0_ldp)/b
end function
pure function div_rd(a,b) result(z)
real(kind=ldp),intent(in)::a
type(dual),intent(in)::b
type(dual)::z
z%v=a/b%v
z%d=(0.0_ldp-z%v*b%d)/b%v
end function
pure function neg(a) result(z)
type(dual),intent(in)::a
type(dual)::z
z=q(-a%v,-a%d)
end function
pure function dexp(a) result(z)
type(dual),intent(in)::a
type(dual)::z
z%v=exp(a%v);z%d=z%v*a%d
end function
pure function dlog(a) result(z)
type(dual),intent(in)::a
type(dual)::z
z=q(log(a%v),a%d/a%v)
end function
pure function dabs(a) result(z)
type(dual),intent(in)::a
type(dual)::z
z=q(abs(a%v),sign(1.0_ldp,a%v)*a%d)
end function
pure function dmin(a,b) result(z)
type(dual),intent(in)::a,b
type(dual)::z
z=a
if(b%v<a%v)z=b
end function
end module

module qs_formal_tangent
use set_kind_module
use qs_dual
use mod_space_grid_v2
use ext_rel_grid_v4, only: r_ext,log_r_ext,nd_ext,nd_add,idmin,idmax
use mod_ray_mom_store, only: qs_lstar
implicit none
private
public :: tangent_map,tangent_scattering,tangent_rsqh,mon_coeff
contains
! Same Steffen branches as MON_INT_FUNS_V2/MON_INTERP, differentiated.
subroutine mon_coeff(y,x,c)
type(dual),intent(in)::y(:)
real(kind=ldp),intent(in)::x(:)
type(dual),intent(out)::c(size(y),4)
type(dual)::s(size(y)),d(size(y)),lim
real(kind=ldp)::h(size(y)),sg
integer::i,n
n=size(y)
do i=1,n-1
 h(i)=x(i+1)-x(i);s(i)=(y(i+1)-y(i))/h(i)
enddo
d(1)=s(1)+(s(1)-s(2))*h(1)/(h(1)+h(2))
do i=2,n-1
 d(i)=(s(i-1)*h(i)+s(i)*h(i-1))/(h(i-1)+h(i))
enddo
d(n)=s(n-1)+(s(n-1)-s(n-2))*h(n-1)/(h(n-2)+h(n-1))
sg=sign(1.0_ldp,s(1)%v)+sign(1.0_ldp,d(1)%v)
d(1)=sg*dmin(abs(s(1)),0.5_ldp*abs(d(1)))
do i=2,n-1
 sg=sign(1.0_ldp,s(i-1)%v)+sign(1.0_ldp,s(i)%v)
 lim=dmin(abs(s(i-1)),abs(s(i)))
 d(i)=sg*dmin(lim,0.5_ldp*abs(d(i)))
enddo
sg=sign(1.0_ldp,s(n-1)%v)+sign(1.0_ldp,d(n)%v)
d(n)=sg*dmin(abs(s(n-1)),0.5_ldp*abs(d(n)))
do i=1,n-1
 c(i,1)=(d(i)+d(i+1)-2.0_ldp*s(i))/h(i)/h(i)
 c(i,2)=(3.0_ldp*s(i)-2.0_ldp*d(i)-d(i+1))/h(i)
 c(i,3)=d(i);c(i,4)=y(i)
enddo
c(n,1:2)=q(0.0_ldp,0.0_ldp);c(n,3)=d(n);c(n,4)=y(n)
end subroutine

subroutine weights(t,ee,wu,wd)
type(dual),intent(in)::t
type(dual),intent(out)::ee,wu,wd
type(dual)::e0,e1,e2
! Differentiate the actual small-t polynomial, including its truncation.
if(t%v>0.01_ldp)then
 ee=exp(-t);e0=1.0_ldp-ee;e1=t-e0
else
 e2=t*t*t*(1.0_ldp-0.25_ldp*t*(1.0_ldp-0.2_ldp*t*(1.0_ldp-t/6.0_ldp*&
    (1.0_ldp-t/7.0_ldp))))/3.0_ldp
 e1=0.5_ldp*(t*t-e2);e0=t-e1;ee=1.0_ldp-e0
endif
wd=e1/t;wu=e0-wd
end subroutine

! Tangent of native MON_INTERP on r^2 H: TVX must use these face values,
! rather than simply averaging H or recycling moment-method face factors.
subroutine tangent_rsqh(r,h,dh,fh,dfh)
real(kind=ldp),intent(in)::r(:),h(:),dh(:)
real(kind=ldp),intent(out)::fh(size(r)-1),dfh(size(r)-1)
type(dual)::y(size(r)),c(size(r),4),z
real(kind=ldp)::x
integer::j,nd
nd=size(r)
if(size(h)/=nd.or.size(dh)/=nd)error stop 'tangent: face shape mismatch'
do j=1,nd
 y(j)=q(r(j)*r(j)*h(j),r(j)*r(j)*dh(j))
enddo
call mon_coeff(y,r,c)
fh=0;dfh=0
do j=1,nd-1
 x=0.5_ldp*(r(j+1)-r(j))
 z=((c(j,1)*x+c(j,2))*x+c(j,3))*x+c(j,4)
 fh(j)=z%v;dfh(j)=z%d
enddo
end subroutine

! eta is total emissivity at the converged scattering fixed point.
! prev_p/m contain derivatives at the PREVIOUS frequency, never current ALI
! iterates. Outputs out_p/m become that state only after scattering closure.
subroutine tangent_map(eta,chi,scat,deta,dchi,dscat,b,dbb,db,ddbb,dlnu,nc,boundary,&
                       prev_p,prev_m,moms,dmoms,out_p,out_m)
real(kind=ldp),intent(in)::eta(:),chi(:),scat(:),deta(:),dchi(:),dscat(:)
real(kind=ldp),intent(in)::b,dbb,db,ddbb,dlnu,prev_p(:,:),prev_m(:,:)
integer,intent(in)::nc
character(len=*),intent(in)::boundary
real(kind=ldp),intent(out)::moms(4,size(eta)),dmoms(4,size(eta)),out_p(:,:),out_m(:,:)
type(dual)::ce(nd_ext),et(nd_ext),es(size(eta)),ch(size(eta)),em(size(eta))
type(dual)::lc(nd_ext),le(nd_ext),cc(nd_ext,4),ec(nd_ext,4),alpha,pow,t1,t2
! The maximum ray size can exceed the original radial grid size.
type(dual),allocatable::cr(:),er(:),ct(:),src(:),dt(:),im(:),ip(:)
type(dual)::ee,wu,wd,planck,grad,m(4,size(eta))
real(kind=ldp)::nuon,lr,z,ds
integer::i,k,j,n,ir,np,nd
nd=size(eta);np=size(prev_p,2)
if(nd<3.or.nd/=nd_sav.or.np/=np_sav)error stop 'tangent: native geometry mismatch'
if(size(chi)/=nd.or.size(scat)/=nd.or.size(deta)/=nd.or.size(dchi)/=nd.or.&
   size(dscat)/=nd)error stop 'tangent: coefficient shape mismatch'
if(any(shape(prev_p)/=shape(prev_m)).or.any(shape(prev_p)/=shape(out_p)).or.&
   any(shape(prev_p)/=shape(out_m)))error stop 'tangent: ray state shape mismatch'
do i=1,np
 if(size(prev_p,1)<ray(i)%nz)error stop 'tangent: insufficient ray capacity'
enddo
if(boundary/='DIFFUSION'.and.boundary/='ZERO_FLUX')error stop 'tangent: unsupported inner boundary'
if(nd_ext/=nd+nd_add)error stop 'tangent: inconsistent native grid'
nuon=0;if(dlnu>0)nuon=1.0_ldp/dlnu
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
   z=log(r_ext(nd_add+idmin)/r_ext(i))
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
endif
m=q(0.0_ldp,0.0_ldp);out_p=0;out_m=0
planck=q(b,db);grad=q(dbb,ddbb)/ch(nd)
! A thick outer boundary integrates all rays; otherwise exclude tangent ray.
j=np-1;if(nd_add>0)j=np
do ir=1,j
 n=ray(ir)%nz
 allocate(cr(n),er(n),ct(n),src(n),dt(n),im(n),ip(n))
 if(ray_points_inserted)then
  k=1
  do i=1,n
   do while(k<nd_ext-1.and.ray(ir)%r_ray(i)<r_ext(k+1))
    k=k+1
   enddo
   if(ray(ir)%r_ray(i)==r_ext(k))then
    cr(i)=ce(k);er(i)=et(k)
   else if(ray(ir)%r_ray(i)==r_ext(k+1))then
    cr(i)=ce(k+1);er(i)=et(k+1)
   else
    z=log(ray(ir)%r_ray(i)/r_ext(k))
    cr(i)=exp(((cc(k,1)*z+cc(k,2))*z+cc(k,3))*z+cc(k,4))
    er(i)=exp(((ec(k,1)*z+ec(k,2))*z+ec(k,3))*z+ec(k,4))
   endif
  enddo
 else
  cr=ce(:n);er=et(:n)
 endif
 do i=1,n
  ct(i)=cr(i)+(nuon+3.0_ldp)*ray(ir)%b_m(i)
  src(i)=(er(i)+nuon*ray(ir)%b_m(i)*q(ray(ir)%i_m_prev(i),prev_m(i,ir)))/ct(i)
 enddo
 do i=1,n-1
  ds=0.5_ldp*(ray(ir)%s_m(i+1)-ray(ir)%s_m(i))
  dt(i)=ds*(ct(i)+ct(i+1))
 enddo
 im(1)=q(0.0_ldp,0.0_ldp)
 do i=2,n
  call weights(dt(i-1),ee,wu,wd)
  im(i)=im(i-1)*ee+wu*src(i-1)+wd*src(i)
 enddo
 if(ir>nc.or.boundary=='ZERO_FLUX')then
  ip(n)=im(n)
 else
  ip(n)=planck+ray(ir)%mu_p(n)*grad
 endif
 do i=1,n
  ct(i)=cr(i)+(nuon+3.0_ldp)*ray(ir)%b_p(i)
  src(i)=(er(i)+nuon*ray(ir)%b_p(i)*q(ray(ir)%i_p_prev(i),prev_p(i,ir)))/ct(i)
 enddo
 do i=1,n-1
  ds=0.5_ldp*(ray(ir)%s_p(i)-ray(ir)%s_p(i+1))
  dt(i)=ds*(ct(i)+ct(i+1))
 enddo
 do i=n-1,1,-1
  call weights(dt(i),ee,wu,wd)
  ip(i)=ip(i+1)*ee+wu*src(i+1)+wd*src(i)
 enddo
 do i=1,n
  out_p(i,ir)=ip(i)%d;out_m(i,ir)=im(i)%d
 enddo
 do i=1,min(nd,np-ir+1)
  k=ray(ir)%lnk(i)
  m(1,i)=m(1,i)+ip(k)*jqw_p(i,ir)+im(k)*jqw_m(i,ir)
  m(2,i)=m(2,i)+ip(k)*hqw_p(i,ir)+im(k)*hqw_m(i,ir)
  m(3,i)=m(3,i)+ip(k)*kqw_p(i,ir)+im(k)*kqw_m(i,ir)
  m(4,i)=m(4,i)+ip(k)*nqw_p(i,ir)+im(k)*nqw_m(i,ir)
 enddo
 deallocate(cr,er,ct,src,dt,im,ip)
enddo
do i=1,nd
 do k=1,4
  moms(k,i)=m(k,i)%v;dmoms(k,i)=m(k,i)%d
 enddo
enddo
end subroutine

! Exact tangent of the converged coherent-scattering map. Dense construction
! is a bounded-fixture reference, not a production preconditioner. No finite
! differences enter the product. dscat affects BOTH extinction extrapolation
! and emissivity; caller includes dscat in dchi if total chi changes with it.
subroutine tangent_scattering(eta,chi,scat,jbase,deta,dchi,dscat,b,dbb,db,ddbb,dlnu,nc,boundary,&
                              prev_p,prev_m,moms,dmoms,out_p,out_m,closure_error,use_krylov,products)
real(kind=ldp),intent(in)::eta(:),chi(:),scat(:),jbase(:),deta(:),dchi(:),dscat(:)
real(kind=ldp),intent(in)::b,dbb,db,ddbb,dlnu,prev_p(:,:),prev_m(:,:)
integer,intent(in)::nc
character(len=*),intent(in)::boundary
real(kind=ldp),intent(out)::moms(4,size(eta)),dmoms(4,size(eta)),out_p(:,:),out_m(:,:),closure_error
real(kind=ldp)::a(size(eta),size(eta)),mat(size(eta),size(eta)),rhs(size(eta)),seed(size(eta)),zero(size(eta))
real(kind=ldp)::p0(size(prev_p,1),size(prev_p,2)),m0(size(prev_m,1),size(prev_m,2))
integer::k,nd,info,piv(size(eta)),nprod
logical,intent(in),optional::use_krylov
integer,intent(out),optional::products
logical::krylov
nd=size(eta);zero=0;p0=0;m0=0
krylov=.false.;if(present(use_krylov))krylov=use_krylov
nprod=0
! eta passed to the map is eta_thermal + sigma*J at the fixed point.
call tangent_map(eta+scat*jbase,chi,scat,deta+dscat*jbase,dchi,dscat,b,dbb,db,ddbb,dlnu,nc,boundary,&
                 prev_p,prev_m,moms,dmoms,out_p,out_m)
rhs=dmoms(1,:)
if(krylov)then
 call scattering_gmres(rhs,nprod)
else
 do k=1,nd
  seed=0;seed(k)=1
  call tangent_map(eta+scat*jbase,chi,scat,seed,zero,zero,b,dbb,0.0_ldp,0.0_ldp,dlnu,nc,boundary,&
                   p0,m0,moms,dmoms,out_p,out_m)
  a(:,k)=dmoms(1,:)
 enddo
 do k=1,nd
  mat(:,k)=-a(:,k)*scat(k);mat(k,k)=mat(k,k)+1
 enddo
 call dgesv(nd,1,mat,nd,piv,rhs,nd,info)
 if(info/=0)error stop 'tangent scattering singular'
 nprod=nd
endif
call tangent_map(eta+scat*jbase,chi,scat,deta+dscat*jbase+scat*rhs,dchi,dscat,b,dbb,db,ddbb,&
                 dlnu,nc,boundary,prev_p,prev_m,moms,dmoms,out_p,out_m)
closure_error=maxval(abs(dmoms(1,:)-rhs))/max(maxval(abs(rhs)),tiny(1.0_ldp))
if(present(products))products=nprod
contains
! Signed scattering response solved without constructing the dense Lambda
! matrix. Local L* is used solely as a diagonal GMRES preconditioner.
! Each Arnoldi product is ONE analytic ray sweep, not a full model iteration.
subroutine scattering_gmres(solution,count)
real(kind=ldp),intent(inout)::solution(nd)
integer,intent(out)::count
real(kind=ldp)::basis(nd,nd+1),h(nd+1,nd),w(nd),diag(nd),g(nd+1),cs(nd),sn(nd),y(nd)
real(kind=ldp)::beta,hn,den,tmp,proj,tolerance
integer::ii,jj,pass,nused
! Do not clip the transfer or change its Jacobian if the approximate diagonal
! is unsuitable: use unity for that preconditioner entry.
diag=1.0_ldp-scat*qs_lstar(1:nd)
where(diag<=0.0_ldp)diag=1.0_ldp
solution=solution/diag
beta=sqrt(sum(solution*solution))
count=0
if(beta==0.0_ldp)then
 solution=0;return
endif
basis=0;h=0;g=0;g(1)=beta;basis(:,1)=solution/beta
tolerance=1.0e-12_ldp*beta
nused=nd
do jj=1,nd
 seed=scat*basis(:,jj)
 call tangent_map(eta+scat*jbase,chi,scat,seed,zero,zero,b,dbb,0.0_ldp,0.0_ldp,dlnu,nc,boundary,&
                  p0,m0,moms,dmoms,out_p,out_m)
 count=count+1
 w=(basis(:,jj)-dmoms(1,:))/diag
 ! Reorthogonalization keeps a cancellation-dominated system accurate.
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
end module
