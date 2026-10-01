! Native finite-difference oracle for the analytic formal-transfer tangent.
! The oracle calls the unmodified build-17 CMF_FORMAL_REL_V4 + scattering ALI.
program formal_tangent_probe
use set_kind_module
use mod_ray_mom_store
use mod_space_grid_v2
use qs_dual
use qs_formal_tangent
implicit none
#ifndef NDV
#define NDV 20
#endif
integer,parameter::nd=NDV,nc=6,np=nd+nc,nf=24
real(kind=ldp)::r(nd),v(nd),sigma(nd),p(np),dop(nd),t(nd),ne(nd),nu(nf)
real(kind=ldp)::eta(nd),chi(nd),esec(nd),ta(nd),tc(nd),fedd(nd),rj(nd),iplus(np)
real(kind=ldp)::de(nd),dc(nd),ds(nd),db,ddbb,bnue,dbb,dlnu,hib,hob,eps,vinf,dopval,frontfac,kaphi,kaplo
real(kind=ldp)::m(4,nd),dm(4,nd),rsqh(nd),df(nd),rf(nd),cerr,maxcerr,maxreplay,maxfp,res,t1,t2,det
real(kind=ldp)::ng(nd,4),q1,q2,q3,a11,a12,a22,cc1,cc2,na,nb,tmp(nd)
real(kind=ldp)::base(5,nd,nf),deriv(5,nd,nf),jp(5,nd,nf),jm(5,nd,nf),fd(5,nd,nf),dropped(5,nd,nf),krylov_deriv(5,nd,nf),history_deriv(5,nd,nf)
real(kind=ldp),allocatable::prevp(:,:),prevm(:,:),outp(:,:),outm(:,:)
integer::i,k,ifr,it,nhist,mode,lu,nray,d0,kk,nproducts,maxproducts,totalproducts
logical::first,newf,conv,thick,krylov
character(len=12)::boundary
read(*,*)vinf,dopval,frontfac,kaphi,kaplo,thick,mode,eps,boundary
r(1)=30.0_ldp
do i=1,nd
 r(i)=30.0_ldp**(real(nd-i,ldp)/(nd-1))
 t1=1.0_ldp-0.999_ldp/r(i)
 v(i)=10.0_ldp+(vinf-10.0_ldp)*t1
 sigma(i)=(vinf-10.0_ldp)*0.999_ldp/r(i)/v(i)-1.0_ldp
 t(i)=2.0_ldp/sqrt(r(i))+0.3_ldp
 ne(i)=1.0_ldp/r(i)**2/(v(i)/10.0_ldp)
enddo
r(nd)=1.0_ldp;dop=dopval
call impar(p,r,r(nd),nc,nd,np)
do k=1,nf
 nu(k)=3.0_ldp*0.5_ldp**(real(k-1,ldp)/(nf-1))
enddo
qs_frm_ali=.true.;maxcerr=0;maxreplay=0;maxfp=0
d0=max(2,nint(0.85_ldp*nd));krylov=.false.;maxproducts=0;totalproducts=0
call sweep(0.0_ldp,base,.true.,.false.)
call sweep(eps,jp,.false.,.false.)
call sweep(-eps,jm,.false.,.false.)
fd=(jp-jm)/(2*eps)
call sweep(0.0_ldp,jp,.true.,.true.)
history_deriv=dropped
krylov=.true.
call sweep(0.0_ldp,jp,.true.,.false.)
krylov_deriv=deriv
krylov=.false.
call sweep(0.0_ldp,jp,.true.,.false.)
open(newunit=lu,file='response.txt',status='replace')
write(lu,'(A)')'# frequency depth base[J H K N r2H_face] tangent[5] finite_difference[5]'
do k=1,nf
 do i=1,nd
  write(lu,'(ES24.15,I5,15ES24.15)')nu(k),i,base(:,i,k),deriv(:,i,k),fd(:,i,k)
 enddo
enddo
close(lu)
do i=1,5
 t1=maxval(abs(fd(i,:,:)))
 if(t1==0)t1=1.0_ldp ! zero-flux boundary direction is exactly null
 write(*,'(A,I1,ES24.15)')'RESULT error_moment_',i,maxval(abs(deriv(i,:,:)-fd(i,:,:)))/t1
 write(*,'(A,I1,ES24.15)')'RESULT krylov_vs_dense_',i,maxval(abs(krylov_deriv(i,:,:)-deriv(i,:,:)))/max(maxval(abs(deriv(i,:,:))),tiny(1.0_ldp))
 write(*,'(A,I1,ES24.15)')'RESULT dropped_history_error_',i,maxval(abs(history_deriv(i,:,:)-fd(i,:,:)))/t1
enddo
write(*,'(A,ES24.15)')'RESULT scattering_tangent_closure ',maxcerr
write(*,'(A,ES24.15)')'RESULT primal_replay_error ',maxreplay
write(*,'(A,ES24.15)')'RESULT scattering_fixed_point_error ',maxfp
write(*,'(A,I8)')'RESULT max_ray_points ',nray
write(*,'(A,I8)')'RESULT max_scattering_krylov_products ',maxproducts
write(*,'(A,I8)')'RESULT total_scattering_krylov_products ',totalproducts
write(*,'(A,L1)')'RESULT ray_points_inserted ',ray_points_inserted
contains
subroutine model(f,ff)
real(kind=ldp),intent(in)::f,ff
real(kind=ldp)::kap,x,b1,profile
integer::j
kap=kaplo;if(f>2.2_ldp)kap=kaphi
do j=1,nd
 chi(j)=kap*ne(j)**1.5_ldp+2.0_ldp*ne(j)
 if(r(j)<6.0_ldp)chi(j)=chi(j)+frontfac*kap*ne(j)**1.5_ldp
 esec(j)=2.0_ldp*ne(j)
 eta(j)=(chi(j)-esec(j))*f**3/(exp(4.8_ldp*f/t(j))-1)
enddo
bnue=f**3/(exp(4.8_ldp*f/t(nd))-1)
b1=f**3/(exp(4.8_ldp*f/t(nd-1))-1)
dbb=(bnue-b1)/(r(nd-1)-r(nd))
de=0;dc=0;ds=0;db=0;ddbb=0
! Local blue-only source/opacity directions explicitly exercise frequency memory.
if(mode==1.and.ifr==9)de(d0)=eta(d0)
if(mode==2.and.ifr==9)dc(d0)=chi(d0)-esec(d0)
if(mode==3)then
 ds=esec*0.3_ldp;dc=ds
endif
if(mode==4)then
 profile=exp(-((real(ifr,ldp)-9.0_ldp)/3.0_ldp)**2)
 do j=1,nd
  de(j)=eta(j)*profile*(0.5_ldp+sin(real(j,ldp)))
  dc(j)=(chi(j)-esec(j))*profile*(0.5_ldp+cos(real(j,ldp)))
 enddo
endif
if(mode==5)then
 db=bnue;ddbb=0.7_ldp*dbb
endif
if(mode==6)then
 de=eta*0.5_ldp;dc=(chi-esec)*0.5_ldp
endif
eta=eta+ff*de;chi=chi+ff*dc;esec=esec+ff*ds
bnue=bnue+ff*db;dbb=dbb+ff*ddbb
end subroutine

subroutine sweep(ff,jout,with_tangent,drop)
real(kind=ldp),intent(in)::ff
real(kind=ldp),intent(out)::jout(5,nd,nf)
logical,intent(in)::with_tangent,drop
rj=0
if(allocated(prevp))deallocate(prevp,prevm,outp,outm)
do ifr=1,nf
 first=(ifr==1);dlnu=0;if(.not.first)dlnu=log(nu(ifr-1)/nu(ifr))
 call model(nu(ifr),ff)
 newf=.true.;nhist=0;conv=.false.;it=0
 do while(.not.conv)
  ta=eta+esec*rj
  call cmf_formal_rel_v4(ta,chi,esec,v,sigma,r,p,tc,fedd,hib,hob,iplus,nu(ifr),dlnu,bnue,dbb,&
                       trim(boundary),thick,dop,0.5_ldp,0.0_ldp,'LOGMON',first,newf,nc,np,nd)
  ! INITIALIZE=T resets the first-frequency blue rays, which are all zero.
  newf=.false.;it=it+1
  res=maxval(abs(tc-rj)/max(tc,rj,tiny(1.0_ldp)))
  if(res<1.0e-13_ldp)then
   rj=tc;conv=.true.
  else
   if(it>=5000)error stop 'native scattering did not converge'
   do i=1,nd
    t1=esec(i)*qs_lstar(i);t2=tc(i)
    if(t1>0.and.t1<1)then
     t2=(tc(i)-t1*rj(i))/(1-t1)
     if(t2<0)t2=tc(i)
    endif
    rj(i)=t2
   enddo
   nhist=nhist+1;ng(:,1:3)=ng(:,2:4);ng(:,4)=rj
   if(nhist>=4.and.it>=5)then
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
 maxfp=max(maxfp,res)
 call get_rsqh_rel(rsqh,r,v,nu(ifr),nd)
 jout(1,:,ifr)=jnu_store;jout(2,:,ifr)=hnu_store
 jout(3,:,ifr)=knu_store;jout(4,:,ifr)=nnu_store
 jout(5,:,ifr)=0;jout(5,1:nd-1,ifr)=rsqh(1:nd-1)
 if(with_tangent)then
  if(first)then
   nray=0
   do kk=1,np
    nray=max(nray,ray(kk)%nz)
   enddo
   allocate(prevp(nray,np),prevm(nray,np),outp(nray,np),outm(nray,np))
   prevp=0;prevm=0
  endif
  if(drop)then
   prevp=0;prevm=0
  endif
  call tangent_scattering(eta,chi,esec,rj,de,dc,ds,bnue,dbb,db,ddbb,dlnu,nc,trim(boundary),&
                          prevp,prevm,m,dm,outp,outm,cerr,krylov,nproducts)
  if(krylov)then
   maxproducts=max(maxproducts,nproducts);totalproducts=totalproducts+nproducts
  endif
  maxcerr=max(maxcerr,cerr)
  do kk=1,4
   maxreplay=max(maxreplay,maxval(abs(m(kk,:)-jout(kk,:,ifr)))/maxval(abs(jout(kk,:,ifr))))
  enddo
  rf=0;df=0
  call tangent_rsqh(r,m(2,:),dm(2,:),rf(1:nd-1),df(1:nd-1))
  maxreplay=max(maxreplay,maxval(abs(rf(1:nd-1)-rsqh(1:nd-1)))/maxval(abs(rsqh(1:nd-1))))
  if(drop)then
   dropped(1:4,:,ifr)=dm;dropped(5,:,ifr)=df
  else
   deriv(1:4,:,ifr)=dm;deriv(5,:,ifr)=df
  endif
  prevp=outp;prevm=outm
 endif
enddo
end subroutine
end program
