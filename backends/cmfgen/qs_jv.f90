! SPDX-License-Identifier: GPL-3.0-only
! QStarS build-18 directional (matrix-free) Jacobian of the native CMFGEN
! rate/radiative-equilibrium residual. Adapts CMFGEN accumulation logic
! (EVALSE_QWVJ_V8, EVALSE_LOWT_V2, ZNET/JBAR line updates, STEQ_T, TWO_PHOT_RATE,
! GENERATE_FULL_MATRIX_V3); upstream GPLv3 provenance: backend README.
!
! For each requested direction s (fractional, dPOPS = POPS*s) this computes the
! RADIATION part of J.dx exactly:  sum_nu (dSTEQ/dJ_nu) dJ_nu, where dJ_nu is the
! analytic nonlocal tangent of the QS_FRM_ALI formal transfer (formal_tangent.f90)
! driven by the exact directional changes of total extinction, thermal emissivity,
! coherent scattering and the diffusion inner boundary, carrying the blue-side ray
! derivative through the frequency sweep. The direct (fixed-J) part comes from the
! native BA assembled with QS_ALI_LIN and QS_LIN_FAC=0 (dumped by QS_JV_BA_DUMP).
! Nothing here changes the primal state, the residual, or the native Newton step.
module qs_jv_state
use set_kind_module
implicit none
logical :: jv_active=.false.
logical :: jv_export=.false.
integer :: nv=0,jv_nt=0,jv_nd=0,nzmax=0,jv_np=0,jv_maxsim=0,jv_ntwo=0,jv_lu_exp=-1
real(kind=ldp),allocatable :: sdir(:,:,:),dpop(:,:,:),jv_pops(:,:)
real(kind=ldp),allocatable :: prevp(:,:,:),prevm(:,:,:)
real(kind=ldp),allocatable :: djrec(:,:),djphot(:,:),dsteq_t(:,:)
real(kind=ldp),allocatable :: dznet(:,:,:),ddown(:,:,:),dup2(:,:,:)
real(kind=ldp),allocatable :: dfull(:,:,:)
logical,allocatable :: lowt_ion(:)
type dse_t
 real(kind=ldp),allocatable :: steq(:,:,:)
end type
type(dse_t),allocatable :: dse(:)
real(kind=ldp) :: st_closure=0,st_replay=0,st_maxdj=0,st_lincheck=0,st_closure_fl=0
integer :: st_products=0,st_maxprod=0,st_nfreq=0,st_negfloor=0
logical :: jv_dual=.false.,jv_checklin=.false.
! Build 30 (Report 22): explicit per-frequency scattering operator for the GMRES of every direction
! (lin_em_build) when nv >= jv_asig_min; QSTARS_JV_NOASIG disables it, QSTARS_JV_ASIG_MIN sets the minimum.
logical :: jv_asig=.true.
integer :: jv_asig_min=4,st_asig=0
! Build 31 (Report 22): with the explicit operator, the directions run in blocks (one per thread); each block
! makes ONE forcing and ONE final ray sweep for all its directions (lin_scattering_multi). QSTARS_JV_NOBATCH:
! per-direction lin_scattering as in build 30. QSTARS_JV_NOSETUP_OMP: serial lin_setup ray loop.
logical :: jv_batch=.true.
integer :: st_batch=0
! Build 25: level-dissolution derivatives of the bound-free opacity/emissivity (occupation
! probabilities W in the LTE populations and the dissolved fraction below each edge), which
! the native VCHI/VETA hold constant. dis_l*: d ln(B,A,X)_LEV_DIS per depth and direction.
logical :: jv_dis=.true.,dis_valid=.false.,jv_vall=.false.
real(kind=ldp),allocatable :: dis_lb(:,:),dis_la(:,:),dis_lx(:,:),dchidis(:,:),detadis(:,:)
! d ln W(row,depth) = dis_cwb*dlnB + dis_cwa*dlnA + dis_cwx*dlnX for levels of dissolving ions
! (one-to-one super levels); used for the line-opacity ratio L_STAR_RATIO = W_up/W_low.
real(kind=ldp),allocatable :: dis_cwb(:,:),dis_cwa(:,:),dis_cwx(:,:)
real(kind=ldp) :: st_dis_wdiff=0,st_dis_rel=0
end module qs_jv_state

! Build 32 (Report 23): conjugate-frequency J of the two-photon USE_RAD method from a QStarS stored-J file
! (QSTARS_TWO_PHOT_J = a QSTARS_RJ.bin on this model's radius and frequency grid). CMFGEN's USE_RAD method takes
! J(nu_B), nu_B = nu_0 - nu, from the lagged J of the previous iteration in EDDFACTOR (GET_J_FOR_TWO_PHOT). A
! single-iteration pass (COMP_F=T, no EDDFACTOR) cannot, and COMP_OPAC then silently replaces USE_RAD by
! OLD_DEFAULT, whose partner-photon occupation is exp(-h nu_B/kT_e). With the variable set, USE_RAD is kept and
! J(nu_B) is interpolated from the file exactly as GET_J_FOR_TWO_PHOT interpolates EDDFACTOR. The file is fixed
! (lagged) data: PHOT_OC_TWO does not depend on the state of the pass. Unset: nothing here is used.
module qs_two_phot_j_mod
use set_kind_module
implicit none
logical :: tpj_first=.true.,tpj_on=.false.,tpj_ready=.false.
integer :: tpj_ncf=0,tpj_nd=0,tpj_plen=0
integer(kind=8) :: tpj_ncall=0
character(len=512) :: tpj_path=' '
real(kind=ldp),allocatable :: tpj_fl(:),tpj_rj(:,:)
end module qs_two_phot_j_mod

! Read directions; allocate accumulators. Called once before the frequency loop.
subroutine qs_jv_begin(pops,nt,nd,compute_ba,lambda_it)
use set_kind_module
use qs_jv_state
use qs_formal_tangent_lin, only: lin_setup_omp
use mod_cmfgen, only: atm,num_ions
use steq_data_mod, only: se
use line_mod, only: znet_sim
use two_phot_mod, only: n_two
implicit none
integer,intent(in)::nt,nd
real(kind=ldp),intent(in)::pops(nt,nd)
logical,intent(in)::compute_ba,lambda_it
character(len=512)::path
integer::plen,ios,lu,hdr(3),id,k,nsim
integer(kind=8)::fsize
jv_active=.false.
call get_environment_variable('QSTARS_JV_DIRS',path,length=plen,status=ios)
if(ios==-1)stop 'QSTARS_JV_DIRS: path longer than 512 characters'
if(ios/=0.or.plen==0)return
if(.not.compute_ba.or.lambda_it)then
 open(newunit=lu,file='QSTARS_JV_FAILURE',status='unknown',position='append')
 write(lu,*)'QSTARS_JV_DIRS set but BA not computed (or lambda iteration); no product'
 close(lu)
 stop 'QSTARS_JV requires COMPUTE_BA and a non-lambda iteration'
endif
inquire(file=path(1:plen),size=fsize)
open(newunit=lu,file=path(1:plen),status='old',access='stream',form='unformatted',iostat=ios)
if(ios/=0)stop 'QSTARS_JV: cannot open direction file'
read(lu,iostat=ios)hdr
if(ios/=0.or.hdr(2)/=nt.or.hdr(3)/=nd.or.hdr(1)<1.or.&
   fsize/=12_8+8_8*int(hdr(1),8)*int(nt,8)*int(nd,8))then
 close(lu)
 open(newunit=lu,file='QSTARS_JV_FAILURE',status='unknown',position='append')
 write(lu,*)'Bad direction header',hdr,' expected NT,ND',nt,nd,' size',fsize
 close(lu)
 stop 'QSTARS_JV: bad direction file'
endif
nv=hdr(1);jv_nt=nt;jv_nd=nd
if(allocated(sdir))deallocate(sdir,dpop,djrec,djphot,dsteq_t,dznet,ddown,dup2,dfull,lowt_ion,jv_pops)
if(allocated(prevp))deallocate(prevp,prevm)
if(allocated(dse))deallocate(dse)
allocate(sdir(nt,nd,nv),dpop(nt,nd,nv),jv_pops(nt,nd))
jv_pops=pops
read(lu,iostat=ios)sdir
close(lu)
if(ios/=0)stop 'QSTARS_JV: cannot read directions'
do k=1,nv
 dpop(:,:,k)=pops*sdir(:,:,k)
enddo
call qs_jv_dis_begin(pops,nt,nd)
nsim=size(znet_sim,2);jv_maxsim=nsim;jv_ntwo=max(n_two,1)
allocate(djrec(nd,nv),djphot(nd,nv),dsteq_t(nd,nv),dznet(nd,nsim,nv))
allocate(ddown(nd,jv_ntwo,nv),dup2(nd,jv_ntwo,nv),dfull(nt,nd,nv))
djrec=0;djphot=0;dsteq_t=0;dznet=0;ddown=0;dup2=0;dfull=0
allocate(dse(size(se)),lowt_ion(size(se)))
lowt_ion=.false.
do id=1,size(se)
 if(.not.se(id)%xzv_pres)cycle
 if(.not.associated(se(id)%steq))cycle
 allocate(dse(id)%steq(size(se(id)%steq,1),nd,nv))
 dse(id)%steq=0
enddo
do id=1,num_ions
 if(atm(id)%xzv_pres)then
  if(minval(atm(id)%xzvlte)==0.0_ldp)lowt_ion(id)=.true.
 endif
enddo
st_closure=0;st_replay=0;st_maxdj=0;st_products=0;st_maxprod=0;st_nfreq=0;st_negfloor=0;st_lincheck=0
call get_environment_variable('QSTARS_JV_EXPORT',path,length=plen,status=ios)
jv_export=((ios==0.or.ios==-1).and.plen>0)
call get_environment_variable('QSTARS_JV_DUAL',path,length=plen,status=ios)
jv_dual=((ios==0.or.ios==-1).and.plen>0)
call get_environment_variable('QSTARS_JV_CHECKLIN',path,length=plen,status=ios)
jv_checklin=((ios==0.or.ios==-1).and.plen>0)
call get_environment_variable('QSTARS_JV_NOASIG',path,length=plen,status=ios)
jv_asig=.not.((ios==0.or.ios==-1).and.plen>0)
call get_environment_variable('QSTARS_JV_ASIG_MIN',path,length=plen,status=ios)
if(ios==0.and.plen>0)read(path(1:plen),*,iostat=ios)jv_asig_min
st_asig=0
call get_environment_variable('QSTARS_JV_NOBATCH',path,length=plen,status=ios)
jv_batch=.not.((ios==0.or.ios==-1).and.plen>0)
call get_environment_variable('QSTARS_JV_NOSETUP_OMP',path,length=plen,status=ios)
lin_setup_omp=.not.((ios==0.or.ios==-1).and.plen>0)
st_batch=0
call get_environment_variable('QSTARS_JV_DIRS',path,length=plen,status=ios)
if(jv_export)then
 open(newunit=jv_lu_exp,file='QSTARS_JV_SPECTRAL.bin',status='replace',access='stream',form='unformatted')
 write(jv_lu_exp)'QSJV0001',int(nv,4),int(nd,4)
endif
jv_active=.true.
open(newunit=lu,file='QSTARS_JV_BEGIN',status='unknown',position='append')
write(lu,*)'directions',nv,' NT',nt,' ND',nd,' max|s|',maxval(abs(sdir)),' file=',path(1:plen)
close(lu)
end subroutine qs_jv_begin

! Build 25: per-direction log-perturbations of the dissolution parameters of COMP_LEV_DIS_BLK,
!   B = (8.589e14*POPION**0.333333/ED)**1.5, A = 0.0009*ED**0.1667/SQRT(T), X = (1+A)**3.15,
! differentiated as written. POPION sums the rows with Z_POP > 0.01 (SET_Z_POP convention).
subroutine qs_jv_dis_begin(pops,nt,nd)
use set_kind_module
use qs_jv_state
use mod_cmfgen, only: atm,num_ions
use mod_lev_dis_blk, only: a_lev_dis,b_lev_dis,x_lev_dis,mod_do_lev_dis
implicit none
integer,intent(in)::nt,nd
real(kind=ldp),intent(in)::pops(nt,nd)
real(kind=ldp)::zp(nt),popion,dpi,dlne,dlnt
integer::id,k,l,ios,plen
character(len=16)::buf
logical,allocatable::ionrow(:)
call get_environment_variable('QSTARS_JV_NODIS',buf,length=plen,status=ios)
jv_dis=.not.((ios==0.or.ios==-1).and.plen>0)
call get_environment_variable('QSTARS_JV_VALL',buf,length=plen,status=ios)
jv_vall=((ios==0.or.ios==-1).and.plen>0)
if(allocated(dis_lb))deallocate(dis_lb,dis_la,dis_lx,dchidis,detadis)
if(allocated(dis_cwb))deallocate(dis_cwb,dis_cwa,dis_cwx)
allocate(dis_lb(nd,nv),dis_la(nd,nv),dis_lx(nd,nv),dchidis(nd,nv),detadis(nd,nv))
allocate(dis_cwb(nt,nd),dis_cwa(nt,nd),dis_cwx(nt,nd))
dis_lb=0;dis_la=0;dis_lx=0;dchidis=0;detadis=0;dis_valid=.false.
dis_cwb=0;dis_cwa=0;dis_cwx=0
st_dis_wdiff=0;st_dis_rel=0
if(.not.mod_do_lev_dis)jv_dis=.false.
if(.not.jv_dis)return
zp=0
do id=1,num_ions-1
 if(.not.atm(id)%xzv_pres)cycle
 zp(atm(id)%eqxzv:atm(id)%eqxzv+atm(id)%nxzv-1)=atm(id)%zxzv-1
 zp(atm(id)%eqxzv+atm(id)%nxzv)=atm(id)%zxzv
enddo
allocate(ionrow(nt))
ionrow=zp>0.01_ldp
do k=1,nv
 do l=1,nd
  popion=sum(pops(:,l),mask=ionrow)
  dpi=sum(dpop(:,l,k),mask=ionrow)
  dlne=dpop(nt-1,l,k)/pops(nt-1,l)
  dlnt=dpop(nt,l,k)/pops(nt,l)
  dis_lb(l,k)=1.5_ldp*(0.333333_ldp*dpi/popion-dlne)
  dis_la(l,k)=0.1667_ldp*dlne-0.5_ldp*dlnt
  dis_lx(l,k)=3.15_ldp*a_lev_dis(l)/(1.0_ldp+a_lev_dis(l))*dis_la(l,k)
 enddo
enddo
! occupation-probability coefficients per level row (OCCUPATION_PROB, as written)
do id=1,num_ions-1
 if(.not.atm(id)%xzv_pres)cycle
 if(atm(id)%zxzv<=0.0_ldp)cycle
 if(atm(id)%nxzv_f/=atm(id)%nxzv)then
  stop 'QSTARS_JV: dissolution line derivatives require one-to-one super levels'
 endif
 call dis_w_coef(id)
enddo
contains
 subroutine dis_w_coef(id)
 integer,intent(in)::id
 real(kind=ldp)::zion,tw,neff,rk,beta,a4,y,f,w,q
 integer::i,l,row
 zion=atm(id)%zxzv
 do i=1,atm(id)%nxzv_f
  row=atm(id)%eqxzv+atm(id)%f_to_s_xzv(i)-1
  tw=3.289395_ldp*zion*zion/atm(id)%edgexzv_f(i)
  neff=0.0_ldp
  if(tw>0.0_ldp.and.tw<961.0_ldp)neff=sqrt(tw)
  if(neff<=2.01_ldp*zion)cycle
  rk=16.0_ldp*neff/(1+neff)/(1+neff)/3.0_ldp
  if(neff<=3.0_ldp)rk=1.0_ldp
  rk=(rk/zion*(zion/neff)**4)**1.5_ldp
  do l=1,nd
   beta=rk*b_lev_dis(l)
   a4=4.0_ldp*(zion-1.0_ldp)*a_lev_dis(l)
   y=1.091_ldp*(x_lev_dis(l)+a4)
   f=y*beta*beta/(7.782_ldp+x_lev_dis(l)*beta)
   w=f/(1.0_ldp+f)
   if(w>0.99999999_ldp)cycle
   q=x_lev_dis(l)*beta/(7.782_ldp+x_lev_dis(l)*beta)
   dis_cwb(row,l)=(1.0_ldp-w)*(2.0_ldp-q)
   dis_cwa(row,l)=(1.0_ldp-w)*a4/(x_lev_dis(l)+a4)
   dis_cwx(row,l)=(1.0_ldp-w)*(x_lev_dis(l)/(x_lev_dis(l)+a4)-q)
  enddo
 enddo
 end subroutine dis_w_coef
end subroutine qs_jv_dis_begin

! Build 25: directional change of the bound-free opacity/emissivity through level dissolution
! at the continuum evaluation frequency NU (= CONT_FREQ), replicating GENOPAETA_V10:
!  nu >= edge : CHI += a*(HN - T1), ETA += 2hnu^3/c^2*a*T1, T1 = exp(LOG_COR_FAC+LOG_HNST)
!  nu <  edge : (route 1, DIS_CONST >= 0, T2 > PHOT_DIS_PARAMETER)
!               CHI += a*T2*(HN - T1'), ETA += 2hnu^3/c^2*a*T2*T1', T1' = exp(LOG_HNST - h nu/kT)
! LOG_HNST = log(W*g) + Saha (LTEPOP_WLD_V2), W from OCCUPATION_PROB: W = F/(1+F),
! F = Y*beta^2/(7.782 + X*beta), beta = K_level*B, Y = 1.091*(X + 4*(Z-1)*A); the dissolved
! fraction T2 = 1/(1+Fc) with beta = DIS_CONST*B. The native VCHI/VETA contain every other
! dependence (populations, n_e and T through the Saha factor) with W and T2 fixed; only
! dW and dT2 are added here.
subroutine qs_jv_dis_freq(nu)
use set_kind_module
use qs_jv_state
use mod_lev_dis_blk, only: b_lev_dis,a_lev_dis,x_lev_dis,mod_do_lev_dis,phot_dis_parameter
use mod_cmfgen, only: atm,num_ions,t
implicit none
real(kind=ldp),intent(in)::nu
real(kind=ldp)::chibf,chiff,hdkt,twohcsq
common /constants/ chibf,chiff,hdkt,twohcsq
real(kind=ldp),allocatable::alpha(:)
real(kind=ldp)::lcf(jv_nd),zion,teta1,tw,neff,rk,beta,y,f,w,q,cy_a,cy_x,cb,ca,cx,dlnw,t1v,t1p
real(kind=ldp)::t1k,disc,fc,t2,dt2,cbc,cac,cxc,emin,hn,a4,cw_b,cw_a,cw_x
integer::id,j,i,l,k,n,ion_lev,nd
logical::wvar
nd=jv_nd
dchidis=0;detadis=0
dis_valid=.true.
if(.not.mod_do_lev_dis)return
teta1=twohcsq*nu**3
do id=1,num_ions-1
 if(.not.atm(id)%xzv_pres)cycle
 zion=atm(id)%zxzv
 if(zion<=0.0_ldp)cycle
 n=atm(id)%nxzv_f
 if(allocated(alpha))deallocate(alpha)
 allocate(alpha(n))
 do j=1,atm(id)%n_xzv_phot
  alpha=0
  call sub_phot_gen(id,alpha,nu,atm(id)%edgexzv_f,n,j,j==1)
  if(count(alpha>0.0_ldp)==0)cycle
  ion_lev=atm(id)%xzv_ion_lev_id(j)
  emin=minval(atm(id)%edgexzv_f(1:n))
  lcf=0
  if(nu>=emin)then
   do l=1,nd
    lcf(l)=log(atm(id+1)%xzv(ion_lev,l)/atm(id+1)%xzv(1,l))+atm(id+1)%log_xzvlte(1,l)-&
           atm(id+1)%log_xzvlte(ion_lev,l)-hdkt*nu/t(l)
   enddo
  endif
  do i=1,n
   if(alpha(i)==0.0_ldp)cycle
   ! occupation-probability constant of the level (OCCUPATION_PROB)
   tw=3.289395_ldp*zion*zion/atm(id)%edgexzv_f(i)
   neff=0.0_ldp
   if(tw>0.0_ldp.and.tw<961.0_ldp)neff=sqrt(tw)
   wvar=neff>2.01_ldp*zion
   rk=0
   if(wvar)then
    rk=16.0_ldp*neff/(1+neff)/(1+neff)/3.0_ldp
    if(neff<=3.0_ldp)rk=1.0_ldp
    rk=(rk/zion*(zion/neff)**4)**1.5_ldp
   endif
   ! dissolution constant below the edge (GENOPAETA, route 1 only)
   disc=-1.0_ldp
   if(j==1.and.nu<atm(id)%edgexzv_f(i))then
    neff=sqrt(3.289395_ldp*zion*zion/(atm(id)%edgexzv_f(i)-nu))
    if(neff>2*zion)then
     t1k=min(1.0_ldp,16.0_ldp*neff/(1+neff)/(1+neff)/3.0_ldp)
     disc=(t1k*zion**3/(neff**4))**1.5_ldp
    endif
   endif
   if(nu<atm(id)%edgexzv_f(i).and.disc<0.0_ldp)cycle
   do l=1,nd
    a4=4.0_ldp*(zion-1.0_ldp)*a_lev_dis(l)
    y=1.091_ldp*(x_lev_dis(l)+a4)
    cy_x=x_lev_dis(l)/(x_lev_dis(l)+a4)
    cy_a=a4/(x_lev_dis(l)+a4)
    ! d ln W = cw_b dlnB + cw_a dlnA + cw_x dlnX
    cw_b=0;cw_a=0;cw_x=0
    if(wvar)then
     beta=rk*b_lev_dis(l)
     f=y*beta*beta/(7.782_ldp+x_lev_dis(l)*beta)
     w=f/(1.0_ldp+f)
     if(w<=0.99999999_ldp)then
      st_dis_wdiff=max(st_dis_wdiff,abs(w-atm(id)%w_xzv_f(i,l)))
      q=x_lev_dis(l)*beta/(7.782_ldp+x_lev_dis(l)*beta)
      cw_b=(1.0_ldp-w)*(2.0_ldp-q)
      cw_a=(1.0_ldp-w)*cy_a
      cw_x=(1.0_ldp-w)*(cy_x-q)
     endif
    endif
    hn=atm(id)%xzv_f(i,l)
    if(nu>=atm(id)%edgexzv_f(i))then
     if(cw_b==0.0_ldp.and.cw_a==0.0_ldp.and.cw_x==0.0_ldp)cycle
     t1v=exp(lcf(l)+atm(id)%log_xzvlte_f(i,l))
     do k=1,nv
      dlnw=cw_b*dis_lb(l,k)+cw_a*dis_la(l,k)+cw_x*dis_lx(l,k)
      dchidis(l,k)=dchidis(l,k)-alpha(i)*t1v*dlnw
      detadis(l,k)=detadis(l,k)+teta1*alpha(i)*t1v*dlnw
     enddo
    else
     beta=disc*b_lev_dis(l)
     fc=y*beta*beta/(7.782_ldp+x_lev_dis(l)*beta)
     t2=1.0_ldp/(1.0_ldp+fc)
     if(t2<=phot_dis_parameter)cycle
     q=x_lev_dis(l)*beta/(7.782_ldp+x_lev_dis(l)*beta)
     cbc=-t2*(1.0_ldp-t2)*(2.0_ldp-q)
     cac=-t2*(1.0_ldp-t2)*cy_a
     cxc=-t2*(1.0_ldp-t2)*(cy_x-q)
     t1p=exp(atm(id)%log_xzvlte_f(i,l)-hdkt*nu/t(l))
     do k=1,nv
      dt2=cbc*dis_lb(l,k)+cac*dis_la(l,k)+cxc*dis_lx(l,k)
      dlnw=cw_b*dis_lb(l,k)+cw_a*dis_la(l,k)+cw_x*dis_lx(l,k)
      dchidis(l,k)=dchidis(l,k)+alpha(i)*(dt2*(hn-t1p)-t2*t1p*dlnw)
      detadis(l,k)=detadis(l,k)+teta1*alpha(i)*(dt2*t1p+t2*t1p*dlnw)
     enddo
    endif
   enddo
  enddo
 enddo
enddo
if(allocated(alpha))deallocate(alpha)
end subroutine qs_jv_dis_freq

! Per-frequency directional products. Call after DO_VAR_CONT (VCHI/VETA valid),
! and before UPDATE_BA_FOR_LINE (QFV_R/QFV_P still as set by EVALSE_QWVJ_V8).
subroutine qs_jv_freq(ml,fl,fqw,first_freq,new_cross,final_cross,nc,np,nd,nt,cont_freq)
use set_kind_module
use qs_jv_state
use qs_formal_tangent, only: tangent_scattering
use qs_formal_tangent_lin, only: lin_setup,lin_scattering,lin_em_build,lin_scattering_multi
use omp_lib, only: omp_get_max_threads
use mod_space_grid_v2, only: ray
use mod_cmfgen, only: atm,num_ions,ion_id,t,ed,pop_atom,clump_fac
use opac_mod, only: veta,vchi,veta_all,vchi_all,eta_clump,chi_clump,chi_scat_clump,chi,eta,chi_cont,eta_cont,&
                    chi_noscat,chi_scat,esec,chi_ray,emhnukt
use radiation_mod, only: rj,dtdr,dbb,ddbbdt,dlog_nu
use var_rad_mod, only: diffw
use line_mod, only: resonance_zone,end_res_zone,line_prof_sim,line_opac_con,line_emis_con,&
                    l_star_ratio,u_star_ratio,gldgu,sim_nl,sim_nup,chil_mat,etal_mat,line_qw_sim,&
                    neg_opac_fac,eina,fl_sim,sim_line_pointer,ave_energy,dl_rat_dt,du_rat_dt
use line_vec_mod, only: vec_id,vec_mnl_f,vec_mnup_f
use steq_data_mod, only: se
use control_variable_mod, only: coherent_es,inner_bnd_meth,neg_opac_option,include_dsldt,&
                                scl_line_cool_rates,scl_line_ht_fac,scl_line_density_limit
use two_phot_mod, only: n_two,two_phot_available,freq_two,type_two,coef_two,g_up_two,g_low_two,&
                        fs_rat_low,fs_rat_up,phot_oc_two,two_method
implicit none
integer,intent(in)::ml,nc,np,nd,nt
real(kind=ldp),intent(in)::fl,fqw,cont_freq
logical,intent(in)::first_freq,new_cross,final_cross
real(kind=ldp)::chibf,chiff,hdkt,twohcsq
common /constants/ chibf,chiff,hdkt,twohcsq
real(kind=ldp)::dchi_abs(nd),deta_c(nd),desec(nd),dray(nd),dchil(nd),detal(nd)
real(kind=ldp)::dchi(nd,nv),deta(nd,nv),dscat(nd,nv),dj(nd,nv),moms(4,nd),dmoms(4,nd),pm(4,nd)
real(kind=ldp)::cerr_k(nv),dmref(4,nd),lindiff(nv)
integer::nprod_k(nv)
real(kind=ldp),allocatable::outp(:,:),outm(:,:),op2(:,:),om2(:,:)
real(kind=ldp)::x,em,den,bnue,dbdt,dtn,db,ddtdr,ddbb,cerr,t1,t2,t3,t4,scl,fb,alpha,y,u,fu,ay
real(kind=ldp)::down_const,up_const,flr,f1,dsum,lbr
integer::k,l,i,sim,nl,nup,ir,id,ipr,ion_lev,ion_eq,j,mnl,mnup,mnl_f,mnup_f,nprod,ihyd
logical::branch1
logical,external::qs_two_phot_j_on
character(len=12)::bnd
if(.not.jv_active)return
if(first_freq)then
 nzmax=0
 do ir=1,np
  nzmax=max(nzmax,ray(ir)%nz)
 enddo
 jv_np=np
 if(allocated(prevp))deallocate(prevp,prevm)
 allocate(prevp(nzmax,np,nv),prevm(nzmax,np,nv))
 prevp=0;prevm=0
endif
do ir=1,np
 if(ray(ir)%nz>nzmax)stop 'QSTARS_JV: ray length grew after first frequency'
enddo
st_nfreq=st_nfreq+1
if(first_freq.and.n_two>0)then
 ! build 32: the tangent's two-photon branch must be the branch of the residual
 if(qs_two_phot_j_on().and.any(two_phot_available(1:n_two)).and.trim(two_method)/='USE_RAD')then
  call qs_two_phot_j_fail('QSTARS_TWO_PHOT_J set but the two-photon method in the continuum loop is '//trim(two_method))
 endif
endif
bnd=trim(inner_bnd_meth)
ihyd=1
if(ion_id(1)=='H0'.or.ion_id(1)=='HMI')ihyd=2
! Diffusion inner boundary primal: identical expressions to COMP_J_BLANK.
x=hdkt*fl/t(nd);em=emhnukt(nd);den=1.0_ldp-em
bnue=twohcsq*(fl**3)*em/den
dbdt=bnue*x/t(nd)/den
! ---------------- phase 1: material directional changes (all directions) ----------
if(jv_dis.and.(new_cross.or.first_freq.or..not.dis_valid))call qs_jv_dis_freq(cont_freq)
do k=1,nv
 do l=1,nd
  if(jv_vall)then
   ! complete variation (important AND not-important variables, COMP_VAR_OPAC VCHI_ALL)
   dchi_abs(l)=dot_product(vchi_all(1:nt,l),dpop(1:nt,l,k))
   deta_c(l)=dot_product(veta_all(1:nt,l),dpop(1:nt,l,k))
  else
   dchi_abs(l)=dot_product(vchi(1:nt,l),dpop(1:nt,l,k))
   deta_c(l)=dot_product(veta(1:nt,l),dpop(1:nt,l,k))
  endif
  if(jv_dis)then
   ! dissolution terms evaluated at CONT_FREQ, moved to FL exactly as COMP_VAR_OPAC moves VCHI/VETA:
   ! dETA(FL) = dETA(CONT)*T4, dCHI(FL) = dCHI(CONT) + dETA(CONT)/T2 - dETA(FL)/T3
   if(cont_freq/=fl)then
    t4=(fl/cont_freq)**3*exp(-hdkt*(fl-cont_freq)/t(l))
    deta_c(l)=deta_c(l)+detadis(l,k)*t4
    dchi_abs(l)=dchi_abs(l)+dchidis(l,k)+detadis(l,k)/(twohcsq*cont_freq**3)-detadis(l,k)*t4/(twohcsq*fl**3)
   else
    dchi_abs(l)=dchi_abs(l)+dchidis(l,k)
    deta_c(l)=deta_c(l)+detadis(l,k)
   endif
  endif
  desec(l)=esec(l)/ed(l)*dpop(nt-1,l,k)
  dray(l)=0
  if(atm(ihyd)%xzv_pres)then
   dray(l)=chi_ray(l)/atm(ihyd)%xzv_f(1,l)*dpop(atm(ihyd)%eqxzv,l,k)
  endif
 enddo
 dscat(:,k)=desec+dray
 dchil=0;detal=0
 do sim=1,jv_maxsim
  if(.not.resonance_zone(sim))cycle
  nl=sim_nl(sim);nup=sim_nup(sim)
  do l=1,nd
   dchil(l)=dchil(l)+line_prof_sim(l,sim)*line_opac_con(sim)*&
        (l_star_ratio(l,sim)*dpop(nl,l,k)-gldgu(sim)*u_star_ratio(l,sim)*dpop(nup,l,k))
   if(jv_dis)then
    ! L_STAR_RATIO = W(up)/W(low) (SET_LINE_OPAC, one-to-one super levels)
    dchil(l)=dchil(l)+line_prof_sim(l,sim)*line_opac_con(sim)*l_star_ratio(l,sim)*jv_pops(nl,l)*&
         ((dis_cwb(nup,l)-dis_cwb(nl,l))*dis_lb(l,k)+(dis_cwa(nup,l)-dis_cwa(nl,l))*dis_la(l,k)+&
          (dis_cwx(nup,l)-dis_cwx(nl,l))*dis_lx(l,k))
   endif
   detal(l)=detal(l)+line_prof_sim(l,sim)*line_emis_con(sim)*u_star_ratio(l,sim)*dpop(nup,l,k)
  enddo
  if(include_dsldt)then
   do l=1,nd
    dchil(l)=dchil(l)+line_prof_sim(l,sim)*line_opac_con(sim)*dpop(nt,l,k)*&
         (dl_rat_dt(l,sim)*jv_pops(nl,l)-gldgu(sim)*du_rat_dt(l,sim)*jv_pops(nup,l))
    detal(l)=detal(l)+line_prof_sim(l,sim)*line_emis_con(sim)*du_rat_dt(l,sim)*jv_pops(nup,l)*dpop(nt,l,k)
   enddo
  endif
 enddo
 deta(:,k)=deta_c+detal
 do l=1,nd
  if(neg_opac_fac(l)==0.0_ldp)then
   ! Native SRCE_CHK/ESEC_CHK floor: CHI=0.1*ETA*CHI_NOSCAT/ETA_CONT or 0.1*CHI_SCAT.
   t1=0.1_ldp*eta(l)*chi_noscat(l)/eta_cont(l)
   t2=0.1_ldp*chi_scat(l)
   branch1=(neg_opac_option=='SRCE_CHK'.and.abs(chi(l)-t1)<=abs(chi(l)-t2))
   if(branch1)then
    dchi(l,k)=0.1_ldp*(deta(l,k)*chi_noscat(l)/eta_cont(l)+eta(l)*dchi_abs(l)/eta_cont(l)-&
                     eta(l)*chi_noscat(l)*deta_c(l)/eta_cont(l)**2)
   else
    dchi(l,k)=0.1_ldp*dscat(l,k)
   endif
   if(k==1)st_negfloor=st_negfloor+1
  else
   dchi(l,k)=dchi_abs(l)+dscat(l,k)+dchil(l)
  endif
 enddo
 dchi(:,k)=dchi(:,k)*clump_fac(1:nd);deta(:,k)=deta(:,k)*clump_fac(1:nd);dscat(:,k)=dscat(:,k)*clump_fac(1:nd)
enddo
! ---------------- phase 2: nonlocal transfer tangent ----------------
pm=0
if(.not.jv_dual)then
 call lin_setup(eta_clump(1:nd)+chi_scat_clump(1:nd)*rj(1:nd),chi_clump(1:nd),chi_scat_clump(1:nd),&
                bnue,dbb,dlog_nu,nc,bnd,pm)
 st_replay=max(st_replay,maxval(abs(pm(1,:)-rj(1:nd)))/max(maxval(abs(rj(1:nd))),tiny(1.0_ldp)))
 if(jv_asig.and.nv>=jv_asig_min)then
  call lin_em_build()
  st_asig=st_asig+1
 endif
endif
lindiff=0
! Directions in parallel only when they can occupy every thread. Otherwise the directions run in
! sequence OUTSIDE any OpenMP construct, so the ray-parallel region of lin_apply is a top-level
! region. Build 28: the former inactive enclosing region (parallel do with a false IF clause) made
! lin_apply's team a nested team; libgomp on macOS then intermittently deadlocked in the team
! barrier or crashed (Report 21). The arithmetic is unchanged.
if(jv_batch.and..not.jv_dual.and..not.jv_checklin.and.jv_asig.and.nv>=jv_asig_min)then
 call dir_batch()
 st_batch=st_batch+1
else if(nv>=omp_get_max_threads())then
!$omp parallel do schedule(dynamic) default(shared) private(k)
 do k=1,nv
  call dir_tangent(k)
 enddo
!$omp end parallel do
else
 do k=1,nv
  call dir_tangent(k)
 enddo
endif
do k=1,nv
 if(maxval(abs(dj(:,k)))>1.0e-200_ldp*maxval(rj(1:nd)).and.cerr_k(k)>st_closure)then
  st_closure=cerr_k(k);st_closure_fl=fl
 endif
enddo
st_products=st_products+sum(nprod_k);st_maxprod=max(st_maxprod,maxval(nprod_k))
st_lincheck=max(st_lincheck,maxval(lindiff))
if(jv_export)then
 write(jv_lu_exp)fl,fqw,rj(1:nd),pm(1,:)
 do k=1,nv
  write(jv_lu_exp)dchi(:,k),deta(:,k),dscat(:,k),dj(:,k)
 enddo
endif
! ---------------- phase 3: radiation part of the residual derivative ---------
do k=1,nv
 ! Bound-free (EVALSE_QWVJ_V8 / BA_UPDATE_V7): band-integrated JREC/JPHOT.
 if(new_cross)then
  djrec(:,k)=0;djphot(:,k)=0
 endif
 djrec(:,k)=djrec(:,k)+dj(:,k)*emhnukt(1:nd)*fqw/fl
 djphot(:,k)=djphot(:,k)+dj(:,k)*fqw/fl
 if(final_cross)then
  do id=1,num_ions-1
   if(.not.atm(id)%xzv_pres)cycle
   if(.not.allocated(dse(id)%steq))cycle
   do l=1,nd
    do i=1,size(se(id)%qfv_r,1)
     dse(id)%steq(i,l,k)=dse(id)%steq(i,l,k)+se(id)%qfv_r(i,l)*djrec(l,k)-se(id)%qfv_p(i,l)*djphot(l,k)
    enddo
   enddo
  enddo
 endif
 ! Low-temperature recombination (EVALSE_LOWT_V2), per frequency.
 do id=1,num_ions
  if(.not.lowt_ion(id))cycle
  do ipr=1,atm(id)%n_xzv_phot
   ion_lev=atm(id)%xzv_ion_lev_id(ipr)
   if(ion_lev==0)exit
   ion_eq=se(id)%ion_lev_to_eq_pnt(ion_lev)
   do l=1,nd
    flr=dj(l,k)*fqw/fl
    lbr=(atm(id+1)%log_xzvlte(1,l)-atm(id+1)%log_xzvlte(ion_lev,l))+&
        log(atm(id+1)%xzv(ion_lev,l)/atm(id+1)%xzv(1,l))
    dsum=0
    do j=1,atm(id)%nxzv
     if(atm(id)%xzvlte(j,l)==0.0_ldp.and.atm(id)%wsxzv(j,l,ipr)>0.0_ldp)then
      t1=flr*atm(id)%wsxzv(j,l,ipr)*exp(atm(id)%log_xzvlte(j,l)+lbr-hdkt*fl/t(l))
      dse(id)%steq(j,l,k)=dse(id)%steq(j,l,k)+t1
      dsum=dsum+t1
     endif
    enddo
    dse(id)%steq(ion_eq,l,k)=dse(id)%steq(ion_eq,l,k)-dsum
   enddo
  enddo
 enddo
 ! Radiative equilibrium row.
 dsteq_t(:,k)=dsteq_t(:,k)+fqw*chi_noscat(1:nd)*dj(:,k)
 if(.not.coherent_es)dsteq_t(:,k)=dsteq_t(:,k)+fqw*esec(1:nd)*dj(:,k)
 ! Lines: ZNET accumulation and S.E./R.E. update at the end of the resonance zone.
 do sim=1,jv_maxsim
  if(resonance_zone(sim))then
   dznet(:,sim,k)=dznet(:,sim,k)-line_qw_sim(1:nd,sim)*dj(:,k)*chil_mat(1:nd,sim)/etal_mat(1:nd,sim)
  endif
  if(end_res_zone(sim))then
   nl=sim_nl(sim);nup=sim_nup(sim)
   i=sim_line_pointer(sim);id=vec_id(i)
   mnl_f=vec_mnl_f(i);mnl=atm(id)%f_to_s_xzv(mnl_f)
   mnup_f=vec_mnup_f(i);mnup=atm(id)%f_to_s_xzv(mnup_f)
   scl=1.0_ldp
   if(scl_line_cool_rates)then
    scl=(ave_energy(nl)-ave_energy(nup))/fl_sim(sim)
    if(abs(scl-1.0_ldp)>scl_line_ht_fac)scl=1.0_ldp
   endif
   do l=1,nd
    t4=scl
    if(pop_atom(l)>=scl_line_density_limit)t4=1.0_ldp
    t2=eina(sim)*atm(id)%xzv_f(mnup_f,l)*dznet(l,sim,k)
    t3=t4*etal_mat(l,sim)*dznet(l,sim,k)
    dse(id)%steq(mnup,l,k)=dse(id)%steq(mnup,l,k)-t2
    dse(id)%steq(mnl,l,k)=dse(id)%steq(mnl,l,k)+t2
    dsteq_t(l,k)=dsteq_t(l,k)-t3
   enddo
   dznet(:,sim,k)=0
  endif
 enddo
 ! Two-photon (TWO_PHOT_RATE), same branch as the residual of this pass (module TWO_METHOD, set by
 ! TWO_PHOT_OPAC_V3 at this frequency's COMP_OPAC call):
 !  OLD_DEFAULT (COMP_F=T without QSTARS_TWO_PHOT_J): partner-photon occupation exp(-h nu_B/kT_e).
 !  otherwise (USE_RAD with QSTARS_TWO_PHOT_J, build 32; LTE; NOSTIM): PHOT_OC_TWO of the residual. For USE_RAD it
 !  is fixed (lagged) data interpolated from the stored-J file and does not depend on the state, so its
 !  derivative is zero and only dJ at the current frequency enters. The opacity/emissivity variation
 !  (VCHI/VETA, TWO_PHOT_VAR_OPAC) and the RE row use the same PHOT_OC_TWO through the native arrays.
 do j=1,n_two
  if(.not.(two_phot_available(j).and.fl<freq_two(j)))cycle
  fb=freq_two(j)-fl
  alpha=twohcsq*(fl**3)
  if(type_two(j)/=1)stop 'QSTARS_JV: unrecognized two-photon type'
  y=fl/freq_two(j);u=y*(1.0_ldp-y);fu=4.0_ldp*u
  ay=24.56_ldp*coef_two(j,1)*(u*(1.0_ldp-fu**0.8_ldp)+0.88_ldp*(u**1.53_ldp)*(fu**0.8_ldp))
  down_const=0.5e-15_ldp*ay/freq_two(j)
  up_const=down_const*g_up_two(j)/g_low_two(j)
  if(two_method=='OLD_DEFAULT')then
   up_const=up_const/alpha
   do l=1,nd
    ddown(l,j,k)=ddown(l,j,k)+down_const*fs_rat_up(l,j)*(dj(l,k)/alpha)*fqw
    dup2(l,j,k)=dup2(l,j,k)+up_const*fs_rat_low(l,j)*exp(-hdkt*fb/t(l))*dj(l,k)*fqw
   enddo
  else
   do l=1,nd
    f1=dj(l,k)/alpha
    ddown(l,j,k)=ddown(l,j,k)+down_const*fs_rat_up(l,j)*(f1+f1*phot_oc_two(l,j))*fqw
    dup2(l,j,k)=dup2(l,j,k)+up_const*fs_rat_low(l,j)*f1*phot_oc_two(l,j)*fqw
   enddo
  endif
 enddo
 st_maxdj=max(st_maxdj,maxval(abs(dj(:,k))/max(rj(1:nd),tiny(1.0_ldp))))
enddo
contains
 ! Transfer tangent of direction k (formerly the body of the direction loop; locals are private).
 subroutine dir_tangent(k)
 integer,intent(in)::k
 real(kind=ldp),allocatable::outp(:,:),outm(:,:),op2(:,:),om2(:,:)
 real(kind=ldp)::moms(4,nd),dmoms(4,nd),dmref(4,nd),dtn,db,ddtdr,ddbb,cerr
 integer::nprod
 allocate(outp(nzmax,np),outm(nzmax,np))
 dtn=dpop(nt,nd,k)
 db=dbdt*dtn
 ddtdr=dot_product(diffw(1:nt),dpop(1:nt,nd,k))
 ddbb=ddbbdt*dtn+dbb/dtdr*ddtdr
 if(jv_dual)then
  call tangent_scattering(eta_clump(1:nd),chi_clump(1:nd),chi_scat_clump(1:nd),rj(1:nd),&
       deta(:,k),dchi(:,k),dscat(:,k),bnue,dbb,db,ddbb,dlog_nu,nc,bnd,&
       prevp(:,:,k),prevm(:,:,k),moms,dmoms,outp,outm,cerr,.true.,nprod)
 else
  call lin_scattering(chi_scat_clump(1:nd),rj(1:nd),deta(:,k),dchi(:,k),dscat(:,k),db,ddbb,&
       prevp(:,:,k),prevm(:,:,k),dmoms,outp,outm,cerr,nprod)
  if(jv_checklin.and.k==1)then
   allocate(op2(nzmax,np),om2(nzmax,np))
   call tangent_scattering(eta_clump(1:nd),chi_clump(1:nd),chi_scat_clump(1:nd),rj(1:nd),&
        deta(:,k),dchi(:,k),dscat(:,k),bnue,dbb,db,ddbb,dlog_nu,nc,bnd,&
        prevp(:,:,k),prevm(:,:,k),moms,dmref,op2,om2,cerr,.true.,nprod)
   lindiff(k)=maxval(abs(dmref(1,:)-dmoms(1,:)))/max(maxval(abs(dmref(1,:))),tiny(1.0_ldp))
   deallocate(op2,om2)
  endif
 endif
 prevp(:,:,k)=outp;prevm(:,:,k)=outm
 dj(:,k)=dmoms(1,:)
 cerr_k(k)=cerr;nprod_k(k)=nprod
 deallocate(outp,outm)
 end subroutine dir_tangent
 ! Build 31: all directions in min(nv,threads) contiguous blocks, blocks in parallel (top-level region; nothing
 ! inside is parallel). Per direction exactly the inputs of dir_tangent/lin_scattering; ray arrays transposed
 ! to (direction,point,ray) for the block and back.
 subroutine dir_batch()
 integer::nbk,ib,k0,k1,m,c,kk
 real(kind=ldp),allocatable::bde(:,:),bdc(:,:),bds(:,:),bdb(:),bddbb(:),bpp(:,:,:),bpm(:,:,:),bop(:,:,:),&
   bom(:,:,:),bmo(:,:,:),bce(:)
 integer,allocatable::bnp(:)
 real(kind=ldp)::dtn_,ddtdr_
 nbk=max(1,min(nv,omp_get_max_threads()))
!$omp parallel do schedule(dynamic) default(shared) private(ib,k0,k1,m,c,kk,bde,bdc,bds,bdb,bddbb,bpp,bpm,&
!$omp bop,bom,bmo,bce,bnp,dtn_,ddtdr_)
 do ib=1,nbk
  k0=((ib-1)*nv)/nbk+1;k1=(ib*nv)/nbk;m=k1-k0+1
  allocate(bde(nd,m),bdc(nd,m),bds(nd,m),bdb(m),bddbb(m),bpp(m,nzmax,np),bpm(m,nzmax,np),&
           bop(m,nzmax,np),bom(m,nzmax,np),bmo(4,nd,m),bce(m),bnp(m))
  do c=1,m
   kk=k0+c-1
   dtn_=dpop(nt,nd,kk)
   bdb(c)=dbdt*dtn_
   ddtdr_=dot_product(diffw(1:nt),dpop(1:nt,nd,kk))
   bddbb(c)=ddbbdt*dtn_+dbb/dtdr*ddtdr_
   bde(:,c)=deta(:,kk);bdc(:,c)=dchi(:,kk);bds(:,c)=dscat(:,kk)
   bpp(c,:,:)=prevp(:,:,kk);bpm(c,:,:)=prevm(:,:,kk)
  enddo
  call lin_scattering_multi(m,chi_scat_clump(1:nd),rj(1:nd),bde,bdc,bds,bdb,bddbb,bpp,bpm,bmo,bop,bom,bce,bnp)
  do c=1,m
   kk=k0+c-1
   prevp(:,:,kk)=bop(c,:,:);prevm(:,:,kk)=bom(c,:,:)
   dj(:,kk)=bmo(1,:,c)
   cerr_k(kk)=bce(c);nprod_k(kk)=bnp(c)
  enddo
  deallocate(bde,bdc,bds,bdb,bddbb,bpp,bpm,bop,bom,bmo,bce,bnp)
 enddo
!$omp end parallel do
 end subroutine dir_batch
end subroutine qs_jv_freq

! Two-photon S.E. terms after the frequency loop (STEQ_BA_TWO_PHOT_RATE_V3).
subroutine qs_jv_two_post(pops,nt,nd)
use set_kind_module
use qs_jv_state
use omp_lib, only: omp_get_max_threads
use steq_data_mod, only: se
use two_phot_mod, only: n_two,two_phot_available,ion_id_two,ion_low_lev_two,ion_up_lev_two,&
                        low_lev_two,up_lev_two,two_method
use qs_two_phot_j_mod, only: tpj_on,tpj_ncall
implicit none
integer,intent(in)::nt,nd
real(kind=ldp),intent(in)::pops(nt,nd)
integer::k,l,j,id
real(kind=ldp)::t1
character(len=512)::path
integer::lu
if(.not.jv_active)return
do k=1,nv
 do l=1,nd
  do j=1,n_two
   if(.not.two_phot_available(j))cycle
   id=ion_id_two(j)
   t1=pops(up_lev_two(j),l)*ddown(l,j,k)-pops(low_lev_two(j),l)*dup2(l,j,k)
   dse(id)%steq(ion_low_lev_two(j),l,k)=dse(id)%steq(ion_low_lev_two(j),l,k)+t1
   dse(id)%steq(ion_up_lev_two(j),l,k)=dse(id)%steq(ion_up_lev_two(j),l,k)-t1
  enddo
 enddo
enddo
if(jv_export)then
 close(jv_lu_exp);jv_lu_exp=-1
endif
open(newunit=lu,file='QSTARS_JV_STATS',status='unknown',position='append')
write(lu,'(A,I8,A,I4,A,ES12.4,A,ES12.4,A,I12,A,I6,A,ES12.4,A,I8,A,ES12.4,A,ES16.8)')'freqs',st_nfreq,' dirs',nv,&
  ' max_tangent_closure',st_closure,' max_primal_replay_rel',st_replay,' scattering_products',st_products,&
  ' max_products',st_maxprod,' max_rel_dJ',st_maxdj,' neg_floor_cells',st_negfloor,' lin_vs_dual',st_lincheck,' closure_freq',st_closure_fl
write(lu,'(A,L2,A,ES12.4,A,L2)')'dissolution_derivatives',jv_dis,' max|W_formula-W_native|',st_dis_wdiff,' VCHI_ALL',jv_vall
write(lu,'(A,L2,A,I6,A,I8)')'explicit_scattering_operator',jv_asig,' min_dirs',jv_asig_min,' freqs_built',st_asig
write(lu,'(A,L2,A,I8,A,I4)')'batched_direction_sweeps',jv_batch,' freqs_batched',st_batch,' blocks',&
  max(1,min(nv,omp_get_max_threads()))
write(lu,'(A,A,A,L2,A,I12)')'two_photon_method ',trim(two_method),' stored_J',tpj_on,' conj_J_evaluations',tpj_ncall
close(lu)
end subroutine qs_jv_two_post

! Map the directional S.E. derivatives to the full-matrix row layout exactly as
! GENERATE_FULL_MATRIX_V3 maps STEQ (same REPLACE and ZERO_STEQ decisions).
subroutine qs_jv_map(kd,nt,nion,replace,zero_steq)
use set_kind_module
use qs_jv_state
use mod_cmfgen, only: atm,num_species,species_beg_id,species_end_id,species_pres,eq_species
use steq_data_mod, only: se
implicit none
integer,intent(in)::kd,nt,nion
logical,intent(in)::replace(nion),zero_steq(nt)
real(kind=ldp)::vec(nt),vion(nion),vnc(num_species)
integer::k,ispec,id,i,eq,lu
if(.not.jv_active)return
do k=1,nv
 vec=0;vion=0;vnc=0
 do ispec=1,num_species
  do id=species_beg_id(ispec),species_end_id(ispec)-1
   if(.not.allocated(dse(id)%steq))cycle
   do i=1,se(id)%n_se-1
    eq=se(id)%eq_in_ba(i)
    vec(eq)=vec(eq)+dse(id)%steq(i,kd,k)
   enddo
   do i=atm(id)%nxzv+1,se(id)%n_se-1
    vion(id)=vion(id)-dse(id)%steq(i,kd,k)
   enddo
   vnc(ispec)=vnc(ispec)+dse(id)%steq(se(id)%n_se,kd,k)
  enddo
 enddo
 vec(nt-1)=0
 vec(nt)=dsteq_t(kd,k)
 do ispec=1,num_species
  if(species_pres(ispec))vec(eq_species(ispec))=vnc(ispec)
 enddo
 do id=1,nion
  if(replace(id).and.atm(id)%xzv_pres)vec(atm(id)%eqxzv)=vion(id)
 enddo
 do i=1,nt
  if(zero_steq(i))vec(i)=0
 enddo
 dfull(:,kd,k)=vec
enddo
if(kd==jv_nd)then
 open(newunit=lu,file='QSTARS_JV_RAD.bin',status='replace',access='stream',form='unformatted')
 write(lu)int(nv,4),int(nt,4),int(jv_nd,4),sdir,dfull
 close(lu)
endif
end subroutine qs_jv_map

! Optional dump of every assembled full block (after REPLACE/FIXPOP and the native
! column scaling by POPS) so that the direct Jacobian action can be applied offline.
subroutine qs_jv_ba_dump(cmat,nt,kd,band)
use set_kind_module
implicit none
integer,intent(in)::nt,kd,band
real(kind=ldp),intent(in)::cmat(nt,nt)
character(len=8)::flag
integer::plen,ios,lu
logical,save::first=.true.,on=.false.
if(first)then
 first=.false.
 call get_environment_variable('QSTARS_BA_DUMP',flag,length=plen,status=ios)
 on=((ios==0.or.ios==-1).and.plen>0)
 if(on)then
  open(newunit=lu,file='QSTARS_BA_TRACE.bin',status='replace',access='stream',form='unformatted')
  close(lu)
 endif
endif
if(.not.on)return
open(newunit=lu,file='QSTARS_BA_TRACE.bin',status='old',access='stream',form='unformatted',position='append')
write(lu)int(kd,4),int(band,4),int(nt,4),cmat
close(lu)
end subroutine qs_jv_ba_dump

! Record (and optionally force) the ground-state/ionization equation selection.
! A forced table fixes the residual definition across a nonlinear solve.
subroutine qs_jv_replace(replace,nion,kd,nd)
use set_kind_module
implicit none
integer,intent(in)::nion,kd,nd
logical,intent(inout)::replace(nion)
character(len=512)::path
integer::plen,ios,lu,kk,i
integer,allocatable,save::forced(:,:)
logical,save::first=.true.,force=.false.
if(first)then
 first=.false.
 call get_environment_variable('QSTARS_REPLACE_FORCE',path,length=plen,status=ios)
 force=((ios==0.or.ios==-1).and.plen>0)
 if(force)then
  allocate(forced(nion,nd))
  open(newunit=lu,file=path(1:plen),status='old',iostat=ios)
  if(ios/=0)stop 'QSTARS_REPLACE_FORCE: cannot open'
  do i=1,nd
   read(lu,*,iostat=ios)kk,forced(:,i)
   if(ios/=0.or.kk/=i)stop 'QSTARS_REPLACE_FORCE: bad table'
  enddo
  close(lu)
 endif
 open(newunit=lu,file='QSTARS_REPLACE_TRACE',status='replace')
 close(lu)
endif
if(force)replace(:)=(forced(:,kd)/=0)
open(newunit=lu,file='QSTARS_REPLACE_TRACE',status='old',position='append')
write(lu,'(I5,*(I2))')kd,merge(1,0,replace)
close(lu)
end subroutine qs_jv_replace

logical function qs_jv_stop_on()
implicit none
character(len=8)::flag
integer::plen,ios
call get_environment_variable('QSTARS_JV_STOP',flag,length=plen,status=ios)
qs_jv_stop_on=((ios==0.or.ios==-1).and.plen>0)
end function qs_jv_stop_on

! Stored converged J per frequency (base pass) and its replay.
subroutine qs_rj_write(ml,fl,rj,nd)
use set_kind_module
implicit none
integer,intent(in)::ml,nd
real(kind=ldp),intent(in)::fl,rj(nd)
character(len=512)::path
integer::plen,ios
integer,save::lu=0
logical,save::first=.true.,on=.false.,opened=.false.
if(first)then
 first=.false.
 call get_environment_variable('QSTARS_RJ_STORE',path,length=plen,status=ios)
 on=((ios==0.or.ios==-1).and.plen>0)
endif
if(.not.on)return
if(ml==1)then
 if(opened)close(lu)
 open(newunit=lu,file='QSTARS_RJ.bin',status='replace',access='stream',form='unformatted')
 write(lu)'QSRJ0001',int(nd,4)
 opened=.true.
endif
if(.not.opened)stop 'QSTARS_RJ_STORE: frequency sweep did not start at ML=1'
write(lu)int(ml,4),fl,rj
flush(lu)
end subroutine qs_rj_write

subroutine qs_rj_read(ml,fl,rj,nd)
use set_kind_module
implicit none
integer,intent(in)::ml,nd
real(kind=ldp),intent(in)::fl
real(kind=ldp),intent(out)::rj(nd)
character(len=512)::path
character(len=8)::tag
integer::plen,ios,mm,n4
real(kind=ldp)::ff
integer,save::lu=0
logical,save::opened=.false.
if(ml==1)then
 call get_environment_variable('QSTARS_JV_REPLAY',path,length=plen,status=ios)
 if(ios/=0)stop 'QSTARS_JV_REPLAY: path missing or longer than 512 characters'
 if(opened)close(lu)
 open(newunit=lu,file=path(1:plen),status='old',access='stream',form='unformatted',iostat=ios)
 opened=(ios==0)
 if(ios/=0)stop 'QSTARS_JV_REPLAY: cannot open stored J'
 read(lu,iostat=ios)tag,n4
 if(ios/=0.or.tag/='QSRJ0001'.or.n4/=nd)stop 'QSTARS_JV_REPLAY: bad stored-J header'
endif
if(.not.opened)stop 'QSTARS_JV_REPLAY: sweep did not start at ML=1'
read(lu,iostat=ios)mm,ff,rj
if(ios/=0.or.mm/=ml.or.abs(ff/fl-1.0_ldp)>1.0e-14_ldp)stop 'QSTARS_JV_REPLAY: stored J does not match frequency'
if(any(.not.(rj>0.0_ldp)))stop 'QSTARS_JV_REPLAY: stored J not positive'
end subroutine qs_rj_read

logical function qs_jv_replay_on()
implicit none
character(len=8)::flag
integer::plen,ios
logical,save::first=.true.,on=.false.
if(first)then
 first=.false.
 call get_environment_variable('QSTARS_JV_REPLAY',flag,length=plen,status=ios)
 on=((ios==0.or.ios==-1).and.plen>0)
endif
qs_jv_replay_on=on
end function qs_jv_replay_on

! Build 29 (Report 22): QSTARS_NO_OBS_FORMAL skips the observer-frame formal solution (CMF_FORM_SOL_V2) in
! solver passes. Its only product is the observed spectrum (IPLUS -> COMP_OBS_V2 -> OBSFLUX); the residual,
! BA, stored J, the tangent and the luminosity guard (comoving RLUMST) never read it. Unset (default) it runs.
logical function qs_no_obs_on()
implicit none
character(len=8)::flag
integer::plen,ios
logical,save::first=.true.,on=.false.
if(first)then
 first=.false.
 call get_environment_variable('QSTARS_NO_OBS_FORMAL',flag,length=plen,status=ios)
 on=((ios==0.or.ios==-1).and.plen>0)
endif
qs_no_obs_on=on
end function qs_no_obs_on

! Optional ALI warm start (build 23): stored J of a nearby state, same grid.
subroutine qs_warm_read(ml,fl,rj,nd)
use set_kind_module
implicit none
integer,intent(in)::ml,nd
real(kind=ldp),intent(in)::fl
real(kind=ldp),intent(inout)::rj(nd)
character(len=512)::path
character(len=8)::tag
integer::plen,ios,mm,n4
real(kind=ldp)::ff,tmp(nd)
integer,save::lu=0
logical,save::first=.true.,on=.false.,opened=.false.
if(first)then
 first=.false.
 call get_environment_variable('QSTARS_ALI_WARM',path,length=plen,status=ios)
 if(ios==-1)stop 'QSTARS_ALI_WARM: path longer than 512 characters'
 on=(ios==0.and.plen>0)
 if(on)then
  open(newunit=lu,file=path(1:plen),status='old',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 'QSTARS_ALI_WARM: cannot open'
  read(lu,iostat=ios)tag,n4
  if(ios/=0.or.tag/='QSRJ0001'.or.n4/=nd)stop 'QSTARS_ALI_WARM: bad header'
  opened=.true.
 endif
endif
if(.not.on)return
if(ml==1)then
 rewind(lu)
 read(lu)tag,n4
endif
read(lu,iostat=ios)mm,ff,tmp
if(ios/=0.or.mm/=ml.or.abs(ff/fl-1.0_ldp)>1.0e-14_ldp)stop 'QSTARS_ALI_WARM: stored J does not match frequency'
if(all(tmp>0.0_ldp))rj=tmp
end subroutine qs_warm_read

! ---------------------------------------------------------------------------------------------------------
! Build 32 (Report 23): QSTARS_TWO_PHOT_J, see module qs_two_phot_j_mod.
logical function qs_two_phot_j_on()
use qs_two_phot_j_mod
implicit none
integer::ios
if(tpj_first)then
 tpj_first=.false.
 call get_environment_variable('QSTARS_TWO_PHOT_J',tpj_path,length=tpj_plen,status=ios)
 if(ios==-1)then
  tpj_plen=len(tpj_path)
  call qs_two_phot_j_fail('path longer than 512 characters')
 endif
 tpj_on=(ios==0.and.tpj_plen>0)
endif
qs_two_phot_j_on=tpj_on
end function qs_two_phot_j_on

subroutine qs_two_phot_j_fail(msg)
use qs_two_phot_j_mod
implicit none
character(len=*),intent(in)::msg
integer,external::error_lu
integer::lu,luer
open(newunit=lu,file='QSTARS_TWO_PHOT_J_FAILURE',status='unknown',position='append')
write(lu,'(A)')trim(msg)
write(lu,'(A)')'file='//tpj_path(1:max(tpj_plen,1))
close(lu)
luer=error_lu()
write(luer,'(1X,A)')'Error: QSTARS_TWO_PHOT_J: '//trim(msg)
write(luer,'(1X,A)')'       file='//tpj_path(1:max(tpj_plen,1))
flush(luer)
stop 'QSTARS_TWO_PHOT_J failure (see QSTARS_TWO_PHOT_J_FAILURE)'
end subroutine qs_two_phot_j_fail

! Load and verify the stored J once per run, after the frequency grid is set and before any opacity is
! evaluated. The file must hold exactly this model's frequency grid (every FL to 1e-14, as the replay reader
! requires) and depth count, with positive J. Any mismatch stops the run.
subroutine qs_two_phot_j_init(nu,ncf,nd)
use set_kind_module
use qs_two_phot_j_mod
use control_variable_mod, only: two_photon_method,incl_two_phot
implicit none
integer,intent(in)::ncf,nd
real(kind=ldp),intent(in)::nu(ncf)
logical,external::qs_two_phot_j_on
integer,external::error_lu
character(len=8)::tag
character(len=160)::msg
integer::lu,luer,ios,n4,mm,ml
integer(kind=8)::fsize,reclen
logical::ex
real(kind=ldp)::ff
if(.not.qs_two_phot_j_on())return
if(tpj_ready)return
if(.not.incl_two_phot)call qs_two_phot_j_fail('set, but two-photon transitions are not included (INC_TWO=F)')
if(trim(two_photon_method)/='USE_RAD')call qs_two_phot_j_fail('set, but TWO_METH is '//trim(two_photon_method)//&
   ', not USE_RAD')
inquire(file=tpj_path(1:tpj_plen),exist=ex,size=fsize)
if(.not.ex)call qs_two_phot_j_fail('stored-J file does not exist')
reclen=12_8+8_8*int(nd,8)
if(fsize/=12_8+reclen*int(ncf,8))then
 write(msg,'(A,I14,A,I8,A,I6)')'stored-J file size ',fsize,' bytes is not the QSRJ0001 layout for NCF=',ncf,' ND=',nd
 call qs_two_phot_j_fail(trim(msg))
endif
open(newunit=lu,file=tpj_path(1:tpj_plen),status='old',action='read',access='stream',form='unformatted',iostat=ios)
if(ios/=0)call qs_two_phot_j_fail('cannot open stored-J file')
read(lu,iostat=ios)tag,n4
if(ios/=0.or.tag/='QSRJ0001'.or.n4/=nd)then
 write(msg,'(A,A,A,I8,A,I6)')'bad stored-J header: tag ',tag,' ND(file)',n4,' ND(model)',nd
 call qs_two_phot_j_fail(trim(msg))
endif
if(allocated(tpj_fl))deallocate(tpj_fl,tpj_rj)
allocate(tpj_fl(ncf),tpj_rj(nd,ncf))
do ml=1,ncf
 read(lu,iostat=ios)mm,ff,tpj_rj(:,ml)
 if(ios/=0.or.mm/=ml)then
  write(msg,'(A,I8,A,I8,A,I6)')'cannot read stored-J record ',ml,' (index in file',mm,', iostat',ios
  call qs_two_phot_j_fail(trim(msg)//')')
 endif
 if(abs(ff/nu(ml)-1.0_ldp)>1.0e-14_ldp)then
  write(msg,'(A,I8,A,ES24.16,A,ES24.16)')'stored-J frequency grid differs at ML=',ml,': file',ff,' model',nu(ml)
  call qs_two_phot_j_fail(trim(msg))
 endif
 if(any(.not.(tpj_rj(:,ml)>0.0_ldp)))then
  write(msg,'(A,I8)')'stored J not positive at ML=',ml
  call qs_two_phot_j_fail(trim(msg))
 endif
 if(ml>1)then
  if(.not.(ff<tpj_fl(ml-1)))then
   write(msg,'(A,I8)')'stored-J frequencies not strictly decreasing at ML=',ml
   call qs_two_phot_j_fail(trim(msg))
  endif
 endif
 tpj_fl(ml)=ff
enddo
close(lu)
tpj_ncf=ncf;tpj_nd=nd;tpj_ncall=0
tpj_ready=.true.
luer=error_lu()
write(luer,'(1X,A)')'QSTARS_TWO_PHOT_J: two-photon method USE_RAD is kept in COMP_F=T passes; the conjugate-'//&
   'frequency J is interpolated from the stored-J file (fixed, lagged data; not EDDFACTOR, not OLD_DEFAULT):'
write(luer,'(1X,A)')'   '//tpj_path(1:tpj_plen)
write(luer,'(1X,A,I8,A,I6,A,ES22.14,A,ES22.14)')'   NCF=',ncf,' ND=',nd,' NU_MAX=',tpj_fl(1),' NU_MIN=',tpj_fl(ncf)
end subroutine qs_two_phot_j_init

! J at the conjugate frequency FREQ from the stored file: the statements of GET_J_FOR_TWO_PHOT (ND = ND_RD
! case) with the EDDFACTOR records replaced by the stored records.
!   FREQ >= FL(1)   : J(1)                       (first record)
!   FREQ <= FL(NCF) : J(NCF)*(FREQ/FL(NCF))**2   (Rayleigh-Jeans extrapolation)
!   otherwise       : T1=(FL_LOW-FREQ)/(FL_LOW-FL_HIGH); J=T1*J_HIGH+(1-T1)*J_LOW  on the bracketing records.
! LST_REC is the native search hint (accepted if it brackets FREQ, otherwise bisection, which then sets it).
subroutine qs_two_phot_get_j(rj,freq,lst_rec,nd)
use set_kind_module
use qs_two_phot_j_mod
implicit none
integer,intent(in)::nd
integer,intent(inout)::lst_rec
real(kind=ldp),intent(out)::rj(nd)
real(kind=ldp),intent(in)::freq
integer::ilow,ihigh,j
real(kind=ldp)::t1
if(.not.tpj_ready)call qs_two_phot_j_fail('conjugate-frequency J requested before the stored-J file was loaded')
if(nd/=tpj_nd)call qs_two_phot_j_fail('conjugate-frequency J requested on a different depth grid')
if(.not.(freq>0.0_ldp))call qs_two_phot_j_fail('conjugate frequency not positive')
tpj_ncall=tpj_ncall+1
if(freq>=tpj_fl(1))then
 rj=tpj_rj(:,1)
 return
endif
if(freq<=tpj_fl(tpj_ncf))then
 t1=(freq/tpj_fl(tpj_ncf))**2
 rj=tpj_rj(:,tpj_ncf)*t1
 return
endif
ilow=0
if(lst_rec>0.and.lst_rec<tpj_ncf)then
 if(freq<=tpj_fl(lst_rec).and.freq>=tpj_fl(lst_rec+1))ilow=lst_rec
endif
if(ilow==0)then
 ilow=1;ihigh=tpj_ncf
 do while(ihigh-ilow>1)
  j=(ilow+ihigh)/2
  if(freq<tpj_fl(j))then
   ilow=j
  else
   ihigh=j
  endif
 enddo
 lst_rec=ilow
endif
t1=(tpj_fl(ilow)-freq)/(tpj_fl(ilow)-tpj_fl(ilow+1))
rj=t1*tpj_rj(:,ilow+1)+(1.0_ldp-t1)*tpj_rj(:,ilow)
end subroutine qs_two_phot_get_j
