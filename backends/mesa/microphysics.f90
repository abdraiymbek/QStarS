! QStarS interface to the public MESA 24.08.1 EOS and opacity APIs.
! No MESA source is copied; link against the separately installed libraries.
module qstars_mesa_microphysics
 use iso_c_binding
 use const_def, only: dp, crad, standard_cgrav, pi, boltz_sigma
 use const_lib, only: const_init
 use math_lib, only: math_init
 use chem_def, only: num_chem_isos, ih1, ihe4
 use chem_lib, only: chem_init
 use eos_def
 use eos_lib
 use kap_def, only: num_kap_fracs, Kap_General_Info
 use kap_lib
 implicit none
 integer :: eh, kh, kh_conduction
 integer, pointer :: chem_id(:),net_iso(:)
 logical :: initialized=.false.
 real(dp) :: atmosphere_x, atmosphere_fixed_kappa
contains
 integer(c_int) function qstars_api_version() bind(C)
  qstars_api_version=5
 end function
 subroutine qstars_constants(arad) bind(C)
  real(c_double), intent(out) :: arad
  arad=crad
 end subroutine

 ! Opacity-only reference: the supplied electrons and degeneracy come from
 ! the same physical state as the caller's EOS. Derivatives are not exported.
 subroutine qstars_opacity_external(n,rho,temp,hydrogen,free_e,eta,kappa,ierr) bind(C)
  integer(c_int), value :: n
  real(c_double), intent(in) :: rho(n),temp(n),free_e(n),eta(n)
  real(c_double), value :: hydrogen
  real(c_double), intent(out) :: kappa(n)
  integer(c_int), intent(out) :: ierr
  real(dp) :: xa(2),fracs(num_kap_fracs),dkxa(2),dkd,dkt
  integer :: i
  ierr=0;xa=[hydrogen,1d0-hydrogen]
  do i=1,n
   if(min(rho(i),temp(i),free_e(i))<=0d0) then
    ierr=-104;return
   endif
   call kap_get(kh,2,chem_id,net_iso,xa,log10(rho(i)),log10(temp(i)), &
       log(free_e(i)),0d0,0d0,eta(i),0d0,0d0,fracs,kappa(i),dkd,dkt,dkxa,ierr)
   if(ierr/=0)return
  enddo
 end subroutine

 subroutine qstars_density_pt(pgas,temp,hydrogen,rho,ierr) bind(C)
  real(c_double), value :: pgas,temp,hydrogen
  real(c_double), intent(out) :: rho
  integer(c_int), intent(out) :: ierr
  real(dp) :: logr,ddp,ddt,xa(2),res(num_eos_basic_results),dd(num_eos_basic_results),tt(num_eos_basic_results)
  real(dp) :: dxa(num_eos_d_dxa_results,2)
  ierr=0
  if (pgas<=0d0.or.temp<=0d0) then
   ierr=-101; return
  endif
  xa=[hydrogen,1d0-hydrogen]
  call eosPT_get(eh,2,chem_id,net_iso,xa,pgas,log10(pgas),temp,log10(temp), &
       rho,logr,ddp,ddt,res,dd,tt,dxa,ierr)
 end subroutine

 ! Reference to MESA's public atmosphere routines, including their unmodified
 ! starting-pressure prescription. mode=0 is fixed opacity; mode=1 varying
 ! integration with a constant-opacity callback; mode=2 native EOS/opacity.
 subroutine qstars_atmosphere(mode,tau,teff,grav,hydrogen,kappa,tolerance,out,ierr) bind(C)
  use atm_lib, only: atm_eval_T_tau_uniform, atm_eval_T_tau_varying
  use atm_def, only: ATM_T_TAU_EDDINGTON
  integer(c_int), value :: mode
  real(c_double), value :: tau,teff,grav,hydrogen,kappa,tolerance
  real(c_double), intent(out) :: out(2)
  integer(c_int), intent(out) :: ierr
  real(dp) :: lum,rad,mass,kap,dtl,dtr,dtm,dtk,dpl,dpr,dpm,dpk
  ierr=0;out=0d0
  atmosphere_x=hydrogen;atmosphere_fixed_kappa=kappa
  if(mode==2) atmosphere_fixed_kappa=-1d0
  rad=1d14;mass=grav*rad**2/standard_cgrav
  lum=4d0*pi*rad**2*boltz_sigma*teff**4
  if(mode==0) then
   call atm_eval_T_tau_uniform(tau,lum,rad,mass,standard_cgrav,kappa,1d0, &
       ATM_T_TAU_EDDINGTON,atmosphere_eos,atmosphere_kap,tolerance,0,.true.,teff,kap, &
       out(1),dtl,dtr,dtm,dtk,out(2),dpl,dpr,dpm,dpk,ierr)
  else
   call atm_eval_T_tau_varying(tau,lum,rad,mass,standard_cgrav, &
       ATM_T_TAU_EDDINGTON,atmosphere_eos,atmosphere_kap,tolerance,10000,.true.,teff, &
       out(1),dtl,dtr,dtm,dtk,out(2),dpl,dpr,dpm,dpk,ierr)
  endif
 end subroutine

 subroutine atmosphere_eos(lnp,lnt,lnrho,res,dd,tt,ierr)
  real(dp), intent(in) :: lnp,lnt
  real(dp), intent(out) :: lnrho,res(:),dd(:),tt(:)
  integer, intent(out) :: ierr
  real(dp) :: pgas,temp,rho,logr,drdp,drdt,xa(2),dxa(num_eos_d_dxa_results,2)
  temp=exp(lnt);pgas=exp(lnp)-crad*temp**4/3d0
  ierr=0
  ! Unlike the example callback, do not turn negative gas pressure into a
  ! tiny positive number. No hydrostatic positive-density state then exists.
  if(pgas<=0d0) then
   ierr=-102;return
  endif
  xa=[atmosphere_x,1d0-atmosphere_x]
  call eosPT_get(eh,2,chem_id,net_iso,xa,pgas,log10(pgas),temp,log10(temp), &
       rho,logr,drdp,drdt,res,dd,tt,dxa,ierr)
  lnrho=logr*log(10d0)
 end subroutine

 subroutine atmosphere_kap(lnrho,lnt,res,dd,tt,kappa,dkd,dkt,ierr)
  real(dp), intent(in) :: lnrho,lnt,res(:),dd(:),tt(:)
  real(dp), intent(out) :: kappa,dkd,dkt
  integer, intent(out) :: ierr
  real(dp) :: xa(2),fracs(num_kap_fracs),dkxa(2)
  ierr=0
  if(atmosphere_fixed_kappa>0d0) then
   kappa=atmosphere_fixed_kappa;dkd=0d0;dkt=0d0;return
  endif
  xa=[atmosphere_x,1d0-atmosphere_x]
  call kap_get(kh,2,chem_id,net_iso,xa,lnrho/log(10d0),lnt/log(10d0), &
       res(i_lnfree_e),dd(i_lnfree_e),tt(i_lnfree_e),res(i_eta),dd(i_eta),tt(i_eta), &
       fracs,kappa,dkd,dkt,dkxa,ierr)
 end subroutine

 subroutine qstars_mlt(mode,n,alpha,inputs,out,jac,ierr) bind(C)
  use auto_diff
  use turb, only: set_mlt
  integer(c_int), value :: mode,n
  real(c_double), value :: alpha
  ! chiT, chiRho, cp, gravity, rho, pressure, T, opacity, grad_rad, grad_ad.
  real(c_double), intent(in) :: inputs(10,n)
  real(c_double), intent(out) :: out(4,n)
  real(c_double), intent(out) :: jac(4,10,n)
  integer(c_int), intent(out) :: ierr
  type(auto_diff_real_star_order1) :: a(10),length,eff,grad,y,velocity,diffusion
  character(len=8) :: option
  integer :: i,j,mixing_type
  ierr=0;out=0d0
  if(mode==0) then
   option='ML1'
  else if(mode==1) then
   option='Henyey'
  else
   ierr=-103;return
  endif
  do i=1,n
   do j=1,10
    a(j)=inputs(j,i)
    a(j)%d1Array(j)=1d0
   enddo
   length=alpha*a(6)/(a(5)*a(4))
   eff=-1d99;grad=a(9);y=0d0;velocity=0d0;diffusion=0d0
   ! The Henyey parameters 8 and 1/3 are MESA's documented defaults.
   call set_mlt(trim(option),alpha,8d0,1d0/3d0,a(1),a(2),a(3),a(4),length, &
       a(5),a(6),a(7),a(8),a(9),a(10),a(10),eff,grad,y,velocity,diffusion,mixing_type,ierr)
   if(ierr/=0) return
   out(:,i)=[grad%val,velocity%val,eff%val,diffusion%val]
   do j=1,10
    jac(:,j,i)=[grad%d1Array(j),velocity%d1Array(j),eff%d1Array(j),diffusion%d1Array(j)]
   enddo
  enddo
 end subroutine
 subroutine qstars_init(ierr) bind(C)
  integer(c_int), intent(out) :: ierr
  type(Kap_General_Info), pointer :: kp
  ierr=0
  if (initialized) return
  call const_init('',ierr)
  if (ierr/=0) return
  call math_init()
  call chem_init('isotopes.data',ierr)
  if (ierr/=0) return
  call eos_init('',.true.,ierr)
  if (ierr/=0) return
  eh=alloc_eos_handle(ierr)
  if (ierr/=0) return
  call kap_init(.true.,'',ierr)
  if (ierr/=0) return
  kh=alloc_kap_handle(ierr)
  if (ierr/=0) return
  call kap_ptr(kh,kp,ierr)
  if (ierr/=0) return
  kp% Zbase=0d0
  kp% use_Type2_opacities=.true.
  kp% include_electron_conduction=.false.
  kh_conduction=alloc_kap_handle(ierr)
  if (ierr/=0) return
  call kap_ptr(kh_conduction,kp,ierr)
  if (ierr/=0) return
  kp% Zbase=0d0
  kp% use_Type2_opacities=.true.
  kp% include_electron_conduction=.true.
  allocate(chem_id(2),net_iso(num_chem_isos))
  chem_id=[ih1,ihe4];net_iso=0
  net_iso(ih1)=1;net_iso(ihe4)=2
  initialized=.true.
 end subroutine

 subroutine qstars_eval(n,rho,temp,hydrogen,out,dr,dt,ierr) bind(C)
  integer(c_int), value :: n
  real(c_double), intent(in) :: rho(n),temp(n)
  real(c_double), value :: hydrogen
  ! EOS first 19 entries, radiative opacity, then radiation+conduction opacity.
  real(c_double), intent(out) :: out(21,n),dr(21,n),dt(21,n)
  integer(c_int), intent(out) :: ierr
  real(dp) :: xa(2),res(num_eos_basic_results),dd(num_eos_basic_results),tt(num_eos_basic_results)
  real(dp) :: dxa(num_eos_d_dxa_results,2),fracs(num_kap_fracs),dkxa(2),dkd,dkt,kappa
  integer :: i,j
  integer, parameter :: ids(19)=[i_lnPgas,i_lnE,i_lnS,i_Cv,i_Cp,i_gamma1,i_grad_ad, &
                                i_chiRho,i_chiT,i_lnfree_e,i_eta,i_mu, &
                                i_frac_HELM,i_frac_OPAL_SCVH,i_frac_FreeEOS,i_frac_PC,i_frac_Skye,i_frac_CMS,i_frac_ideal]
  ierr=0
  xa=[hydrogen,1d0-hydrogen]
  do i=1,n
   call eosDT_get(eh,2,chem_id,net_iso,xa,rho(i),log10(rho(i)),temp(i),log10(temp(i)),res,dd,tt,dxa,ierr)
   if (ierr/=0) return
   do j=1,19
    out(j,i)=res(ids(j));dr(j,i)=dd(ids(j));dt(j,i)=tt(ids(j))
   enddo
   call kap_get(kh,2,chem_id,net_iso,xa,log10(rho(i)),log10(temp(i)), &
       res(i_lnfree_e),dd(i_lnfree_e),tt(i_lnfree_e),res(i_eta),dd(i_eta),tt(i_eta), &
       fracs,kappa,dkd,dkt,dkxa,ierr)
   if (ierr/=0) return
   out(20,i)=kappa;dr(20,i)=kappa*dkd;dt(20,i)=kappa*dkt
   call kap_get(kh_conduction,2,chem_id,net_iso,xa,log10(rho(i)),log10(temp(i)), &
       res(i_lnfree_e),dd(i_lnfree_e),tt(i_lnfree_e),res(i_eta),dd(i_eta),tt(i_eta), &
       fracs,kappa,dkd,dkt,dkxa,ierr)
   if (ierr/=0) return
   out(21,i)=kappa;dr(21,i)=kappa*dkd;dt(21,i)=kappa*dkt
  enddo
 end subroutine
end module
