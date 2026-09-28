! QStarS interface to the public MESA 24.08.1 EOS and opacity APIs.
! No MESA source is copied; link against the separately installed libraries.
module qstars_mesa_microphysics
 use iso_c_binding
 use const_def, only: dp, crad
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
contains
 integer(c_int) function qstars_api_version() bind(C)
  qstars_api_version=3
 end function
 subroutine qstars_constants(arad) bind(C)
  real(c_double), intent(out) :: arad
  arad=crad
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
