! Thin CGS interface to separately installed, unmodified FreeEOS 3.0.0 physics.
! EOS1 without radiation: (ifoption,ifmodified,ifion)=(3,11,-2).
! No source from the GPL library is copied into this interface.
module qstars_freeeos
 use iso_c_binding
 use mod_free_eos, only: free_eos
 use mod_isotopic_mass_data, only: isotopic_mass
 use mod_free_eos_constants, only: prad_const
 implicit none
contains
 integer(c_int) function qstars_freeeos_version() bind(C)
  qstars_freeeos_version=1
 end function
 subroutine qstars_freeeos_constants(arad) bind(C)
  real(c_double), intent(out) :: arad
  arad=3d0*prad_const
 end subroutine
 subroutine qstars_freeeos_eval(mode,n,values,temps,hydrogen,out,ierr) bind(C)
  integer(c_int), value :: mode,n
  real(c_double), intent(in) :: values(n),temps(n)
  real(c_double), value :: hydrogen
  real(c_double), intent(out) :: out(24,n)
  integer(c_int), intent(out) :: ierr
  real(c_double) :: eps(20),fl,t,rho,rl,p,pl,cf,cp,s,sf,st,grada,rtp
  real(c_double) :: qe,qv,rmue,fh2,fhe2,fhe3,xmu1,xmu3,eta,g1,g2,g3,h2,h2plus,lambda,gamma_e,sound2
  real(c_double) :: pressure(3),energy(3),entropy(3),density(3)
  integer :: i,iterations,info
  ierr=0;out=0d0
  if ((mode/=1.and.mode/=2).or.hydrogen<0d0.or.hydrogen>1d0) then
   ierr=-1;return
  endif
  eps=0d0;eps(1)=hydrogen/isotopic_mass(1);eps(2)=(1d0-hydrogen)/isotopic_mass(2)
  do i=1,n
   if(values(i)<=0d0.or.temps(i)<=0d0) then
    ierr=-2;return
   endif
   call free_eos(0,3,11,-2,mode,eps,log(values(i)),log(temps(i)),fl, &
       t,rho,rl,p,pl,cf,cp,s,sf,st,grada,rtp,qe,qv,rmue,fh2,fhe2,fhe3,xmu1,xmu3,eta, &
       g1,g2,g3,h2,h2plus,lambda,gamma_e,sound2,iterations,info, &
       pressure=pressure,energy=energy,entropy=entropy,density=density)
   if(info/=0) then
    ierr=info;return
   endif
   ! Derivatives are w.r.t. ln rho for mode=2, ln Pgas for mode=1, and ln T.
   out(:,i)=[pressure,energy,entropy,density,cp,grada,g1,xmu3,eta,1d0/xmu1, &
       fh2,fhe2,fhe3,h2,h2plus,real(iterations,c_double)]
  enddo
 end subroutine
end module
