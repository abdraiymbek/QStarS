program check_microphysics
 use iso_c_binding
 use qstars_mesa_microphysics
 implicit none
 real(c_double) :: rho(1),temp(1),x,arr(21,1),dd(21,1),tt(21,1)
 integer(c_int) :: ierr
 integer :: ios
 call qstars_init(ierr)
 if (ierr/=0) stop 1
 do
  read(*,*,iostat=ios) rho(1),temp(1),x
  if (ios/=0) exit
  call qstars_eval(1,rho,temp,x,arr,dd,tt,ierr)
  if (ierr/=0) stop 2
  write(*,'(21es26.17)') arr(:,1)
 enddo
end program
