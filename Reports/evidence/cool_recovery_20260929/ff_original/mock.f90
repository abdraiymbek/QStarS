module mod_lev_dis_blk
 use set_kind_module
 logical :: MOD_DO_LEV_DIS=.false.
 real(kind=LDP) :: PHOT_DIS_PARAMETER=0.0001_LDP
 real(kind=LDP) :: X_LEV_DIS(3),A_LEV_DIS(3),B_LEV_DIS(3)
end module
module mod_cmfgen
 use set_kind_module
 type ion
  logical :: XzV_PRES
  real(kind=LDP),allocatable :: XzV(:,:),GIONXzV_F(:)
  integer :: EQXzV,NXzV
  real(kind=LDP) :: ZXzV
 end type
 integer,parameter :: NUM_IONS=2
 type(ion) :: ATM(3)
 character(3) :: ION_ID(2)
 real(kind=LDP) :: ED(3),T(3)
end module
