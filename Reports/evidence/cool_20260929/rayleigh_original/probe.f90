program probe
 use set_kind_module
 implicit none
 integer, parameter :: ND=3
 type ion
   real(kind=LDP),allocatable :: XzV_F(:,:), AXzV_F(:,:), EDGEXZV_F(:)
   integer :: NXzV_F, EQXzV
 end type
 type(ion) :: ATM(2)
 character(3) :: ION_ID(2)
 logical :: SPECIES_PRES(1),INCL_RAY_SCAT
 real(kind=LDP) :: CHI_RAY(ND), CONT_FREQ, VJ(40,1,ND), dj_LOC(5,1,ND)
 integer :: ID,I,J,K,L
 ION_ID=[character(3)::'H0','HI']; SPECIES_PRES=.true.; INCL_RAY_SCAT=.true.
 do ID=1,2
   ATM(ID)%NXzV_F=1
   if(ID==2)ATM(ID)%NXzV_F=30
   I=ATM(ID)%NXzV_F
   allocate(ATM(ID)%XzV_F(I,ND),ATM(ID)%AXzV_F(I,I),ATM(ID)%EDGEXZV_F(I))
   ATM(ID)%XzV_F=0.7_LDP; ATM(ID)%AXzV_F=0; ATM(ID)%EDGEXZV_F=0.182365_LDP
   ATM(ID)%EQXzV=ID
 end do
 ATM(2)%XzV_F(1,:)=[2.0e10_LDP,5.0e11_LDP,9.0e12_LDP]
 ATM(2)%EDGEXZV_F(1)=3.288087_LDP
 CONT_FREQ=0.6_LDP
	  CHI_RAY(1:ND)=0.0_LDP
	  ID=1; IF(ION_ID(1) .EQ. 'HMI')ID=2
	  IF(SPECIES_PRES(1) .AND. INCL_RAY_SCAT)THEN
	    CALL RAYLEIGH_SCAT(CHI_RAY,ATM(ID)%XzV_F,ATM(ID)%AXzV_F,ATM(ID)%EDGEXZV_F, &
     &             ATM(1)%NXzV_F,CONT_FREQ,ND)
	  END IF

 write(*,'(3ES25.16)') CHI_RAY
 VJ=0; dj_LOC=1; J=1
 do K=1,ND
   L=K
	      IF(SPECIES_PRES(1))THEN
	        IF(ION_ID(1) .EQ. 'HMI')THEN
	          VJ(2,J,K)=VJ(2,J,K) + CHI_RAY(L)*dj_LOC(5,J,K)/ATM(2)%XzV_F(1,L)
	        ELSE
	          VJ(1,J,K)=VJ(1,J,K) + CHI_RAY(L)*dj_LOC(5,J,K)/ATM(1)%XzV_F(1,L)
	        END IF
	      END IF
 end do
 write(*,'(3ES25.16)') VJ(2,1,:)
 write(*,'(3ES25.16)') VJ(1,1,:)
end program
