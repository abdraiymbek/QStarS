program probe
 use set_kind_module
 use mod_cmfgen
 implicit none
 integer,parameter :: ND=3,NT=40
 real(kind=LDP) :: chi(ND),eta(ND),cool(ND),bp(ND),jp(ND),em(ND),nu
 real(kind=LDP) :: vc(NT,ND),ve(NT,ND),vcn(NT,ND),ven(NT,ND)
 real(kind=LDP) :: hs(1,ND),di(30,ND),dl(30,ND),edge(1),g(1),x
 real(kind=LDP) :: chiE(ND),etaE(ND),vcE(NT,ND),veE(NT,ND)
 real(kind=LDP) :: CHIBF,CHIFF,HDKT,TWOHCSQ
 common/constants/CHIBF,CHIFF,HDKT,TWOHCSQ
 integer :: map(1),i,k,mode,eq
 logical :: imp(NT)
 read(*,*)mode,nu,T,ED,di(1,:)
 CHIBF=1;CHIFF=1;HDKT=4.79924307336622_LDP;TWOHCSQ=0.014745003464_LDP
 chi=0;eta=0;cool=0;hs=1;edge=3.288;g=2;dl=0;map=1;imp=.true.;eq=2
 do k=1,ND
  x=HDKT*nu/T(k);em(k)=exp(-x)
  bp(k)=TWOHCSQ*nu**3*em(k)/(1-em(k));jp(k)=0.37_LDP*bp(k)
 end do
 if(mode==1)call DO_H0_FF_COOL(cool,di(1,:),ED,T,bp,jp,nu,ND)
 call DO_H0_FF(eta,chi,di(1,:),ED,T,em,nu,7,ND)
 if(mode/=1)call DO_H0_FF_COOL(cool,di(1,:),ED,T,bp,jp,nu,ND)
 do k=1,ND
  write(*,'(A,3ES25.16)')'DATA ',chi(k),eta(k),cool(k)
 end do
 if(mode/=0)return
 di(2:,:)=0
 vc=0;ve=0;vcn=0;ven=0
 call VAR_OP_V10(vc,ve,vcn,ven,hs,hs,hs,hs,1,hs,edge,1,map, &
  di,dl,dl,30,1,1,ED,T,em,imp,nu,0.0_LDP,1,.true.,1,eq,NT,ND,.false.)
 do k=1,ND
  write(*,'(A,6ES25.16)')'DERIV ',vc(eq,k),ve(eq,k),vc(NT-1,k),ve(NT-1,k),vc(NT,k),ve(NT,k)
 end do
 ION_ID=[character(3)::'H0','HI']
 do i=1,3
  ATM(i)%XzV_PRES=.false.;ATM(i)%NXzV=30;ATM(i)%EQXzV=eq;ATM(i)%ZXzV=1
  allocate(ATM(i)%XzV(30,ND),ATM(i)%GIONXzV_F(30))
  ATM(i)%XzV=0;ATM(i)%GIONXzV_F=2
 end do
 ATM(1)%XzV_PRES=.true.;ATM(1)%ZXzV=0
 ATM(2)%XzV(1,:)=di(1,:);ATM(2)%XzV(2,:)=0.2_LDP*di(1,:)
 call COMP_FREE_FREE_V2(chiE,etaE,vcE,veE,nu,nu,.true.,.true.,.true.,ND,NT)
 do k=1,ND
  write(*,'(A,8ES25.16)')'EHB ',chiE(k),etaE(k),vcE(eq,k),veE(eq,k), &
    vcE(NT-1,k),veE(NT-1,k),vcE(NT,k),veE(NT,k)
 end do
end program
subroutine SUB_PHOT_GEN(ID,a,nu,edge,n,phot,dis)
 use set_kind_module
 integer :: ID,n,phot
 real(kind=LDP) :: a(n),nu,edge(n)
 logical :: dis
 a=0
end subroutine
subroutine GFF_VEC(a,nu,t,z,n)
 use set_kind_module
 integer :: n
 real(kind=LDP) :: a(n),nu,t(n),z
 a=1
end subroutine
subroutine FF_RES_GAUNT(a,nu,t,id,g,z,n)
 use set_kind_module
 integer :: n,id
 real(kind=LDP) :: a(n),nu,t(n),g(*),z
end subroutine
