program probe
use set_kind_module
use mod_ray_mom_store
implicit none
integer,parameter::ndp=5,nm=2
real(kind=ldp)::eta(ndp),chi(ndp),esec(ndp),theta(ndp),v(ndp),sigma(ndp),r(ndp),j(ndp),h(ndp),dop(ndp)
real(kind=ldp)::tx(ndp,ndp,nm),tvx(ndp-1,ndp,nm),ki(ndp,ndp,nm),work(ndp,ndp),rhsc(ndp-1,ndp)
real(kind=ldp)::txt(ndp),txdt(ndp),tvxt(ndp),tvxdt(ndp),dbb,hob,hib,fr,amp,z
logical::init,do_tx(nm)
integer::a,b,mode
read(*,*)mode,amp
nd_store=ndp
allocate(r_store(ndp),rmid_store(ndp),ext_rmid_store(ndp+1),jnu_store(ndp),hnu_store(ndp),nnu_store(ndp),knu_store(ndp),gam_rel_store(ndp))
r=[5.0_ldp,4.0_ldp,3.0_ldp,2.0_ldp,1.0_ldp];v=50; sigma=-1;dop=100
r_store=r;gam_rel_store=1/sqrt(1-(v*3.33564e-6_ldp)**2);jnu_store=1
! Actual nonnegative two-angle distributions. All four stored moments agree.
do a=1,ndp
 z=real(a,ldp)/(ndp+1)
 hnu_store(a)=z*.3_ldp+(1-z)*.8_ldp
 knu_store(a)=z*.3_ldp**2+(1-z)*.8_ldp**2
 nnu_store(a)=z*.3_ldp**3+(1-z)*.8_ldp**3
enddo
hn_def_on_nodes=.true.
hnu_at_ob=hnu_store(1);nnu_at_ob=nnu_store(1);hnu_at_ib=hnu_store(ndp);nnu_at_ib=nnu_store(ndp)
hbc=hnu_store(1);nbc=nnu_store(1);in_hbc=hnu_store(ndp)
tx=0;tvx=0;ki=0;work=0;rhsc=0;txt=0;txdt=0;tvxt=0;tvxdt=0;do_tx=.true.;dbb=1
eta=[.1_ldp,.2_ldp,.3_ldp,.4_ldp,.5_ldp];chi=[1._ldp,1.5_ldp,2._ldp,2.5_ldp,3._ldp];esec=.1_ldp
if(mode==1)eta(3)=eta(3)+amp
if(mode==2)chi(3)=chi(3)+amp
! Preserve coherent scattering at fixed scattering coefficient in opacity FD.
theta=esec/chi
 do b=1,2
 init=b==1;fr=1._ldp-.01_ldp*(b-1)
 call mom_jrel_v9(eta,chi,esec,v,sigma,r,j,h,hib,hob,dop,.5_ldp,fr,.01_ldp,dbb,&
 'NONE','NONE','NONE',0._ldp,'DIFFUSION','HONJ','LINEAR',.true.,'N_ON_J',.true.,.true.,init,ndp)
 call var_jrel_v5(eta,chi,esec,theta,v,sigma,r,tx,tvx,txt,txdt,tvxt,tvxdt,ki,work,rhsc,&
 init,fr,.01_ldp,'NONE','NONE','NONE','DIFFUSION','HONJ',0._ldp,1._ldp,dbb,1._ldp,1._ldp,&
 .true.,.true.,do_tx,'LINEAR',ndp,nm,nm)
 enddo
 do a=1,ndp
 write(*,'(A,3ES25.16)')'DATA ',j(a),tx(a,3,2),tx(a,3,1)
 enddo
end program
