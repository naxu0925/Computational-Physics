program MCintegration

implicit none
real(8):: r1,r2,rho1,rho2,rho,eps,I_x,I_z,ave_x=0,err_x=0,ave_z=0,err_z=0 
!I_x means moment of inertia with x as the axis, similar to I_z
real:: x,y,z
real(8),external:: ran
integer:: n=1000000,nbin=10,i,j
open(1,file='read.in')
read(1,*) r1,r2,rho1,rho2,eps
close(1)
call initran(1)

do j=1,nbin
 I_x=0
 I_z=0
 do i=1,n
 x=(dble(ran())-0.5)*2.*r1
 y=(dble(ran())-0.5)*2.*r1
 z=(dble(ran())-0.5)*2.*r1
 if((x**2+y**2+z**2)>(r1**2)) then
  rho=0 
 else if((x**2+y**2)>(r2**2)) then
 rho=rho1
 else if((x**2+y**2)<(r2**2)) then
 rho=rho2
 end if
 I_z=I_z+rho*(x**2+y**2)
 I_x=I_x+rho*(z**2+y**2)
 end do
 I_z=I_z*8.d0*r1**3/dble(n)
 I_x=I_x*8.d0*r1**3/dble(n)
 ave_z=ave_z+I_z
 ave_x=ave_x+I_x
 err_z=err_z+I_z**2
 err_x=err_x+I_x**2
end do
ave_z=ave_z/dble(nbin)
ave_x=ave_x/dble(nbin)
err_z=err_z/dble(nbin)
err_x=err_x/dble(nbin)
err_z=sqrt((err_z-ave_z**2)/dble(nbin-1))
err_x=sqrt((err_x-ave_x**2)/dble(nbin-1))
do while ( ((err_z/ave_z) .GE. eps) .or.( (err_x/ave_x) .GE. eps) )
 call binup(nbin,ave_z,err_z,ave_x,err_x,r1,r2,rho1,rho2)
end do
print*,"The bin number is: ", nbin-1
print*,"The moment of inertia about z axis Iz is: ", ave_z
print*,"The relative standard deviation of Iz is: ", err_z/ave_z
print*,"The moment of inertia about x axis Ix is: ", ave_x
print*,"The relative standard deviation of Ix is: ", err_x/ave_x

end program MCintegration

!It takes about 40 seconds to get the result, the result is as blow
!The bin number is 181~185 :
!The moment of inertia around z axis is: 4.691(6)*10^-3 
!The relative standard deviation of Iz is : 9.0(6)*10^-5
!The moment of inertia around x axis is: 4.947(7)*10^-3
!The relative standard deviation of Ix is:9.(8)*10^-5

 
!subroutine binup can add one bin number for each time,and update the new average and error of Iz and Ix 

subroutine binup(nbin,ave_z,err_z,ave_x,err_x,r1,r2,rho1,rho2)

integer:: nbin, m=1000000,k
real(8):: ave_z,err_z,ave_x,err_x,x,y,z,r1,r2,rho1,rho2,rho
real(8):: Ix=0,Iz=0
do k=0,m
 x=(dble(ran())-0.5)*2.*r1
 y=(dble(ran())-0.5)*2.*r1
 z=(dble(ran())-0.5)*2.*r1
 if((x**2+y**2+z**2)>(r1**2)) then
 rho=0 
 else if((x**2+y**2)>(r2**2)) then
 rho=rho1
 else if((x**2+y**2)<(r2**2)) then
 rho=rho2
 end if
 Iz=Iz+rho*(x**2+y**2)
 Ix=Ix+rho*(z**2+y**2)
 end do
 Iz=Iz*8.d0*r1**3/dble(m)
 Ix=Ix*8.d0*r1**3/dble(m)
 err_z=(err_z**2*(nbin-1)+ave_z**2)*nbin+Iz**2
 err_x=(err_x**2*(nbin-1)+ave_x**2)*nbin+Ix**2   
 Iz=ave_z*nbin+Iz
 Ix=ave_x*nbin+Ix
 nbin=nbin+1
 ave_z=Iz/dble(nbin)
 ave_x=Ix/dble(nbin)
 err_z=sqrt((err_z/dble(nbin)-ave_z**2)/dble(nbin-1))
 err_x=sqrt((err_x/dble(nbin)-ave_x**2)/dble(nbin-1))
 return
 end

 

 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common  dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif
 return
 end subroutine initran

