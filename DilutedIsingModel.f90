!This program simulates 2D diluted Ising model with Monte Carlo method
! Ennergy, magnetization are measured

!-----------------------------!
 module systemvariables
!-----------------------------!
 implicit none

 integer :: l                        ! system length
 integer :: n                        ! number of spins (n=l*l)
 real(8) :: pflip(-1:1,-4:4)         ! flip probabilities
 integer, allocatable :: spin(:)     ! spin array

 end module systemvariables
!-----------------------------!



!-----------------------------!
program dilute_ising
!------------------------------!
use systemvariables; implicit none

 integer :: i,j,k,bins,binsteps,seed,nt
 real(8) :: tmax,t,dt
 print *, "input the size of the lattice, l:"
 read *,l
 n=l*l
 open(10,file='read.in',status='old')
 read(10,*) nt,tmax,dt
 read(10,*) bins, binsteps
 close(10)
 t=tmax
 do k=1,nt
    call initialize(t)
    do i=1,binsteps
       call mcstep
    enddo
    do j=1,bins
       call resetdatasums
       do i=1,binsteps
          call mcstep
          call measure
       enddo
      call writebindata(binsteps)
    enddo
    deallocate(spin)
    t=t-dt
  enddo

end program dilute_ising
!--------------------------------!


!-----------------------------!
 subroutine writebindata(steps)
!------------------------------!
 use systemvariables; implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,num1,num2
 common/measurments/enrg1,enrg2,magn1,magn2,num1,num2

 integer :: steps
 
 open(1,file='bindata.dat',status='unknown',position='append')
 enrg1=enrg1/(dble(steps)*dble(n))
 enrg2=enrg2/(dble(steps)*dble(n)**2)
 magn1=magn1/(dble(steps)*dble(n))
 magn2=magn2/(dble(steps)*dble(n)**2)
 num1=num1/(dble(steps)*dble(n))
 num2=num2/(dble(steps)*dble(n)**2)
 write(1,*) enrg1,enrg2,magn1,magn2,num1,num2
 close(1)

end subroutine writebindata
!----------------------------------!



!----------------------------------!
subroutine measure
!----------------------------------!
 use systemvariables; implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,num1,num2
 common/measurments/enrg1,enrg2,magn1,magn2,num1,num2

 integer :: s,x,y,e,m
 real(8) :: h

 e=0
 do s=0,n-1
    x=mod(s,l); y=s/l
    e=e-spin(s)*(spin(mod(x+1,l)+y*l)+spin(x+mod(y+1,l)*l))
    num1=num1+abs(spin(s))
 enddo
 m=sum(spin)
 enrg1=enrg1+dble(e)
 enrg2=enrg2+dble(e)**2

 magn1=magn1+dble(abs(m))
 magn2=magn2+dble(m)**2
 num2=num2+num1**2

 end subroutine measure
!-------------------------!


!------------------------!
subroutine resetdatasums
!------------------------!
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,num1,num2
 common/measurments/enrg1,enrg2,magn1,magn2,num1,num2

 enrg1=0.d0
 enrg2=0.d0
 magn1=0.d0
 magn2=0.d0
 num1=0.d0
 num2=0.d0


 end subroutine resetdatasums
!------------------------------!


!----------------------------------!
subroutine mcstep
!----------------------------------!
 use systemvariables; implicit none

 integer :: i,s,x,y,s1,s2,s3,s4
 real(8), external :: ran
 real(8) :: r

 do i=1,n
    s=int(ran()*n)
    x=mod(s,l); y=s/l
    s1=spin(mod(x+1,l)+y*l)
    s2=spin(x+mod(y+1,l)*l)
    s3=spin(mod(x-1+l,l)+y*l)
    s4=spin(x+mod(y-1+l,l)*l)
    r=ran()
  if(spin(s)==0) then
    if(r<0.5) then
      if(ran()<(pflip(0,s1+s2+s3+s4)/pflip(1,s1+s2+s3+s4))) then
      spin(s)=1
      end if
    else if(r>0.5) then
      if(ran()<(pflip(0,s1+s2+s3+s4)/pflip(-1,s1+s2+s3+s4))) then
      spin(s)=-1
      end if
    end if
  else if(spin(s)==1) then
    if(r<0.5) then
      if(ran()<(pflip(1,s1+s2+s3+s4)/pflip(0,s1+s2+s3+s4))) then
      spin(s)=0
      end if
    else if(r>0.5) then
      if(ran()<(pflip(1,s1+s2+s3+s4)/pflip(-1,s1+s2+s3+s4))) then
      spin(s)=-1
      end if
    end if
   else if(spin(s)==-1) then
    if(r<0.5) then
      if(ran()<(pflip(-1,s1+s2+s3+s4)/pflip(1,s1+s2+s3+s4))) then
      spin(s)=1
      end if
    else if(r>0.5) then
      if(ran()<(pflip(-1,s1+s2+s3+s4)/pflip(0,s1+s2+s3+s4))) then
      spin(s)=0
      end if
    end if
   end if  
 end do

 end subroutine mcstep
!-------------------------------!



!-------------------------------!
 subroutine initialize(t)
!-------------------------------!
 use systemvariables
 implicit none

 integer :: i,j,ns
 real(8) :: t,r

 integer, allocatable :: seed(:)
 real(8), external :: ran

 n=l*l
 do i=-4,4
    pflip(-1,i)=exp(1.d0*i/t)
    pflip(+1,i)=exp(-1.d0*i/t)
    pflip(0,i)=1
 enddo

 call initran(1)

 allocate (spin(0:n-1))
 do i=0,n-1
    r=ran()
    if(r<1.0/3.0)then
    spin(i)=0
    else if(r<2.0/3.0) then
    spin(i)=1
    else 
    spin(i)=-1
    end if
 enddo

 end subroutine initialize
!-------------------------------!



!----------------------------------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

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
 common/bran64/dmu64,ran64,mul64,add64
      
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

end subroutine initran
