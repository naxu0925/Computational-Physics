!This program uses variational method to solve Schrodinger equation

!-----------------------------------!
 module system
!-----------------------------------!
 implicit none
 real(8), parameter:: lx=5.d0, ly=10.d0
 real(8), parameter:: x0=2.d0, y0=2.d0
 real(8), parameter:: v=0.1
 real(8), parameter:: pi=3.141592653589793238d0, beta=13.11d0
 end module system
!-----------------------------------!




!-----------------------------------!
 program variational
!-----------------------------------!
 use system; implicit none

 real(8):: wavefunction,energy,function,velement
 integer:: nbasis,nx,ny
 integer:: i=0,j=0,k=0,l=0,m=0,n=0
 real(8), allocatable :: eig(:)
 real(8), allocatable :: ham(:,:)
 integer:: kx,ky,px,py

 print*,'Number of basis states,nx and ny'
 read(*,*) nx,ny
 nbasis=nx*ny
 allocate(eig(nbasis))
 allocate(ham(nbasis,nbasis))
 ham=0.d0
 open(10, file='mat.txt')
 write(10,*), nbasis
 do j=1,ny
  do i=1,nx
   do n=1,ny
    do m=1,nx
       k=(i-1)*nx+j
       l=(m-1)*nx+n
       if((i.eq.m).and.(j.eq.n)) then
         ham(k,k)=energy(i,j)+velement(j,i,n,m)
       else 
         ham(l,k)=velement(j,i,n,m)
       endif
     enddo
    enddo
  enddo
 enddo
 do i=1,nbasis
    write(10,*) ' '
    do j=1,nbasis
       write(10,*) ham(i,j)/beta
    enddo
 enddo
 close(10)
 deallocate (eig)
 deallocate (ham)

 end program variational
!--------------------------------------!




!--------------------------------------!
 real(8) function energy(kx,ky)
!---------------------------------------!
 use system;  implicit none
 integer:: kx,ky
 energy=pi**2/2*((kx/lx)**2+(ky/ly)**2)
 
 end function energy
!----------------------------------------!




!----------------------------------------!
 real(8) function  wavefunction(kx,ky,x,y)
!-----------------------------------------!
 use system; implicit none
 integer:: kx,ky
 real(8):: x,y
 wavefunction=2/sqrt(lx*ly)*sin(kx*pi*x/lx)*sin(ky*pi*y/ly)
 
 end function wavefunction
!----------------------------------------------!



!---------------------------------------!
real(8) function potential(x,y)
!----------------------------------------!
 use system; implicit none
 real(8)::x,y
 
 if(((x.gt.2.0).and.(x.lt.3.0)).and.(((y.gt.2.0).and. (y.lt.2.5)).or.((y.gt.7.5).and.(y.lt.8.0)))) then
    potential=v*beta
 else
   potential=0.0d0
 endif

end function potential
!---------------------------------------------!





!---------------------------------------------!
 real(8) function velement(ky,kx,py,px)
!---------------------------------------------!
 use system; implicit none
 integer:: kx,ky,px,py
 real(8):: a,b,c,d,a_a,c_c,c_d,a_b,x,y
 a=kx*pi/real(lx)
 b=px*pi/real(lx)
 c=ky*pi/real(ly)
 d=py*pi/real(ly)
 x=2*(c-d)
 y=2*(c+d)
 a_a=0.5+(1/(4*a))*(sin(4*a)-sin(6*a))
 c_c=0.5+(1/(4*c))*(sin(15*c)+sin(4*c)-sin(16*c)-sin(5*c));
 c_d=(sin((c-d)*2.5)-sin((c-d)*2)+sin((c-d)*8)-sin((c-d)*7.5))/x-(sin((c+d)*8)-sin((c+d)*7.5)+sin((c+d)*2.5)-sin((c+d)*2))/y
 a_b=(sin((a-b)*3.0)-sin((a-b)*2.0))/(2*a-2*b)-(sin((a+b)*3.0)-sin((a+b)*2.0))/(2*a+2*b)
 if((a.eq.b).and.(c.eq.d)) then
    velement=(4/(lx*ly))*v*beta*(a_a*c_c)
 else if((a.eq.b).and.(c.ne.d)) then
    velement=(4/(lx*ly))*v*beta*a_a*c_d
 else if((a.ne.b).and.(c.eq.d)) then
    velement=(4/(lx*ly))*v*beta*c_c*a_b
 else
    velement=(4/(lx*ly))*v*beta*a_b*c_d
 endif

end function velement
!-----------------------------------------------------!


