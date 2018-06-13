

program randomwalk

integer:: b,Nw,i,j,k,x,n
real, allocatable :: c(:,:)
integer, dimension(0:63):: s
real, dimension(0:63):: d
real, dimension(-100:100) :: p
integer*8::r,r0
character(len=2):: bit
real(8)::f

open (1,file="read.in")
read(1,*) Nw,n,r0
close(1)
allocate(c(-n:n,0:63))
p=0
c=0
s=0
d=0
do j=1,Nw
 call walk(s,n,r0)
 do k=0,63
    x=s(k)
    c(x,k)=c(x,k)+1
 end do
end do


do i=-n+2,n-2,2
   call lnfac(n,f)
   p(i)=f
   call lnfac((n+i)/2,f)
   p(i)=p(i)-f
   call lnfac((n-i)/2,f)
   p(i)=p(i)-f-n*log(2.)
   p(i)=exp(p(i))
end do
p(n)=-n*log(2.0)
p(n)=exp(p(n))
p(-n)=p(n)
open(2,file="d.dat") 
do b=0,63
  d(b)=0
  write(bit,'(i2)')b
  open(3,file='p'//trim(adjustl(bit))//'.dat')
  do i=-n,n,2
  d(b)=d(b)+((c(i,b)/Nw)-p(i))**2
  write(3,*), i, c(i,b)/Nw, p(i)
  end do
  write(2,*), b, sqrt(Nw*d(b))
end do
close(2)
close(3)

deallocate (c)

end program randomwalk


subroutine ran(r0)

integer*8::a=2862933555777941757_8, c=1013904243_8,r0,r
r=a*r0+c
r0=r
return
end

subroutine walk(s,n,r0)
integer:: n,b,i
integer*8::r0,r
integer, dimension(0:63):: s
s=0
 do i=1,n
   call ran(r0)
    do b=0,63
       if( btest(r0,b) ) then
          s(b) = s(b)+1
       else
       s(b) = s(b)-1
       end if
    end do
  end do
return
end 


subroutine lnfac(n,f)

integer::n,i
real(8)::f
f=0
do i=1,n
 f=f+log(real(i))
end do
return
end
