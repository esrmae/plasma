Module solvers
Use shared_data
Implicit None
!
Contains
!
!******************************************************************************
Subroutine trip3(a,b,c,f,beg,dir)
!******************************************************************************
Implicit None
Integer :: i,j,k,dir
Integer, Dimension(3) :: beg
Real(mytype), Dimension(zsize(1),zsize(2),zsize(3)) :: a,b,c,f 
Real(mytype), Dimension(zsize(1),zsize(2)) :: fn 
Real(mytype), Dimension(:,:,:), Allocatable :: s,q
Real(mytype) :: p
!
!       Forward elimination sweep
!
!*******************
If (dir==3) then
!*******************
!
  Allocate(s(zsize(1),zsize(2),zsize(3)),q(zsize(1),zsize(2),zsize(3)))
  s=0.0e0; q=0.0e0
!
  Do j=beg(2),zsize(2); Do i=beg(1),zsize(1)
   fn(i,j) = f(i,j,nz)
   q(i,j,1) = -c(i,j,1)/b(i,j,1)
   f(i,j,1) = f(i,j,1)/b(i,j,1)
   s(i,j,1) = -a(i,j,1)/b(i,j,1)
  End Do; End Do
!
  Do k = 2,nz-1; Do j=beg(2),zsize(2); Do i=beg(1),zsize(1)
    p = 1.0e0/(b(i,j,k)+a(i,j,k)*q(i,j,k-1))
    q(i,j,k) = -c(i,j,k)*p
    f(i,j,k) = (f(i,j,k)-a(i,j,k)*f(i,j,k-1))*p
    s(i,j,k) = -a(i,j,k)*s(i,j,k-1)*p
   End Do; End Do; End Do
!
!       Backward pass
!
  Do j=beg(2),zsize(2); Do i=beg(1),zsize(1)
   f(i,j,nz) = 0.0e0
   s(i,j,nz) = 1.0e0
  End Do; End Do
!
  Do k = nz-1,1,-1; Do j=beg(2),zsize(2); Do i=beg(1),zsize(1)
    s(i,j,k) = s(i,j,k)+q(i,j,k)*s(i,j,k+1)
    f(i,j,k) = f(i,j,k)+q(i,j,k)*f(i,j,k+1)
   End Do; End Do; End Do
!
  Do j=beg(2),zsize(2); Do i=beg(1),zsize(1)
   f(i,j,nz) = (fn(i,j)-c(i,j,nz)*f(i,j,1)-a(i,j,nz)*                                  &
   f(i,j,nz-1))/(c(i,j,nz)*s(i,j,1)+a(i,j,nz)*s(i,j,nz-1)+b(i,j,nz))
  End Do; End Do
!
!       Backward elimination pass
!
  Do k = nz-1,1,-1; Do j=beg(2),zsize(2); Do i=beg(1),zsize(1)
    f(i,j,k) = f(i,j,nz)*s(i,j,k)+f(i,j,k)
   End Do; End Do; End Do
!

  Deallocate(s,q)
!
Else
 Write(*,*) 'Direction not accounted for'
!*******************
End If
!
!
Return
End Subroutine trip3
!
!******************************************************************************
        Subroutine pentay(a,b,c,d,e,r,beg,dir)
!******************************************************************************
Implicit None
Integer :: i,j,k,dir
Integer, Dimension(3) :: beg
Real(mytype), Dimension(:,:,:) :: a,b,c,d,e,r
Real(mytype) :: t1,t2,t3
!
!*******************
If (dir==2) then
!*******************
!
  Do k=beg(3),ysize(3); Do j = 1,ny-2; Do i=beg(1),ysize(1)
    t1 = b(i,j+1,k)/c(i,j,k)
    c(i,j+1,k) = c(i,j+1,k)-t1*d(i,j,k)
    d(i,j+1,k) = d(i,j+1,k)-t1*e(i,j,k)
    r(i,j+1,k) = r(i,j+1,k)-t1*r(i,j,k)
    t2 = a(i,j+2,k)/c(i,j,k)
    b(i,j+2,k) = b(i,j+2,k)-t2*d(i,j,k)
    c(i,j+2,k) = c(i,j+2,k)-t2*e(i,j,k)
    r(i,j+2,k) = r(i,j+2,k)-t2*r(i,j,k)
  End Do; End Do; End Do
!
  Do k=beg(3),ysize(3); Do i=beg(1),ysize(1)
   t3 = b(i,ny,k)/c(i,ny-1,k)
   c(i,ny,k) = c(i,ny,k)-t3*d(i,ny-1,k)
   r(i,ny,k) = r(i,ny,k)-t3*r(i,ny-1,k)
   r(i,ny,k) = r(i,ny,k)/c(i,ny,k)
   r(i,ny-1,k) = (r(i,ny-1,k)-d(i,ny-1,k)*r(i,ny,k))/c(i,ny-1,k)
  End Do; End Do
!
  
  Do k=beg(3),ysize(3); Do j = ny-2,1,-1; Do i=beg(1),ysize(1)
    r(i,j,k) = (r(i,j,k)-d(i,j,k)*r(i,j+1,k)-e(i,j,k)*r(i,j+2,k))/c(i,j,k)
  End Do; End Do; End Do
!
Else
 Write(*,*) 'Direction not accounted for'
!*******************
End If
!
Return
End Subroutine pentay
!
!******************************************************************************
Subroutine dtridy(a,b,c,r,beg,dir)
!******************************************************************************
Implicit None
Integer :: i,j,k,dir
Integer, Dimension(3) :: beg
Real(mytype), Dimension(:,:,:) :: a,b,c,r
Real(mytype) :: t
!
!*******************
If (dir==2) then
!*******************
!
  Do k=beg(3),ysize(3); Do j = beg(2)+1,ny; Do i=beg(1),ysize(1)
    t = a(i,j,k)/b(i,j-1,k)
    b(i,j,k) = b(i,j,k)-c(i,j-1,k)*t
    r(i,j,k) = r(i,j,k)-r(i,j-1,k)*t
  End Do; End Do; End Do
!
  Do k=beg(3),ysize(3); Do i=beg(1),ysize(1)
    r(i,ny,k) = r(i,ny,k)/b(i,ny,k)
  End Do; End Do
!
  Do k=beg(3),ysize(3); Do j = ny-1,beg(2),-1; Do i=beg(1),ysize(1)
    r(i,j,k) = (r(i,j,k)-c(i,j,k)*r(i,j+1,k))/b(i,j,k)
  End Do; End Do; End Do
!
Else
 Write(*,*) 'Direction not accounted for'

!*******************
End If
!
Return
End Subroutine dtridy
!
!******************************************************************************
Subroutine trip1_new(a,b,c,f,beg,dir)
!******************************************************************************
Implicit None
Integer :: i,j,k,dir,ibeg
Integer, Dimension(3) :: beg
Real(mytype), Dimension(:,:,:) :: a,b,c,f
Real(mytype), Dimension(:,:,:), Allocatable  :: s,q
Real(mytype), Dimension(:,:), Allocatable  :: fn
Real(mytype) :: p
!
!*******************
If (dir==1) then
!*******************
!
  ibeg=beg(1)
  Allocate(s(xsize(1),xsize(2),xsize(3)),q(xsize(1),xsize(2),xsize(3)))
  Allocate(fn(xsize(2),xsize(3)))
  s=0.0e0; q=0.0e0; fn=0.0e0
!
!       Forward elimination sweep
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2)
   fn(j,k) = f(nx,j,k)
   q(ibeg,j,k) = -c(ibeg,j,k)/b(ibeg,j,k)
   f(ibeg,j,k) =  f(ibeg,j,k)/b(ibeg,j,k)
   s(ibeg,j,k) = -a(ibeg,j,k)/b(ibeg,j,k)
  End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2); Do i = beg(1)+1,nx-1 
    p = 1.0e0/(b(i,j,k)+a(i,j,k)*q(i-1,j,k))
    q(i,j,k) = -c(i,j,k)*p
    f(i,j,k) = (f(i,j,k)-a(i,j,k)*f(i-1,j,k))*p
    s(i,j,k) = -a(i,j,k)*s(i-1,j,k)*p
  End Do; End Do; End Do
!
!       Backward pass
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2)
    f(nx,j,k) = 0.0e0
    s(nx,j,k) = 1.0e0
  End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2); Do i = nx-1,beg(1),-1
    s(i,j,k) = s(i,j,k)+q(i,j,k)*s(i+1,j,k)
    f(i,j,k) = f(i,j,k)+q(i,j,k)*f(i+1,j,k)
  End Do; End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2)
    f(nx,j,k) = (fn(j,k)-c(nx,j,k)*f(ibeg,j,k)-a(nx,j,k)*f(nx-1,j,k))/  &
                (c(nx,j,k)*s(ibeg,j,k)+a(nx,j,k)*s(nx-1,j,k)+b(nx,j,k))
  End Do; End Do
!
!       Backward elimination pass
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2); Do i = nx-1,beg(1),-1
    f(i,j,k) = f(nx,j,k)*s(i,j,k)+f(i,j,k)
  End Do; End Do; End Do
!
  Deallocate(s,q)
!
Else
 Write(*,*) 'Direction not accounted for'

!*******************
End If
!
!
Return
End Subroutine trip1_new
!
!******************************************************************************
Subroutine penta1(a,b,c,d,e,f,beg,dir)
!******************************************************************************
Implicit None
Integer :: i,j,k,dir
Integer, Dimension(3) :: beg
Real(mytype), Dimension(:,:,:) :: a,b,c,d,e,f
Real(mytype), Dimension(:,:,:), Allocatable  :: r,s,q
Real(mytype) :: p
!
!       Forward elimination sweep
!
!*******************
If (dir==1) then
!*******************
!
  Allocate(r(xsize(1),xsize(2),xsize(3)),s(xsize(1),xsize(2),xsize(3)), &
           q(xsize(1),xsize(2),xsize(3)))
!
  r=0.0e0; s=0.0e0; q=0.0e0
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2)
  q(1,j,k) = -d(1,j,k)/c(1,j,k)
  r(1,j,k) = -e(1,j,k)/c(1,j,k)
  f(1,j,k) =  f(1,j,k)/c(1,j,k)
  s(1,j,k) = -b(1,j,k)/c(1,j,k)
  End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2); Do i = 2,nx-2
    p = 1.0e0/(c(i,j,k)+b(i,j,k)*q(i-1,j,k))
    q(i,j,k) = -(d(i,j,k)+b(i,j,k)*r(i-1,j,k))*p
    r(i,j,k) = -e(i,j,k)*p
    f(i,j,k) = (f(i,j,k)-b(i,j,k)*f(i-1,j,k))*p
    s(i,j,k) = (-b(i,j,k)*s(i-1,j,k)+s(i,j,k))*p
    b(i+1,j,k) = b(i+1,j,k)+a(i+1,j,k)*q(i-1,j,k)
    c(i+1,j,k) = c(i+1,j,k)+a(i+1,j,k)*r(i-1,j,k)
    s(i+1,j,k) = -a(i+1,j,k)*s(i-1,j,k)
    f(i+1,j,k) = f(i+1,j,k)-a(i+1,j,k)*f(i-1,j,k)
  End Do; End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2)
  d(nx-1,j,k) = d(nx-1,j,k)-s(nx-1,j,k)
  s(nx-2,j,k) = s(nx-2,j,k)+r(nx-2,j,k)
  p = 1.0e0/(c(nx-1,j,k)+b(nx-1,j,k)*q(nx-2,j,k))
  s(nx-1,j,k) = -(d(nx-1,j,k)+b(nx-1,j,k)*s(nx-2,j,k))*p
  f(nx-1,j,k) = (f(nx-1,j,k)-b(nx-1,j,k)*f(nx-2,j,k))*p
  b(nx,j,k) = b(nx,j,k)+a(nx,j,k)*q(nx-2,j,k)
  c(nx,j,k) = c(nx,j,k)+a(nx,j,k)*s(nx-2,j,k)
  f(nx,j,k) = f(nx,j,k)-a(nx,j,k)*f(nx-2,j,k)
!
!       Backward pass
!
  s(nx-2,j,k) = s(nx-2,j,k)+q(nx-2,j,k)*s(nx-1,j,k)
  f(nx-2,j,k) = f(nx-2,j,k)+q(nx-2,j,k)*f(nx-1,j,k)
  End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2); Do i = nx-3,1,-1
    s(i,j,k) = q(i,j,k)*s(i+1,j,k)+s(i,j,k)+r(i,j,k)*s(i+2,j,k)
    f(i,j,k) = f(i,j,k)+q(i,j,k)*f(i+1,j,k)+r(i,j,k)*f(i+2,j,k)
  End Do; End Do; End Do
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2)
    f(nx,j,k) = (f(nx,j,k)-d(nx,j,k)*f(1,j,k)-b(nx,j,k)*f(nx-1,j,k))/                         &
                (d(nx,j,k)*s(1,j,k)+b(nx,j,k)*s(nx-1,j,k)+c(nx,j,k))
  End Do; End Do
!
!       Backward elimination pass
!
  Do k=beg(3),xsize(3); Do j=beg(2),xsize(2); Do i = nx-1,1,-1
    f(i,j,k) = f(nx,j,k)*s(i,j,k)+f(i,j,k)
  End Do; End Do; End Do
!
  Deallocate(r,s,q)
!
Else
 Write(*,*) 'Direction not accounted for'

!*******************
End If
!
Return
End Subroutine penta1
!
!******************************************************************************
Subroutine trip1(a,b,c,f,j2)
!******************************************************************************
Implicit None
Integer :: i,j,k,j2,j1,ja,jj
Real(mytype), Dimension(j2) :: a,b,c,f,s,q
Real(mytype) :: fn,p
!
j1=1
ja=j1+1
!
fn = f(j2)
!
!       Forward elimination sweep
!
q(1) = -c(1)/b(1)
f(1) = f(1)/b(1)
s(1) = -a(1)/b(1)
!
Do i = 2,j2-1
p = 1.0e0/(b(i)+a(i)*q(i-1))
q(i) = -c(i)*p
f(i) = (f(i)-a(i)*f(i-1))*p
s(i) = -a(i)*s(i-1)*p
End Do
!
!       Backward pass
!
jj=1+j2
f(j2) = 0.0e0
s(j2) = 1.0e0
!
Do k = 2,j2
i=jj-k
s(i) = s(i)+q(i)*s(i+1)
f(i) = f(i)+q(i)*f(i+1)
End Do
!
f(j2) = (fn-c(j2)*f(1)-a(j2)*f(j2-1))/                    &
(c(j2)*s(1)+a(j2)*s(j2-1)+b(j2))
!
!       Backward elimination pass
!
Do k = 2,j2
i=jj-k
f(i) = f(j2)*s(i)+f(i)
End Do
!
Return
End Subroutine trip1
!
!******************************************************************************
Subroutine ctriby(a,b,c,r,n)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,n2,i1d,k1d
Real(mytype), Dimension(:,:,:) :: a,b,c
Real(mytype) :: t
Complex(mytype), Dimension(:,:,:) :: r
!
  Do k = 1,fft_ysize(3); Do j = 2,n; Do i = 1,fft_ysize(1)
      t = a(i,j,k)/b(i,j-1,k)
      b(i,j,k) = b(i,j,k)-c(i,j-1,k)*t
      r(i,j,k) = r(i,j,k)-r(i,j-1,k)*t
  End Do; End Do; End Do
!
  Do k = 1,fft_ysize(3); Do i = 1,fft_ysize(1)
    r(i,n,k) = r(i,n,k)/b(i,n,k)
  End Do; End Do
!
  Do k = 1,fft_ysize(3); Do j = n-1,1,-1; Do i = 1,fft_ysize(1)
      r(i,j,k) = (r(i,j,k)-c(i,j,k)*r(i,j+1,k))/b(i,j,k)
  End Do; End Do; End Do
!
Return
End Subroutine ctriby
!
!******************************************************************************
Subroutine ctriby_mins(a,b,c,r,n)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,n2,i1d,k1d
Real(mytype), Dimension(:,:,:) :: a,b,c
Real(mytype) :: t
Complex(mytype), Dimension(:,:,:) :: r
!
! For i=1, k=1
!
  i=1; k=1
  n2=n/2
!
  Do j = 2,n2
      t = a(i,j,k)/b(i,j-1,k)
      b(i,j,k) = b(i,j,k)-c(i,j-1,k)*t
      r(i,j,k) = r(i,j,k)-r(i,j-1,k)*t
  End Do
!
  Do j = n-1,n2+1,-1
      t = c(i,j,k)/b(i,j+1,k)
      b(i,j,k) = b(i,j,k)-a(i,j+1,k)*t
      r(i,j,k) = r(i,j,k)-r(i,j+1,k)*t
  End Do

!
  If(n2*2 == n) Then	! n == even
      r(i,n2,k) = (r(i,n2,k)-r(i,n2+1,k))/(b(i,n2,k)-c(i,n2,k)-a(i,n2+1,k)+b(i,n2+1,k))
      r(i,n2+1,k) = -r(i,n2,k)
  Else			! n == odd
      r(i,n2,k) = (r(i,n2,k)-r(i,n2+1,k))/(b(i,n2,k)-a(i,n2+1,k))
      r(i,n2+1,k) = 0.0e0
  End If
!
  Do j = n2-1,1,-1
      r(i,j,k) = (r(i,j,k)-c(i,j,k)*r(i,j+1,k))/b(i,j,k)
  End Do
!
  Do j = n2+2,n
      r(i,j,k) = (r(i,j,k)-a(i,j,k)*r(i,j-1,k))/b(i,j,k)
  End Do
!
! For rest of k=1
! 
  k = 1
!
  Do j = 2,n; Do i = 2,fft_ysize(1)
      t = a(i,j,k)/b(i,j-1,k)
      b(i,j,k) = b(i,j,k)-c(i,j-1,k)*t
      r(i,j,k) = r(i,j,k)-r(i,j-1,k)*t
  End Do; End Do
!
  Do i = 2,fft_ysize(1)
    r(i,n,k) = r(i,n,k)/b(i,n,k)
  End Do
!
  Do j = n-1,1,-1; Do i = 2,fft_ysize(1)
      r(i,j,k) = (r(i,j,k)-c(i,j,k)*r(i,j+1,k))/b(i,j,k)
  End Do; End Do
!
! For rest of k>1
!
  Do k = 2,fft_ysize(3); Do j = 2,n; Do i = 1,fft_ysize(1)
      t = a(i,j,k)/b(i,j-1,k)
      b(i,j,k) = b(i,j,k)-c(i,j-1,k)*t
      r(i,j,k) = r(i,j,k)-r(i,j-1,k)*t
  End Do; End Do; End Do
!
  Do k = 2,fft_ysize(3); Do i = 1,fft_ysize(1)
    r(i,n,k) = r(i,n,k)/b(i,n,k)
  End Do; End Do
!
  Do k = 2,fft_ysize(3); Do j = n-1,1,-1; Do i = 1,fft_ysize(1)
      r(i,j,k) = (r(i,j,k)-c(i,j,k)*r(i,j+1,k))/b(i,j,k)
  End Do; End Do; End Do
!
Return
End Subroutine ctriby_mins
!
!******************************************************************************
Subroutine ctribl(a,b,c,r,n) 
! for boundary layer
!******************************************************************************
Implicit None
Integer :: i,j,k,n,n2
Real(mytype), Dimension(ny) :: a,b,c
Real(mytype) :: t
Complex(mytype), Dimension(ny) :: r
!
!
n2=n-1
Do j = 2,n2
t = a(j)/b(j-1)
b(j) = b(j)-c(j-1)*t
r(j) = r(j)-r(j-1)*t
End Do
!
r(n2) = (r(n2)-r(n2+1))/(b(n2)-a(n2+1))
r(n2+1) = 0.0e0
!
Do k = 1,n2-1
j=n2-k
r(j) = (r(j)-c(j)*r(j+1))/b(j)
End Do
!
Return
End Subroutine ctribl
!
!******************************************************************************
Subroutine ctri2j(a,b,c,r,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m
Real(mytype), Dimension(2*nx,ny) :: a,b,c
Real(mytype), Dimension(2*nx) :: t
Complex(mytype), Dimension(2*nx,ny) :: r
!
!------ to solve tri-diagonal system of eq. from Joel's' book
!
Do j = 2,n
Do i = 1,m
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
End Do
End Do
!
Do i = 1,m
r(i,n) = r(i,n)/b(i,n)
End Do
!
Do k = 1,n-1
Do i = 1,m
j=n-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
End Do
!
Return
End Subroutine ctri2j
!
!******************************************************************************
Subroutine ctri2y(a,b,c,r,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m,n2,j2
Real(mytype), Dimension(2*nx,ny) :: a,b,c
Real(mytype), Dimension(2*nx) :: t
Complex(mytype), Dimension(2*nx,ny) :: r
!
Do j = 2,n
Do i = 2,m
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
End Do
End Do
!
Do i = 2,m
r(i,n) = r(i,n)/b(i,n)
End Do
!
Do k = 1,n-1
Do i = 2,m
j=n-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
End Do
!
i=1
n2=n/2
!
Do j = 2,n2
j2=n+1-j
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
r(i,j2) = r(i,j2)-r(i,j2+1)*t(i)
b(i,j2) = b(i,j)
End Do
!
r(i,n2) = 0.5e0*(r(i,n2)-r(i,n2+1))/(b(i,n2)-c(i,n2))
r(i,n2+1) = -r(i,n2)
!
Do k = 1,n/2-1
j=n2-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
!
Do j = n/2+2,n
r(i,j) = (r(i,j)-a(i,j)*r(i,j-1))/b(i,j)
End Do
!
Return
End Subroutine ctri2y
!
!******************************************************************************
Subroutine ctrib0(a,b,c,r,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m
Real(mytype), Dimension(nx,ny) :: a,b,c
Real(mytype), Dimension(nx) :: t
Complex(mytype), Dimension(nx,ny) :: r
!
Do j = 2,n
Do i = 1,m
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
End Do
End Do
!
Do i = 1,m
r(i,n) = r(i,n)/b(i,n)
End Do
!
Do k = 1,n-1
Do i = 1,m
j=n-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
End Do
!
i=1
!
Do j = n,1
r(i,j) = r(i,j)-r(i,1)
End Do
!
Return
End Subroutine ctrib0
!
!******************************************************************************
Subroutine ctrib2(a,b,c,r,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m
Real(mytype), Dimension(2*nx,ny) :: a,b,c
Real(mytype), Dimension(2*nx) :: t
Complex(mytype), Dimension(2*nx,ny) :: r
!
Do j = 2,n
Do i = 1,m
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
End Do
End Do
!
Do i = 1,m
r(i,n) = r(i,n)/b(i,n)
End Do
!
Do k = 1,n-1
Do i = 1,m
j=n-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
End Do
!
Return
End Subroutine ctrib2
!
!******************************************************************************
Subroutine ctri2d(a,b,c,r,u,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m
Real(mytype), Dimension(2*nx,ny) :: a,b,c,gam
Real(mytype), Dimension(2*nx) :: t,bet
Complex(mytype), Dimension(2*nx,ny) :: r,u
!
Do i = 1,m
bet(i) = b(i,1)
u(i,1) = r(i,1)/bet(i)
End Do
!
Do j = 2,n
Do i = 1,m
gam(i,j) = c(i,j-1)/bet(i)
bet(i) = b(i,j)-a(i,j)*gam(i,j)
u(i,j) = (r(i,j)-a(i,j)*u(i,j-1))/bet(i)
End Do
End Do
!
Do j = n-1,1,-1
Do i = 1,m
u(i,j) = u(i,j)-gam(i,j+1)*u(i,j+1)
End Do
End Do
!
Return
End Subroutine ctri2d
!
!******************************************************************************
Subroutine ctridx(a,b,c,r,ibeg,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m,ibeg
Real(mytype), Dimension(nx,ny) :: a,b,c
Real(mytype), Dimension(ny) :: t
Complex(mytype), Dimension(nx,ny) :: r
!
Do i = ibeg+1,n
Do j = 1,m
t(j) = a(i,j)/b(i-1,j)
b(i,j) = b(i,j)-c(i-1,j)*t(j)
r(i,j) = r(i,j)-r(i-1,j)*t(j)
End Do
End Do
!
Do j = 1,m
r(n,j) = r(n,j)/b(n,j)
End Do
!
Do k = 1,n-ibeg
i=n-k
Do j = 1,m
r(i,j) = (r(i,j)-c(i,j)*r(i+1,j))/b(i,j)
End Do
End Do
!
Return
End Subroutine ctridx
!
!******************************************************************************
Subroutine ctrip1(a,b,c,f,j2,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m,j2,j1,jj,ja
Real(mytype), Dimension(nx,ny) :: a,b,c,s,q
Real(mytype), Dimension(ny) :: p
Complex(mytype), Dimension(nx,ny) :: f
Complex(mytype), Dimension(ny) :: fn
!
j1=1
ja=j1+1
!
Do j = 1,m
fn(j) = f(j2,j)
End Do
!
!       Forward elimination sweep
!
Do j = 1,m
q(1,j) = -c(1,j)/b(1,j)
f(1,j) = f(1,j)/b(1,j)
s(1,j) = -a(1,j)/b(1,j)
End Do
!
Do i = 2,j2-1
Do j = 1,m
p(j) = 1.0e0/(b(i,j)+a(i,j)*q(i-1,j))
q(i,j) = -c(i,j)*p(j)
f(i,j) = (f(i,j)-a(i,j)*f(i-1,j))*p(j)
s(i,j) = -a(i,j)*s(i-1,j)*p(j)
End Do
End Do
!
!       Backward pass
!
jj=1+j2
Do j = 1,m
f(j2,j) = 0.0e0
s(j2,j) = 1.0e0
End Do
!
Do k = 2,j2
i=jj-k
Do j = 1,m
s(i,j) = s(i,j)+q(i,j)*s(i+1,j)
f(i,j) = f(i,j)+q(i,j)*f(i+1,j)
End Do
End Do
!
Do j = 1,m
f(j2,j) = (fn(j)-c(j2,j)*f(1,j)-a(j2,j)*f(j2-1,j))/                    &
(c(j2,j)*s(1,j)+a(j2,j)*s(j2-1,j)+b(j2,j))
End Do
!
!       Backward elimination pass
!
Do k = 2,j2
i=jj-k
Do j = 1,m
f(i,j) = f(j2,j)*s(i,j)+f(i,j)
End Do
End Do
!
Return
End Subroutine ctrip1
!
!******************************************************************************
Subroutine ctrip2(a,b,c,r,n,m)
!******************************************************************************
Implicit None
Integer :: i,j,k,n,m,j2,j1,jj,ja
Real(mytype), Dimension(nx,ny) :: a,b,c
Real(mytype), Dimension(ny) :: t
Complex(mytype), Dimension(nx,ny) :: r
!
Do j = 2,n
Do i = 1,m
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
End Do
End Do
!
Do i = 1,m
r(i,n) = r(i,n)/b(i,n)
End Do
!
Do k = 1,n-1
Do i = 1,m
j=n-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
End Do
!
Return
End Subroutine ctrip2
!
!******************************************************************************
Subroutine tridy(a,b,c,r,n,m)
!******************************************************************************
!	This subroutine solvers triadiagonal system of equations in case of 2-point 
!	formulation of gcdy_f in u and w mom equatino in SOLVEY in  COMPUTE_LEFT.f90
Implicit None
Integer :: i,j,k,n,m
Real(mytype), Dimension(m,n) :: a,b,c
Real(mytype), Dimension(m) :: t
Real(mytype), Dimension(m,n) :: r
!
Do j = 2,n; Do i = 1,m
t(i) = a(i,j)/b(i,j-1)
b(i,j) = b(i,j)-c(i,j-1)*t(i)
r(i,j) = r(i,j)-r(i,j-1)*t(i)
End Do; End Do
!
Do i = 1,m
r(i,n) = r(i,n)/b(i,n)
End Do
!
Do k = 1,n-1; Do i = 1,m
j=n-k
r(i,j) = (r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do; End Do
!
Return
End Subroutine tridy
!
!******************************************************************************
Subroutine penta_dir(a,b,c,d,e,r,n,m,ibeg)
!******************************************************************************
!       (1-A2)U*=U*: N=NY(length of array), M=NX(number of arrays), IBEG=2(U) or 1(W)
!       (1-A3)U**=RHS
Implicit none
Integer :: i,j,k,ibeg,n,m
Real(mytype), Dimension(m,n) :: a,b,c,d,e,r
Real(mytype), Dimension(m) :: t1,t2,t3
!
Do j=1,n-2
Do i=ibeg,m
t1(i)=b(i,j+1)/c(i,j)
c(i,j+1)=c(i,j+1)-t1(i)*d(i,j)
d(i,j+1)=d(i,j+1)-t1(i)*e(i,j)
r(i,j+1)=r(i,j+1)-t1(i)*r(i,j)
t2(i)=a(i,j+2)/c(i,j)
b(i,j+2)=b(i,j+2)-t2(i)*d(i,j)
c(i,j+2)=c(i,j+2)-t2(i)*e(i,j)
r(i,j+2)=r(i,j+2)-t2(i)*r(i,j)
End Do
End Do
!
Do i=ibeg,m
t3(i)=b(i,n)/c(i,n-1)
c(i,n)=c(i,n)-t3(i)*d(i,n-1)
r(i,n)=r(i,n)-t3(i)*r(i,n-1)
r(i,n)=r(i,n)/c(i,n)
r(i,n-1)=(r(i,n-1)-d(i,n-1)*r(i,n))/c(i,n-1)
End Do
!
Do k=1,n-2
j=n-k-1
Do i=ibeg,m
r(i,j)=(r(i,j)-d(i,j)*r(i,j+1)-e(i,j)*r(i,j+2))/c(i,j)
End Do
End Do
!
Return
End Subroutine penta_dir
!
!******************************************************************************
        Subroutine penta3(a,b,c,d,e,f,n,m,jbeg)
!******************************************************************************
Implicit none
Integer :: i,j,k,jbeg,n,m,j1
Real(mytype), Dimension(nx,nz) :: a,b,c,d,e,f,r,s,q
Real(mytype), Dimension(nx) :: p
!
j1=1
Do i=1,m
Do j=1,n
s(i,j)=0.0e0
End Do
End Do
!
!       Forward elimination sweep
!
Do i=jbeg,m
q(i,1)=-d(i,1)/c(i,1)
r(i,1)=-e(i,1)/c(i,1)
f(i,1)=f(i,1)/c(i,1)
s(i,1)=-b(i,1)/c(i,1)
End Do
!
Do j=j1+1,n-2
Do i=jbeg,m
p(i)=1.0e0/(c(i,j)+b(i,j)*q(i,j-1))
q(i,j)=-(d(i,j)+b(i,j)*r(i,j-1))*p(i)
r(i,j)=-e(i,j)*p(i)
f(i,j)=(f(i,j)-b(i,j)*f(i,j-1))*p(i)
s(i,j)=(-b(i,j)*s(i,j-1)+s(i,j))*p(i)
b(i,j+1)=b(i,j+1)+a(i,j+1)*q(i,j-1)
c(i,j+1)=c(i,j+1)+a(i,j+1)*r(i,j-1)
s(i,j+1)=-a(i,j+1)*s(i,j-1)
f(i,j+1)=f(i,j+1)-a(i,j+1)*f(i,j-1)
End Do
End Do
!
Do i=jbeg,m
d(i,n-1)=d(i,n-1)-s(i,n-1)
s(i,n-2)=s(i,n-2)+r(i,n-2)
p(i)=1.0e0/(c(i,n-1)+b(i,n-1)*q(i,n-2))
s(i,n-1)=-(d(i,n-1)+b(i,n-1)*s(i,n-2))*p(i)
f(i,n-1)=(f(i,n-1)-b(i,n-1)*f(i,n-2))*p(i)
b(i,n)=b(i,n)+a(i,n)*q(i,n-2)
c(i,n)=c(i,n)+a(i,n)*s(i,n-2)
f(i,n)=f(i,n)-a(i,n)*f(i,n-2)
!
!       Backward pass
!
s(i,n-2)=s(i,n-2)+q(i,n-2)*s(i,n-1)
f(i,n-2)=f(i,n-2)+q(i,n-2)*f(i,n-1)
End Do
!
Do j=n-3,1,-1
Do i=jbeg,m
s(i,j)=q(i,j)*s(i,j+1)+s(i,j)+r(i,j)*s(i,j+2)
f(i,j)=f(i,j)+q(i,j)*f(i,j+1)+r(i,j)*f(i,j+2)
End Do
End Do
!
Do i=jbeg,m
f(i,n)=(f(i,n)-d(i,n)*f(i,1)-b(i,n)*f(i,n-1))/                         &
(d(i,n)*s(i,1)+b(i,n)*s(i,n-1)+c(i,n))
End Do
!
!       Backward elimination pass
!
Do j=n-1,1,-1
Do i=jbeg,m
f(i,j)=f(n,j)*s(i,j)+f(i,j)
End Do
End Do
!
Return
End Subroutine penta3
!
!******************************************************************************
Subroutine dtridyz(a,b,c,r,jbeg,n,m)
!******************************************************************************
Implicit none
Integer :: i,j,k,jbeg,n,m
Real(mytype), Dimension(m,n) :: a,b,c,r
Real(mytype), Dimension(m) :: t
!
Do j=jbeg+1,n
Do i=1,m
t(i)=a(i,j)/b(i,j-1)
b(i,j)=b(i,j)-c(i,j-1)*t(i)
r(i,j)=r(i,j)-r(i,j-1)*t(i)
End Do
End Do
!
Do i=1,m
r(i,n)=r(i,n)/b(i,n)
End Do
!
Do k=1,n-jbeg
Do i=1,m
j=n-k
r(i,j)=(r(i,j)-c(i,j)*r(i,j+1))/b(i,j)
End Do
End Do
!
Return
End Subroutine dtridyz
!
End Module solvers


