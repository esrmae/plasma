Module profiles
Use shared_data
Implicit None
Contains
!
!******************************************************************************
Subroutine parabolic
!******************************************************************************
!       This subroutine determines the initial velcoity profiles given to the flow field at start.
Implicit None
Integer :: i,j,k,j1d,i1d,k1d
Real :: utemp
!
u1 = 0.0e0; v1 = 0.0e0; w1 = 0.0e0
!
!   Parabolic laminar solution for u 
!
utemp=1.5e0
Do k = 1,xsize(3); Do i = -1,xsize(1)+2; DO  j=1,xsize(2)+1
  j1d = xstart(2)-1+j
  u1(i,j,k) = utemp*yc(j1d)*(height-yc(j1d))
!  k1d = xstart(3)-1+k
!  u1(i,j,k) = utemp*yc(j1d)*(height-yc(j1d))*sin(2.0e0*pi/width*z(k1d))
!  u1(i,j,k) = utemp*yc(j1d)*(height-yc(j1d))*sin(2.0e0*pi/length*x(i))
End Do; End Do; End Do
!
Return
End Subroutine parabolic
!
End Module profiles
