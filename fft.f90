Module fft_solve
Use shared_data
Use solvers
Use decomp_2d
Use decomp_2d_fft
Use boundary_conditions
Implicit None
!
Contains
!
!******************************************************************************
Subroutine dft
!******************************************************************************
!       FFT For periodic B.c
Implicit None
Integer :: i,j,k,i1d,k1d
Real(mytype) :: dxx,dzz
Real(mytype), Dimension(:), Allocatable :: ak1,ak3,ajp,ajc,ajm
Real(mytype), Dimension(:,:,:), Allocatable :: ajjc,ajjm,ajjp
Complex(mytype), Dimension(:,:,:), Allocatable :: ccapz,ccapy
!
Allocate (ak1(nx+1),ak3(nz+1),ajp(ny),ajc(ny),ajm(ny), &
          ajjc(fft_ysize(1),fft_ysize(2),fft_ysize(3)),ajjm(fft_ysize(1),fft_ysize(2),fft_ysize(3)), &
          ajjp(fft_ysize(1),fft_ysize(2),fft_ysize(3)),ccapy(fft_ysize(1),fft_ysize(2),fft_ysize(3)), &
          ccapz(fft_zsize(1),fft_zsize(2),fft_zsize(3)))
ak1=0.0e0; ak3=0.0e0; ajp=0.0e0; ajc=0.0e0; ajm=0.0e0
ajjc=0.0e0; ajjm=0.0e0; ajjp=0.0e0; ccapy=0.0e0; ccapz=0.0e0
!
dxx=dx(2)*dx(2)
dzz=dz(2)*dz(2)
!
Do j = 1,ny
   ajm(j) = 2.0e0/dy(j)/(dy(j-1)+dy(j))
   ajp(j) = 2.0e0/dy(j)/(dy(j+1)+dy(j))
End Do
   ajm(1) = 0.0e0
   ajp(ny) = 0.0e0
   ajc = -ajp-ajm
!
!       MODIFIED WAVE NUMBER DEFINITION
Do k = 1,nz/2+1
   ak3(k) = 2.0e0*(1.0e0-cos(real(k-1,kind=mytype)*2.0e0*pi/real(nz,kind=mytype)))/dzz
End Do
!
Do i = 1,nx/2+1
   ak1(i) = 2.0e0*(1.0e0-cos(real(i-1,kind=mytype)*2.0e0*pi/real(nx,kind=mytype)))/dxx
End Do
!
Do k = nz,nz/2+2,-1
   ak3(k) = ak3(nz+2-k)
End Do
!
!	Set Rhs for fft
!Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
!   rcap(i,j,k) = farray(i,j,k)
!End Do; End Do; End Do
!
!********************
!       FORWARD FFT
!********************
!	FFTW OPTION
 Call decomp_2d_fft_3d(farray,ccapz)
!
 Call transpose_z_to_y(ccapz,ccapy,fft)
!
  ccapy=ccapy/real(nx*nz,kind=mytype)
!
!	TDMA in wall normal direction
Do k = 1,fft_ysize(3)
k1d = fft_ystart(3)-1+k
Do i = 1,fft_ysize(1)
i1d = fft_ystart(1)-1+i
Do j = 1,fft_ysize(2)
  ajjm(i,j,k) = ajm(j)
  ajjp(i,j,k) = ajp(j)
  ajjc(i,j,k) = ajc(j)-ak3(k1d)-ak1(i1d)
End Do; End Do; End Do

!
If (if2d == 0) Then
!
  If ((fft_ystart(1)==1).and.(fft_ystart(3)==1)) then
    Call ctriby_mins(ajjm,ajjc,ajjp,ccapy,ny)
  Else
    Call ctriby(ajjm,ajjc,ajjp,ccapy,ny)
  End If
!
Else If (if2d_noy == 1) Then

  ajjm = 0.0e0; ajjp = 0.0e0
  Do k = 1,fft_ysize(3); Do j = 1,fft_ysize(2); Do i = 1,fft_ysize(1)
    ajjc(i,j,k) = -ak3(k)-ak1(i)
    If (i == 1 .And. k == 1) Then
      ccapy(i,j,k) = 0.0e0
    Else 
      ccapy(i,j,k) = ccapy(i,j,k)/ajjc(i,j,k)
    End If
  End Do; End Do; End Do
!   
End If
!
 Call transpose_y_to_z(ccapy,ccapz,fft)
!
!**************************************
!      INVERSE FFT.
!************************************
!	FFTW OPTION
!
  Call decomp_2d_fft_3d(ccapz,farray)
!
!
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
   varray(i,j,k) = farray(i,j,k)
End Do; End Do; End Do
!
!
!	B.C.s
  Call update_halo(varray,2)
  If (xmin(2)==1) varray(:,0,:)=0.0e0
  If (xmax(2)==1) varray(:,xsize(2)+1,:)=0.0e0

Deallocate (ak1,ak3,ajp,ajc,ajm,ajjc, &
		ajjm,ajjp,ccapy,ccapz)

Return
End Subroutine dft
!
End Module fft_solve
