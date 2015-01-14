Module poisson_solver
Use shared_data
Use fft_solve
Use boundary_conditions
Implicit None
Contains

!******************************************************************************
Subroutine compute_poisson 
!******************************************************************************
Implicit None

Call poisson_rhs
Call cmpspe1 

Return
End Subroutine compute_poisson
!
!******************************************************************************
Subroutine poisson_rhs
!******************************************************************************
!------ This procedure computes the righthand side vector.
Implicit None
Integer :: i,j,k
Real(mytype) :: myfsum,fsum
!
if (if2d_noy == 1) vhat1=0.0e0
 Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
 farray(i,j,k) =                                 &
 ((uhat1(i+1,j,k)-uhat1(i,j,k))/dx(xstart(1)-1+i)              &
 +(vhat1(i,j+1,k)-vhat1(i,j,k))/dy(xstart(2)-1+j)              &
 +(what1(i,j,k+1)-what1(i,j,k))/dz(xstart(3)-1+k))/dt        
 End Do; End Do; End Do
!
!       Compatibility condition test.
!
myfsum = 0.0e0; fsum = 0.0e0
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
myfsum = myfsum + farray(i,j,k)
End Do; End Do; End Do
!
  Call MPI_Allreduce(myfsum,fsum,1,real_type,MPI_SUM,MPI_COMM_CART,ierror)
!
If (nrank == 0) Write(21,1048) ii,fsum
 1048   Format(' II=',I6,' FSUM=',E20.9)
!
Return
End Subroutine poisson_rhs
!
!******************************************************************************
Subroutine cmpspe1
!******************************************************************************
Implicit None
Integer :: i,j,k
!
!	FFT solver
!
   Call dft 
!
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
p1(i,j,k) = p1(i,j,k)+varray(i,j,k)
End Do; End Do; End Do
Call bc_pressure
!
Return
End Subroutine cmpspe1
!
!******************************************************************************
Subroutine bc_pressure
!******************************************************************************
Implicit none
Integer :: i,j,k
!
Call update_halo(p1,2)
!If (jconbc==1) then
!  If(xmin(2)==1) then
!    Do k=-1,xsize(3)+2; Do i=-1,xsize(1)+2
!      p1(i,0,k)=p1(i,1,k)
!      p1(i,-1,k)=p1(i,1,k)
!    End Do; End Do
!  End If
!  If(xmax(2)==1) then
!    Do k=-1,xsize(3)+2; Do i=-1,xsize(1)+2
!      p1(i,xsize(2)+1,k)=p1(i,xsize(2),k)
!      p1(i,xsize(2)+2,k)=p1(i,xsize(2),k)
!    End Do; End Do
!  End If
!End If
!
Return
End Subroutine bc_pressure
!
End Module poisson_solver
