Module mass_calc
Use shared_data
Use decomp_2d
Use mpi
Implicit None
Contains
!
!******************************************************************************
Subroutine rmsfl
!******************************************************************************
!-------This procedure calculates mass flux at inlet and outlet
!       and their difference.
Implicit None
Integer :: i,j,k,i1d,j1d,k1d
Real(mytype) :: rmsfli_dum, rmsflm_dum, rmsflo_dum, rmsflc_dum, rmsflb_dum
!
rmsfli=0.0e0; rmsflm=0.0e0; rmsflo=0.0e0; rmsflc=0.0e0; rmsflb=0.0e0
rmsfli_dum=0.0e0; rmsflm_dum=0.0e0; rmsflo_dum=0.0e0; rmsflc_dum=0.0e0
rmsflb_dum=0.0e0
!
! Finds flux at inlet, middle and outlet 
Do k = 1,xsize(3); Do j = 1,xsize(2)
  j1d = xstart(2)-1+j
  k1d = xstart(3)-1+k
!
  rmsfli_dum=rmsfli_dum+u1(1,j,k)*dy(j1d)*dz(k1d)*roh
  rmsflm_dum=rmsflm_dum+u1(nx/2+1,j,k)*dy(j1d)*dz(k1d)*roh
  rmsflo_dum=rmsflo_dum+u1(nx+1,j,k)*dy(j1d)*dz(k1d)*roh
  rmsflc_dum=rmsflc_dum+xinlet(j,k,1,2,2)*dy(j1d)*dz(k1d)*roh
End Do; End Do
!
! Sums over results on all processes
  Call MPI_Allreduce(rmsfli_dum,rmsfli,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  Call MPI_Allreduce(rmsflm_dum,rmsflm,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  Call MPI_Allreduce(rmsflo_dum,rmsflo,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  Call MPI_Allreduce(rmsflc_dum,rmsflc,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
!
! Finds difference in flux  MUST BE SORTED OUT
Do k = 1,xsize(3); Do i = 1,xsize(1)
  j1d = xstart(2)-1+j
  If (coord(1) == 0) rmsflb_dum=rmsflb_dum+yinlet1(i,k,2,1,2)*dx(i)*dz(k1d)*roh
  If (coord(1) == p_row-1) rmsflb_dum=rmsflb_dum-yinlet1(i,k,2,2,2)*dx(i)*dz(k1d)*roh
End Do; End Do
!
! Sums over results on all processes
  Call MPI_Allreduce(rmsflb_dum,rmsflb,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
!
dmsfl=rmsfli+rmsflb-rmsflo
!
Return
End Subroutine rmsfl
!
!******************************************************************************
Subroutine compute_dpdx_mean
!******************************************************************************
!	This procedure ensures the mass consevation of the simulation by calculating
!       pressure gradient for next time step which conserve mass
Implicit None
Integer i,j,k,iblock,ib,j1d,k1d
Real(mytype) :: sum1,dpdx_old
!
iblock = 0; ib=iblock+1
sum1 = 0.0e0
dpdx_old=dpdx_mean
!**************** must be in X pencil *******************
!
If (xmax(2) == 1) then
  j=xsize(2)+1
  j1d=xstart(2)-1+j
  Do k = 1,xsize(3); Do i = 1,xsize(1)
    k1d=xstart(3)-1+k
    sum1=sum1 + roh*nu*(gcdy_f(j1d,3,1,ib)*(u1(i,j-1,k)+u1(i+1,j-1,k))*0.5e0 +     &
		        gcdy_f(j1d,4,1,ib)*(u1(i,j-2,k)+u1(i+1,j-2,k))*0.5e0)*dx(i)*dz(k1d)
  End Do; End Do
End If
!
If (xmin(2) == 1) then
  j=1
  j1d=xstart(2)-1+j
  Do k = 1,xsize(3); Do i = 1,xsize(1)
    k1d=xstart(3)-1+k
    sum1=sum1 - roh*nu*(gcdy_f(j1d,1,1,ib)*(u1(i,j+1,k)+u1(i+1,j+1,k))*0.5e0 +     &
		        gcdy_f(j1d,2,1,ib)*(u1(i,j,k)  +u1(i+1,j,k)  )*0.5e0)*dx(i)*dz(k1d)  
  End Do; End Do
End If
!
!	Sum on all processes
  Call MPI_Allreduce(sum1,dpmean,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
!
dpmean=dpmean/2.0e0/width
dpdx_mean=dpmean/length
!
!
If (nrank == 0) Write(21,'(I6,E15.8,A,E15.8,A,E15.8)') ii,time,'  dpdx_old=',dpdx_old,'  dpdx_mean=',dpdx_mean
!
Return
End Subroutine compute_dpdx_mean
!
!******************************************************************************
Subroutine res_mass(resms)
!*****************************************************************************
!	This procedure calculates the residual mass 
Implicit None
Integer :: i,j,k,j1d,k1d
Real(mytype) :: resms,myres
!
  myres = 0.0e0; resms = 0.0e0
!
!**************** must be in X pencil *******************
!
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
  j1d = xstart(2)-1+j
  k1d = xstart(3)-1+k
  myres = myres + abs(((u1(i,j,k)-u1(i+1,j,k))*dy(j1d)*dz(k1d)+     &
		       (v1(i,j,k)-v1(i,j+1,k))*dx(i)*dz(k1d)+     &
		       (w1(i,j,k)-w1(i,j,k+1))*dx(i)*dy(j1d))*roh)
End Do; End Do; End Do
!
  Call MPI_Allreduce(myres,resms,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
!
Return
End Subroutine res_mass
!
End Module mass_calc
