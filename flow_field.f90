Module flow_field
Use shared_data
Use mass_calc
Use boundary_conditions
Use profiles
Use decomp_2d_io
Use decomp_2d
!
Implicit None
Integer :: nxo,nyo,nzo
Real(mytype) :: leng,hei,wid,reo,conyo
!
Contains
!
!******************************************************************************
Subroutine read_restart
!******************************************************************************
! 	This procedure reads the velocit flow field to restart the simulation from 
!	previously calculated flow field
Implicit None
!
  Call read_restart_bin
  If (nrank == 0) Call read_restart_print
  If (nrank == 0) Call check_parameter
!
If (ifdt_fixed == 1) dt = dt_fixed
If (ifdt_variable == 1) dt  = dt_initial
If (if_dpdxin == 1) dpdx_mean = -dpdx_mean_in
!
Return
End Subroutine read_restart
!******************************************************************************
Subroutine read_restart_bin
!******************************************************************************
!-------This procedure reads velocity and pressure fields for future runs.
Implicit None
Integer :: samples_les,fh_rres,i,j,k
Integer(KIND=MPI_OFFSET_KIND) :: disp_res
Real(mytype), dimension(:,:,:), allocatable :: dum
disp_res = 0_MPI_OFFSET_KIND
!
 If (nrank == 0) write(*,*) trim(folder)//trim(unit22)
 Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//trim(unit22),MPI_MODE_RDONLY, &
			MPI_INFO_NULL,fh_rres,ierror)
!
 If (restrt_old==1) then
   Call decomp_2d_read_head3dold(fh_rres,disp_res,nxo,nyo,nzo,leng,hei,wid,reo)
 Else
   Call decomp_2d_read_head3d(fh_rres,disp_res,nxo,nyo,nzo,leng,hei,wid,reo,conyo)
 End If
 allocate(dum(xwrite(1),xwrite(2),xwrite(3)))
 Call decomp_2d_read_var(fh_rres,disp_res,nx,ny,nz,1,dum)
 u1(lo(1):up(1),lo(2):up(2),lo(3):up(3))=dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))
 Call decomp_2d_read_var(fh_rres,disp_res,nx,ny,nz,1,dum)
 v1(lo(1):up(1),lo(2):up(2),lo(3):up(3))=dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))
 Call decomp_2d_read_var(fh_rres,disp_res,nx,ny,nz,1,dum)
 w1(lo(1):up(1),lo(2):up(2),lo(3):up(3))=dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))
 Call decomp_2d_read_var(fh_rres,disp_res,nx,ny,nz,1,dum)
 p1(lo(1):up(1),lo(2):up(2),lo(3):up(3))=dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))
 deallocate(dum)
!
!       Reading varray (Only required in case of multigrid)
!
!If (ic_varray_rd == 1) Read(22) varray
!
  Call MPI_FILE_CLOSE(fh_rres,ierror) 
!
Return
End Subroutine read_restart_bin

!******************************************************************************
Subroutine read_restart_print
!******************************************************************************
! This procedure reads velocity and pressure fields for future runs.
Implicit None
Integer :: i,j,k,m,n

Write(21,*)leng,hei,wid,reo,dt,time,nxo,nyo,nzo,dpdx_mean,conyo
Write(*,*)leng,hei,wid,reo,dt,time,nxo,nyo,nzo,dpdx_mean,conyo
Write(*,*)'READVF:dpdx_mean=',dpdx_mean

Return
End Subroutine read_restart_print

!******************************************************************************
Subroutine check_parameter
!*****************************************************************************
! This procedure checks the number of grid points,domain and renold number
! from restart file and compare these values with those given by input file.
Implicit None

If (nxo /= nx .Or. nyo /= ny .Or. nzo /= nz) Then
   Write(21,155)
   Write(21,*)'nxo=',nxo,' nx=',nx
   Write(21,*)'nyo=',nyo,' ny=',ny
   Write(21,*)'nzo=',nzo,' nz=',nz
   Print*,'WRONG GRID SIZE'
   Stop
EndIf
!
If (leng /= length) Then
   Write(21,155)
   Write(21,*)' Old length =',leng,' New length =',length
   Print*,'WRONG Length'
   Stop
EndIf
!
If (hei /= height) Then
   Write(21,155)
   Write(21,*)' Old height =',hei,' New height =',height
   Print*,'WRONG height'
   Stop
EndIf
!
If (wid /= width) Then
   Write(21,155)
   Write(21,*)' Old width =',wid,' New width =',width
   Print*,'WRONG width'
   Stop
EndIf
!
If (reo /= re) Then
   Write(21,155)
   Write(21,*)'Old Re =',reo,' New Re =',re
   Print*,'Reynolds number changed'
End If
!
If (conyo /= cony .and. restrt_old /= 1) Then
   Write(21,155)
   Write(21,*)' Old cony =',conyo,'New cony =',cony
   Print*,'Cony changed'
End If
!
   155 Format(' Input parameters different. Abort...')
!
Return 
End Subroutine check_parameter

!******************************************************************************
Subroutine initial_flow
!******************************************************************************
!	This procedure determines the initial velocity and pressure field in case of new simulation.
Implicit None
Integer :: i,j,k
!
time = 0.0e0
!
If (ibmpfg == 0) Call parabolic          !Initial parabolic profiles
!
Call setup_pressure                           
!
Return
End Subroutine initial_flow
!
!******************************************************************************
Subroutine setup_pressure
!******************************************************************************
!	This procedure determines initial pressure field
Implicit None
Integer :: i,j,k
Real :: utemp,wtemp,ptemp
Real, Dimension(:), Allocatable :: temp_arr
!
Allocate (temp_arr(nx))
temp_arr = 0.0e0; dpz = 0.0e0; dpdx_mean = 0.0e0; ptemp = 0.0e0; p1 = 0.0e0
!
If (ichann == 1) Then
   utemp=1.5e0
   dpdx_mean = -2.0e0*utemp*nu
!       Turbulent flow
 If (irandm.eq.1)  dpdx_mean = -roh*retau*retau*nu*nu
 If (iconbc == 1) dpdx_mean = 0.0e0
   wtemp = 0.0e0
   dpdz_mean = -2.0e0*wtemp*nu
Else If (ichann == 0) Then
   dpdx_mean = 0.0e0
   dpdz_mean = 0.0e0
End If
!
If (Abs(dpdx_mean) > very_small) Then
 Do i = 1,nx
  temp_arr(i) = dpdx_mean*xc(i)
 End Do
 Else If (Abs(dpdz_mean) > very_small) Then
 Do k = 1,nz
  temp_arr(k) = dpdz_mean*zc(k)
 End Do
End If
!
If (Abs(dpdx_mean) > very_small) then
 If (iconbc == 1) then
  Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
  p1(i,j,k) = temp_arr(i)
 End Do; End Do; End Do
  Else If (iconbc == 0) then
  p1 = 0.0e0
 End If
End If
!
dpmean = dpdx_mean*length
dpz = dpdz_mean*width
!
Deallocate (temp_arr)

Return
End Subroutine setup_pressure
!
!******************************************************************************
Subroutine write_restart 
!******************************************************************************
! This procedure saves velocity and pressure fields for future runs.
Implicit None
!
  Call write_restart_bin

Return
End Subroutine write_restart
!
!******************************************************************************
Subroutine write_restart_bin
!******************************************************************************
!-------This procedure saves velocity and pressure fields for future runs.
Implicit None
integer :: fh_wres
integer(KIND=MPI_OFFSET_KIND) :: disp_res
real(mytype), dimension(:,:,:), allocatable :: dum
disp_res = 0_MPI_OFFSET_KIND
!
 Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//trim(unit23),MPI_MODE_CREATE+MPI_MODE_WRONLY, &
			MPI_INFO_NULL,fh_wres,ierror)
!
!If (irite_les == 1) Then
!  If (nrank == 0) Write(23) sgc
!End If
 Call decomp_2d_write_head3d(fh_wres,disp_res)

 allocate(dum(xwrite(1),xwrite(2),xwrite(3)))
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=u1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_wres,disp_res,nx,ny,nz,1,dum)
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=v1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_wres,disp_res,nx,ny,nz,1,dum)
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=w1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_wres,disp_res,nx,ny,nz,1,dum)
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=p1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_wres,disp_res,nx,ny,nz,1,dum)
 deallocate(dum)
!
!If (ic_varray_rite == 1) Write(23) varray
!
  Call MPI_FILE_CLOSE(fh_wres,ierror) 
!
Return
End Subroutine write_restart_bin
!
!******************************************************************************
Subroutine write_lamb2_par      !This subroutine writes output for lambda2
!******************************************************************************
Implicit None
Integer :: i,j,k,i1d,j1d,k1d
Integer, Dimension(3) :: xs,xe,xd
Character*10 :: ranknm
!
xs(:)=lo(:);xe(:)=up(:)
Do i=1,3
  If(xmin(i).Eq.1) Then
    xs(i)=xs(i)+2
  Else
    xs(i)=xs(i)-1
  End If
  If(xmax(i).Eq.1) Then
    xe(i)=xe(i)-2
  Else
    xe(i)=xe(i)+1
  End If
End Do
xd(:)=xe(:)-xs(:)+1
Write(ranknm,'(I10)') nrank+1000*ii
Open(Unit=403,File=trim(folder)//'lambda2-ins'//'_'//trim(adjustl(ranknm))//'.vtk')
Write(403,'(A)') '# vtk DataFile Version 2.0'
Write(403,'(A)') 'Volume example'
Write(403,'(A)') 'ASCII'
Write(403,'(A)') 'DATASET RECTILINEAR_GRID'
Write(403,'(A,3I10)') 'DIMENSIONS',xd(1),xd(2),xd(3)
Write(403,'(A,I10,A)') 'X_COORDINATES', xd(1), 'float'
Do i=xs(1),xe(1)
  i1d=xstart(1)+i-1
  Write(403,'(E15.7)') xc(i1d)
End Do
Write(403,'(A,I10,A)') 'Y_COORDINATES', xd(2), 'float'
Do j=xs(2),xe(2)
  j1d=xstart(2)+j-1
  Write(403,'(E15.7)') yc(j1d)
End Do
Write(403,'(A,I10,A)') 'Z_COORDINATES', xd(3), 'float'
Do k=xs(3),xe(3)
  k1d=xstart(3)+k-1
  Write(403,'(E15.7)') zc(k1d)
End Do
Write(403,'(A,I10)') 'POINT_DATA',xd(1)*xd(2)*xd(3)
Write(403,'(A)') 'SCALARS lam2 float 1'
Write(403,'(A)') 'LOOKUP_TABLE default'
DO k=xs(3),xe(3); DO j=xs(2),xe(2); DO i=xs(1),xe(1)
  Write(403,'(F20.10)') l2_array(i,j,k)
End Do; End Do; End Do
Close(403)
Return
End Subroutine write_lamb2_par
!
!******************************************************************************
Subroutine write_3dfield_par      !This subroutine writes 3d flow fields
!******************************************************************************
Implicit None
Integer :: i,j,k,i1d,j1d,k1d
Integer, Dimension(3) :: xs,xe,xd
Character*10 :: ranknm
!
xs(:)=lo(:);xe(:)=up(:)
Do i=1,3
  If(xmin(i).Eq.1) Then
    xs(i)=xs(i)+2
  Else
    xs(i)=xs(i)-1
  End If
  If(xmax(i).Eq.1) Then
    xe(i)=xe(i)-2
  Else
    xe(i)=xe(i)+1
  End If
End Do
xd(:)=xe(:)-xs(:)+1
Write(ranknm,'(I10)') nrank+1000*ii
Open(Unit=404,File=trim(folder)//'uvwp'//'_'//trim(adjustl(ranknm))//'.vtk')
WRITE(404,'(A)') '# vtk DataFile Version 2.0'
WRITE(404,'(A)') 'Volume example'
WRITE(404,'(A)') 'ASCII'
WRITE(404,'(A)') 'DATASET RECTILINEAR_GRID'
WRITE(404,'(A,3I10)') 'DIMENSIONS',xd(1),xd(2),xd(3)
Write(404,'(A,I10,A)') 'X_COORDINATES', xd(1), 'float'
Do i=xs(1),xe(1)
  i1d=xstart(1)+i-1
  Write(404,'(E15.7)') xc(i1d)
End Do
Write(404,'(A,I10,A)') 'Y_COORDINATES', xd(2), 'float'
Do j=xs(2),xe(2)
  j1d=xstart(2)+j-1
  Write(404,'(E15.7)') yc(j1d)
End Do
Write(404,'(A,I10,A)') 'Z_COORDINATES', xd(3), 'float'
Do k=xs(3),xe(3)
  k1d=xstart(3)+k-1
  Write(404,'(E15.7)') zc(k1d)
End Do
Write(404,'(A,I10)') 'POINT_DATA',xd(1)*xd(2)*xd(3)
!Write(404,'(A)') 'SCALARS p float 1'
!Write(404,'(A)') 'LOOKUP_TABLE default'
!DO k=xs(3),xe(3); DO j=xs(2),xe(2); DO i=xs(1),xe(1)
!  write(404,'(f20.10)') p1(i,j,k)
!End Do; End Do; End Do
If(ibcjet==1.And.jet_type==2) Then
Write(404,'(A)') 'SCALARS Fz float 1'
Write(404,'(A)') 'LOOKUP_TABLE default'
DO k=xs(3),xe(3); DO j=xs(2),xe(2); DO i=xs(1),xe(1)
  write(404,'(f20.10)') fbody(i,j,k)
End Do; End Do; End Do
End If
!Write(404,'(A)') 'VECTORS U float'
!DO k=xs(3),xe(3); DO j=xs(2),xe(2); DO i=xs(1),xe(1)
!  write(404,'(3f20.10)') u1(i,j,k),v1(i,j,k),w1(i,j,k)
!End Do; End Do; End Do
Close(404)
Return
End Subroutine write_3dfield_par
!
!******************************************************************************
Subroutine write_flucs_par      !This subroutine writes 3d streamwise u'
!******************************************************************************
Implicit None
Integer :: i,j,k,i1d,j1d,k1d
Integer, Dimension(3) :: xs,xe,xd
Character*10 :: ranknm
!
xs(:)=lo(:);xe(:)=up(:)
Do i=1,3
  If(xmin(i).Eq.1) Then
    xs(i)=xs(i)+2
  Else
    xs(i)=xs(i)-1
  End If
  If(xmax(i).Eq.1) Then
    xe(i)=xe(i)-2
  Else
    xe(i)=xe(i)+1
  End If
End Do
xd(:)=xe(:)-xs(:)+1
Write(ranknm,'(I10)') nrank+1000*ii
Open(Unit=405,File=trim(folder)//'uflucs'//'_'//trim(adjustl(ranknm))//'.vtk')
WRITE(405,'(A)') '# vtk DataFile Version 2.0'
WRITE(405,'(A)') 'Volume example'
WRITE(405,'(A)') 'ASCII'
WRITE(405,'(A)') 'DATASET RECTILINEAR_GRID'
WRITE(405,'(A,3I10)') 'DIMENSIONS',xd(1),xd(2),xd(3)
Write(405,'(A,I10,A)') 'X_COORDINATES', xd(1), 'float'
Do i=xs(1),xe(1)
  i1d=xstart(1)+i-1
  Write(405,'(E15.7)') xc(i1d)
End Do
Write(405,'(A,I10,A)') 'Y_COORDINATES', xd(2), 'float'
Do j=xs(2),xe(2)
  j1d=xstart(2)+j-1
  Write(405,'(E15.7)') yc(j1d)
End Do
Write(405,'(A,I10,A)') 'Z_COORDINATES', xd(3), 'float'
Do k=xs(3),xe(3)
  k1d=xstart(3)+k-1
  Write(405,'(E15.7)') zc(k1d)
End Do
Write(405,'(A,I10)') 'POINT_DATA',xd(1)*xd(2)*xd(3)
Write(405,'(A)') 'SCALARS u float 1'
Write(405,'(A)') 'LOOKUP_TABLE default'
DO k=xs(3),xe(3); DO j=xs(2),xe(2); DO i=xs(1),xe(1)
  write(405,'(f20.10)') u1(i,j,k)-umean(j)
End Do; End Do; End Do
Close(405)
Return
End Subroutine write_flucs_par
!
End Module flow_field
