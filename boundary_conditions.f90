Module boundary_conditions
Use shared_data
Use mass_calc
!Use cal_coefficient
!Use profiles
Use decomp_2d
Use mpi
Implicit None
Contains
!
!******************************************************************************
Subroutine initia_bc_update
!******************************************************************************
!	This procedure sets the boundary conditions for flow field in "initia"
!	subroutine.
Integer :: i,j,k,iuvw,is,it,i1d,j1d,k1d
Real(mytype), dimension(:,:,:), allocatable :: buf1,buf2,buf3
!
!**************** must be in X pencil *******************
!
  call MPI_Cart_shift(MPI_COMM_CART,0,1,left_row,right_row,ierror)
  call MPI_Cart_shift(MPI_COMM_CART,1,1,left_col,right_col,ierror)
!
  call MPI_COMM_RANK(MPI_COMM_ROW,row_rank,ierror)
  call MPI_COMM_RANK(MPI_COMM_COL,col_rank,ierror)
!
  xmin=0; ymin=0; zmin=0; xmax=0; ymax=0; zmax=0
  xmin(1)=1; xmax(1)=1; ymin(2)=1; ymax(2)=1; zmin(3)=1; zmax(3)=1
!
  If (col_rank == 0) then
    xmin(2)=1; ymin(1)=1; zmin(1)=1
  End If
  If (col_rank == p_row-1) then
    xmax(2)=1; ymax(1)=1; zmax(1)=1
  End If
!
  If (row_rank == 0) then
    xmin(3)=1; ymin(3)=1; zmin(2)=1
  End If
  If (row_rank == p_col-1) then
    xmax(3)=1; ymax(3)=1; zmax(2)=1
  End If
!
 lo=1
 up=xsize
 Do i=1,3
   If (xmin(i)==1) lo(i)=lo(i)-2
   If (xmax(i)==1) up(i)=up(i)+2
   xwrite(i)=up(i)-lo(i)+1
 End Do
!
Return
End Subroutine initia_bc_update
!
!******************************************************************************
Subroutine  inlet_initial
!******************************************************************************
Implicit None
!
!       Initialises xinlet & yinlet
!
  call xinlet_initial
  call yinlet_initial
!
Return
End Subroutine inlet_initial
!
!******************************************************************************
Subroutine  xinlet_initial
!******************************************************************************
Implicit None
Integer :: i,j,k,iuvw,is,it
!
 If (iconbc == 1) Then
   Do k = -1,xsize(3)+2; Do j = -1,xsize(2)+2
     xinlet(j,k,1,1,1) = u1(1,j,k)
     xinlet(j,k,2,1,1) = v1(0,j,k)
     xinlet(j,k,3,1,1) = w1(0,j,k)
     xinlet(j,k,1,2,1) = u1(nxt,j,k)
     xinlet(j,k,2,2,1) = v1(nxt,j,k)
     xinlet(j,k,3,2,1) = w1(nxt,j,k)
   End Do; End Do
   xinlet(:,:,:,:,2) = xinlet(:,:,:,:,1)
 End If
!
Return
End Subroutine xinlet_initial
!
!******************************************************************************
Subroutine  yinlet_initial
!******************************************************************************
Implicit None
Integer :: i,j,k,iuvw,is,it
!
  Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
    If (xmin(2) == 1) then
      yinlet1(i,k,1,1,1) = u1(i,0,k)
      yinlet1(i,k,2,1,1) = v1(i,1,k)
      yinlet1(i,k,3,1,1) = w1(i,0,k)
    End If
    If (xmax(2) == 1) then
      yinlet1(i,k,1,2,1) = u1(i,xsize(2)+1,k)
      yinlet1(i,k,2,2,1) = v1(i,xsize(2)+1,k)
      yinlet1(i,k,3,2,1) = w1(i,xsize(2)+1,k)
    End If
  End Do; End Do
!
  Do is=1,2; Do iuvw=1,3; Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
    yinlet1(i,k,iuvw,is,2) = yinlet1(i,k,iuvw,is,1)
  End Do; End Do; End Do; End Do
!
Return
End Subroutine yinlet_initial
!
!******************************************************************************
Subroutine  xbc_update
!******************************************************************************
Implicit None
Integer :: i,j,k,iuvw,is,it
!
 If (iconbc == 1) Then
!
   Do i = -1,0; Do k = -1,xsize(3)+2; Do j = -1,xsize(2)+2
     u1(i,j,k) = xinlet(j,k,1,1,2)
     u1(1,j,k) = xinlet(j,k,1,1,2)
     v1(i,j,k) = xinlet(j,k,2,1,2)
     w1(i,j,k) = xinlet(j,k,3,1,2)
   End Do; End Do; End Do
   Do i = 1,2; Do k = -1,xsize(3)+2; Do j = -1,xsize(2)+2
     u1(i+xsize(1),j,k) = xinlet(j,k,1,2,2)
     v1(i+xsize(1),j,k) = xinlet(j,k,2,2,2)
     w1(i+xsize(1),j,k) = xinlet(j,k,3,2,2)
   End Do; End Do; End Do
!
   Do is=1,2; Do iuvw=1,3; Do k = -1,xsize(3)+2; Do j = -1,xsize(2)+2
     xinlet(j,k,iuvw,is,1) = xinlet(j,k,iuvw,is,2)
   End Do; End Do; End Do; End Do
!
 End If
!
Return
End Subroutine xbc_update
!
!******************************************************************************
Subroutine  pbc_update
!******************************************************************************
Implicit None
Integer :: i,j,k

If (jconbc==1) then
  If(xmin(2)==1) then
    Do k=-1,xsize(3)+2; Do i=-1,xsize(1)+2
      p1(i,0,k)=p1(i,1,k)
      p1(i,-1,k)=p1(i,1,k)
    End Do; End Do
  End If
  If(xmax(2)==1) then
    Do k=-1,xsize(3)+2; Do i=-1,xsize(1)+2
      p1(i,xsize(2)+1,k)=p1(i,xsize(2),k)
      p1(i,xsize(2)+2,k)=p1(i,xsize(2),k)
    End Do; End Do
  End If
End If

!
Return
End Subroutine pbc_update
!
!******************************************************************************
Subroutine  ybc_update
!******************************************************************************
Implicit None
Integer :: i,j,k,iuvw,is,it
!
  If (xmin(2) == 1) then
    Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
      u1(i,0,k) = yinlet1(i,k,1,1,2)
      v1(i,0,k) = yinlet1(i,k,2,1,2)
      v1(i,1,k) = yinlet1(i,k,2,1,2)
      w1(i,0,k) = yinlet1(i,k,3,1,2)
      u1(i,-1,k) = yinlet1(i,k,1,1,2)
      v1(i,-1,k) = yinlet1(i,k,2,1,2)
      w1(i,-1,k) = yinlet1(i,k,3,1,2)
    End Do; End Do
  End If
!
  If (xmax(2) == 1) then
    Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
      u1(i,xsize(2)+1,k) = yinlet1(i,k,1,2,2)
      v1(i,xsize(2)+1,k) = yinlet1(i,k,2,2,2)
      w1(i,xsize(2)+1,k) = yinlet1(i,k,3,2,2)
      u1(i,xsize(2)+2,k) = yinlet1(i,k,1,2,2)
      v1(i,xsize(2)+2,k) = yinlet1(i,k,2,2,2)
      w1(i,xsize(2)+2,k) = yinlet1(i,k,3,2,2)
    End Do; End Do
  End If
!
  Do is=1,2; Do iuvw=1,3; Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
    yinlet1(i,k,iuvw,is,1) = yinlet1(i,k,iuvw,is,2)
  End Do; End Do; End Do; End Do
!
Return
End Subroutine ybc_update
!
!******************************************************************************
Subroutine  bc_left_update
!******************************************************************************
!	This procedure is to update boundary conditions in COMPUTE_LEFT at  icmpnt == 3
Implicit None
Integer :: ibeg,iuvw,i,j,k
!
!	Stream wise Boundary conditons
!
!**************** must be in X pencil *******************
!
  call update_halo(uhat1,2)
  call update_halo(vhat1,2)
  call update_halo(what1,2)
!
! Use next step boundary conditions for u,v,w
! Boundary conditions in y
!
  If (xmin(2) == 1) then
    Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
      vhat1(i,1,k) = yinlet1(i,k,2,1,2)
    End Do; End Do
  End If
!
  If (xmax(2) == 1) then
    Do k = -1,xsize(3)+2; Do i = -1,xsize(1)+2
      vhat1(i,xsize(2)+1,k) = yinlet1(i,k,2,2,2)
    End Do; End Do
  End If
!
! Boundary conditions in x
!
  If (iconbc == 1) Then
    Do k = -1,xsize(3)+1; Do j = -1,xsize(2)+2
      uhat1(1,j,k) = xinlet(j,k,1,1,2)
      uhat1(nxt,j,k) = xinlet(j,k,1,2,2)
    End Do; End Do
  End If
!
Return
End Subroutine bc_left_update
!
!******************************************************************************
Subroutine bc_update
!******************************************************************************
!	This procedure is to update boundary conditions and transposes velocities
!       in update_field.f90 module
Implicit None
Real(mytype), dimension(:,:,:), allocatable :: buf1,buf2,buf3
!
!********* Transpose u,v,w  ********************
!
Allocate(buf1(1:xsize(1),1:xsize(2),1:xsize(3)), &
	 buf2(1:ysize(1),1:ysize(2),1:ysize(3)), &
	 buf3(1:zsize(1),1:zsize(2),1:zsize(3)) )
!
  buf1(1:xsize(1),1:xsize(2),1:xsize(3))=v1(1:xsize(1),1:xsize(2),1:xsize(3))
  Call transpose_x_to_y(buf1,buf2)
  v2(1:ysize(1),1:ysize(2),1:ysize(3))=buf2(1:ysize(1),1:ysize(2),1:ysize(3))
!
  buf1(1:xsize(1),1:xsize(2),1:xsize(3))=w1(1:xsize(1),1:xsize(2),1:xsize(3))
  Call transpose_x_to_y(buf1,buf2)
  w2(1:ysize(1),1:ysize(2),1:ysize(3))=buf2(1:ysize(1),1:ysize(2),1:ysize(3))
  Call transpose_y_to_z(buf2,buf3)
  w3(1:zsize(1),1:zsize(2),1:zsize(3))=buf3(1:zsize(1),1:zsize(2),1:zsize(3))
!
Deallocate(buf1,buf2,buf3)
!
!	Boundary Conditions
!
  call update_halo(u1,2)
  call update_halo(v1,2)
  call update_halo(w1,2)
  call update_halo(p1,2)
!
  call xbc_update
  call ybc_update
  call pbc_update
!
  call update_halo(v2,2)
  call update_halo(w2,2)
!
  call update_halo(w3,2)
!
Return 
End Subroutine bc_update
!
!*****************************************************************************
Subroutine bc_outflow_right
!*****************************************************************************
!       This procedure implements the convective boundary conditions in compute_right_implicit module 
!	and determines the outflow boundary velocities
Implicit None
!
If (icmpnt == 1) Call bc_outflow_u
If (icmpnt == 2) Call bc_outflow_v
If (icmpnt == 3) Call bc_outflow_w
!
Return
End Subroutine bc_outflow_right
!
!******************************************************************************
Subroutine bc_outflow_u
!*****************************************************************************
!	This procedure calculates the outflow bounday conditions xinlet(i,k,1,2,2) based on du/dt+Udu/dx = 0
!	in the u-momentum equation.
Implicit None
Integer :: j,k
Real(mytype) :: corate
!
  Call rmsfl
!
uconv=dt*rmsflo/height/width/roh/dx(nx)

Do k = 1,xsize(3); Do j = -1,xsize(2)+2
xinlet(j,k,1,2,2) = u1(nxt,j,k) - uconv*(u1(nxt,j,k)-u1(nx,j,k))
End Do; End Do
!
  Call rmsfl
!
corate=(rmsfli+rmsflb)/rmsflc
Write(21,*)'corate=',corate
!
Do k = 1,xsize(3); Do j = -1,xsize(2)+2
xinlet(j,k,1,2,2) = xinlet(j,k,1,2,2)*corate
End Do; End Do
!
Do k = 1,xsize(3); Do j = -1,0
xinlet(j,k,1,2,2) = yinlet1(nxt,k,1,1,2)
End Do; End Do
!
Do k = 1,xsize(3); Do j = xsize(2)+1,xsize(2)+2
xinlet(j,k,1,2,2) = xinlet(ny,k,1,2,2)
End Do; End Do
!
Return
End Subroutine bc_outflow_u
!
!******************************************************************************
Subroutine bc_outflow_v
!******************************************************************************
!	This procedure calculate the outflow bounday conditions xinlet(i,k,2,2,2) based on dv/dt+Udv/dx = 0
!	in the v-momentum equation.
Implicit None
Integer :: j,k
!
Do k = 1,xsize(3); Do j = -1,xsize(2)+2
xinlet(j,k,2,2,2) = v1(nxt,j,k)-uconv*(v1(nxt,j,k)-v1(nx,j,k))
End Do; End Do
!
Do k = 1,xsize(3); Do j = -1,0
xinlet(j,k,2,2,2) = yinlet1(nxt,k,2,1,2)
End Do; End Do
!
Do k = 1,xsize(3); Do j = xsize(2)+1,xsize(2)+2
xinlet(j,k,2,2,2) = yinlet1(nxt,k,2,2,2)
End Do; End Do
!
Return
End Subroutine bc_outflow_v
!
!******************************************************************************
Subroutine bc_outflow_w
!******************************************************************************
!	This procedure calculate the outflow bounday conditions xinlet(i,k,3,2,2) based on dw/dt+Udw/dx = 0
!	in the w-momentum equation.
Implicit None
Integer :: j,k
!
Do k = 1,xsize(3); Do j = -1,xsize(2)+2
xinlet(j,k,3,2,2) = w1(nxt,j,k)-uconv*(w1(nxt,j,k)-w1(nx,j,k))
End Do; End Do
!
Do k = 1,xsize(3); Do j = -1,0
xinlet(j,k,3,2,2) = yinlet1(nxt,k,3,1,2)
End Do; End Do
!
Do k = 1,xsize(3); Do j = xsize(2)+1,xsize(2)+2
xinlet(j,k,3,2,2) = yinlet1(nxt,k,3,2,2)
End Do; End Do
!
Return
End Subroutine bc_outflow_w
!	
End Module boundary_conditions
