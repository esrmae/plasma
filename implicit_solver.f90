Module implicit_solver
Use shared_data
Use compute_right_implicit
Use compute_left_implicit
!
Contains
!
!******************************************************************************
Subroutine implicit_setup
!******************************************************************************
Implicit none
!
Allocate(iux_vis(0:nx+1,2),iuy_vis(0:ny+1,2,2),iuy_con(0:ny+1,2),ivx_vis(0:nx+1,2,2),&
ivy_vis(0:ny+1,2),ivx_con(0:nx+1,2),iux_grdvis(0:nxt,2),ivy_grdvis(nyt,2),             &
iuy_grdvis(ny+1,2,2),iuy_grdcon(ny+1,2),ivx_grdvis(nx+1,2,2),ivx_grdcon(nx+1,2),     &
iux_grdcon(0:nxt,2),ivy_grdcon(nyt,2),iux_con(0:nxt,2),ivy_con(nyt,2))
iux_vis=0.0e0;iuy_vis=0.0e0;iuy_con=0.0e0;iux_grdvis=0.0e0;ivy_grdvis=0.0e0;iuy_grdvis=0.0e0;ivx_grdvis=0.0e0
ivx_vis=0.0e0;ivy_vis=0.0e0;ivx_con=0.0e0;iuy_grdcon=0.0e0;ivx_grdcon=0.0e0;iux_grdcon=0.0e0;ivy_grdcon=0.0e0
iux_con=0.0e0;ivy_con=0.0e0
!
! Both du1 and du2 must be defined in x pencil
Allocate (du1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),du2(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2))
du1=0.0e0;du2=0.0e0
!
If(jconbc.eq.1) then
ivy_vis(2,1)=1.0e0;ivy_vis(ny,2)=1.0e0;
ivy_grdvis(1,1)=1.0e0;ivy_grdvis(ny,2)=1.0e0
iuy_vis(1,1,1)=1.0e0;iuy_vis(2,2,1)=1.0e0; iuy_vis(ny-1,2,2)=1.0e0;iuy_vis(ny,1,2)=1.0e0
iuy_grdvis(1,1,1)=1.0e0;iuy_grdvis(2,2,1)=1.0e0; iuy_grdvis(ny,2,2)=1.0e0;iuy_grdvis(nyt,1,2)=1.0e0
iuy_con(1,1)=1.0e0; iuy_con(ny,2)=1.0e0
iuy_grdcon(1,1)=1.0e0; iuy_grdcon(nyt,2)=1.0e0
ivy_con(2,1)=1.0e0;ivy_con(ny,2)=1.0e0
ivy_grdcon(1,1)=1.0e0;ivy_grdcon(ny,2)=1.0e0
End If
!
Return
End Subroutine implicit_setup
!
!******************************************************************************
Subroutine implicit_solve
!******************************************************************************
Implicit None
!
Do icmpnt=1,ncmpnt,if2d_noy+1
  Call compute_rhs_implicit
  Call compute_lhs_implicit
 End Do
!
Return
End Subroutine implicit_solve
!
End Module implicit_solver
