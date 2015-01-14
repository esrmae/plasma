Module compute_left_implicit
Use shared_data
Use solvers
Use boundary_conditions
Implicit none
Integer, private :: i,j,k,iblock,ib,is
Real(mytype), private :: bt,nutxz1,nutxz2,nutxy1,nutxy2,nutc1,nutc2,  &
nutyz1, nutyz2,rcouxu1,rcouxu2,a1yvu,b1yvu, c1yvu,rcouyv1,rcouyv2,gt,    &
rcouyu1,rcouyu2,rcovxu1,rcovxu2,a2xuv,b2xuv,c2xuv,rcovyv1,rcovyv2,rcovxv1, &
rcovxv2,rcowxu1,rcowxu2,a3xuw,b3xuw,c3xuw,rcowyv1,rcowyv2,a3yvw,b3yvw,     &
c3yvw,rcowxw1,rcowxw2,rcowyw1,rcowyw2,rcouzw1,a1zwu,b1zwu,c1zwu,rcovzw1,   &
rcovzw2,a2zwv,b2zwv,c2zwv,rcowzw1,rcowzw2,rcovzv1,rcovzv2,rcouzu1,rcouzu2, &
rcouzw2
Integer, Dimension(3) :: beg
!
Real(mytype), dimension(:,:,:),allocatable :: a,b,c,d,e
contains
!
!******************************************************************************
Subroutine compute_lhs_implicit
!******************************************************************************
Implicit none
!
bt = dt/2.0e0; gt = dt/2.0e0
iblock = 0; ib = iblock+1
!
 Call solvez_implicit       	     !(1-a3)v* = RHS
!
 Call solvey_implicit                !(i-a2)v** = v*
!
 Call solvex_implicit                !(1-a1)v = v**
!
!       Boundary Conditions
!
If (icmpnt == 3) Call bc_left_update
!
Return
End Subroutine compute_lhs_implicit
!
!******************************************************************************
Subroutine beg_setup(pencil)
!******************************************************************************
Implicit none
Integer :: pencil
!
beg = 1 !=(ibeg,jbeg,kbeg)
!
If (pencil == 1) then
  If (icmpnt == 1) Then
    If ((iconbc == 1).and.(xmin(1) == 1)) beg(1) = 2
  Else If (icmpnt == 2) Then
    If ((jconbc == 1).and.(xmin(2) == 1)) beg(2) = 2
  Else If (icmpnt == 3) Then
    If ((kconbc == 1).and.(xmin(3) == 1)) beg(3) = 2
  End If 
Else If (pencil == 2) then
  If (icmpnt == 1) Then
    If ((iconbc == 1).and.(ymin(1) == 1)) beg(1) = 2
  Else If (icmpnt == 2) Then
    If ((jconbc == 1).and.(ymin(2) == 1)) beg(2) = 2
  Else If (icmpnt == 3) Then
    If ((kconbc == 1).and.(ymin(3) == 1)) beg(3) = 2
  End If 
Else If (pencil == 3) then
  If (icmpnt == 1) Then
    If ((iconbc == 1).and.(zmin(1) == 1)) beg(1) = 2
  Else If (icmpnt == 2) Then
    If ((jconbc == 1).and.(zmin(2) == 1)) beg(2) = 2
  Else If (icmpnt == 3) Then
    If ((kconbc == 1).and.(zmin(3) == 1)) beg(3) = 2
  End If 
End If 
!
Return
End Subroutine beg_setup
!******************************************************************************
Subroutine solvez_implicit
!******************************************************************************
!       This procedure solves the left hand side in Z-direction
Implicit none
Integer :: i1d,j1d
!
!**************** must be in X pencil *******************
 Call beg_setup(1)
!
 If (icmpnt == 2) then; Call rhs2_calc;
 Else If (icmpnt == 3) then; Call rhs3_calc; End If
!
!********* Transpose intf  ********************
 Call transpose_x_to_y(intf1,intf2)
 Call transpose_y_to_z(intf2,intf3)
!
 Call beg_setup(3)
!
!**************** must be in Z pencil *******************
!
Allocate(a(zsize(1),zsize(2),zsize(3)),b(zsize(1),zsize(2),zsize(3)), &
         c(zsize(1),zsize(2),zsize(3)),d(zsize(1),zsize(2),zsize(3)), &
         e(zsize(1),zsize(2),zsize(3)))
a=0.0e0; b=0.0e0; c=0.0e0; d=0.0e0; e=0.0e0
!
!*******************************
	If(icmpnt == 1) Then
!*******************************
!
Do k = beg(3),zsize(3)
Do j = beg(2),zsize(2)
 j1d=zstart(2)-1+j
Do i = beg(1),zsize(1)
 i1d=zstart(1)-1+i
!
!       Cofficient of d2u/dz2
 nutxz1 = bre(i1d)*(nu+rnutxz3(i,j,k))/dz(k)
 nutxz2 = bre(i1d)*(nu+rnutxz3(i,j,k+1))/dz(k)
!       Cofficient of dwu/dz
 rcouzw1 = gt/dz(k)/(dx(i1d)+dx(i1d-1))*(w3(i,j,k)*dx(i1d-1)+w3(i-1,j,k)*dx(i1d))
 rcouzw2 = gt/dz(k)/(dx(i1d)+dx(i1d-1))*(w3(i,j,k+1)*dx(i1d-1)+ w3(i-1,j,k+1)*dx(i1d))
 a1zwu = -rcouzw1*dz(k)/(dz(k)+dz(k-1))
 b1zwu = rcouzw2*dz(k+1)/(dz(k)+dz(k+1)) -rcouzw1*dz(k-1)/(dz(k)+dz(k-1))
 c1zwu = rcouzw2*dz(k)/(dz(k)+dz(k+1))
 !a(i,j,k) = bt*nutxz1*gcdz_f(k,4,1,1)
 b(i,j,k) = -bt *(nutxz2*gcdz_f(k+1,4,1,1)-nutxz1*gcdz_f(k,3,1,1))+a1zwu
 c(i,j,k) = 1.0e0-bt *(nutxz2*gcdz_f(k+1,3,1,1)-nutxz1*gcdz_f(k,2,1,1))+b1zwu
 d(i,j,k) = -bt *(nutxz2*gcdz_f(k+1,2,1,1)-nutxz1*gcdz_f(k,1,1,1))+c1zwu
 !e(i,j,k) = -bt*nutxz2*gcdz_f(k+1,1,1,1)
End Do; End Do; End Do
!
 Call trip3(b,c,d,intf3,beg,3)
!
!*******************************
        Else If(icmpnt == 2) Then
!*******************************
!
Do k = beg(3),zsize(3)
Do j = beg(2),zsize(2)
 j1d=zstart(2)-1+j
Do i = beg(1),zsize(1)
 i1d=zstart(1)-1+i
!       Cofficients for d2v/dz2
nutyz1 = brec(i1d)*(nu+rnutyz3(i,j,k))/dz(k)
nutyz2 = brec(i1d)*(nu+rnutyz3(i,j,k+1))/dz(k)
!       Cofficients for d(wv)/dz
rcovzw1 = gt/dz(k)/(dy(j1d)+dy(j1d-1))*(w3(i,j,k)*dy(j1d-1)+w3(i,j-1,k)*dy(j1d))
rcovzw2 = gt/dz(k)/(dy(j1d)+dy(j1d-1))*(w3(i,j,k+1)*dy(j1d-1)+ w3(i,j-1,k+1)*dy(j1d))
a2zwv = -rcovzw1*dz(k)/(dz(k)+dz(k-1))
b2zwv = rcovzw2*dz(k+1)/(dz(k)+dz(k+1))-rcovzw1*dz(k-1)/(dz(k)+dz(k-1))
c2zwv = rcovzw2*dz(k)/(dz(k)+dz(k+1))
!a(i,j,k) = bt*nutyz1*gcdz_f(k,4,1,1)
b(i,j,k) = bt*(nutyz1*gcdz_f(k,3,1,1)-nutyz2*gcdz_f(k+1,4,1,1))+a2zwv
c(i,j,k) = 1.0e0+bt*(nutyz1*gcdz_f(k,2,1,1)-nutyz2*gcdz_f(k+1,3,1,1))+b2zwv
d(i,j,k) = bt*(nutyz1*gcdz_f(k,1,1,1)-nutyz2*gcdz_f(k+1,2,1,1))+c2zwv
!e(i,j,k) = -bt*nutyz2*gcdz_f(k+1,1,1,1)
End Do; End Do; End Do
!
 Call trip3(b,c,d,intf3,beg,3)
!
!*******************************
        Else If(icmpnt == 3) Then
!*******************************
!
Do k = beg(3),zsize(3)
Do j = beg(2),zsize(2)
 j1d=zstart(2)-1+j
Do i = beg(1),zsize(1)
 i1d=zstart(1)-1+i
!       Cofficients for d2w/dz2
nutc1 = brec(i1d)*(nu+2.0e0*pitero3(i,j,k-1))/dz(k-1)/dzs(k)
nutc2 = brec(i1d)*(nu+2.0e0*pitero3(i,j,k))/dz(k)/dzs(k)
!       Cofficints for d(ww)/dz
rcowzw1 = gt/2.0e0/dzs(k)*(w3(i,j,k)+w3(i,j,k-1))
rcowzw2 = gt/2.0e0/dzs(k)*(w3(i,j,k)+w3(i,j,k+1))
!
b(i,j,k) = -bt*nutc1-rcowzw1
c(i,j,k) = 1.0e0+bt*(nutc2+nutc1)+(rcowzw2-rcowzw1)
d(i,j,k) = -bt*nutc2+rcowzw2
End Do; End Do; End Do
!
 Call trip3(b,c,d,intf3,beg,3)
!
End If
!
Deallocate(a,b,c,d,e)
!
 Return
End Subroutine solvez_implicit
!
!******************************************************************************
Subroutine solvey_implicit   
!******************************************************************************
!	This procedure solves the left hand side in Y-direction 
Implicit none
Integer :: i1d,k1d
!
!********* Transpose intf *************
 Call transpose_z_to_y(intf3,intf2)

!**************** must be in Y pencil *******************
!
 Call beg_setup(2)
!
Allocate(a(ysize(1),ysize(2),ysize(3)),b(ysize(1),ysize(2),ysize(3)), &
         c(ysize(1),ysize(2),ysize(3)),d(ysize(1),ysize(2),ysize(3)), &
         e(ysize(1),ysize(2),ysize(3)))
a=0.0e0; b=0.0e0; c=0.0e0; d=0.0e0; e=0.0e0
!
!*****************************
	If(icmpnt == 1) Then
!*****************************
Do k = beg(3),ysize(3)
 k1d=ystart(3)-1+k
Do j = beg(2),ysize(2)
Do i = beg(1),ysize(1)
 i1d=ystart(1)-1+i
!
nutxy1 = bre(i1d)*(nu+rnutxy2(i,j,k))/dy(j)
nutxy2 = bre(i1d)*(nu+rnutxy2(i,j+1,k))/dy(j)
!	Cofficients for d(vu)/dy
rcouyv1 = gt/dy(j)*(v2(i,j,k)*dx(i1d-1)+v2(i-1,j,k)*dx(i1d))/ (dx(i1d)+dx(i1d-1))
rcouyv2 = gt/dy(j)*(v2(i,j+1,k)*dx(i1d-1)+v2(i-1,j+1,k)*dx(i1d))/ (dx(i1d)+dx(i1d-1))
a1yvu = -rcouyv1*dy(j)/(dy(j)+dy(j-1))
b1yvu =  rcouyv2*dy(j+1)/(dy(j)+dy(j+1))-rcouyv1*dy(j-1)/(dy(j)+dy(j-1))
c1yvu = rcouyv2*dy(j)/(dy(j)+dy(j+1))
If (j == 1.and.ibsoth == 1) Then
a1yvu = -rcouyv1
b1yvu =  rcouyv2*dy(j+1)/(dy(j)+dy(j+1))
End If
If (j == ny.and.ibnoth == 1) Then
b1yvu = -rcouyv1*dy(j-1)/(dy(j)+dy(j-1))
c1yvu =  rcouyv2
End If
!
a(i,j,k) = bt*nutxy1*gcdy_f(j,4,1,ib)
b(i,j,k) = -bt*(nutxy2*gcdy_f(j+1,4,1,ib)-nutxy1*gcdy_f(j,3,1,ib))+(1.0e0-iuy_con(j,1))*a1yvu
c(i,j,k) = 1.0e0-bt*(nutxy2*gcdy_f(j+1,3,1,ib)-nutxy1*gcdy_f(j,2,1,ib)) + b1yvu
d(i,j,k) = -bt*(nutxy2*gcdy_f(j+1,2,1,ib)-nutxy1* gcdy_f(j,1,1,ib))+(1.0e0-iuy_con(j,2))*c1yvu
e(i,j,k) = -bt*nutxy2*gcdy_f(j+1,1,1,ib)
End Do; End Do; End Do
!
Do k = beg(3),ysize(3);Do i = beg(1),ysize(1)
b(i,1,k) = 0.0e0; a(i,1,k) = 0.0e0; a(i,2,k) = 0.0e0
e(i,ny-1,k) = 0.0e0; d(i,ny,k) = 0.0e0; e(i,ny,k) = 0.0e0
End Do; End Do
!
 Call pentay(a,b,c,d,e,intf2,beg,2)
!
!*****************************
	Else If(icmpnt == 2) Then
!*****************************
!
Do k = beg(3),ysize(3)
Do j = beg(2),ysize(2)
Do i = beg(1),ysize(1)
 i1d=ystart(1)-1+i
!
nutc1 = brec(i1d)*(nu+2.0e0*pitero2(i,j-1,k))/dy(j-1)/dys(j)
nutc2 = brec(i1d)*(nu+2.0e0*pitero2(i,j,k))/dy(j)/dys(j)
!	Cofficients for d(vv)/dy
rcovyv1 = gt/2.0e0/dys(j)*(v2(i,j,k)+v2(i,j-1,k))
rcovyv2 = gt/2.0e0/dys(j)*(v2(i,j,k)+v2(i,j+1,k))
!
b(i,j,k) = -bt*nutc1-(1.0e0-ivy_vis(j,1))*rcovyv1
c(i,j,k) = 1.0e0+bt*(nutc1+nutc2)+(rcovyv2-rcovyv1)
d(i,j,k) = -bt*nutc2+(1.0e0-ivy_vis(j,2))*rcovyv2
End Do; End Do; End Do
!
Do k = beg(3),ysize(3);Do i = beg(1),ysize(1)
b(i,beg(2),k) = 0.0e0
d(i,ny,k) = 0.0e0
End Do; End Do
!
 Call dtridy(b,c,d,intf2,beg,2)
!
!*****************************
	Else If(icmpnt == 3) Then
!*****************************
!
Do k = beg(3),ysize(3)
 k1d=ystart(3)-1+k
Do j = beg(2),ysize(2)
Do i = beg(1),ysize(1)
 i1d=ystart(1)-1+i
!
!	Cofficients for d2w/dy2
nutyz1 = brec(i1d)*(nu+rnutyz2(i,j,k))/dy(j)
nutyz2 = brec(i1d)*(nu+rnutyz2(i,j+1,k))/dy(j)
!	Cofficients for d(vw)/du
rcowyv1 = gt/dy(j)/(dz(k1d)+dz(k1d-1))*(v2(i,j,k)*dz(k1d-1)+v2(i,j,k-1)*dz(k1d))
rcowyv2 = gt/dy(j)/(dz(k1d)+dz(k1d-1))*(v2(i,j+1,k)*dz(k1d-1)+v2(i,j+1,k-1)*dz(k1d))
a3yvw = -rcowyv1*dy(j)/(dy(j)+dy(j-1))
b3yvw =  rcowyv2*dy(j+1)/(dy(j)+dy(j+1))-rcowyv1*dy(j-1)/(dy(j)+dy(j-1))
c3yvw =  rcowyv2*dy(j)/(dy(j)+dy(j+1))
If (j == 1.and.ibsoth == 1) Then
a3yvw = -rcowyv1
b3yvw =  rcowyv2*dy(j+1)/(dy(j)+dy(j+1))
End If
If (j == ny.and.ibnoth == 1) Then
b3yvw =  -rcowyv1*dy(j-1)/(dy(j)+dy(j-1))
c3yvw = rcowyv2
End If
!
a(i,j,k) = bt*nutyz1*gcdy_f(j,4,1,ib)
b(i,j,k) = bt*(nutyz1*gcdy_f(j,3,1,ib)-nutyz2* gcdy_f(j+1,4,1,ib))+a3yvw*(1.0e0-iuy_con(j,1))
c(i,j,k) = 1.0e0+bt *(nutyz1*gcdy_f(j,2,1,ib)-nutyz2* gcdy_f(j+1,3,1,ib))+b3yvw
d(i,j,k) = bt*(nutyz1*gcdy_f(j,1,1,ib)-nutyz2* gcdy_f(j+1,2,1,ib))+c3yvw*(1.0e0-iuy_con(j,2))
e(i,j,k) = -bt*nutyz2*gcdy_f(j+1,1,1,ib)
End Do; End Do; End Do
!
Do k = beg(3),ysize(3);Do i = beg(1),ysize(1)
a(i,1,k) = 0.0e0; b(i,1,k) = 0.0e0; a(i,2,k) = 0.0e0
e(i,ny-1,k) = 0.0e0; d(i,ny,k) = 0.0e0; e(i,ny,k) = 0.0e0
End Do; End Do
!
 Call pentay(a,b,c,d,e,intf2,beg,2)
!
End If
!
Deallocate(a,b,c,d,e)
!
Return
End Subroutine solvey_implicit    
!
!******************************************************************************
Subroutine solvex_implicit   
!******************************************************************************
!	This procedure solves the left hand side in X-direction 
Implicit none
Integer :: j1d,k1d
!
 Call transpose_y_to_x(intf2,intf1)
!
 Call beg_setup(1)
!
Allocate(a(xsize(1),xsize(2),xsize(3)),b(xsize(1),xsize(2),xsize(3)), &
         c(xsize(1),xsize(2),xsize(3)),d(xsize(1),xsize(2),xsize(3)), &
         e(xsize(1),xsize(2),xsize(3)))
a=0.0e0; b=0.0e0; c=0.0e0; d=0.0e0; e=0.0e0
!
!**************** must be in X pencil *******************
!
!*****************************
If(icmpnt == 1) Then
!*****************************
Do k = beg(3),xsize(3)
 k1d=xstart(3)-1+k
Do j = beg(2),xsize(2)
 j1d=xstart(2)-1+j
Do i = beg(1),xsize(1)
!
nutc1 = buffc(i-1)*brec(i-1)*(nu+2.0e0*pitero1(i-1,j,k))/dx(i-1)/dxs(i)
nutc2 = buffc(i)*brec(i)*(nu+2.0e0*pitero1(i,j,k))/dx(i)/dxs(i)
!	Cofficient for duu/dx
rcouxu1 = gt/2.0e0/dxs(i)*(u1(i-1,j,k)+u1(i,j,k))
rcouxu2 = gt/2.0e0/dxs(i)*(u1(i+1,j,k)+u1(i,j,k))
b(i,j,k) = -bt*nutc1-(1.0e0-iux_vis(i,1))*rcouxu1
c(i,j,k) = 1.0e0+bt*(nutc1+nutc2)+(rcouxu2-rcouxu1)
d(i,j,k) = -bt*nutc2+(1.0e0-iux_vis(i,2))*rcouxu2
End Do; End Do; End Do
!
!	Periodic B.C.
!
 Call trip1_new(b,c,d,intf1,beg,1)
!
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = beg(1),xsize(1)
 du1(i,j,k) = intf1(i,j,k)
End Do; End Do; End Do
!
!	Periodic boundary conditons
  call update_halo(du1,2)
!
Do k = 0,xsize(3)+1; Do i = 0,xsize(1)+1
If (xmin(2) == 1) du1(i,0,k) = 0.0e0
If (xmax(2) == 1) du1(i,xsize(2)+1,k) = 0.0e0
End Do; End Do
!	Gets missing zeros
If (xmin(3)==1) then
  du1(0,:,0) = 0.0e0
  du1(nx+1,:,0) = 0.0e0
End If
If (xmax(3)==1) then
  du1(0,:,xsize(3)+1) = 0.0e0
  du1(xsize(1)+1,:,xsize(3)+1) = 0.0e0
End If
!
!*****************************
Else If(icmpnt == 2) Then
!*****************************
!
Do k = beg(3),xsize(3)
Do j = beg(2),xsize(2)
 j1d=xstart(2)-1+j
Do i = beg(1),xsize(1)
nutxy1 = buff(i)*bre(i)*(nu+rnutxy1(i,j,k))/dx(i)
nutxy2 = buff(i+1)*bre(i+1)*(nu+rnutxy1(i+1,j,k))/dx(i)
!	Cofficients for d(uv/dx)
rcovxu1 = gt/dx(i)*(u1(i,j,k)*dy(j1d-1)+u1(i,j-1,k)*dy(j1d)) /(dy(j1d)+dy(j1d-1))
rcovxu2 = gt/dx(i)*(u1(i+1,j,k)*dy(j1d-1)+u1(i+1,j-1,k)*dy(j1d))/(dy(j1d)+dy(j1d-1))
a2xuv = -rcovxu1*dx(i)/(dx(i)+dx(i-1))
b2xuv =  rcovxu2*dx(i+1)/(dx(i)+dx(i+1)) -rcovxu1*dx(i-1)/(dx(i)+dx(i-1))
c2xuv = rcovxu2*dx(i)/(dx(i)+dx(i+1))
If (i == 1.and.ibwest == 1) Then
a2xuv = -rcovxu1
b2xuv =  rcovxu2*dx(i+1)/(dx(i)+dx(i+1))
End If
If (i == nx.and.ibeast == 1) Then
b2xuv = -rcovxu1*dx(i-1)/(dx(i)+dx(i-1))
c2xuv =  rcovxu2
End If
!
a(i,j,k) = bt*nutxy1*gcdx_f(i,4,1,ib)
b(i,j,k) = bt*(nutxy1*gcdx_f(i,3,1,ib)-nutxy2* gcdx_f(i+1,4,1,ib))+(1.0e0-ivx_con(i,1))*a2xuv
c(i,j,k) = 1.0e0+bt*(nutxy1*gcdx_f(i,2,1,ib)-nutxy2* gcdx_f(i+1,3,1,ib))+b2xuv
d(i,j,k) = bt*(nutxy1*gcdx_f(i,1,1,ib)-nutxy2* gcdx_f(i+1,2,1,ib))+(1.0e0-ivx_con(i,2))*c2xuv
e(i,j,k) = -bt*nutxy2*gcdx_f(i+1,1,1,ib)
End Do; End Do; End Do
!
 Call penta1(a,b,c,d,e,intf1,beg,1)
!
Do k = 1,xsize(3); Do j = beg(2),xsize(2); Do i = 1,xsize(1)
du2(i,j,k) = intf1(i,j,k)
End Do; End Do; End Do
!
! Periodic BCs
  call update_halo(du2,2)
!
! Wall BCs
Do k = 0,xsize(3)+1; Do i = 0,xsize(1)+1
If (xmin(2) == 1) then; du2(i,1,k) = 0.0e0; du2(i,0,k) = 0.0e0; End If
If (xmax(2) == 1) du2(i,xsize(2)+1,k) = 0.0e0
End Do; End Do
!	Gets missing zeros
If (xmin(3)==1) then
  du2(0,:,0) = 0.0e0
  du2(nx+1,:,0) = 0.0e0
End If
If (xmax(3)==1) then
  du2(0,:,xsize(3)+1) = 0.0e0
  du2(xsize(1)+1,:,xsize(3)+1) = 0.0e0
End If
!
!*****************************
Else If(icmpnt == 3) Then
!*****************************
Do k = beg(3),xsize(3)
 k1d=xstart(3)-1+k
Do j = beg(2),xsize(2)
Do i = beg(1),xsize(1)
!
!	Cofficints for  (d2w/dx2)
nutxz1 = buff(i)*bre(i)*(nu+rnutxz1(i,j,k))/dx(i)
nutxz2 = buff(i+1)*bre(i+1)*(nu+rnutxz1(i+1,j,k))/dx(i)
!	Cofficients for d(uw)/dx
rcowxu1 = gt/dx(i)/(dz(k1d)+dz(k1d-1))*(u1(i,j,k)*dz(k1d-1)+u1(i,j,k-1)*dz(k1d))
rcowxu2 = gt/dx(i)/(dz(k1d)+dz(k1d-1))*(u1(i+1,j,k)*dz(k1d-1)+u1(i+1,j,k-1)*dz(k1d))
a3xuw = -rcowxu1*dx(i)/(dx(i)+dx(i-1))
b3xuw =  rcowxu2*dx(i+1)/(dx(i)+dx(i+1))-rcowxu1*dx(i-1)/(dx(i)+dx(i-1))
c3xuw =  rcowxu2*dx(i)/(dx(i)+dx(i+1))
If (i == 1.and.ibwest == 1) Then
a3xuw = -rcowxu1
b3xuw =  rcowxu2*dx(i+1)/(dx(i)+dx(i+1))
End If
If (i == nx.and.ibeast == 1) Then
b3xuw = -rcowxu1*dx(i-1)/(dx(i)+dx(i-1))
c3xuw =  rcowxu2
End If
!
a(i,j,k) = bt*nutxz1*gcdx_f(i,4,1,ib)
b(i,j,k) = bt*(nutxz1*gcdx_f(i,3,1,ib)-nutxz2*gcdx_f(i+1,4,1,ib)) +a3xuw*(1.0e0-ivx_con(i,1))
c(i,j,k) = 1.0e0+bt*(nutxz1*gcdx_f(i,2,1,ib)-nutxz2*gcdx_f(i+1,3,1,ib))+b3xuw
d(i,j,k) = bt*(nutxz1*gcdx_f(i,1,1,ib)-nutxz2*gcdx_f(i+1,2,1,ib)) +c3xuw*(1.0e0-ivx_con(i,2))
e(i,j,k) = -bt*nutxz2*gcdx_f(i+1,1,1,ib)
End Do; End Do; End Do
!
!       Periodic B.C.
!
 Call penta1(a,b,c,d,e,intf1,beg,1)
!
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
 what1(i,j,k) = intf1(i,j,k)+w1(i,j,k)
End Do; End Do; End Do
!
!	Caculation of uhat and vhat from du1 and du2
 Call uvhat_calc
!
End If
!
Deallocate(a,b,c,d,e)
!
Return
End Subroutine solvex_implicit    
!
!******************************************************************************
Subroutine rhs2_calc
!******************************************************************************
!	This procedure calculates the right hand side of V-Momentum equation by
!	the decoupling procedure by kim et.al. 2002
Implicit none
Integer :: j1d
!
!**************** must be in X pencil *******************
Do k = beg(3),xsize(3)
Do j = beg(2),xsize(2)
 j1d=xstart(2)-1+j
Do i = beg(1),xsize(1)
!
!	Cofficients for d2(du1**)/dx2
nutc1 = buffc(i-1)*brec(i-1)*(nu+2.0e0*pitero1(i-1,j,k))/dx(i-1)/dxs(i)
nutc2 = buffc(i)*brec(i)*(nu+2.0e0*pitero1(i,j,k))/dx(i)/dxs(i)
!	Cofficients for d(vdu1**)/dx
rcovxv1 = gt/dx(i)/(dy(j1d)+dy(j1d-1))/(dx(i)+dx(i-1))                &
*(v1(i,j,k)*dx(i-1)+v1(i-1,j,k)*dx(i))
rcovxv2 = gt/dx(i)/(dy(j1d)+dy(j1d-1))/(dx(i)+dx(i+1))                &
*(v1(i+1,j,k)*dx(i)+v1(i,j,k)*dx(i+1))
If (i == 1.and.ibwest == 1)&
rcovxv1 = gt/dx(i)/(dy(j1d)+dy(j1d-1))*v1(0,j,k)
If (i == nx.and.ibeast == 1)& 
rcovxv2 = gt/dx(i)/(dy(j1d)+dy(j1d-1))*v1(nx+1,j,k)
!
 intf1(i,j,k) = intf1(i,j,k)                                        &
+(1.0e0-ivx_con(i,1))*rcovxv1*(du1(i,j,k)*dy(j1d-1)+                  &
du1(i,j-1,k)*dy(j1d))                                             &
-(1.0e0-ivx_con(i,2))*rcovxv2*(du1(i+1,j,k)*dy(j1d-1)+                &
du1(i+1,j-1,k)*dy(j1d))
End Do;End Do;End Do
!
Return
End Subroutine rhs2_calc
!
!******************************************************************************
Subroutine rhs3_calc
!******************************************************************************
!	This procedure calculates the Righ hand side of W-Momentum equation by
!	the decoupling procedure by kim et.al. 2002
Implicit none
Integer :: j1d,k1d
!
!**************** must be in X pencil *******************
Do k = beg(3),xsize(3)
 k1d=xstart(3)-1+k
Do j = beg(2),xsize(2)
 j1d=xstart(2)-1+j
Do i = beg(1),xsize(1)
!
!	Cofficients for d(wdu1**)/dx
rcowxw1 = gt/dx(i)/(dz(k1d)+dz(k1d-1))/(dx(i)+dx(i-1))                  &
*(w1(i,j,k)*dx(i-1)+w1(i-1,j,k)*dx(i))
rcowxw2 = gt/dx(i)/(dz(k1d)+dz(k1d-1))/(dx(i)+dx(i+1))                  &
*(w1(i+1,j,k)*dx(i)+w1(i,j,k)*dx(i+1))
If (i == 1.and.ibwest == 1)&
rcowxw1 = gt/dx(i)/(dz(k1d)+dz(k1d-1))*w1(0,j,k)
If (i == nx.and.ibeast == 1)& 
rcowxw2 = gt/dx(i)/(dz(k1d)+dz(k1d-1))*w1(nx+1,j,k)
!	Cofficients for d(wdu2**)/dy
rcowyw1 = gt/dy(j1d)/(dz(k1d)+dz(k1d-1))/(dy(j1d)+dy(j1d-1))*                &
(w1(i,j,k)*dy(j1d-1)+w1(i,j-1,k)*dy(j1d))
rcowyw2 = gt/dy(j1d)/(dz(k1d)+dz(k1d-1))/(dy(j1d)+dy(j1d+1))*                &
(w1(i,j+1,k)*dy(j1d)+w1(i,j,k)*dy(j1d+1))
  If (j == 1.and.ibsoth == 1.and.xmin(2) == 1)& 
rcowyw1 = gt/dy(j1d)/(dz(k1d)+dz(k1d-1))*w1(i,0,k)
If (j == xsize(2).and.ibnoth == 1.and.xmax(2) == 1)& 
rcowyw2 = gt/dy(j1d)/(dz(k1d)+dz(k1d-1))*w1(i,xsize(2)+1,k)
!
intf1(i,j,k) = intf1(i,j,k)                                           &
+(1.0e0-ivx_con(i,1))*rcowxw1*(du1(i,j,k)*dz(k1d-1)+                    &
du1(i,j,k-1)*dz(k1d))                                               &
-(1.0e0-ivx_con(i,2))*rcowxw2*(du1(i+1,j,k)*dz(k1d-1)+                  &
du1(i+1,j,k-1)*dz(k1d))                                             &
!
+(1.0e0-iuy_con(j1d,1))*rcowyw1*(du2(i,j,k)*dz(k1d-1)+                    &
du2(i,j,k-1)*dz(k1d))                                               &
-(1.0e0-iuy_con(j1d,2))*rcowyw2*(du2(i,j+1,k)*dz(k1d-1)+                  &
du2(i,j+1,k-1)*dz(k1d))   
End Do;End Do;End Do
!
Return
End Subroutine rhs3_calc
!
!******************************************************************************
        Subroutine uvhat_calc
!******************************************************************************
!	This subroutines determines the uhat and what from du3*
Implicit none
Integer :: j1d,k1d
Real(mytype), dimension(:,:,:), allocatable :: dvhat,duhat,intf1_halo
!  
Allocate(dvhat(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2), &
	 duhat(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2), &
         intf1_halo(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2))
!
 intf1_halo(1:xsize(1),1:xsize(2),1:xsize(3))=intf1(1:xsize(1),1:xsize(2),1:xsize(3))
 Call update_halo(intf1_halo,2)
!
!	vhat calculation
!
!dvhat = 0.0e0;duhat = 0.0e0
!ibeg = 1; jbeg = 1
!If (xmin(2) == 1) jbeg = 2
!
Do k = beg(3),xsize(3)
 k1d=xstart(3)-1+k
Do j = beg(2),xsize(2)
 j1d=xstart(2)-1+j
Do i = beg(1),xsize(1)
!	Cofficients for d(vdu3*)/dz
rcovzv1 = gt/dz(k1d)/(dy(j1d)+dy(j1d-1))/(dz(k1d)+dz(k1d-1))*               &
(v1(i,j,k)*dz(k1d-1)+v1(i,j,k-1)*dz(k1d))
rcovzv2 = gt/dz(k1d)/(dy(j1d)+dy(j1d-1))/(dz(k1d)+dz(k1d+1))*               &
(v1(i,j,k+1)*dz(k1d)+v1(i,j,k)*dz(k1d+1))
dvhat(i,j,k) = du2(i,j,k)                                                   &
+rcovzv1*(intf1_halo(i,j,k)*dy(j1d-1)+intf1_halo(i,j-1,k)*dy(j1d))          &
-rcovzv2*(intf1_halo(i,j,k+1)*dy(j1d-1)+intf1_halo(i,j-1,k+1)*dy(j1d))
End Do;End Do;End Do
!
Do k = beg(3),xsize(3); Do j = beg(2),xsize(2); Do i = beg(1),xsize(1)
vhat1(i,j,k) = dvhat(i,j,k)+v1(i,j,k)
End Do; End Do; End Do
!
!	Periodic Boundary conditon
!
  Call update_halo(dvhat,2)
!
!	uhat calculation
!
!ibeg = 1; jbeg = 1
!If (xmin(1) == 1) ibeg = 2
!If (iconbc == 0) ibeg = 1
!
Do k = beg(3),xsize(3)
 k1d=xstart(3)-1+k
Do j = beg(2),xsize(2)
 j1d=xstart(2)-1+j
Do i = beg(1),xsize(1)
!
!	Cofficients  for d(udu2*)/dy
rcouyu1 = gt/dy(j1d)*(u1(i,j,k)*dy(j1d-1)+u1(i,j-1,k)*dy(j1d))              &
/(dy(j1d)+dy(j1d-1))/(dx(i)+dx(i-1))
rcouyu2 = gt/dy(j1d)*(u1(i,j+1,k)*dy(j1d)+u1(i,j,k)*dy(j1d+1))              &
/(dy(j1d)+dy(j1d+1))/(dx(i)+dx(i-1))
!	Cofficients  for d(udu3*)/dz
rcouzu1 = gt/dz(k1d)/(dx(i)+dx(i-1))/(dz(k1d)+dz(k1d-1))*                   &
(u1(i,j,k)*dz(k1d-1)+u1(i,j,k-1)*dz(k1d))
rcouzu2 = gt/dz(k1d)/(dx(i)+dx(i-1))/(dz(k1d)+dz(k1d+1))*                   &
(u1(i,j,k+1)*dz(k1d)+u1(i,j,k)*dz(k1d+1))
duhat(i,j,k) = du1(i,j,k)                                                   &
+(1.0e0-iuy_con(j1d,1))*rcouyu1*(dvhat(i,j,k)*dx(i-1)+dvhat(i-1,j,k)*dx(i))     &
-(1.0e0-iuy_con(j1d,2))*rcouyu2*(dvhat(i,j+1,k)*dx(i-1)+dvhat(i-1,j+1,k)*dx(i)) &
+rcouzu1*(intf1_halo(i,j,k)*dx(i-1)+intf1_halo(i-1,j,k)*dx(i))              &
-rcouzu2*(intf1_halo(i,j,k+1)*dx(i-1)+intf1_halo(i-1,j,k+1)*dx(i))
End Do;End Do;End Do
!
Do k = beg(3),xsize(3); Do j = beg(2),xsize(2); Do i = beg(1),xsize(1)
uhat1(i,j,k) = duhat(i,j,k)+u1(i,j,k)
End Do; End Do; End Do
Deallocate(dvhat,duhat,intf1_halo)
!
Return
End Subroutine uvhat_calc
!
End Module compute_left_implicit
