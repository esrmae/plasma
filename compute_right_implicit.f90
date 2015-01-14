MODULE compute_right_implicit
!No LES or Block/IBM at the moment 01/03/2010
Use shared_data
Use boundary_conditions
!
Implicit none
integer, private :: i,j,k,ibeg,ib,jbeg,iblock,it,is,iuvwc,kbeg,j1d,k1d
Real(mytype), dimension (:,:,:),allocatable :: uh1
Real(mytype),dimension(:,:),allocatable :: dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, &
dudz2,dvdz2,dwdz2,dwdx2,dwdy2,vv,ww,uc,vc,wc,uc2,vc2,wc2,uv,vw,uw,uw2,      &
vw2,ww2,uu,dudx,hij1,hij2,hij3,hij4,hij5
!
Contains
!	
!******************************************************************************
Subroutine compute_rhs_implicit
!*****************************************************************************
!------ This procedure computes convection, diffusion and pressure integrals
!       using known velocity field.
Implicit none
!
Allocate(hij1(xsize(1),xsize(2)),hij2(xsize(1),xsize(2)),hij3(xsize(1),xsize(2)),hij4(xsize(1),xsize(2)),hij5(xsize(1),xsize(2)),&
uh1(xsize(1),xsize(2),xsize(3)),dudy(xsize(1)+1,xsize(2)+1),dudz(xsize(1)+1,xsize(2)+1),dvdx(xsize(1)+1,xsize(2)+1),     &
dvdy(xsize(1)+1,0:xsize(2)+1),dvdz(xsize(1)+1,xsize(2)+1),dwdx(xsize(1)+1,xsize(2)+1),dwdy(xsize(1)+1,xsize(2)+1),     &
dwdz(xsize(1)+1,xsize(2)+1),dudz2(xsize(1)+1,xsize(2)+1),dvdz2(xsize(1)+1,xsize(2)+1),dwdz2(xsize(1)+1,xsize(2)+1),  &
dwdx2(xsize(1)+1,xsize(2)+1),dwdy2(xsize(1)+1,xsize(2)+1),vv(xsize(1)+1,0:xsize(2)+1),ww(xsize(1)+1,xsize(2)+1),       &
uc(xsize(1)+1,xsize(2)+1),vc(xsize(1)+1,xsize(2)+1),wc(xsize(1)+1,xsize(2)+1),uc2(xsize(1)+1,xsize(2)+1),            &
vc2(xsize(1)+1,xsize(2)+1),wc2(xsize(1)+1,xsize(2)+1),uv(xsize(1)+1,xsize(2)+1),vw(xsize(1)+1,xsize(2)+1),           &
uw(xsize(1)+1,xsize(2)+1),uw2(xsize(1)+1,xsize(2)+1), vw2(xsize(1)+1,xsize(2)+1),ww2(xsize(1)+1,xsize(2)+1),         &
dudx(0:xsize(1),xsize(2)),uu(0:xsize(1)+1,xsize(2)+1))
!	
intf1= 0.0e0
iblock= 0;  ib= iblock+1
!If (iconbc ==  1 .and. if2d /=  1) Call bc_right
!
If (icmpnt ==  1)  Then
  Call compute_rhs_u1 
Elseif (icmpnt ==  2) Then
  Call compute_rhs_v1 
Elseif (icmpnt ==  3)  Then
  Call compute_rhs_w1 
End If
!
Deallocate(hij1,hij2,hij3,hij4,hij5,uh1,dudy,dudz,dvdx,dvdy,dvdz,  &
dwdx,dwdy,dwdz,dudz2,dvdz2,dwdz2,dwdx2,dwdy2,vv,ww,uc,vc,wc,uc2,vc2, &
wc2,uv,vw,uw,uw2,vw2,ww2,dudx,uu)
!
Return
End subroutine compute_rhs_implicit
!
!******************************************************************************
Subroutine compute_rhs_u1
!******************************************************************************
!	This procedure calculates the right hand side of NS equation for u1 momentum equation
Implicit none
!
iuvwc= 1;  ibeg= 1; jbeg= 1; kbeg= 1
If ((iconbc ==  1).and.(xmin(1) == 1)) ibeg= 2
!
!       Viscous contribution
!
!
Do k= kbeg,xsize(3)
 k1d=xstart(3)-1+k
 Call gradient_x(dudx,gfdx_c,u1,ibeg-1,xsize(1),jbeg,xsize(2),0,xsize(1),1,xsize(2),0,nx+1,2,iuvwc)
 Call gradient_y(dudy,gcdy_f,u1,ibeg,xsize(1),jbeg,xsize(2)+1,1,xsize(1)+1,1,xsize(2)+1,1,ny+1,1,iuvwc)
 Call gradient_z(dudz,gcdz_f,u1,ibeg,xsize(1),jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nz+1,1,iuvwc,0)
 Call gradient_z(dudz2,gcdz_f2,u1,ibeg,xsize(1),jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nz+1,2,iuvwc,1)
 Call visc_cont(ibeg,xsize(1),jbeg,xsize(2))  
!
!       Convective Contribution
!
 Call convective_ii
 Call convective_ij
 Call convective_jk
 Call conv_cont(ibeg,xsize(1),jbeg,xsize(2))               
 Call intfo_calc(ibeg,xsize(1),jbeg,xsize(2))                 
!
!       Total right hand side for u1 momentum equation
!
  Call intf_calc(ibeg,xsize(1),jbeg,xsize(2))         
!
!	MBC contribution
!
  If (if_mbc ==  1) then
    Call mbc_dir_visc
    Call mbc_dir_conv
  End If
End Do
!
Return
End Subroutine compute_rhs_u1
!
!******************************************************************************
subroutine compute_rhs_v1
!*****************************************************************************
!	This procedure calculates the RHS of NS for v1 momentum equation
implicit none
!
iuvwc= 2; ibeg= 1; jbeg= 1; kbeg= 1
If ((jconbc == 1).and.(xmin(2) == 1)) jbeg= 2 
!
!       Viscous contribution
!
Do k= kbeg,xsize(3)
 k1d=xstart(3)-1+k
 Call gradient_x(dvdx,gcdx_f,v1,ibeg,xsize(1)+1,jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nx+1,1,iuvwc)
 Call gradient_y(dvdy,gfdy_c,v1,ibeg,xsize(1),jbeg-1,xsize(2),1,xsize(1)+1,0,xsize(2)+1,1,ny+1,2,iuvwc)
 Call gradient_z(dvdz,gcdz_f,v1,ibeg,xsize(1),jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nz+1,1,iuvwc,0)
 Call gradient_z(dvdz2,gcdz_f2,v1,ibeg,xsize(1),jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nzt,2,iuvwc,1)
 Call visc_cont(ibeg,xsize(1),jbeg,xsize(2))               
!
!       Convective contribution
!
 Call convective_ij
 Call convective_ii
 Call convective_jk
 Call conv_cont(ibeg,xsize(1),jbeg,xsize(2))               
 Call intfo_calc(ibeg,xsize(1),jbeg,xsize(2))               
!
!       Total right hand side for v1 momentum equation
!
 Call intf_calc(ibeg,xsize(1),jbeg,xsize(2))               
!
!	MBC
!
  If (if_mbc ==  1) then
    Call mbc_dir_visc
    Call mbc_dir_conv
  End If
End Do
!
Return
End Subroutine compute_rhs_v1
!
!******************************************************************************
Subroutine compute_rhs_w1
!*****************************************************************************
!	This procedure calculates the RHS of NS for w1 momentum equation
implicit none
!
iuvwc= 3;  ibeg= 1;  jbeg= 1;  kbeg= 1
If ((kconbc ==  1).and.(xmin(3) == 1)) kbeg= 2
!
!       Viscous contribution
!
Do k= kbeg,xsize(3)
 k1d=xstart(3)-1+k
  Call gradient_x(dwdx,gcdx_f,w1,ibeg,xsize(1)+1,jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nx+1,1,iuvwc)
  Call gradient_y(dwdy,gcdy_f,w1,ibeg,xsize(1),jbeg,xsize(2)+1,1,xsize(1)+1,1,xsize(2)+1,1,ny+1,1,iuvwc)
  Call gradient_z(dwdz,gfdz_c,w1,ibeg,xsize(1),jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nz+1,2,iuvwc,0)
  Call gradient_z(dwdz2,gfdz_c,w1,ibeg,xsize(1),jbeg,xsize(2),1,xsize(1)+1,1,xsize(2)+1,1,nz+1,1,iuvwc,-1)
  Call visc_cont(ibeg,xsize(1),jbeg,xsize(2))               
!
!       Convective contribution
!
  Call convective_ij
  Call convective_jk
  Call convective_ii
  Call conv_cont(ibeg,xsize(1),jbeg,xsize(2))               
  Call intfo_calc(ibeg,xsize(1),jbeg,xsize(2))               
!
!       Total right hand side for v1 momentum equation
!
  Call intf_calc(ibeg,xsize(1),jbeg,xsize(2))               
!
!	MBC
  If (if_mbc ==  1) then
    Call mbc_dir_visc
    Call mbc_dir_conv
  End If
End Do
!
Return
End subroutine compute_rhs_w1
!
!******************************************************************************
Subroutine gradient_x(temp1,temp2,temp3,mm,nn,i1,jj,m1,m2,n1,n2,k1,k2,ind,iuvw)
!******************************************************************************
!       This subroutine calculates gradient in x_direction at cell center and cell face in z-plane
Implicit none
Integer, Intent(In) :: mm,nn,i1,jj,m1,m2,n1,n2,k1,k2,ind,iuvw
Real(mytype),dimension(m1:m2,n1:n2),Intent(In Out) :: temp1
Real(mytype),dimension(k1:k2,4,ngr,nbl),Intent(In) :: temp2
Real(mytype),dimension(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),intent(in) :: temp3
!
If (iuvw ==  1) Then
 Do j= i1,jj;  Do i= mm,nn
 temp1(i,j)= temp2(i,1,1,1)*temp3(i+ind,j,k)+             &
 (1.0e0-iux_grdvis(i,2))*temp2(i,2,1,1)*temp3(i+ind-1,j,k)+    &
 (1.0e0-iux_grdvis(i,1))*temp2(i,3,1,1)*temp3(i+ind-2,j,k)+    &
 temp2(i,4,1,1)*temp3(i+ind-3,j,k)
 temp1(i,j)= (2.0e0*pitero1(i,j,k)+nu)*temp1(i,j)
End Do;  End Do
!
!*************************************
        Else If (iuvw ==  2) Then
!************************************
!
Do j= i1,jj;  Do i= mm,nn
temp1(i,j)= (1.0e0-ivx_grdvis(i,2,2))*temp2(i,1,1,1)*temp3(i+ind,j,k)+      &
(1.0e0-ivx_grdvis(i,1,2))*temp2(i,2,1,1)*temp3(i+ind-1,j,k)+               &
(1.0e0-ivx_grdvis(i,1,1))*temp2(i,3,1,1)*temp3(i+ind-2,j,k)+               &
(1.0e0-ivx_grdvis(i,2,1))*temp2(i,4,1,1)*temp3(i+ind-3,j,k)
temp1(i,j)= (rnutxy1(i,j,k)+nu)*temp1(i,j)
End Do;  End Do
!
!*************************************
        Else If (iuvw ==  3) Then
!*************************************
!
Do j= i1,jj;  Do i= mm,nn
temp1(i,j)= (1.0e0-ivx_grdvis(i,2,2))*temp2(i,1,1,1)*temp3(i+ind,j,k)+ &
(1.0e0-ivx_grdvis(i,1,2))*temp2(i,2,1,1)*temp3(i+ind-1,j,k)+&
(1.0e0-ivx_grdvis(i,1,1))*temp2(i,3,1,1)*temp3(i+ind-2,j,k)+           &
(1.0e0-ivx_grdvis(i,2,1))*temp2(i,4,1,1)*temp3(i+ind-3,j,k)
temp1(i,j)= (rnutxz1(i,j,k)+nu)*temp1(i,j)
End Do;  End Do
!
End If
!
Return
End Subroutine gradient_x
!                                                                                                                                
!******************************************************************************
Subroutine gradient_y(temp1,temp2,temp3,mm,nn,i1,jj,m1,m2,n1,n2,k1,k2,ind,iuvw)
!******************************************************************************
Implicit none
Integer, Intent(In) :: mm,nn,i1,jj,m1,m2,n1,n2,k1,k2,ind,iuvw
Real(mytype),dimension(m1:m2,n1:n2),Intent(In Out) :: temp1
Real(mytype),dimension(k1:k2,4,ngr,nbl),Intent(In) :: temp2
Real(mytype),dimension(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),Intent(In) :: temp3
!
If (iuvw ==  1) Then
Do j= i1,jj; j1d=xstart(2)-1+j;  Do i= mm,nn
temp1(i,j)= (1.0e0-iuy_grdvis(j1d,2,2))*temp2(j1d,1,1,1)*temp3(i,j+ind,k)&
+(1.0e0-iuy_grdvis(j1d,1,2))*temp2(j1d,2,1,1)* temp3(i,j+ind-1,k)+      &
(1.0e0-iuy_grdvis(j1d,1,1))* temp2(j1d,3,1,1)*temp3(i,j+ind-2,k)+       &
(1.0e0-iuy_grdvis(j1d,2,1))*temp2(j1d,4,1,1)*temp3(i,j+ind-3,k)
temp1(i,j)= (rnutxy1(i,j,k)+nu)*temp1(i,j)
End Do;  End Do
!
!*****************************************
        Else If (iuvw ==  2) Then
!*****************************************
!
Do j= i1,jj; j1d=xstart(2)-1+j;  Do i= mm,nn
temp1(i,j)= temp2(j1d,1,1,1)*temp3(i,j+ind,k)+                &
(1.0e0-ivy_grdvis(j1d,2))*temp2(j1d,2,1,1)*temp3(i,j+ind-1,k)+       &
(1.0e0-ivy_grdvis(j1d,1))*temp2(j1d,3,1,1)*temp3(i,j+ind-2,k)+       &
temp2(j1d,4,1,1)*temp3(i,j+ind-3,k)
temp1(i,j)= (2.0e0*pitero1(i,j,k)+nu)*temp1(i,j)
End Do;  End Do
!
!*********************************************
        Else If (iuvw ==  3) Then
!*********************************************
!
Do j= i1,jj; j1d=xstart(2)-1+j;  Do i= mm,nn
temp1(i,j)= (1.0e0-iuy_grdvis(j1d,2,2))*temp2(j1d,1,1,1)*temp3(i,j+ind,k)+&
(1.0e0-iuy_grdvis(j1d,1,2))*temp2(j1d,2,1,1)*temp3(i,j+ind-1,k)+         &
(1.0e0-iuy_grdvis(j1d,1,1))*temp2(j1d,3,1,1)*temp3(i,j+ind-2,k)+             &
(1.0e0-iuy_grdvis(j1d,2,1))*temp2(j1d,4,1,1)*temp3(i,j+ind-3,k)
temp1(i,j)= (rnutyz1(i,j,k)+nu)*temp1(i,j)
End Do;  End Do
!
End If
!
Return
End Subroutine gradient_y
!
!******************************************************************************
Subroutine gradient_z(temp1,temp2,temp3,i1,i2,j1,j2,m1,m2,n1,n2,k1,k2,ind,iuvw,indx)
!******************************************************************************
Implicit none
Integer, Intent(In) :: i1,i2,j1,j2,m1,m2,n1,n2,k1,k2,ind,iuvw,indx
Real(mytype),dimension(m1:m2,n1:n2),Intent(In Out) :: temp1
Real(mytype),dimension(k1:k2,4,ngr,nbl),Intent(In) :: temp2
Real(mytype),dimension(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),Intent(In) :: temp3
!
outer : If (iuvw ==  1) Then
Do j= j1,j2; Do i= i1,i2
temp1(i,j)= temp2(k1d,1,1,1)*temp3(i,j,k+ind)+temp2(k1d,2,1,1)*temp3(i,j,k+ind-1)+&
temp2(k1d,3,1,1)*temp3(i,j,k+ind-2)+temp2(k1d,4,1,1)*temp3(i,j,k+ind-3)
temp1(i,j)= (rnutxz1(i,j,k+indx)+nu)*temp1(i,j)
End Do;  End Do
!
!*******************************
Else If (iuvw ==  2) Then
!*******************************
!
Do j= j1,j2; Do i= i1,i2
temp1(i,j)= temp2(k1d,1,1,1)*temp3(i,j,k+ind)+temp2(k1d,2,1,1)*temp3(i,j,k+ind-1)+&
temp2(k1d,3,1,1)*temp3(i,j,k+ind-2)+ temp2(k1d,4,1,1)*temp3(i,j,k+ind-3)
temp1(i,j)= (rnutyz1(i,j,k+indx)+nu)*temp1(i,j)
End Do;  End Do
!
!*******************************
Else If (iuvw ==  3) Then
!*******************************
!
Do j= j1,j2; Do i= i1,i2
temp1(i,j)= temp2(k1d,1,1,1)*temp3(i,j,k+ind)+temp2(k1d,2,1,1)*temp3(i,j,k+ind-1)+&
temp2(k1d,3,1,1)*temp3(i,j,k+ind-2)+temp2(k1d,4,1,1)*temp3(i,j,k+ind-3)
temp1(i,j)= (2.0e0*pitero1(i,j,k+indx)+nu)*temp1(i,j)
End Do;  End Do
!
End If outer
!
Return
End Subroutine gradient_z
! 
!******************************************************************************
Subroutine visc_cont(ii1,ii2,jj1,jj2)
!*****************************************************************************
!	This procedure calculates the viscous contribution on RHS of NS equation.
Implicit none
Integer, Intent(In) :: ii1,ii2,jj1,jj2
Integer :: iii
!
outer : If (iuvwc ==  1) Then
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
hij1(i,j)= (dudx(i,j)-dudx(i-1,j))
hij2(i,j)= (dudy(i,j+1)-dudy(i,j))*dxs(i)/dy(j1d)
hij3(i,j)= (dudz2(i,j)-dudz(i,j))*dxs(i)/dz(k1d)
End Do;  End Do
!
Do j= jj1,jj2;  Do i= ii1,ii2
  uh1(i,j,k)= (hij1(i,j)*buff(i)+hij2(i,j)+hij3(i,j))*bre(i)
End Do;  End Do
!
!******************************
Else If (iuvwc ==  2) Then
!******************************
!
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
hij1(i,j)= (dvdx(i+1,j)-dvdx(i,j))
hij2(i,j)= (dvdy(i,j)-dvdy(i,j-1))*dx(i)/dys(j1d)
hij3(i,j)= (dvdz2(i,j)-dvdz(i,j))*dx(i)/dz(k1d)
End Do;  End Do
!
Do j= jj1,jj2;  Do i= ii1,ii2
  uh1(i,j,k)= (hij1(i,j)*buff(i)+hij2(i,j)+hij3(i,j))*bre(i)
End Do;  End Do
!
!******************************
Else If (iuvwc ==  3) Then
!******************************
!
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
hij1(i,j)= (dwdx(i+1,j)-dwdx(i,j))
hij2(i,j)= (dwdy(i,j+1)-dwdy(i,j))*dx(i)/dy(j1d)
hij3(i,j)= (dwdz(i,j)-dwdz2(i,j))*dx(i)/dzs(k1d)
End Do;  End Do
!
Do j= jj1,jj2;  Do i= ii1,ii2
uh1(i,j,k)= (hij1(i,j)*buffc(i)+hij2(i,j)+hij3(i,j))*brec(i)
End Do;  End Do
!
End If outer
!
Return
End Subroutine visc_cont
!
!******************************************************************************
Subroutine convective_ii
!******************************************************************************
!       This procedure calculates the integral of  d(UU)/dx d(VV)/dy and d(WW)/dz
!       contribution from viscous part on RHS of the NS equation
Implicit none
!
If (iuvwc ==  1) Then
Do j= jbeg,xsize(2);  Do i= ibeg-1,xsize(1)
uu(i,j)= 0.25e0*(u1(i,j,k)+u1(i+1,j,k))&
 *((1.0e0-iux_grdcon(i,1))*u1(i,j,k)+(1.0e0-iux_grdcon(i,2))*u1(i+1,j,k))
End Do;  End Do
!
!
!********************************
        Else If (iuvwc ==  2) Then
!*******************************
!
Do j= jbeg-1,xsize(2); j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)
vv(i,j)= 0.25e0*(v1(i,j,k)+v1(i,j+1,k))*&
((1.0e0-ivy_grdcon(j1d,1))*v1(i,j,k)+(1.0e0-ivy_grdcon(j1d,2))*v1(i,j+1,k))
End Do;  End Do
!
!*********************************
        Else If (iuvwc ==  3) Then
!********************************
!
Do j= jbeg,xsize(2);  Do i= ibeg,xsize(1)
ww(i,j)= 0.25e0*(w1(i,j,k)+w1(i,j,k+1))*(w1(i,j,k)+w1(i,j,k+1))
ww2(i,j)= 0.25e0*(w1(i,j,k-1)+w1(i,j,k))*(w1(i,j,k-1)+w1(i,j,k))
End Do;  End Do
!
End If
!
Return
End Subroutine convective_ii
!
!******************************************************************************
Subroutine convective_ij
!******************************************************************************
!       This procedure calculates the integral of  d(UV)/dy d(VU)/dx and d(UW)/dx
!       contribution from viscous part on RHS of the NS equation
Implicit none
Integer :: iii
!
If (iuvwc ==  1) Then
Do j= jbeg,xsize(2)+1; j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)
uc(i,j)= ((1.0e0-iuy_grdcon(j1d,2))*u1(i,j,k)*dy(j1d-1)+&
(1.0e0-iuy_grdcon(j1d,1))*u1(i,j-1,k)*dy(j1d))/(dy(j1d-1)+dy(j1d))
vc(i,j)= (v1(i,j,k)*dx(i-1)+v1(i-1,j,k)*dx(i))/(dx(i-1)+dx(i))
uv(i,j)= uc(i,j)*vc(i,j)
End Do;  End Do
!
If (jconbc ==  1) Then
If(ibsoth ==  1.and.xmin(2)==1) Then
j= 1
 j1d=xstart(2)-1+j
Do i= ibeg,xsize(1)
uc(i,j)= (1.0e0-iuy_grdcon(j1d,1))*u1(i,0,k)
uv(i,j)= uc(i,j)*vc(i,j)
End Do
End If
!
If(ibnoth ==  1.and.xmax(2)==1) Then
j= xsize(2)+1
 j1d=xstart(2)-1+j
Do i= ibeg,xsize(1)
uc(i,j)= (1.0e0-iuy_grdcon(j1d,2))*u1(i,xsize(2)+1,k)
uv(i,j)= uc(i,j)*vc(i,j)
End Do
End If
End If
!
!
!************************************
        Else If (iuvwc ==  2) Then
!***********************************
!
Do j= jbeg,xsize(2); j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)+1
uc(i,j)= (u1(i,j,k)*dy(j1d-1)+u1(i,j-1,k)*dy(j1d))/(dy(j1d-1)+dy(j1d))
vc(i,j)= ((1.0e0-ivx_grdcon(i,2))*v1(i,j,k)*dx(i-1)+&
(1.0e0-ivx_grdcon(i,1))*v1(i-1,j,k)*dx(i))/(dx(i-1)+dx(i))
uv(i,j)= uc(i,j)*vc(i,j)
End Do;  End Do
!
!*********************************
        Else If (iuvwc ==  3) Then
!********************************
!
Do j= jbeg,xsize(2);  Do i= ibeg,xsize(1)+1
uc(i,j)= (u1(i,j,k-1)*dz(k1d)+u1(i,j,k)*dz(k1d-1))/(dz(k1d-1)+dz(k1d))
wc(i,j)= ((1.0e0-ivx_grdcon(i,2))*w1(i,j,k)*dx(i-1)+&
(1.0e0-ivx_grdcon(i,1))*w1(i-1,j,k)*dx(i))/(dx(i-1)+dx(i))
uw(i,j)= uc(i,j)*wc(i,j)
End Do;  End Do
!
End If
!
Return
End Subroutine convective_ij
!
!******************************************************************************
Subroutine convective_jk
!******************************************************************************
!       This procedure calculates the integral of  d(UW)/dz d(VW)/dz and d(WV)/dy
!       contribution from viscous part on RHS of the NS equation
Implicit none
!
If (iuvwc ==  1) Then
Do j= jbeg,xsize(2);  Do i= ibeg,xsize(1)
uc(i,j)= (u1(i,j,k)*dz(k1d-1)+u1(i,j,k-1)*dz(k1d))/(dz(k1d-1)+dz(k1d))
uc2(i,j)= (u1(i,j,k+1)*dz(k1d)+u1(i,j,k)*dz(k1d+1))/(dz(k1d)+dz(k1d+1))
wc(i,j)= (w1(i,j,k)*dx(i-1)+w1(i-1,j,k)*dx(i))/(dx(i-1)+dx(i))
wc2(i,j)= (w1(i,j,k+1)*dx(i-1)+w1(i-1,j,k+1)*dx(i))/(dx(i-1)+dx(I))
uw(i,j)= uc(i,j)*wc(i,j)
uw2(i,j)= uc2(i,j)*wc2(i,j)
End Do;  End Do
!
!*********************************
Else If (iuvwc ==  2) Then
!*********************************
!
Do j= jbeg,xsize(2); j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)
vc(i,j)= (v1(i,j,k)*dz(k1d-1)+v1(i,j,k-1)*dz(k1d))/(dz(k1d-1)+dz(k1d))
vc2(i,j)= (v1(i,j,k+1)*dz(k1d)+v1(i,j,k)*dz(k1d+1))/(dz(k1d)+dz(k1d+1))
wc(i,j)= (w1(i,j,k)*dy(j1d-1)+w1(i,j-1,k)*dy(j1d))/(dy(j1d-1)+dy(j1d))
wc2(i,j)= (w1(i,j,k+1)*dy(j1d-1)+w1(i,j-1,k+1)*dy(j1d))/(dy(j1d-1)+dy(j1d))
vw(i,j)= vc(i,j)*wc(i,j)
vw2(i,j)= vc2(i,j)*wc2(i,j)
End Do;  End Do
!
!*********************************
Else If (iuvwc ==  3) Then
!*********************************
!
Do j= jbeg,xsize(2)+1; j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)
vc(i,j)= (v1(i,j,k-1)*dz(k1d)+v1(i,j,k)*dz(k1d-1))/(dz(k1d-1)+dz(k1d))
wc(i,j)= ((1.0e0-iuy_grdcon(j1d,2))*w1(i,j,k)*dy(j1d-1)+&
(1.0e0-iuy_grdcon(j1d,1))*w1(i,j-1,k)*dy(j1d))/(dy(j1d-1)+dy(j1d))
vw(i,j)= vc(i,j)*wc(i,j)
End Do; End Do
!
If (jconbc ==  1) Then
If(ibsoth ==  1.and.xmin(2)==1) Then
j= 1
 j1d=xstart(2)-1+j
Do i= ibeg,xsize(1)
wc(i,j)= (1.0e0-iuy_grdcon(j1d,1))*w1(i,0,k)
vw(i,j)= vc(i,j)*wc(i,j)
End Do
End If
!
If(ibnoth ==  1.and.xmax(2)==1) Then
j= xsize(2)+1
 j1d=xstart(2)-1+j
Do i= ibeg,xsize(1)
wc(i,j)= (1.0e0-iuy_grdcon(j1d,2))*w1(i,xsize(2)+1,k)
vw(i,j)= vc(i,j)*wc(i,j)
End Do
End If
End If
!
End If
!
Return
End Subroutine convective_jk
!
!******************************************************************************
Subroutine conv_cont(ii1,ii2,jj1,jj2)
!******************************************************************************
!	This module calculates the convective contribution on RHS of NS from previous time step
Implicit none
Integer, Intent(In) :: ii1,ii2,jj1,jj2
!
If (iuvwc ==  1) Then
!
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
hij1(i,j)= uu(i,j)-uu(i-1,j)
hij2(i,j)= (uv(i,j+1)-uv(i,j))*dxs(i)/dy(j1d)
hij3(i,j)= (uw2(i,j)-uw(i,j))*dxs(i)/dz(k1d)
End Do;  End Do
!
!********************************
Else If (iuvwc ==  2) Then
!********************************
!
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
hij1(i,j)= uv(i+1,j)-uv(i,j)
hij2(i,j)= (vv(i,j)-vv(i,j-1))*dx(i)/dys(j1d)
hij3(i,j)= (vw2(i,j)-vw(i,j))*dx(i)/dz(k1d)
End Do;  End Do
!
!********************************
Else If (iuvwc ==  3) Then
!********************************
!
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
hij1(i,j)= uw(i+1,j)-uw(i,j)
hij2(i,j)= (vw(i,j+1)-vw(i,j))*dx(i)/dy(j1d)
hij3(i,j)= (ww(i,j)-ww2(i,j))*dx(i)/dzs(k1d)
End Do;  End Do
!
End If
!
Return
End Subroutine conv_cont
!
!******************************************************************************
Subroutine intfo_calc(ii1,ii2,jj1,jj2)
!******************************************************************************
!	This module detemines the convective contribution to be used in the next time step
!	iuvwc is to choose between u1,v1 and w1 Momentum equation
!	iuvwc= 1 u1 momentum equation, iuvwc= 2 v1 Momentum equation,iuvwc= 3,w1 Momentum equation
Implicit none
Integer, Intent(In) :: ii1,ii2,jj1,jj2
!
If (iuvwc ==  1) Then
!
Do j= jj1,jj2;  Do i= ii1,ii2
intfuo(i,j,k)= hij1(i,j)+hij2(i,j)+hij3(i,j)
End Do;  End Do
!
!**********************************
Else If (iuvwc ==  2) Then
!**********************************
!
Do j= jj1,jj2;  Do i= ii1,ii2
intfvo(i,j,k)= hij1(i,j)+hij2(i,j)+hij3(i,j)
End Do;  End Do
!
!**********************************
Else If (iuvwc ==  3) Then
!**********************************
!
Do j= jj1,jj2;  Do i= ii1,ii2
intfwo(i,j,k)= hij1(i,j)+hij2(i,j)+hij3(i,j)
End Do;  End Do
!
End If
!
Return
End Subroutine intfo_calc
!
!******************************************************************************
Subroutine intf_calc(ii1,ii2,jj1,jj2)
!******************************************************************************
!	This module detemines the  total RHS of NS equation
Implicit none
Integer, Intent(In) :: ii1,ii2,jj1,jj2
!
If (iuvwc ==  1) Then
!
Do j= jj1,jj2;  Do i= ii1,ii2
  intf1(i,j,k)= intf1(i,j,k)+dt*uh1(i,j,k)-dt*intfuo(i,j,k)        &
  -dt*(p1(i,j,k)-p1(i-1,j,k)+dpdx_mean*dxs(i))/roh
  intf1(i,j,k)= intf1(i,j,k)/dxs(i)
End Do;  End Do
!
!***************************************
Else If (iuvwc ==  2) Then
!***************************************
!
Do j= jj1,jj2; j1d=xstart(2)-1+j;  Do i= ii1,ii2
  intf1(i,j,k)= intf1(i,j,k)+dt*uh1(i,j,k)-dt*intfvo(i,j,k)        &
 -dt*(p1(i,j,k)-p1(i,j-1,k))*dx(i)/dys(j1d)/roh
  intf1(i,j,k)= intf1(i,j,k)/dx(i)
End Do;  End Do
!
!***************************************
Else If (iuvwc ==  3) Then
!***************************************
!
Do j= jj1,jj2;  Do i= ii1,ii2
  intf1(i,j,k)= intf1(i,j,k)+dt*uh1(i,j,k)-dt*intfwo(i,j,k)        &
  -dt*(p1(i,j,k)-p1(i,j,k-1))*dx(i)/dzs(k1d)/roh
  intf1(i,j,k)= intf1(i,j,k)/dx(i)
End Do;  End Do
If(ibcjet==1.And.jet_type==2) Then
  Do j= jj1,jj2;  Do i= ii1,ii2
    intf1(i,j,k)= intf1(i,j,k)+fbody(i,j,k)*dt
  End Do; End Do
End If
!
End If
!
Return
End Subroutine intf_calc
!
!******************************************************************************
Subroutine mbc_dir_visc
!******************************************************************************
!       This subroutine adds the MBC matrix on the right hand side fo the Navier-stokes
!       equations for the dirichlet boundary condtions.
Implicit none
Real(mytype) :: dt2,nutc1,nutc2,nutxy1,nutxy2,nutxz1,nutxz2,nutyz1,nutyz2
!
!       Real variables asre defined by following topology
!       rnu : Viscous term, rco : convective term
!       rco(u/v/w :control volume)(x,y,z : gradient direction)(u/v/w : variable for k-1 step)
!       a(1/2/3 : Control volume)(xyx : gradient direction) (u/v/w : variable for k-1 step)*&
!       (u/v/w : variable for k time step)
!
iblock= 0;  ib= iblock+1
dt2= dt/2.0e0
!
If (icmpnt ==  1) Then
Do j= jbeg,xsize(2); j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)
!       Cofficients for d2u/dx2
nutc1= buffc(i-1)*brec(i-1)*(nu+2.0e0*pitero1(i-1,j,k))/dx(i-1)/dxs(i)
nutc2= buffc(i)*brec(i)*(nu+2.0e0*pitero1(i,j,k))/dx(i)/dxs(i)
!       Cofficients for d2u/dy2
nutxy1= bre(i)*(nu+rnutxy1(i,j,k))/dy(j1d)
nutxy2= bre(i)*(nu+rnutxy1(i,j+1,k))/dy(j1d)
!
intf1(i,j,k)= intf1(i,j,k)                                           &
!       Viscous Term Contribution
!       d2u/dx2 Contribution
 +iux_vis(i,1)*nutc1*(dt2*xinlet(j,k,1,1,2)+dt2*xinlet(j,k,1,1,1))&
 +iux_vis(i,2)*nutc2*(dt2*xinlet(j,k,1,2,2)+dt2*xinlet(j,k,1,2,1))&
!       d2u/dy2 Contribution
 +iuy_vis(j1d,1,1)*(nutxy2*gcdy_f(j1d+1,4,1,ib)-nutxy1*gcdy_f(j1d,3,1,ib))&
  *(dt2*yinlet1(i,k,1,1,2)+dt2*yinlet1(i,k,1,1,1))           &
 -iuy_vis(j1d,2,1)*nutxy1*gcdy_f(j1d,4,1,ib)*(dt2*yinlet1(i,k,1,1,2)+    &
  dt2*yinlet1(i,k,1,1,1))                                  &
 +iuy_vis(j1d,2,2)*nutxy2*gcdy_f(j1d+1,1,1,ib)*(dt2*yinlet1(i,k,1,2,2)+  &
  dt2*yinlet1(i,k,1,2,1))&
 +iuy_vis(j1d,1,2)*(nutxy2*gcdy_f(j1d+1,2,1,ib)-nutxy1*gcdy_f(j1d,1, &
  1,ib))*(dt2*yinlet1(i,k,1,2,2)+dt2*yinlet1(i,k,1,2,1))
End Do;  End Do
!
!
!*****************************************
        Else If (icmpnt ==  2) Then
!*****************************************
!
Do j= jbeg,xsize(2); j1d=xstart(2)-1+j;  Do i= ibeg,xsize(1)
!       Cofficients for d2v/dx2
nutxy1= buff(i)*bre(i)*(nu+rnutxy1(i,j,k))/dx(i)
nutxy2= buff(i+1)*bre(i+1)*(nu+rnutxy1(i+1,j,k))/dx(i)
!       Cofficients for d2v/dy2
nutc1= brec(i)*(nu+2.0e0*pitero1(i,j-1,k))/dy(j1d-1)/dys(j1d)
nutc2= brec(i)*(nu+2.0e0*pitero1(i,j,k))/dy(j1d)/dys(j1d)
!
intf1(i,j,k)= intf1(i,j,k)                                          &
!       Visocous Term contribution
!       d2v/dx2 Contribution
 -ivx_vis(i,1,1)*(nutxy1*gcdx_f(i,3,1,ib)-nutxy2*gcdx_f(i+1,  &
  4,1,ib))*(dt2*xinlet(j,k,2,1,2)+dt2*xinlet(j,k,2,1,1))       &
 -ivx_vis(i,2,1)*nutxy1*gcdx_f(i,4,1,ib)*(dt2*xinlet(j,k,2,1,2)   &
 +dt2*xinlet(j,k,2,1,1))                                      &
 +ivx_vis(i,2,2)*nutxy2*gcdx_f(i+1,1,1,ib)*(dt2*xinlet(j,k,2,2,2) &
 +dt2*xinlet(j,k,2,2,1))                                      &
 -ivx_vis(i,1,2)*(nutxy1*gcdx_f(i,2,1,ib)-nutxy2*gcdx_f(i+1,  &
  3,1,ib))*(dt2*xinlet(j,k,2,2,2)+dt2*xinlet(j,k,2,2,1))       &
!       d2v/dy2 Contribution
 +ivy_vis(j1d,1)*nutc1*(dt2*yinlet1(i,k,2,1,2)+dt2*yinlet1(i,k,2,1,1))  &
 +ivy_vis(j1d,2)*nutc2*(dt2*yinlet1(i,k,2,2,2)+dt2*yinlet1(i,k,2,2,1))
End Do;  End Do
!
!*************************************
        Else If (icmpnt ==  3) Then
!*************************************
!
Do j= ibeg,xsize(2); j1d=xstart(2)-1+j;  Do i= jbeg,xsize(1)
!       Cofficients for d2w/dx2
  nutxz1= buff(i)*bre(i)*(nu+rnutxz1(i,j,k))/dx(i)
  nutxz2= buff(i+1)*bre(i+1)*(nu+rnutxz1(i+1,j,k))/dx(i)
!       Cofficients for d2w/dy2
  nutyz1= brec(i)*(nu+rnutyz1(i,j,k))/dy(j1d)
  nutyz2= brec(i)*(nu+rnutyz1(i,j+1,k))/dy(j1d)
!
          intf1(i,j,k)= intf1(i,j,k)                                          &
!       Visocous Term contribution
!       d2w/dx2 Contribution
  -ivx_vis(i,1,1)*(nutxz1*gcdx_f(i,3,1,ib)-nutxz2*gcdx_f(i+1, &
   4,1,ib))*(dt2*xinlet(j,k,3,1,2)+dt2*xinlet(j,k,3,1,1))       &
  -ivx_vis(i,2,1)*nutxz1*gcdx_f(i,4,1,ib)*(dt2*xinlet(j,k,3,1,2)  &
  +dt2*xinlet(j,k,3,1,1))                                       &
  +ivx_vis(i,2,2)*nutxz2*gcdx_f(i+1,1,1,ib)*(dt2*xinlet(j,k,3,2,2)  &
  +dt2*xinlet(j,k,3,2,1)) &
  -ivx_vis(i,1,2)*(nutxz1*gcdx_f(i,1,1,ib)-nutxz2*gcdx_f(i+1, &
   2,1,ib))*(dt2*xinlet(j,k,3,2,2)+dt2*xinlet(j,k,3,2,1))      &
!       d2w/dy2 Contribution
 -iuy_vis(j1d,1,1)*(nutyz1*gcdy_f(j1d,3,1,ib)-nutyz2*gcdy_f(j1d+1,  &
  4,1,ib))*(dt2*yinlet1(i,k,3,1,2)+dt2*yinlet1(i,k,3,1,1))       &
 -iuy_vis(j1d,2,1)*nutyz1*gcdy_f(j1d,4,1,ib)*(dt2*yinlet1(i,k,3,1,2)     &
 +dt2*yinlet1(i,k,3,1,1))                                      &
 +iuy_vis(j1d,2,2)*nutyz2*gcdy_f(j1d+1,1,1,ib)*(dt2*yinlet1(i,k,3,2,2)   &
 +dt2*yinlet1(i,k,3,2,1)) &
 -iuy_vis(j1d,1,2)*(nutyz1*gcdy_f(j1d,1,1,ib)-nutyz2*gcdy_f(j1d+1,  &
  2,1,ib))*(dt2*yinlet1(i,k,3,2,2)+dt2*yinlet1(i,k,3,2,1))
End Do;  End Do
End If
!
Return
End Subroutine mbc_dir_visc
!
!******************************************************************************
Subroutine mbc_dir_conv
!******************************************************************************
!       This subroutine adds the MBC matrix on the right hand side fo the Navier-stokes
!       equations for the dirichlet boundary condtions.
Implicit none
Integer :: i,j,iblock,ib
Real(mytype) :: dt2,rcouxu1,rcouxu2,a1yvu,b1yvu,c1yvu,rcouyv1,rcouyv2,rcouyu1, &
rcouyu2,rcovxu1,rcovxu2,a2xuv,b2xuv,c2xuv,rcovyv1,rcovyv2,rcovxv1,    &
rcovxv2,rcowxu1,rcowxu2,a3xuw,b3xuw,c3xuw,rcowyv1,rcowyv2,a3yvw,b3yvw,&
c3yvw,rcowxw1,rcowxw2,rcowyw1,rcowyw2
!
!       Real variables asre defined by following topology
!       rnu : Viscous term, rco : convective term
!       rco(u/v/w :control volume)(x,y,z : gradient direction)(u/v/w : variable for k-1 step)
!       a(1/2/3 : Control volume)(xyx : gradient direction) (u/v/w : variable for k-1 step)*&
!       (u/v/w : variable for k time step)
!
iblock= 0;  ib= iblock+1
dt2= dt/2.0e0
!
If (icmpnt ==  1) Then
Do j= jbeg,xsize(2); j1d=xstart(2)-1+j; Do i= ibeg,xsize(1)
!       Cofficient for duu/dx
rcouxu1= dt2/2.0e0/dxs(i)*(u1(i-1,j,k)+u1(i,j,k))
rcouxu2= dt2/2.0e0/dxs(i)*(u1(i+1,j,k)+u1(i,j,k))
!       Cofficients for dvu/dy
rcouyv1= dt2/dy(j1d)*(v1(i,j,k)*dx(i-1) + v1(i-1,j,k)*dx(i))/(dx(i)+dx(i-1))
rcouyv2= dt2/dy(j1d)*(v1(i,j+1,k)*dx(i-1) + v1(i-1,j+1,k)*dx(i))/(dx(i)+dx(i-1))
a1yvu=  -rcouyv1*dy(j1d)/(dy(j1d)+dy(j1d-1))
b1yvu=   rcouyv2*dy(j1d+1)/(dy(j1d)+dy(j1d+1))-rcouyv1*dy(j1d-1)/(dy(j1d)+dy(j1d-1))
c1yvu=  rcouyv2*dy(j1d)/(dy(j1d)+dy(j1d+1))
If (j1d ==  1.and.ibsoth ==  1) Then
a1yvu=  -rcouyv1
b1yvu=   rcouyv2*dy(j1d+1)/(dy(j1d)+dy(j1d+1))
End If
If (j1d ==  ny.and.ibnoth ==  1) Then
b1yvu=   -rcouyv1*dy(j1d-1) /(dy(j1d)+dy(j1d-1))
c1yvu=  rcouyv2
End If
!       Cofficient for d(uv)/dy
rcouyu1=  dt2/dy(j1d)/(dx(i)+dx(i-1))/(dy(j1d)+dy(j1d-1))*(u1(i,j,k)*dy(j1d-1)+u1(i,j-1,k)*dy(j1d))
rcouyu2=  dt2/dy(j1d)/(dx(i)+dx(i-1))/(dy(j1d)+dy(j1d+1))*(u1(i,j+1,k)*dy(j1d)+u1(i,j,k)*dy(j1d+1))
If (j1d ==  1.and.ibsoth ==  1)&
rcouyu1=  dt2/dy(j1d)/(dx(i)+dx(i-1))* u1(i,0,k)
If (j1d ==  ny.and.ibnoth ==  1)&
rcouyu2=  dt2/dy(j1d)/(dx(i)+dx(i-1))* u1(i,xsize(2)+1,k)
!
intf1(i,j,k)= intf1(i,j,k)                                           &
!       d(uu)/dx Contribution
+iux_con(i,1)*rcouxu1*xinlet(j,k,1,1,2)                            &
-iux_con(i,2)*rcouxu2*xinlet(j,k,1,2,2)                            &
!       d(vu)/dy Contribution
-iuy_con(j1d,1)*a1yvu*yinlet1(i,k,1,1,2)                              &
-iuy_con(j1d,2)*c1yvu*yinlet1(i,k,1,2,2)                              &
!       d(uv)/dy  Contribution
+iuy_con(j1d,1)*rcouyu1*(yinlet1(i,k,2,1,2)*dx(i-1)+                  &
 yinlet1(i-1,k,2,1,2)*dx(i))                                        &
-iuy_con(j1d,2)*rcouyu2*(yinlet1(i,k,2,2,2)*dx(i-1)+                  &
 yinlet1(i-1,k,2,2,2)*dx(i))
End Do;  End Do
!
!*****************************************
        Else If (icmpnt ==  2) Then
!*****************************************
!
Do j= jbeg,xsize(2);j1d=xstart(2)-1+j; Do i= ibeg,xsize(1)
!       Cofficients for d(uv/dx)
rcovxu1= dt2/dx(i)*(u1(i,j,k)*dy(j1d-1)+u1(i,j-1,k)*dy(j1d))/(dy(j1d)+dy(j1d-1))
rcovxu2= dt2/dx(i)*(u1(i+1,j,k)*dy(j1d-1)+u1(i+1,j-1,k)*dy(j1d))/(dy(j1d)+dy(j1d-1))
a2xuv=  -rcovxu1*dx(i)/(dx(i)+dx(i-1))
b2xuv=   rcovxu2*dx(i+1)/(dx(i)+dx(i+1)) -rcovxu1*dx(i-1)/(dx(i)+dx(i-1))
c2xuv=  rcovxu2*dx(i)/(dx(i)+dx(i+1))
If (i ==  1.and.ibwest ==  1) Then
a2xuv=  -rcovxu1
b2xuv=   rcovxu2*dx(i+1)/(dx(i)+dx(i+1))
End If
If (i ==  nx.and.ibeast ==  1) Then
b2xuv=  -rcovxu1*dx(i-1)/(dx(i)+dx(i-1))
c2xuv=   rcovxu2
End If
!       Cofficients for d(vv)/dy
rcovyv1= dt2/2.0e0/dys(j1d)*(v1(i,j,k)+v1(i,j-1,k))
rcovyv2= dt2/2.0e0/dys(j1d)*(v1(i,j,k)+v1(i,j+1,k))
!       Cofficients for d(vu)/dx
rcovxv1= dt2/dx(i)/(dy(j1d)+dy(j1d-1))/(dx(i)+dx(i-1))*(v1(i,j,k)*dx(i-1)+v1(i-1,j,k)*dx(i))
rcovxv2= dt2/dx(i)/(dy(j1d)+dy(j1d-1))/(dx(i)+dx(i+1))*(v1(i+1,j,k)*dx(i)+v1(i,j,k)*dx(i+1))
If (i ==  1.and.ibwest ==  1)&
rcovxv1= dt2/dx(i)/(dy(j1d)+dy(j1d-1))*v1(0,j,k)
If (i ==  nx.and.ibeast ==  1)&
rcovxv2= dt2/dx(i)/(dy(j1d)+dy(j1d-1))*v1(nxt,j,k)
!
intf1(i,j,k)= intf1(i,j,k)                                          &
!       d(uv)/dx Contribution
-ivx_con(i,1)*a2xuv*xinlet(j,k,2,1,2)                             &
-ivx_con(i,2)*c2xuv*xinlet(j,k,2,2,2)                             &
!       Contribution from d(vv)/dy
+ivy_con(j1d,1)*rcovyv1*yinlet1(i,k,2,1,2)                           &
-ivy_con(j1d,2)*rcovyv2*yinlet1(i,k,2,2,2)                           &
!       Contribution from d(vu)/dy
+ivx_con(i,1)*rcovxv1*(xinlet(j,k,1,1,2)*dy(j1d-1)+                 &
 xinlet(j-1,k,1,1,2)*dy(j1d))                                       &
-ivx_con(i,2)*rcovxv2*(xinlet(j,k,1,2,2)*dy(j1d-1)+                 &
 xinlet(j-1,k,1,2,2)*dy(j1d))
End Do;  End Do
!*************************************
        Else If (icmpnt ==  3) Then
!*************************************
!
        Do j= ibeg,xsize(2);j1d=xstart(2)-1+j; Do i= jbeg,xsize(1)
!       Cofficients for d(uw)/dx 
rcowxu1= dt2/dx(i)/(dz(k1d)+dz(k1d-1))*(u1(i,j,k)*dz(k1d-1)+u1(i,j,k-1)*dz(k1d))
rcowxu2= dt2/dx(i)/(dz(k1d)+dz(k1d-1))*(u1(i+1,j,k)*dz(k1d-1)+u1(i+1,j,k-1)*dz(k1d))
a3xuw=  -rcowxu1*dx(i)/(dx(i)+dx(i-1))
b3xuw=   rcowxu2*dx(i+1)/(dx(i)+dx(i+1))-rcowxu1*dx(i-1)/(dx(i)+dx(i-1))
c3xuw=   rcowxu2*dx(i)/(dx(i)+dx(i+1))
If (i ==  1.and.ibwest ==  1) Then
a3xuw=  -rcowxu1
b3xuw=   rcowxu2*dx(i+1)/(dx(i)+dx(i+1))
End If
If (i ==  nx.and.ibeast ==  1) Then
b3xuw=  -rcowxu1*dx(i-1)/(dx(i)+dx(i-1))
c3xuw=   rcowxu2
End If
!       Cofficients for d(vw)/dy
rcowyv1= dt2/dy(j1d)/(dz(k1d)+dz(k1d-1))*(v1(i,j,k)*dz(k1d-1)+v1(i,j,k-1)*dz(k1d))
rcowyv2= dt2/dy(j1d)/(dz(k1d)+dz(k1d-1))*(v1(i,j+1,k)*dz(k1d-1)+v1(i,j+1,k-1)*dz(k1d))
a3yvw=  -rcowyv1*dy(j1d)/(dy(j1d)+dy(j1d-1))
b3yvw=   rcowyv2*dy(j1d+1)/(dy(j1d)+dy(j1d+1))-rcowyv1*dy(j1d-1)/(dy(j1d)+dy(j1d-1))
c3yvw=   rcowyv2*dy(j1d)/(dy(j1d)+dy(j1d+1))
If (j1d ==  1.and.ibsoth ==  1) Then
a3yvw=  -rcowyv1
b3yvw=   rcowyv2*dy(j1d+1)/(dy(j1d)+dy(j1d+1))
End If
If (j1d ==  ny.and.ibnoth ==  1) Then
b3yvw=   -rcowyv1*dy(j1d-1)/(dy(j1d)+dy(j1d-1))
c3yvw=  rcowyv2
End If
!       Cofficients for d(wu)/dx
rcowxw1= dt2/dx(i)/(dz(k1d)+dz(k1d-1))/(dx(i)+dx(i-1))                  &
        *(w1(i,j,k)*dx(i-1)+w1(i-1,j,k)*dx(i))
rcowxw2= dt2/dx(i)/(dz(k1d)+dz(k1d-1))/(dx(i)+dx(i+1))                  &
        *(w1(i+1,j,k)*dx(i)+w1(i,j,k)*dx(i+1))
If (i ==  1.and.ibwest ==  1)&
rcowxw1= dt2/dx(i)/(dz(k1d)+dz(k1d-1))*w1(0,j,k)
If (i ==  nx.and.ibeast ==  1)&
rcowxw2= dt2/dx(i)/(dz(k1d)+dz(k1d-1))*w1(nxt,j,k)
!       Cofficients for d(wv)/dx
rcowyw1=  dt2/dy(j1d)/(dz(k1d)+dz(k1d-1))/(dy(j1d)+dy(j1d-1))*               &
         (w1(i,j,k)*dy(j1d-1)+w1(i,j-1,k)*dy(j1d))
rcowyw2=  dt2/dy(j1d)/(dz(k1d)+dz(k1d-1))/(dy(j1d)+dy(j1d+1))*               &
         (w1(i,j+1,k)*dy(j1d)+w1(i,j,k)*dy(j1d+1))
If (j1d ==  1.and.ibsoth ==  1)&
rcowyw1=  dt2/dy(j1d)/(dz(k1d)+dz(k1d-1))*w1(i,0,k)
If (j1d ==  ny.and.ibnoth ==  1)&
rcowyw2=  dt2/dy(j1d)/(dz(k1d)+dz(k1d-1))*w1(i,xsize(2)+1,k)
!
          intf1(i,j,k)= intf1(i,j,k)                                          &
!       d(uw)/dx Contribution
-ivx_con(i,1)*a3xuw*xinlet(j,k,3,1,2)                             &
-ivx_con(i,2)*c3xuw*xinlet(j,k,3,2,2)                             &
!       d(vw)/dy Contribution
!-iuy_con(j1d,1)*a3yvw*yinlet1(i,k,3,1,2)                             &
!-iuy_con(j1d,2)*c3yvw*yinlet1(i,k,3,2,2)                             &
!       d(wu)/dx Contribution
+ivx_con(i,1)*rcowxw1*(xinlet(j,k,1,1,2)*dz(k1d-1)+                 &
 xinlet(j,k-1,1,1,2)*dz(k1d))                                       &
-ivx_con(i,2)*rcowxw2*(xinlet(j,k,1,2,2)*dz(k1d-1)+                 &
 xinlet(j,k-1,1,2,2)*dz(k1d))                                       &
!       d(wv)/dy Contribution
+iuy_con(j1d,1)*rcowyw1*(yinlet1(i,k,2,1,2)*dz(k1d-1)+                 &
 yinlet1(i,k-1,2,1,2)*dz(k1d))                                       &
-iuy_con(j1d,2)*rcowyw2*(yinlet1(i,k,2,2,2)*dz(k1d-1)+                 &
 yinlet1(i,k-1,2,2,2)*dz(k1d))
End Do;  End Do
End If
!
Return
        End Subroutine mbc_dir_conv

!******************************************************************************
!
End MODULE compute_right_implicit

