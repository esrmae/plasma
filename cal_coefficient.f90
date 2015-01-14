Module cal_coefficient
Use shared_data
Use boundary_conditions
Implicit None
Integer , Private :: isec
Contains
!
!******************************************************************************
Subroutine compute_fact
!******************************************************************************
!      This procedure computes geometrical factors to calculate gradients
!      using 3rd order interpolation.
Implicit None
!
Call gfds_c_comp
Call gcds_f_comp
Call gcds_c_comp
!Call coff_calc_out
!
Return
End Subroutine compute_fact
!
!******************************************************************************
Subroutine gfds_c_comp
!******************************************************************************
!	This procedure calculates gradient at cel center with variable laocated
!	at cell face using two points
Implicit None
Integer :: ib,iblock
!
iblock = 0
ib=iblock+1
!
!	Gradient at Cell Center When Variable is at Cell Face (Two-Point)
!
Call scheme1(gfdx_c,dx,0,xsize(1),1,ib)
Call scheme1(gfdy_c,dy,1,ysize(2),1,ib)
Call scheme1(gfdz_c,dz,1,zsize(3),1,ib)
Call scheme1(gfdz_c2,dz,1,zsize(3),0,ib)
!
Return
End Subroutine gfds_c_comp 
!
!******************************************************************************
Subroutine gcds_f_comp
!*****************************************************************************
!	This procedure calculates the gradient at cell face with variable located
!	at cell center using four points
Implicit None
Integer :: i,ib,iblock,j,m,k,low,up
!
iblock = 0
ib=iblock+1
!
!	Gradients at Cell Face when Variable is at Cell Center(Four_Point)
!
!**************** must be in X pencil *******************
If (igr_2pt.eq.0) Then
!
If (iconbc == 0) then
  Call scheme2(dx,1,xsize(1),1,1,ib,xsize(1))
Else 
  Call scheme2(dx,3,xsize(1)-2,1,1,ib,xsize(1))
End If
!
If (jconbc == 0) then
  Call scheme2(dy,1,ysize(2),2,1,ib,ysize(2))
Else 
  Call scheme2(dy,3,ysize(2)-2,2,1,ib,ysize(2))
End If
!
Else If (igr_2pt.eq.1) Then
Do i = 1,xsize(1)+1
gcdx_f(i,1,1,ib) = 0.0e0
gcdx_f(i,2,1,ib) = 1.0/dxs(i)
gcdx_f(i,3,1,ib) = -1.0/dxs(i)
gcdx_f(i,4,1,ib) = 0.0e0
End Do
!**************** must be in Y pencil *******************
Do i = 1,ysize(2)+1
gcdy_f(i,1,1,ib) = 0.0e0
gcdy_f(i,2,1,ib) = 1.0/dys(i)
gcdy_f(i,3,1,ib) = -1.0/dys(i)
gcdy_f(i,4,1,ib) = 0.0e0
End Do
End If
!
!**************** must be in Z pencil *******************
If (ifgrz_2pt.eq.0) Then
!
Call scheme2(dz,1,zsize(3),3,1,ib,zsize(3))
!
gcdz_f2(:,:,1,ib)=gcdz_f(:,:,1,ib)
Else If (ifgrz_2pt == 1) Then
Do i = 1,zsize(3)+1
gcdz_f(i,1,1,ib) = 0.0e0
gcdz_f(i,2,1,ib) = 1.0/dzs(i)
gcdz_f(i,3,1,ib) = -1.0/dzs(i)
gcdz_f(i,4,1,ib) = 0.0e0
End Do
Do i = 1,zsize(3)
gcdz_f2(i,1,1,ib) = 0.0e0
gcdz_f2(i,2,1,ib) = 1.0/dzs(i+1)
gcdz_f2(i,3,1,ib) = -1.0/dzs(i+1)
gcdz_f2(i,4,1,ib) = 0.0e0
End Do
End If
!
!	Pressure Coeffiients
!
!**************** must be in X pencil *******************
Do i = 1,xsize(1)+1
Do m = 1,4
  gcdx_f(i,m,2,ib) = gcdx_f(i,m,1,ib)
End Do
End Do
!**************** must be in Y pencil *******************
Do i = 1,ysize(2)+1
Do m = 1,4
 gcdy_f(i,m,2,ib) = gcdy_f(i,m,1,ib)
End Do
End Do
!**************** must be in Z pencil *******************
Do i = 1,zsize(3)+1
Do m = 1,4
 gcdz_f(i,m,2,ib) = gcdz_f(i,m,1,ib)
 gcdz_f2(i,m,2,ib) = gcdz_f2(i,m,1,ib)
End Do
End Do
!
!	Boundary Gradeints
!
If (igr_2pt.eq.0) Then
!
If (iconbc == 1)Call scheme2_bc(dx,1,xsize(1),ib,1)
If (jconbc /= 0)Call scheme2_bc(dy,1,ysize(2),ib,2)
!
Else If (igr_2pt.eq.1) Then
!**************** must be in X pencil *******************
If (iconbc == 1) then
i=1
gcdx_f(i,1,1,ib) = 0.0e0
gcdx_f(i,2,1,ib) = 2.0/dx(i)
gcdx_f(i,3,1,ib) = -2.0/dx(i)
gcdx_f(i,4,1,ib) = 0.0e0
gcdx_f(i,1:4,2,ib) = 0.0e0
i=xsize(1)+1
gcdx_f(i,1,1,ib) = 0.0e0
gcdx_f(i,2,1,ib) = 2.0/dx(xsize(1))
gcdx_f(i,3,1,ib) = -2.0/dx(xsize(1))
gcdx_f(i,4,1,ib) = 0.0e0
gcdx_f(i,1:4,2,ib) = 0.0e0
End If
!**************** must be in Y pencil *******************
If (jconbc /= 0) Then
j=1
gcdy_f(j,1,1,ib) = 0.0e0
gcdy_f(j,2,1,ib) = 2.0/dy(j)
gcdy_f(j,3,1,ib) = -2.0/dy(j)
gcdy_f(j,4,1,ib) = 0.0e0
gcdy_f(j,1:4,2,ib) = 0.0e0
j=ysize(2)+1
gcdy_f(j,1,1,ib) = 0.0e0
gcdy_f(j,2,1,ib) = 2.0/dy(ysize(2))
gcdy_f(j,3,1,ib) = -2.0/dy(ysize(2))
gcdy_f(j,4,1,ib) = 0.0e0
gcdy_f(j,1:4,2,ib) = 0.0e0
End If
End If
!
Return
End Subroutine gcds_f_comp 
!
!******************************************************************************
Subroutine gcds_c_comp
!*****************************************************************************
!	This subroutine calculates the gradient at cell center with variable located
!	cell center using five point formulation.
Implicit None
Integer :: i,ib,iblock
!
iblock = 0
ib=iblock+1
!
!	Gradients at Cell Center when Variable is at Cell Center(Five_Point)
Call scheme3(gcdx_c,dx,1,xsize(1),1,ib)
Call scheme3(gcdy_c,dy,1,ysize(2),2,ib)
Call scheme3(gcdz_c,dz,1,zsize(3),3,ib)
!Call scheme3_LES(grdsxx,dx,1,xsize(1)/2,1,ib)
!Call scheme3_les(grdsyy,dy,1,ysize(2)/2,1,ib)
!Call scheme3_les(grdszz,dz,1,zsize(3)/2,2,ib)
If (iconbc == 1)Call scheme3_bc(gcdx_c,grdsxx,dx,1,xsize(1),1,ib) 
If (jconbc /= 0)Call scheme3_bc(gcdy_c,grdsyy,dy,1,ysize(2),1,ib) 
!If (kconbc.ne.0) Call scheme3_bc(gcdz_c,grdszz,dz,1,zsize(3),2,ib) 
!
Return
End Subroutine gcds_c_comp 
!
!******************************************************************************
Subroutine third(ij,a,b,c,d,iuv,ik,ib)
!*****************************************************************************
!------ This subroutine computes gradient coefficients using 3rd order
Implicit None
Integer,Intent( In ) :: ij,iuv,ik,ib
Real(mytype),Intent( In ) :: a,b,c,d
Real(mytype), Dimension(4,4) :: mat
Real(mytype), Dimension(4) :: vec
Integer :: i,j
Real(mytype) :: div
!
mat(1,:)=1.0e0
Do i=2,4
  mat(i,1)=a**(i-1)
  mat(i,2)=b**(i-1)
  mat(i,3)=c**(i-1)
  mat(i,4)=d**(i-1)
End Do
vec(:)=0.0e0
vec(2)=1.0e0
Call gauss(mat,vec)
!
If (iuv == 1) Then
  gcdx_f(ij,1,ik,ib) = vec(1)
  gcdx_f(ij,2,ik,ib) = vec(2)
  gcdx_f(ij,3,ik,ib) = vec(3)
  gcdx_f(ij,4,ik,ib) = vec(4)
Else If (iuv == 2) Then
  gcdy_f(ij,1,ik,ib) = vec(1)
  gcdy_f(ij,2,ik,ib) = vec(2)
  gcdy_f(ij,3,ik,ib) = vec(3)
  gcdy_f(ij,4,ik,ib) = vec(4)
Else If (iuv == 3) Then
  gcdz_f(ij,1,ik,ib) = vec(1)
  gcdz_f(ij,2,ik,ib) = vec(2)
  gcdz_f(ij,3,ik,ib) = vec(3)
  gcdz_f(ij,4,ik,ib) = vec(4)
End If
!
Return
End Subroutine third 
!
!******************************************************************************
Subroutine gauss(mat,vec)
!*****************************************************************************
!------ This subroutine inverts 4x4 matrix: 
! On entry - mat = matrix, vec = rhs
! On exit  - vec = lhs = (mat^-1)*rhs
Implicit None
Real(mytype), Dimension(4,4) :: mat
Real(mytype), Dimension(4) :: vec
Real(mytype) :: dum
Integer :: i,j
!
Do i=2,4
dum=mat(i,1)/mat(1,1)
Do j=1,4
  mat(i,j)=mat(i,j)-mat(1,j)*dum
End Do
vec(i)=vec(i)-vec(1)*dum
End Do
!
Do i=3,4
dum=mat(i,2)/mat(2,2)
Do j=2,4
  mat(i,j)=mat(i,j)-mat(2,j)*dum
End Do
vec(i)=vec(i)-vec(2)*dum
End Do 
!
dum=mat(4,3)/mat(3,3)
Do j=3,4
  mat(4,j)=mat(4,j)-mat(3,j)*dum
End Do
vec(4)=vec(4)-vec(3)*dum
!
vec(4)=vec(4)/mat(4,4)
vec(3)=(vec(3)-mat(3,4)*vec(4))/mat(3,3)
vec(2)=(vec(2)-mat(2,4)*vec(4)-mat(2,3)*vec(3))/mat(2,2)
vec(1)=(vec(1)-mat(1,4)*vec(4)-mat(1,3)*vec(3)-mat(1,2)*vec(2))/mat(1,1)
!
Return
End Subroutine gauss 
!
!******************************************************************************
	Subroutine third3 (ij,a,b,c,d,iuv,nabv,ik,ib)
!*****************************************************************************
!       This subroutine computes gradient coefficients using 3rd order
!       NABV is the number of points above the face under consideration.
Implicit None
Integer, Intent( In ) :: ij,nabv,iuv,ik,ib
Real(mytype), Intent( In ) :: a,b,c,d
Real(mytype) :: deno
!
If (nabv == 1) Then
  deno=a*(b+c)*(b-c)+b*(c+a)*(c-a)+c*(a+b)*(a-b)
  If (iuv == 1) Then
    gcdx_f(ij,1,ik,ib) = (b+c)*(b-c)/deno
    gcdx_f(ij,2,ik,ib) = (c+a)*(c-a)/deno
    gcdx_f(ij,3,ik,ib) = (a+b)*(a-b)/deno
    gcdx_f(ij,4,ik,ib) = 0.0e0
  Else If (iuv == 2) Then
    gcdy_f(ij,1,ik,ib) = (b+c)*(b-c)/deno
    gcdy_f(ij,2,ik,ib) = (c+a)*(c-a)/deno
    gcdy_f(ij,3,ik,ib) = (a+b)*(a-b)/deno
    gcdy_f(ij,4,ik,ib) = 0.0e0
  Else If (iuv == 3) Then
    gcdz_f(ij,1,ik,ib) = (b+c)*(b-c)/deno
    gcdz_f(ij,2,ik,ib) = (c+a)*(c-a)/deno
    gcdz_f(ij,3,ik,ib) = (a+b)*(a-b)/deno
    gcdz_f(ij,4,ik,ib) = 0.0e0
  End If
Else If (nabv == 2) Then
  deno=b*(c+d)*(c-d)+c*(d+b)*(d-b)+d*(b+c)*(b-c)
  If (iuv == 1) Then
    gcdx_f(ij,1,ik,ib) = 0.0e0
    gcdx_f(ij,2,ik,ib) = (c+d)*(c-d)/deno
    gcdx_f(ij,3,ik,ib) = (d+b)*(d-b)/deno
    gcdx_f(ij,4,ik,ib) = (b+c)*(b-c)/deno
  Else If (iuv == 2) Then
    gcdy_f(ij,1,ik,ib) = 0.0e0
    gcdy_f(ij,2,ik,ib) = (c+d)*(c-d)/deno
    gcdy_f(ij,3,ik,ib) = (d+b)*(d-b)/deno
    gcdy_f(ij,4,ik,ib) = (b+c)*(b-c)/deno
  Else If (iuv == 3) Then
    gcdz_f(ij,1,ik,ib) = 0.0e0
    gcdz_f(ij,2,ik,ib) = (c+d)*(c-d)/deno
    gcdz_f(ij,3,ik,ib) = (d+b)*(d-b)/deno
    gcdz_f(ij,4,ik,ib) = (b+c)*(b-c)/deno
  End If
End If
!
Return
End Subroutine third3
!
!******************************************************************************
Subroutine comp_pres_coff
!*****************************************************************************
!       This procedure computes coefficients a's before solving Poisson eq.'
!       for pressure.
Implicit None
Integer :: i,j,k,j1d,k1d
!
!**************** must be in X pencil *******************
Do k = 1,xsize(3)+1; Do j = 1,xsize(2)+1; Do i = 1,xsize(1)+1
  j1d = xstart(2)-1+j
  k1d = xstart(3)-1+k
  ae(i,j,k) = dy(j1d)*dz(k1d)/dxs(i)
  as(i,j,k) = dx(i)*dz(k1d)/dys(j1d)
  at(i,j,k) = dx(i)*dy(j1d)/dzs(k1d)
  ap(i,j,k) = 1.0e0
End Do; End Do; End Do

If (if2d_noy == 1) as = 0.0e0
!
  Call update_halo(ae,2)
  Call update_halo(as,2)
  Call update_halo(at,2)

If (xmin(2) == 1) as(:,1,:) = 0.0e0
If (xmax(2) == 1) then
  Do k = 1,xsize(3)+1; Do i = 1,xsize(1)+1
    ae(i,xsize(2)+1,k) = 0.0e0
    as(i,xsize(2)+1,k) = 0.0e0
    at(i,xsize(2)+1,k) = 0.0e0
  End Do; End Do
End If
!
!       Periodic B.C.
!       AE, AS, AT are O.K at the inlet plane.
!
If (iconbc == 0) Then
   ae(xsize(1)+1,:,:) = ae(1,:,:)
!
!       Convective B.C.
!
Else If (iconbc == 1) Then
   ae(1,:,:) = 0.0e0
   ae(xsize(1)+1,:,:) = 0.0e0
End If
!
!       For cutting the links
!
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
   ap(i,j,k) = ae(i,j,k)+ae(i+1,j,k)+as(i,j,k)+as(i,j+1,k)+at(i,j,k)+at(i,j,k+1)
   ap(i,j,k) = 1.0e0/ap(i,j,k)
End Do; End Do; End Do
!
!       Periodic B.C.
!
If (iconbc == 0) ap(xsize(1)+1,:,:) = ap(1,:,:)
!
Return
End Subroutine comp_pres_coff
!
!******************************************************************************
Subroutine scheme1(temp,ds,mm,nn,jj,ib)
!******************************************************************************
!	This subroutine calculates the gradient at cell center using two points 
!	when the varables are located at cell face.
!       jj is an argument which is needed to calculate gradient at next cell or
!	the previous one
Implicit None
Integer, Intent( In ) :: nn,mm,jj,ib
Real(mytype), Dimension(mm:nn+1,4,ngr,nbl), Intent( InOut ) :: temp
Real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
Integer :: i,j
!
If (jj == 1) Then
  Do i = mm,nn
    temp(i,1,1,ib) = 0.0e0
    temp(i,2,1,ib) = 1.0e0/ds(i)
    temp(i,3,1,ib) = -1.0e0/ds(i)
    temp(i,4,1,ib) = 0.0e0
  End Do
Else If (jj == 0) Then
  Do i = mm,nn
    temp(i,1,1,ib) = 0.0e0
    temp(i,2,1,ib) = 1.0e0/ds(i-1)
    temp(i,3,1,ib) = -1.0e0/ds(i-1)
    temp(i,4,1,ib) = 0.0e0
  End Do
Else If (jj == 2) Then
  Do i = mm,nn
    temp(i,1,1,ib) = 0.0e0
    temp(i,2,1,ib) = 1.0e0/ds(i+1)
    temp(i,3,1,ib) = -1.0e0/ds(i+1)
    temp(i,4,1,ib) = 0.0e0
  End Do
Else
Print*,'jj must be between 0 and 2 in scheme1 subroutine'
Stop
End If
!
If (ibmpfg.eq.1) Then
Do j=1,4
Do i = mm,nn
temp(i,j,1,2) = temp(i,j,1,1)
End Do
End Do
End If
!
Return
End Subroutine scheme1
!
!******************************************************************************
Subroutine scheme2(ds,mm,nn,iuvw,ik,ib,nn1)
!******************************************************************************
!	This subroutine calculates the gradient at cell face with variables at the
!	cell center(4-point)
Implicit None
Integer, Intent( In ) :: nn,mm,ib,iuvw,ik,nn1
Real(mytype), Dimension(-1:nn1+2), Intent( In ) :: ds
Integer :: i,m
Real(mytype) :: a,b,c,d
!
Do i = mm,nn+1
a = ds(i)+0.5e0*ds(i+1)
b = 0.5e0*ds(i)
c = -0.5e0*ds(i-1)
d = -ds(i-1)-0.5e0*ds(i-2)
Call third (i,a,b,c,d,iuvw,ik,ib)
End Do
!
Return
End Subroutine scheme2
!
!******************************************************************************
Subroutine scheme2_bc(ds,mm,nn,ib,ij)
!******************************************************************************
! 	This subroutine determines the gradient at cell face at boundaries with
!	varibles at cell center(IBWEST=1=variable at boundary, ibwest = 0=at cell center)
!	ij=1 x_direction, ij=2 y-direction,ij=3 z direction
Implicit None
Integer, Intent( In ) :: nn,mm,ib,ij
Real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
Integer :: i,m
Real(mytype) :: a,b,c,d
!
Do i = mm,mm+1
a=ds(i)+0.5e0*ds(i+1)
b = 0.5e0*ds(i)
c = -0.5e0*ds(i-1)
d = -ds(i-1)-0.5e0*ds(i-2)
If (ij == 1.and.i == 1) Then
  If (ibwest == 1) Then
    c = 0.0e0
    d = 0.0e0
    Call third3 (i,a,b,c,d,ij,1,1,ib)
  Else If (ibwest == 0) Then
    d = 0.0e0
    Call third3 (i,a,b,c,d,ij,1,1,ib)
  End If
  Do m = 1,4
    gcdx_f(i,m,2,ib) = 0.0e0
  End Do
Else If (ij == 1.and.i == 2) Then
  If (ibwest == 1) Then
    d = -ds(i-1)
    Call third (i,a,b,c,d,ij,1,ib)
  End If
  Call third3 (i,a,b,c,d,ij,1,2,ib) 
Else If (ij == 2.and.i == 1) Then
  If (ibsoth == 1) Then
    c = 0.0e0
    d = 0.0e0
    Call third3 (i,a,b,c,d,ij,1,1,ib)
  Else If (ibsoth == 0) Then
    d = 0.0e0
    Call third3 (i,a,b,c,d,ij,1,1,ib)
  End If
  Do m = 1,4
    gcdy_f(i,m,2,ib) = 0.0e0
  End Do
Else If (ij == 2.and.i == 2) Then
  If (ibsoth == 1) Then
    d = -ds(i-1)
   Call third (i,a,b,c,d,ij,1,ib)
  End If
  Call third3 (i,a,b,c,d,ij,1,2,ib)
Else If (ij == 3.and.i == 1) Then
  If (kbsoth == 1) Then
    c = 0.0e0
    d = 0.0e0
    Call third3 (i,a,b,c,d,ij,1,1,ib)
  Else If (kbsoth == 0) Then
    d = 0.0e0
    Call third3 (i,a,b,c,d,ij,1,1,ib)
  End If
  Do m = 1,4
    gcdz_f(i,m,2,ib) = 0.0e0
  End Do
Else If (ij == 3.and.i == 2) Then
  If (kbsoth == 1) Then
    d = -ds(i-1)
   Call third (i,a,b,c,d,ij,1,ib)
  End If
  Call third3 (i,a,b,c,d,ij,1,2,ib)
Else
  Print*,'ij value must be between 1 and 3 in schem2_bc subroutine'
  stop
End If
End Do
!
Do i = nn,nn+1
 a=ds(i)+0.5e0*ds(i+1)
 b = 0.5e0*ds(i)
 c = -0.5e0*ds(i-1)
 d = -ds(i-1)-0.5e0*ds(i-2)
!
If (ij == 1.and.i == nn) Then 
    If (ibeast == 1) Then
      a=ds(i)
      Call third (i,a,b,c,d,ij,1,ib)
    End If
    Call third3 (i,a,b,c,d,ij,2,2,ib)
  Else If (ij == 1.and.i == (nn+1)) Then
    If (ibeast == 1) Then
      a = 0.0e0
      b = 0.0e0
      Call third3 (i,a,b,c,d,ij,2,1,ib)
    Else If (ibeast == 0) Then
      a = 0.0e0
      Call third3 (i,a,b,c,d,ij,2,1,ib)
    End If
    Do m = 1,4
      gcdx_f(i,m,2,ib) = 0.0e0
    End Do
  Else If (ij == 2.and.i == nn) Then 
   If (ibnoth == 1) Then
     a=ds(i)
     Call third (i,a,b,c,d,ij,1,ib)
   End If
   Call third3 (i,a,b,c,d,ij,2,2,ib)
  Else If (ij == 2 .And. i == nn+1) Then
    If (ibnoth == 1) Then
      a = 0.0e0
      b = 0.0e0
      Call third3 (i,a,b,c,d,ij,2,1,ib)
    Else If (ibnoth == 0) Then
      a = 0.0e0
      Call third3 (i,a,b,c,d,ij,2,1,ib)
    End If
    Do m = 1,4
      gcdy_f(i,m,2,ib) = 0.0e0
    End Do
  Else If (ij == 3.and.i == nn) Then
   If (kbnoth == 1) Then
     a=ds(i)
     Call third (i,a,b,c,d,ij,1,ib)
   End If
   Call third3 (i,a,b,c,d,ij,2,2,ib)
  Else If (ij == 3.and.i == (nn+1)) Then
    If (kbnoth == 1) Then
      a = 0.0e0
      b = 0.0e0
      Call third3 (i,a,b,c,d,ij,2,1,ib)
    Else If (kbnoth == 0) Then
      a = 0.0e0
      Call third3 (i,a,b,c,d,ij,2,1,ib)
    End If
    Do m = 1,4
      gcdz_f(i,m,2,ib) = 0.0e0
    End Do
  Else
    Print*,'IJ VALUE MUST BE BETWEEN 1 and 3 IN SCHEME2_BC'
  End If
End Do
!
Return
End Subroutine scheme2_bc
!
!******************************************************************************
Subroutine scheme3(temp,ds,mm,nn,jj,ib)
!******************************************************************************
!       This subroutine calculates the gradient at cell center using three points
!       when the varables are located at cell center.
!       jj=1 x direction ,jj=2 y direction, jj=3 z direction 
Implicit None
Integer, Intent( In ) :: nn,mm,ib,jj
Real(mytype), Dimension(mm-1:nn+1,5,ngr,nbl), Intent( InOut ) :: temp
Real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
Integer :: i,m
Real(mytype) :: c,d,deno
!
Do i = mm,nn
  c = -(ds(i-1)+ds(i))*0.5e0
  d=(ds(i)+ds(i+1))*0.5e0
  deno=c*d*(c-d)
  temp(i,1,1,ib) = 0.0e0
  temp(i,2,1,ib) = c*c/deno
  temp(i,3,1,ib) = (d+c)*(d-c)/deno
  temp(i,4,1,ib) = -d*d/deno
  temp(i,5,1,ib) = 0.0e0
  Do m = 1,5
    temp(i,m,2,ib) = temp(i,m,1,ib)
  End Do
End Do
!
Return
End Subroutine scheme3
!
!******************************************************************************
Subroutine scheme3_bc(temp1,temp2,ds,mm,nn,jj,ib)
!******************************************************************************
!       This subroutine calculates the gradient at cell center using three points
!       when the varables are located at cell center for convective boundary conditions.
!       jj argument is only used because the dimension of grdszz is different from 
!	othrs two array jj=2 for z direction AND 1 for others
Implicit None
Integer, Intent( In ) :: nn,mm,ib,jj
Real(mytype), Dimension(mm-1:nn+1,5,ngr,nbl), Intent( InOut ) :: temp1
Real(mytype), Dimension(mm:nn/2+jj,5,2), Intent( InOut ) :: temp2
Real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
Integer :: i,m
Real(mytype) :: c,d,deno
!
i=mm
c=(ds(i)+ds(i+1))*0.5e0
d=ds(i+1)+0.5e0*(ds(i)+ds(i+2))
deno=c*d*(c-d)
temp1(i,1,1,ib) = c*c/deno
temp1(i,2,1,ib) = -d*d/deno
temp1(i,3,1,ib) = (d+c)*(d-c)/deno
temp1(i,4,1,ib) = 0.0e0
temp1(i,5,1,ib) = 0.0e0
c = -ds(i)*0.5e0
d=(ds(i)+ds(i+1))*0.5e0
deno=c*d*(c-d)
temp1(i,1,2,ib) = 0.0e0
temp1(i,2,2,ib) = c*c/deno
temp1(i,3,2,ib) = (d+c)*(d-c)/deno
temp1(i,4,2,ib) = -d*d/deno !0.0e0        !Neumann B.c.
temp1(i,5,2,ib) = 0.0e0
!
deno=ds(i)*0.5e0
temp1(i-1,1,2,ib) = 0.0e0
temp1(i-1,2,2,ib) = 1.0e0/deno
temp1(i-1,3,2,ib) = -1.0e0/deno
temp1(i-1,4,2,ib) = 0.0e0 
temp1(i-1,5,2,ib) = 0.0e0
!
i=mm
c=(ds(1)+ds(2)+ds(3)+ds(4))*0.5e0
d=ds(3)+ds(4)+0.5e0*(ds(1)+ds(2)+ds(5)+ds(6))
deno=c*d*(c-d)
temp2(i,1,ib) = c*c/deno
temp2(i,2,ib) = -d*d/deno
temp2(i,3,ib) = (d+c)*(d-c)/deno
temp2(i,4,ib) = 0.0e0
temp2(i,5,ib) = 0.0e0
!
i=nn
c = -(ds(i-1)+ds(i))*0.5e0
d = -ds(i-1)-(ds(i)+ds(i-2))*0.5e0
deno=c*d*(c-d)
temp1(i,1,1,ib) = 0.0e0
temp1(i,2,1,ib) = 0.0e0
temp1(i,3,1,ib) = (d+c)*(d-c)/deno
temp1(i,4,1,ib) = -d*d/deno
temp1(i,5,1,ib) = c*c/deno
c = -(ds(i-1)+ds(i))*0.5e0
d=ds(i)*0.5e0
deno=c*d*(c-d)
temp1(i,1,2,ib) = 0.0e0
temp1(i,2,2,ib) = c*c/deno !0.0e0     !Neumann B.c.
temp1(i,3,2,ib) = (d+c)*(d-c)/deno
temp1(i,4,2,ib) = -d*d/deno
temp1(i,5,2,ib) = 0.0e0
!
deno=ds(i)*0.5
temp1(i+1,1,2,ib) = 0.0e0
temp1(i+1,2,2,ib) = 0.0e0     
temp1(i+1,3,2,ib) = 1.0e0/deno
temp1(i+1,4,2,ib) = -1.0e0/deno
temp1(i+1,5,2,ib) = 0.0e0
!
i=nn/2
c = -(ds(nn-3)+ds(nn-2)+ds(nn-1)+ds(nn))*0.5e0
d = -(ds(nn-2)+ds(nn-3))-(ds(nn)+ds(nn-1)+ds(nn-4)+ds(nn-5))*0.5e0
deno=c*d*(c-d)
temp2(i,1,ib) = 0.0e0
temp2(i,2,ib) = 0.0e0
temp2(i,3,ib) = (d+c)*(d-c)/deno
temp2(i,4,ib) = -d*d/deno
temp2(i,5,ib) = c*c/deno
!
Return
End Subroutine scheme3_bc
!
End Module cal_coefficient
