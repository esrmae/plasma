Module update_field
   Use shared_data
   Use mass_calc
   Use boundary_conditions
   Use mpi
   Use decomp_2d
   Use decomp_2d_io
   Use randgen
   Use flow_field
Implicit None
Integer :: total_step,mpi_comm_used,mpi_comm_other
Real(mytype), Dimension(:), Allocatable :: tauw,rnut,dudy1
Real(mytype), Dimension(:,:), Allocatable :: quad,omega
Real(mytype), Dimension(:,:,:), Allocatable :: wx,temp,arrx,arry,tke,stat
Real(mytype) :: utau,utau0,utau1,utau2,energy,esum,pspent
Character (Len = 40) :: crite
!
! 2D
Real(mytype), Dimension(:,:), Allocatable :: tauw_2d
Real(mytype), Dimension(:,:,:), Allocatable :: quad2d,omega2d
Real(mytype), Dimension(:,:,:,:), Allocatable :: stat2d,tke2d,wx2d
!
Contains
!
!******************************************************************************
Subroutine update 
!******************************************************************************
! This procedure saves current U, V and P and updates U and V after
! solving the Poisson equation for P.
Implicit None
Integer :: i,j,k,iuvw,is,ib,iblock
Real(mytype) :: rms_tmp,rms_new,rms_old,rms_diff,rms_chk,dpdx_old
!
!************ X pencil ************
mpi_comm_used = MPI_COMM_ROW
mpi_comm_other = MPI_COMM_COL
!
iblock = 0; ib=iblock+1
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
   u1(i,j,k) = uhat1(i,j,k)-ae(i,j,k)/dy(xstart(2)-1+j)/dz(xstart(3)-1+k)*dt*(varray(i,j,k)-varray(i-1,j,k))
   v1(i,j,k) = vhat1(i,j,k)-as(i,j,k)/dx(xstart(1)-1+i)/dz(xstart(3)-1+k)*dt*(varray(i,j,k)-varray(i,j-1,k))
   w1(i,j,k) = what1(i,j,k)-at(i,j,k)/dx(xstart(1)-1+i)/dy(xstart(2)-1+j)*dt*(varray(i,j,k)-varray(i,j,k-1))
End Do; End Do; End Do
!
  Call bc_update
  Call div_calc
!
  time=time+dt
  Call statistics 
  Call write_update
!
  If(ipart == 1) Call vis_particle
  If((ilamq == 1.Or.iwfield == 1.Or.iflucs == 1 ).And. Mod(ii,np_snap_3d) == 0) Call calc3d
!
Return
End Subroutine update
!******************************************************************************
Subroutine statistics
!******************************************************************************
! This procedure calculates the turbulent statiscics
Implicit None
!
Allocate (tauw(-1:xsize(2)+2),rnut(-1:xsize(2)+2),stat(-1:xsize(2)+2,4,4),quad(-1:xsize(2)+2,9),         &
        omega(-1:xsize(2)+2,6),tke(-1:xsize(2)+2,4,7),wx(-1:xsize(2)+2,2,8),temp(xsize(1),xsize(3),14),            &
        arrx(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),arry(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),     &
        dudy1(-1:xsize(2)+2))
!
tauw = 0.0e0;rnut = 0.0e0;stat = 0.0e0;quad = 0.0e0;omega = 0.0e0;tke = 0.0e0
wx = 0.0e0;temp = 0.0e0;arrx = 0.0e0;arry = 0.0e0
!
If (lpstat_2d == 1) Then
   If (lpstat_x == 1) Then
      Allocate (tauw_2d(xsize(3),2))
      Allocate (stat2d(-1:xsize(3)+2,-1:xsize(2)+2,4,4),quad2d(-1:xsize(3)+2,-1:xsize(2)+2,9), &
      tke2d(-1:xsize(3)+2,-1:xsize(2)+2,4,7),wx2d(-1:xsize(3)+2,-1:xsize(2)+2,2,8),omega2d(-1:xsize(3)+2,-1:xsize(2)+2,6))
   End If
   If (lpstat_z == 1) Then
      Allocate (tauw_2d(xsize(1),2))
      Allocate (stat2d(-1:xsize(1)+2,-1:xsize(2)+2,4,4),quad2d(-1:xsize(1)+2,-1:xsize(2)+2,9), &
      tke2d(-1:xsize(1)+2,-1:xsize(2)+2,4,7),wx2d(-1:xsize(1)+2,-1:xsize(2)+2,2,8),omega2d(-1:xsize(1)+2,-1:xsize(2)+2,6))
   End If
   stat2d = 0.0e0;quad2d = 0.0e0;tke2d = 0.0e0;wx2d = 0.0e0;omega2d = 0.0e0
   tauw_2d = 0.0e0
!
End If

   Call cal_statistics
   Call write_statistics
!
Deallocate (tauw,rnut,stat,quad,omega,tke,wx,temp,arrx,arry,dudy1)
If (lpstat_2d == 1) Deallocate (stat2d,quad2d,tke2d,wx2d,omega2d,tauw_2d)
!
End Subroutine statistics
!******************************************************************************
Subroutine cal_statistics
!******************************************************************************
!------This subroutine computes turbulence statistics at every nprint time steps.
Implicit None
Integer :: i,j,k,ib,iblock,m,n,jj,j2
!Real(mytype) :: temp_val,ucen,vcen,wcen,uc,vc,dudx,dudy,dudz,dvdx,dvdy,dvdz,       &
!dwdx,dwdy,dwdz,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,omegax,omegay,      &
!omegaz,denox,denoy,denoz,omegx2,omegy2,omegz2,dx2,dy2,dz2,        &
!wwmn,ucent,vcent,wcent,umeany
Real(mytype) :: energy1
!
! For 2d unsteady case
i2dust=0
If(nphase.Ne.1) Then
Do i=1,nphase-1
  If(Mod(Int(time/dt_fixed),Int(tperiod/dt_fixed)).Eq.(i-1)*Int(tperiod/dt_fixed/Real(nphase-1))) Then
    iphase=i; i2dust=1
    If(lpstat_2d == 1) ntmavg2(iphase+1)=ntmavg2(iphase+1)+1;  Exit
  End If
End Do
End If
If(nrank.Eq.0.And.i2dust.Eq.1) Write(*,*) 'ii=',ii,'time=',time,'iphase=',iphase
!
! Caclulate (1d and 2d) statistics
If ((itmavg == 1 .And. (Mod(ii,ntmavg) == 0 .Or. ii .le. 5)) .Or. &
  (lp_ins_1d == 1 .And. (Mod(ii,np_ins_1d) == 0 .Or. ii .le. 5)) .Or. &
   i2dust == 1 .Or. lstop == 1) Then
!
 Call calc_stat
!
 Do j = 0,xsize(2)+1
   umean(j) = stat(j,1,1)
   vmean(j) = stat(j,2,1)
   vmean(j) = 0.0e0
 End Do
!
  Call calc_quad
  Call calc_tke
  Call calc_vorticity
  If(ispec.Eq.1) Call calc_spectra
!
 energy = 0.0e0; energy1 = 0.0e0
 Do j = 1,xsize(2)
   energy1 = energy1 + stat(j,1,2)-stat(j,1,1)**2+stat(j,2,2)+stat(j,3,2)
 End Do
 Call MPI_Allreduce(energy1,energy,1,real_type,MPI_SUM,mpi_comm_other,ierror)
 energy = energy/ny
Else
 crite = 'stat'
 Do j = 1,xsize(2)+1
 temp = 0.0e0
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,1) = u1(i,j,k)
   End Do; End Do
!
   Call avg_stat(temp(1,1,1),crite,j,1,1)
   umean(j) = stat(j,1,1)
 End Do
!
End If
!
! Calculate point statistics
  Call calc_stat_0d
!
Return
End Subroutine cal_statistics
!******************************************************************************
Subroutine avg_stat(temp2d,word,j,m,n)
!******************************************************************************
Implicit none
Integer, Intent (In) :: j,m
Integer, Intent (In), Optional  :: n
Real(mytype), Dimension(xsize(1),xsize(3)), Intent(In) :: temp2d
Real(mytype), Dimension(:), Allocatable :: temp_1d
Real(mytype) :: temp_val
Integer :: dim_val
Character (Len = 40) :: word

temp_val = 0.0e0
If (lpstat_2d == 1) Then
   If (lpstat_x == 1) Then       ! average in the streamwise (x) direction
      Allocate (temp_1d(xsize(3)))   
      dim_val=xsize(3)
   End If
   If (lpstat_z == 1) Then       ! average in the spanwise (z) direction
      Allocate (temp_1d(xsize(1)))   
      dim_val=xsize(1)
   End If
   temp_1d = 0.0e0
End If

If (lpstat_1d == 1) Then
   Call avg_xz(temp2d(1,1),temp_val)
      If (word(1:4) == 'stat') stat(j,m,n) = temp_val
      If (word(1:4) == 'quad') quad(j,m) = temp_val
      If (word(1:3) == 'tke' ) tke(j,m,n) = temp_val
      If (word(1:2) == 'wx'  ) wx(j,m,n) = temp_val
      If (word(1:5) == 'omega') omega(j,m-8) = temp_val
End If
If (lpstat_2d == 1) Then
   If (lpstat_x == 1) Call avg_x(temp2d(1,1),temp_1d(1))
   If (lpstat_z == 1) Call avg_z(temp2d(1,1),temp_1d(1))
      If (word(1:4) == 'stat') stat2d(1:dim_val,j,m,n) = temp_1d(1:dim_val)
      If (word(1:4) == 'quad') quad2d(1:dim_val,j,m) = temp_1d(1:dim_val)
      If (word(1:3) == 'tke' ) tke2d(1:dim_val,j,m,n) = temp_1d(1:dim_val)
      If (word(1:2) == 'wx'  ) wx2d(1:dim_val,j,m,n) = temp_1d(1:dim_val)
      If (word(1:5) == 'omega') omega2d(1:dim_val,j,m-8) = temp_1d(1:dim_val)
End If
 If (lpstat_2d == 1)     Deallocate (temp_1d)   
Return
End Subroutine avg_stat
!
!******************************************************************************
Subroutine avg_x(u2,umnb)
!******************************************************************************
!------ This procedure make a line average in z.
Implicit None
Integer :: i,k
Real(mytype) :: sum1
Real(mytype), Dimension(xsize(1),xsize(3)),Intent(In) :: u2
Real(mytype), Dimension(xsize(3)):: umnb
Do k = 1,xsize(3)
   sum1 = 0.0e0
   Do i = 1,xsize(1)
      sum1=sum1+u2(i,k)
   End Do
   umnb(k) = sum1/real(nx,kind=mytype)
End Do
If(itravz.Eq.1) Call shift_z_fft(umnb)
Return
End Subroutine avg_x
!
!******************************************************************************
Subroutine avg_z(u2,umnb)
!******************************************************************************
!------ This procedure make a line average in z.
Implicit None
Integer :: i,k
Real(mytype) :: sum1, sum1_tot
Real(mytype), Dimension(xsize(1),xsize(3)),Intent(In) :: u2
Real(mytype), Dimension(xsize(1)):: umnb
Do i = 1,xsize(1)
   sum1 = 0.0e0
   Do k = 1,xsize(3)
      sum1=sum1+u2(i,k)
   End Do
   Call MPI_Allreduce(sum1,sum1_tot,1,real_type,MPI_SUM,mpi_comm_used,ierror)
   umnb(i) = sum1_tot/real(nz,kind=mytype)
End Do
If(itravx.Eq.1) Call shift_x_fft(umnb)
Return
End Subroutine avg_z
!******************************************************************************
Subroutine avg_xz(u2,umnb)
!******************************************************************************
!------ This procedure make a plane average in xz plane
Implicit None
Integer :: i,k
Real(mytype) :: sum1,umnb,sum1_tot
Real(mytype), Dimension(xsize(1),xsize(3)),Intent( In ) :: u2
!
sum1 = 0.0e0
Do k = 1,xsize(3)
Do i = 1,xsize(1)
   sum1 = sum1 + u2(i,k)
End Do; End Do
 call MPI_Allreduce(sum1,sum1_tot,1,real_type,MPI_SUM,mpi_comm_used,ierror)
umnb=sum1_tot/real(nx,kind=mytype)/real(nz,kind=mytype)
!
Return
End Subroutine avg_xz
!******************************************************************************
Subroutine shift_x_fft(temp_1d)
!******************************************************************************
Implicit none
Include "fftw3.f"
Real(mytype), Dimension(nx), Intent(inout) :: temp_1d
Real(mytype) :: newx,kappax,cplxr,cplxi,cplxa
Integer :: i
Complex(mytype), Allocatable, Dimension(:) :: outfftxu
Integer*8 :: plan1,plan2
!
Allocate(outfftxu(nx/2+1))
If(iplasma==1.And.iosci==0) Then; newx=cspeed*time;
Else If(iplasma==1.And.iosci==1) Then; newx=cspeed/omeg*sin(omeg*time)
Else; newx=omeg*time/kx; End If
! FORWARD
Call dfftw_plan_dft_r2c_1d(plan1,nx,temp_1d,outfftxu,FFTW_ESTIMATE)
Call dfftw_execute(plan1)
Call dfftw_destroy_plan(plan1)
outfftxu(:)=outfftxu(:)/Real(nx,mytype)
! SHIFT
Do i=1,nx/2+1
  kappax=2.0e0*pi/length*Real(i-1)
  cplxr=Real(outfftxu(i)); cplxi=Aimag(outfftxu(i))
  cplxa=Atan2(cplxi,cplxr)
  outfftxu(i)=Abs(outfftxu(i))*cmplx(cos(cplxa+kappax*newx), &
                  sin(cplxa+kappax*newx))
End Do
! BACKWARD
Call dfftw_plan_dft_c2r_1d(plan2,nx,outfftxu,temp_1d,FFTW_ESTIMATE)
Call dfftw_execute(plan2)
Call dfftw_destroy_plan(plan2)
!
Deallocate(outfftxu)
Return
End Subroutine shift_x_fft
!******************************************************************************
Subroutine shift_z_fft(temp_1d)
!******************************************************************************
Implicit none
Include "fftw3.f"
Real(mytype), Dimension(xsize(3)), Intent(inout) :: temp_1d
Real(mytype), Allocatable, Dimension(:) :: dum_1d
Real(mytype) :: newz,kappaz,cplxr,cplxi,cplxa
Integer :: k,k1d
Complex(mytype), Allocatable, Dimension(:) :: outfftzu
Integer*8 :: plan1,plan2
!
Allocate(outfftzu(nz/2+1),dum_1d(nz))
Call MPI_Allgather(temp_1d,xsize(3),real_type,dum_1d,xsize(3),real_type,mpi_comm_used,ierror)
!
If(iplasma==1.And.iosci==0) Then; newz=cspeed*time;
Else If(iplasma==1.And.iosci==1) Then; newz=cspeed/omeg*sin(omeg*time)
Else; newz=omeg*time/kz; End If
! FORWARD
Call dfftw_plan_dft_r2c_1d(plan1,nz,dum_1d,outfftzu,FFTW_ESTIMATE)
Call dfftw_execute(plan1)
Call dfftw_destroy_plan(plan1)
outfftzu(:)=outfftzu(:)/Real(nz,mytype)
! SHIFT
Do k=1,nz/2+1
  kappaz=2.0e0*pi/width*Real(k-1)
  cplxr=Real(outfftzu(k)); cplxi=Aimag(outfftzu(k))
  cplxa=Atan2(cplxi,cplxr)
  outfftzu(k)=Abs(outfftzu(k))*cmplx(cos(cplxa+kappaz*newz), &
                  sin(cplxa+kappaz*newz))
End Do
! BACKWARD
Call dfftw_plan_dft_c2r_1d(plan2,nz,outfftzu,dum_1d,FFTW_ESTIMATE)
Call dfftw_execute(plan2)
Call dfftw_destroy_plan(plan2)
!
Do k=1,xsize(3); k1d=k+xstart(3)-1
  temp_1d(k)=dum_1d(k1d)
End Do
!
Deallocate(outfftzu,dum_1d)
Return
End Subroutine shift_z_fft
!******************************************************************************
Subroutine shift_x_fd(temp_1d)
!******************************************************************************
Implicit none
Real(mytype), Dimension(nx), Intent(inout) :: temp_1d
Real(mytype), Dimension(nx) :: dum_1d
Real(mytype) :: newx,c,d
Integer :: bel,ab,i
Do i=1,nx
  newx=modulo(xc(i)+omeg*time/kx,length)
  ab=ceiling(newx/dx(1))
  bel=floor(newx/dx(1))
  c=xc(ab)-newx
  d=newx-xc(bel)
  If((c+d)==0.0) then; c=1.0e0; d=1.0e0; End If
  ab=modulo(ab,nx); bel=modulo(bel,nx)
  if (ab==0) ab=nx; if (bel==0) bel=nx 
  dum_1d(i)=(d*temp_1d(ab)+c*temp_1d(bel))/(c+d)
End Do
temp_1d=dum_1d
Return
End Subroutine shift_x_fd
!******************************************************************************
Subroutine shift_z_fd(temp_1d)
!******************************************************************************
Implicit none
Real(mytype), Dimension(xsize(3)), Intent(inout) :: temp_1d
Real(mytype), Dimension(nz) :: dum_1d
Real(mytype) :: newz,c,d
Integer :: bel,ab,k,k1d
dum_1d(:)=0.0e0
Call MPI_Allgather(temp_1d,xsize(3),real_type,dum_1d,xsize(3),real_type,mpi_comm_used,ierror)
Do k=1,xsize(3)
  k1d=k+xstart(3)-1
  newz=modulo(zc(k1d)+omeg*time/kz,width)
  ab=ceiling(newz/dz(1))
  bel=floor(newz/dz(1))
  c=zc(ab)-newz
  d=newz-zc(bel)
  If((c+d)==0.0) then; c=1.0e0; d=1.0e0; End If
  ab=modulo(ab,nz); bel=modulo(bel,nz)
  if (ab==0) ab=nz; if (bel==0) bel=nz
  temp_1d(k)=(d*dum_1d(ab)+c*dum_1d(bel))/(c+d)
End Do
Return
End Subroutine shift_z_fd
!******************************************************************************
Subroutine write_statistics
!******************************************************************************
!------This subroutine computes turbulence statistics at every nprint time steps.
Implicit None
Integer :: i,j,k,ib,iblock,m,n,jj,j2
Real(mytype) :: temp_val,ucen,vcen,wcen,uc,vc,dudx,dudy,dudz,dvdx,dvdy,dvdz,       &
dwdx,dwdy,dwdz,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,omegax,omegay,      &
omegaz,denox,denoy,denoz,omegx2,omegy2,omegz2,dx2,dy2,dz2,        &
wwmn,ucent,vcent,wcent,umn,umeany,taulmo1,taulmo2
!
ib=1
!
!	Instaneous  0d statistics
   Call write_statis_0d
!
!	New binary Instantaneous files
!	SNAP SHOTS (Instaneous Velocity Field Information)
!
   If ((lp_snap_1d == 1 .And. (Mod(ii,np_snap_1d) == 0) .Or. &
       (lp_snap_1d == 1 .And. (lstop == 1 .Or. ii .eq. 1)))) Call snap_1d
   If ((lp_snap_2d == 1 .And. (Mod(ii,np_snap_2d) == 0) .Or. &
       (lp_snap_2d == 1 .And. (lstop == 1 .Or. ii .eq. 1)))) Call snap_2d
!   If (lp_snap_2d == 1 .And. i2dust == 1) Call snap_2d
   If ((lp_snap_3d == 1 .And. (Mod(ii,np_snap_3d) == 0) .Or. &
       (lp_snap_3d == 1 .And. (lstop == 1 .Or. ii .eq. 1)))) Call snap_3d
!
!	Instantaneous statistics (can be space or line averaged)
!
   If (lp_ins_1d == 1 .And. (Mod(ii,np_ins_1d) == 0 .or. lstop == 1 .Or. ii .le. 5)) Call ins_1d
   If (lp_ins_2d == 1 .And. (Mod(ii,np_ins_2d) == 0 .or. lstop == 1 .Or. ii .le. 5)) Call ins_2d
!
!	Time averaged files
!
If ((itmavg == 1 .And. ii == nstep) .Or. (itmavg == 1 .And. lstop == 1)) Then
  taulmo1 = abs(taulmo(1))
  taulmo2 = abs(taulmo(xsize(2)+1))
  Call MPI_Bcast(taulmo1,1,real_type,0,MPI_COMM_CART,ierror)
  Call MPI_Bcast(taulmo2,1,real_type,nproc-1,MPI_COMM_CART,ierror)

  utau0=sqrt(0.5e0*(taulmo1+taulmo2)/samples)
  Call write_output_1d_avg_binary
  If (lpstat_2d == 1)  Call write_output_2d_avg_binary
  If (nrank == 0) Call stat_parameters
  If (ispec == 1) Call write_3d_spectra
  If (ishear == 1) Call write_shear_2d
End If
!
Return
End Subroutine write_statistics
!******************************************************************************
Subroutine write_statis_0d
!******************************************************************************
! This subroutine calculates total energy and creates file for time history of instentaneous
!ity and pressure at selected positions at every time step.
Implicit None
Integer :: i,j,k,j2
Real(mytype) :: temp_val
!
!	Statistics every time step
   Call rmsfl
   If(ibcjet==1.And.jet_type==2) Call pow_spent
   If (nrank==0) Then
     If(ibcjet==1.And.jet_type==2) Then; Write(30,97) ii,time,dt,utau,utau1,utau2,dpdx_mean,rmsflo/height/width,pspent
     Else; Write(30,97) ii,time,dt,utau,utau1,utau2,dpdx_mean,rmsflo/height/width; End If
   End If
!
!       turbulence statistics

If ((lp_ins_1d == 1 .And.(Mod(ii,np_ins_1d) == 0 .Or. ii .le. 5)) .Or. lstop == 1) Then
!	Statistics every 5 time steps
   If (nrank==0) Write(31,97) ii,time,dt,utau,rmsflo/height/width,esum,energy
End If
!
97  Format(I5,16E15.7)
!
Return
End Subroutine write_statis_0d
!******************************************************************************
Subroutine pow_spent
!******************************************************************************
! This subroutine calculates energy spent by flow control
Implicit None
Integer :: i,j,k,j1d,k1d
Real(mytype) :: pspent_dum
!
pspent=0.0e0; pspent_dum=0.0e0
!
! Calculate the power spent
Do i=1,xsize(1); Do k = 1,xsize(3); k1d=xstart(3)-1+k
  Do j = 1,xsize(2); j1d = xstart(2)-1+j
    pspent_dum=pspent_dum+Abs(w1(i,j,k)*fbody(i,j,k)*dx(i)*dy(j1d)*dzs(k1d))
End Do; End Do; End Do
!
! Sums over results on all processes
  Call MPI_Allreduce(pspent_dum,pspent,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
!
Return
End Subroutine pow_spent
!
!******************************************************************************
Subroutine write_update
!*****************************************************************************
! THis subroutine writes the velocity field information to the log file at the 
! end of simulation
Implicit None
Integer :: i,j,k
Real(mytype) :: resms
!
! CALCULATION OF RESIDUAL MASS
  Call res_mass(resms)
!
  If (irandm /= 0 .And. Mod(ii,20) == 0) Call random_flow_generator
!
  Call rmsfl
!
   If (nrank==0) Write(21,1904) rmsfli,rmsflm,rmsflo,dmsfl
1904   Format(' mass flow at_inlet=',E15.8,/' mass flow at_middle=',      &
        E15.8,/' mass flow at_outlet=',                                    &
        E15.8,/' Difference in massflow=',E15.8)                     
   If (nrank==0) Write(21,223) resms/rmsfli
223    Format(' Residual mass/mass flux at inlet=',E15.8)
!
Return
End Subroutine write_update
!******************************************************************************
Subroutine random_flow_generator
!******************************************************************************
Implicit None
Integer :: i,j,k,jmin,jmax
!
If (nrank==0) Write(21,*)'Random flow amplitude=',amp

jmin=1; jmax=xsize(2)
If (xmin(2)==1) jmin=3
If (xmax(2)==1) jmax=xsize(2)-2

Do k = 1,xsize(3)
  Do j = jmin,jmax
    Do i = 1,xsize(1)
      u1(i,j,k) = u1(i,j,k)+amp*(rand_zbql()-0.5e0)
      v1(i,j,k) = v1(i,j,k)+amp*(rand_zbql()-0.5e0)
      w1(i,j,k) = w1(i,j,k)+amp*(rand_zbql()-0.5e0)
    End Do
  End Do
End Do
!
!       Wall B.C.
!
Call bc_update
!
Return
End Subroutine random_flow_generator
!******************************************************************************
Subroutine cal_dt
!******************************************************************************
!------ This procedure adjusts Dt after completing a whole time step according
!       to the given CFL number.
Implicit None
Integer :: i,j,k
Real(mytype) :: temp3,tmpmx,rem,dtold,dtnew,dtmax,tmpmx_tot
!
tmpmx = 0.0e0; tmpmx_tot = 0.0e0
Do k = 1,xsize(3)
Do j = 1,xsize(2)
Do i = 1,xsize(1)
temp3 = (u1(i,j,k)+u1(i+1,j,k))/2.0e0/dx(i)&
+(v1(i,j,k)+v1(i,j+1,k))/2.0e0/dy(j)&
+(w1(i,j,k)+w1(i,j,k+1))/2.0e0/dz(k)
If (tmpmx < temp3) tmpmx = temp3
End Do; End Do; End Do
!
 call MPI_Allreduce(tmpmx,tmpmx_tot,1,real_type,MPI_MAX,MPI_COMM_CART,ierror)
!
dtold = dt
dtmax = 0.4e0*re/(retau*retau)
dtnew = cfl/tmpmx_tot
dt = dtold*0.25e0 + dtnew*0.75e0
If (idt_ret == 1 .Or. idt_ret_un == 1) dt = Min(dt,dtmax)
If (ifdt_fixed == 1) Then; dt =dt_fixed
ElseIf (ifdt_variable == 1) Then ; Call calc_variable_dt
EndIf
!
   If (nrank==0) Write(21,142) tmpmx_tot*dt,dt,dtnew
 142    Format('CFL condition CFL= ',E15.7,' DT= ',2E15.7)
 98     Format(5E15.7)
!
Return
End Subroutine cal_dt
!
!******************************************************************************
Subroutine calc_variable_dt
!******************************************************************************
Implicit None
!
If ((time-time0).le.time_dtcon) Then
dt=dt_initial
If (nrank == 0) Write (21,*) 'Time step constant initially',time,dt
Else If ((time-time0-time_dtcon).lt.deltat_dt) Then
dt =dt*(time-time0-time_dtcon)/deltat_dt
dt=Max(dt,dt_initial)
If (nrank == 0) Write (21,*) 'Time step increases gradually',time,dt
End If
!
Return
End Subroutine calc_variable_dt
!
!******************************************************************************
Subroutine snap_1d
!******************************************************************************
 Implicit None
 Integer :: i,j,k,idec,jdec,kdec,arb_write
 Real(mytype), Dimension(:), Allocatable  :: u_line,v_line,w_line,p_line,dum
!
	idec=xp_snap_1d; jdec=yp_snap_1d; kdec=zp_snap_1d
        nsmpl_snap_1d=nsmpl_snap_1d+1
	arb_write=0
!
!	Snapshot along x
  Call decomp_2d_write_simphead(fh_sn1d,disp_sn1d)
If (lp_snap_x == 1) Then
Allocate (u_line(-1:xsize(1)+2),v_line(-1:xsize(1)+2),w_line(-1:xsize(1)+2),p_line(-1:xsize(1)+2), &
          dum(1:xwrite(1)))
  u_line(:)=0.0e0;v_line(:)=0.0e0;w_line(:)=0.0e0;p_line(:)=0.0e0
  If (((xstart(2).le.jdec).and.(xend(2).ge.jdec)).and.((xstart(3).le.kdec).and.(xend(3).ge.kdec))) then
    jdec=jdec-xstart(2);kdec=kdec-xstart(3)
    Do i=-1,xsize(1)+1
      u_line(i)=0.5e0*(u1(i,jdec,kdec)+u1(i+1,jdec,kdec))
      v_line(i)=0.5e0*(v1(i,jdec,kdec)+v1(i,jdec+1,kdec))
      w_line(i)=0.5e0*(w1(i,jdec,kdec)+w1(i,jdec,kdec+1))
      p_line(i)=p1(i,jdec,kdec)
    End Do
    arb_write=1
  End If
  !Write(218)u_line; Write(218)v_line; Write(218)w_line; Write(218)p_line
  dum(1:xwrite(1))=u_line(lo(1):up(1))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nx,1,dum,arb_write)
  dum(1:xwrite(1))=v_line(lo(1):up(1))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nx,1,dum,arb_write)
  dum(1:xwrite(1))=w_line(lo(1):up(1))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nx,1,dum,arb_write)
  dum(1:xwrite(1))=p_line(lo(1):up(1))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nx,1,dum,arb_write)
Deallocate(u_line,v_line,w_line,p_line,dum)
Endif
!	Snapshot along y
If (lp_snap_y == 1) Then
Allocate (u_line(-1:xsize(2)+2),v_line(-1:xsize(2)+2),w_line(-1:xsize(2)+2),p_line(-1:xsize(2)+2), &
          dum(1:xwrite(2)))
  u_line(:)=0.0e0;v_line(:)=0.0e0;w_line(:)=0.0e0;p_line(:)=0.0e0
  If (((xstart(1).le.idec).and.(xend(1).ge.idec)).and.((xstart(3).le.kdec).and.(xend(3).ge.kdec))) then
    idec=idec-xstart(1);kdec=kdec-xstart(3)
    Do j=-1,xsize(2)+1
      u_line(j)=0.5e0*(u1(idec,j,kdec)+u1(idec+1,j,kdec))
      v_line(j)=0.5e0*(v1(idec,j,kdec)+v1(idec,j+1,kdec))
      w_line(j)=0.5e0*(w1(idec,j,kdec)+w1(idec,j,kdec+1))
      p_line(j)=p1(idec,j,kdec)
    End Do
    arb_write=1
  End If
  !  IS SPLIT, USE MPI WRITE
  !Write(218)u_line; Write(218)v_line; Write(218)w_line; Write(218)p_line
  dum(1:xwrite(2))=u_line(lo(2):up(2))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,ny,2,dum,arb_write)
  dum(1:xwrite(2))=v_line(lo(2):up(2))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,ny,2,dum,arb_write)
  dum(1:xwrite(2))=w_line(lo(2):up(2))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,ny,2,dum,arb_write)
  dum(1:xwrite(2))=p_line(lo(2):up(2))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,ny,2,dum,arb_write)
Deallocate(u_line,v_line,w_line,p_line,dum)
Endif
!	Snapshot along z
If (lp_snap_z == 1) Then
Allocate (u_line(-1:xsize(3)+2),v_line(-1:xsize(3)+2),w_line(-1:xsize(3)+2),p_line(-1:xsize(3)+2), &
          dum(1:xwrite(3)))
  u_line(:)=0.0e0;v_line(:)=0.0e0;w_line(:)=0.0e0;p_line(:)=0.0e0
  If (((xstart(1).le.idec).and.(xend(1).ge.idec)).and.((xstart(2).le.jdec).and.(xend(2).ge.jdec))) then
    idec=idec-xstart(1);jdec=jdec-xstart(2)
    Do k=-1,xsize(3)+1
      u_line(k)=0.5e0*(u1(idec,jdec,k)+u1(idec+1,jdec,k))
      v_line(k)=0.5e0*(v1(idec,jdec,k)+v1(idec,jdec+1,k))
      w_line(k)=0.5e0*(w1(idec,jdec,k)+w1(idec,jdec,k+1))
      p_line(k)=p1(idec,jdec,k)
    End Do
    arb_write=1
  End If
  !  IS SPLIT, USE MPI WRITE
  !Write(218)u_line; Write(218)v_line; Write(218)w_line; Write(218)p_line
  dum(1:xwrite(3))=u_line(lo(3):up(3))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nz,3,dum,arb_write)
  dum(1:xwrite(3))=v_line(lo(3):up(3))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nz,3,dum,arb_write)
  dum(1:xwrite(3))=w_line(lo(3):up(3))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nz,3,dum,arb_write)
  dum(1:xwrite(3))=p_line(lo(3):up(3))
  Call decomp_2d_write_1dx_oli(fh_sn1d,disp_sn1d,nz,3,dum,arb_write)
Deallocate(u_line,v_line,w_line,p_line,dum)
Endif
!
Return
End Subroutine snap_1d
!
!******************************************************************************
Subroutine snap_2d
!******************************************************************************
 Implicit None
 Integer :: i,j,k,idec,jdec,kdec,arb_write
 Integer, Dimension(:), Allocatable :: arbwrite
 Real(mytype), Dimension(:,:), Allocatable :: u_line,v_line,w_line,p_line,dum
 Real(mytype), dimension(:), Allocatable :: vc,vcpx,vcmx,vcpy,vcmy,vcpz,vcmz, &
                       wc,wcpx,wcmx,wcpy,wcmy,wcpz,wcmz
 Real(mytype) :: dwdy,dvdz,omegax
 Integer :: j1d,k1d,ib,isn
!
idec=xp_snap_2d; jdec=yp_snap_2d; kdec=zp_snap_2d
nsmpl_snap_2d=nsmpl_snap_2d+1
arb_write=0
!
!       Snapshot along x
!
Do isn=1,nposxz
  Call decomp_2d_write_simphead(fh_sn2d(isn),disp_sn2d(isn))
End Do
If (lp_snap_xy == 1) Then
Allocate (u_line(-1:xsize(1)+2,-1:xsize(2)+2),v_line(-1:xsize(1)+2,-1:xsize(2)+2), &
          w_line(-1:xsize(1)+2,-1:xsize(2)+2),p_line(-1:xsize(1)+2,-1:xsize(2)+2), &
          dum(xwrite(1),xwrite(2)))
  u_line(:,:)=0.0e0;v_line(:,:)=0.0e0;w_line(:,:)=0.0e0;p_line(:,:)=0.0e0
  If ((xstart(3).le.kdec).and.(xend(3).ge.kdec)) then
    kdec=kdec-xstart(3)
    Do j=-1,xsize(2)+1; Do i=-1,xsize(1)+1
      u_line(i,j)=0.5e0*(u1(i,j,kdec)+u1(i+1,j,kdec))
      v_line(i,j)=0.5e0*(v1(i,j,kdec)+v1(i,j+1,kdec))
      w_line(i,j)=0.5e0*(w1(i,j,kdec)+w1(i,j,kdec+1))
      p_line(i,j)=p1(i,j,kdec)
    End Do; End Do
    arb_write=1
  End If
  !  IS SPLIT, USE MPI WRITE
  !Write(220)u_line; Write(220)v_line; Write(220)w_line; Write(220)p_line
  dum(1:xwrite(1),1:xwrite(2))=u_line(lo(1):up(1),lo(2):up(2))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),nx,ny,1,dum,arb_write)
  dum(1:xwrite(1),1:xwrite(2))=v_line(lo(1):up(1),lo(2):up(2))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),nx,ny,1,dum,arb_write)
  dum(1:xwrite(1),1:xwrite(2))=w_line(lo(1):up(1),lo(2):up(2))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),nx,ny,1,dum,arb_write)
  dum(1:xwrite(1),1:xwrite(2))=p_line(lo(1):up(1),lo(2):up(2))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),nx,ny,1,dum,arb_write)
Deallocate(u_line,v_line,w_line,p_line,dum)
Endif
!
!       Snapshot along y
!
If (lp_snap_yz == 1) Then
Allocate (u_line(-1:xsize(2)+2,-1:xsize(3)+2),v_line(-1:xsize(2)+2,-1:xsize(3)+2), &
          w_line(-1:xsize(2)+2,-1:xsize(3)+2),p_line(-1:xsize(2)+2,-1:xsize(3)+2), &
          dum(xwrite(2),xwrite(3)))
  idec=idec-xstart(1)
  u_line(:,:)=0.0e0;v_line(:,:)=0.0e0;w_line(:,:)=0.0e0;p_line(:,:)=0.0e0
  Do k=-1,xsize(3)+1; Do j=-1,xsize(2)+1
    u_line(j,k)=0.5e0*(u1(idec,j,k)+u1(idec+1,j,k))
    v_line(j,k)=0.5e0*(v1(idec,j,k)+v1(idec,j+1,k))
    w_line(j,k)=0.5e0*(w1(idec,j,k)+w1(idec,j,k+1))
    p_line(j,k)=p1(idec,j,k)
  End Do; End Do
  arb_write=1
  !  IS SPLIT, USE MPI WRITE
  !Write(220)u_line; Write(220)v_line; Write(220)w_line; Write(220)p_line
  dum(1:xwrite(2),1:xwrite(3))=u_line(lo(2):up(2),lo(3):up(3))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),ny,nz,3,dum,arb_write)
  dum(1:xwrite(2),1:xwrite(3))=v_line(lo(2):up(2),lo(3):up(3))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),ny,nz,3,dum,arb_write)
  dum(1:xwrite(2),1:xwrite(3))=w_line(lo(2):up(2),lo(3):up(3))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),ny,nz,3,dum,arb_write)
  !  WRITE THE STREAMWISE VORTICITY WX
!  !------------------------------------------------------------------------------------
!  Allocate(vc(1:nx),vcpx(1:nx),vcmx(1:nx),vcpy(1:nx),vcmy(1:nx),vcpz(1:nx),vcmz(1:nx), &
!	   wc(1:nx),wcpx(1:nx),wcmx(1:nx),wcpy(1:nx),wcmy(1:nx),wcpz(1:nx),wcmz(1:nx))
!  ib=1
!  p_line(:,:)=0.0e0
!  Do j = 0,xsize(2)+1
!    j1d=xstart(2)-1+j
!    Do k = 1,xsize(3)
!     k1d=xstart(3)-1+k
!     Do i = idec,idec
!      vc(i)    = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
!      vcmx(i)  = 0.5e0*(v1(i-1,j,k)+v1(i-1,j+1,k))
!      vcpx(i)  = 0.5e0*(v1(i+1,j,k)+v1(i+1,j+1,k))
!      vcmz(i)  = 0.5e0*(v1(i,j,k-1)+v1(i,j+1,k-1))
!      vcpz(i)  = 0.5e0*(v1(i,j,k+1)+v1(i,j+1,k+1))
!      wc(i)     = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
!      wcmx(i)   = 0.5e0*(w1(i-1,j,k)+w1(i-1,j,k+1))
!      wcpx(i)   = 0.5e0*(w1(i+1,j,k)+w1(i+1,j,k+1))
!      wcmy(i)   = 0.5e0*(w1(i,j-1,k)+w1(i,j-1,k+1))
!      wcpy(i)   = 0.5e0*(w1(i,j+1,k)+w1(i,j+1,k+1))
!      dwdy = gcdy_c(j1d,2,2,ib)*wcpy(i)+gcdy_c(j1d,3,2,ib)*wc(i) &
!            +gcdy_c(j1d,4,2,ib)*wcmy(i)
!      dvdz = gcdz_c(k1d,2,1,ib)*vcpz(i)+gcdz_c(k1d,3,1,ib)*vc(i) &
!            +gcdz_c(k1d,4,1,ib)*vcmz(i)
!      omegax = dwdy - dvdz
!      p_line(j,k) = omegax
! End Do; End Do; End Do
! Deallocate(vc,vcpx,vcmx,vcpy,vcmy,vcpz,vcmz,wc,wcpx,wcmx,wcpy,wcmy,wcpz,wcmz)
! !------------------------------------------------------------------------------------
  dum(1:xwrite(2),1:xwrite(3))=p_line(lo(2):up(2),lo(3):up(3))
  Call decomp_2d_write_2dx_oli(fh_sn2d(1),disp_sn2d(1),ny,nz,3,dum,arb_write)
Deallocate(u_line,v_line,w_line,p_line,dum)
Endif
!
!       Snapshot along z
!
If (lp_snap_xz == 1) Then
Allocate (u_line(-1:xsize(1)+2,-1:xsize(3)+2),v_line(-1:xsize(1)+2,-1:xsize(3)+2), &
          w_line(-1:xsize(1)+2,-1:xsize(3)+2),p_line(-1:xsize(1)+2,-1:xsize(3)+2), &
          dum(xwrite(1),xwrite(3)))
  u_line(:,:)=0.0e0;v_line(:,:)=0.0e0;w_line(:,:)=0.0e0;p_line(:,:)=0.0e0
!
  Allocate(arbwrite(1:nposxz)); arbwrite(1:nposxz)=0
  Do isn=1,nposxz
    If (nposxz .Gt. 1) jdec=posxz(isn)
    If ((xstart(2).le.jdec).and.(xend(2).ge.jdec)) then
      jdec=jdec-xstart(2)
      Do k=-1,xsize(3)+1; Do i=-1,xsize(1)+1
        u_line(i,k)=0.5e0*(u1(i,jdec,k)+u1(i+1,jdec,k))
        v_line(i,k)=0.5e0*(v1(i,jdec,k)+v1(i,jdec+1,k))
        w_line(i,k)=0.5e0*(w1(i,jdec,k)+w1(i,jdec,k+1))
        p_line(i,k)=p1(i,jdec,k)
      End Do; End Do
      arbwrite(isn)=1
    End If
  !  IS SPLIT, USE MPI WRITE
  !Write(220)u_line; Write(220)v_line; Write(220)w_line; Write(220)p_line
    dum(1:xwrite(1),1:xwrite(3))=u_line(lo(1):up(1),lo(3):up(3))
    Call decomp_2d_write_2dx_oli(fh_sn2d(isn),disp_sn2d(isn),nx,nz,2,dum,arbwrite(isn))
    dum(1:xwrite(1),1:xwrite(3))=v_line(lo(1):up(1),lo(3):up(3))
    Call decomp_2d_write_2dx_oli(fh_sn2d(isn),disp_sn2d(isn),nx,nz,2,dum,arbwrite(isn))
    dum(1:xwrite(1),1:xwrite(3))=w_line(lo(1):up(1),lo(3):up(3))
    Call decomp_2d_write_2dx_oli(fh_sn2d(isn),disp_sn2d(isn),nx,nz,2,dum,arbwrite(isn))
    dum(1:xwrite(1),1:xwrite(3))=p_line(lo(1):up(1),lo(3):up(3))
    Call decomp_2d_write_2dx_oli(fh_sn2d(isn),disp_sn2d(isn),nx,nz,2,dum,arbwrite(isn))
  End Do
Deallocate(u_line,v_line,w_line,p_line,dum)
Endif
!
Return
End Subroutine snap_2d
!
!******************************************************************************
Subroutine Snap_3d
!******************************************************************************
Implicit None
Real(mytype), Dimension(:,:,:), Allocatable :: dum
!
nsmpl_snap_3d=nsmpl_snap_3d+1
 Call decomp_2d_write_head3d(fh_sn3d,disp_sn3d)

 allocate(dum(xwrite(1),xwrite(2),xwrite(3)))
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=u1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_sn3d,disp_sn3d,nx,ny,nz,1,dum)
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=v1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_sn3d,disp_sn3d,nx,ny,nz,1,dum)
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=w1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_sn3d,disp_sn3d,nx,ny,nz,1,dum)
 dum(1:xwrite(1),1:xwrite(2),1:xwrite(3))=p1(lo(1):up(1),lo(2):up(2),lo(3):up(3))
 Call decomp_2d_write_var(fh_sn3d,disp_sn3d,nx,ny,nz,1,dum)
 deallocate(dum)
!
Return
End Subroutine Snap_3d
!
!******************************************************************************
Subroutine write_3d_spectra
!******************************************************************************
Implicit None
Integer :: nx1
!
 nx1=nx/2+1
 Call decomp_2d_write_var_fft(fh_spec,disp_spec,nx1,ny,nz,3,spec_sum(:,:,:,1))
 Call decomp_2d_write_var_fft(fh_spec,disp_spec,nx1,ny,nz,3,spec_sum(:,:,:,2))
 Call decomp_2d_write_var_fft(fh_spec,disp_spec,nx1,ny,nz,3,spec_sum(:,:,:,3))
 Call decomp_2d_write_var_fft(fh_spec,disp_spec,nx1,ny,nz,3,spec_sum(:,:,:,4))
!Write(*,*) spec_sum
!
Return
End Subroutine write_3d_spectra
!
!******************************************************************************
Subroutine div_calc
!******************************************************************************
!	This subroutine calculates the divergence for final flow field
Implicit None
Integer :: i,j,k,iloc_m,jloc_m,kloc_m
Real(mytype) :: fsum1,fsum,div_max
Real(mytype), Dimension(:,:,:), Allocatable :: farray_temp
iloc_m = 0.0e0;jloc_m = 0.0e0;kloc_m = 0.0e0

!
Allocate (farray_temp(xsize(1),xsize(2),xsize(3)))
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
  farray_temp(i,j,k) =  ((u1(i+1,j,k)-u1(i,j,k))/dx(xstart(1)-1+i)    &
                        +(v1(i,j+1,k)-v1(i,j,k))/dy(xstart(2)-1+j)    &
                        +(w1(i,j,k+1)-w1(i,j,k))/dz(xstart(3)-1+k))
End Do; End Do; End Do
!
!       Compatibility condition test.
!
fsum = 0.0e0; fsum1 = 0.0e0 
Do k = 1,xsize(3); Do j = 1,xsize(2); Do i = 1,xsize(1)
  fsum1=fsum1+farray_temp(i,j,k)
End Do; End Do; End Do
Call MPI_Allreduce(fsum1,fsum,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
!
If (nrank == 0) Write(21,1048) ii,fsum
 1048   Format(' Sum of divergence for final flowfield at  II=',I6,' FSUM_FINAL=',E20.9)
Deallocate (farray_temp)
!
Return
End Subroutine div_calc
!******************************************************************************
Subroutine write_output_1d_avg_binary
!******************************************************************************
!	This subroutine write 1d time averaged data in binary format
Implicit None
!
  Integer :: i,j,k,arb_write
  Real(mytype), Dimension(:), Allocatable :: dum
!
  Allocate(dum(xwrite(2)))
  arb_write=0
  If (xmin(3)==1) arb_write=1
!
!  do j=0,xsize(2)+1
!    if(nrank==0) write(*,*) nrank,j+xstart(2)-1,omega_sum(j,3)
!  end do

  Do j=1,4; Do i=1,4 
    dum(1:xwrite(2))=stat_sum(lo(2):up(2),i,j)
    Call decomp_2d_write_1dx_oli(fh_av1d,disp_av1d,ny,2,dum,arb_write)
  End Do; End Do
  Do i=1,9 
    dum(1:xwrite(2))=quad_sum(lo(2):up(2),i)
    Call decomp_2d_write_1dx_oli(fh_av1d,disp_av1d,ny,2,dum,arb_write)
  End Do
  Do i=1,6 
    dum(1:xwrite(2))=omega_sum(lo(2):up(2),i)
    Call decomp_2d_write_1dx_oli(fh_av1d,disp_av1d,ny,2,dum,arb_write)
  End Do
  Do j=1,7; Do i=1,4 
    dum(1:xwrite(2))=tke_sum(lo(2):up(2),i,j)
    Call decomp_2d_write_1dx_oli(fh_av1d,disp_av1d,ny,2,dum,arb_write)
  End Do; End Do
  Do j=1,8; Do i=1,2 
    dum(1:xwrite(2))=wx_sum(lo(2):up(2),i,j)
    Call decomp_2d_write_1dx_oli(fh_av1d,disp_av1d,ny,2,dum,arb_write)
  End Do; End Do
Close(211)
!
Return
End Subroutine write_output_1d_avg_binary
!******************************************************************************
 Subroutine write_output_2d_avg_binary
!******************************************************************************
Implicit None
!
  Integer :: i,j,k,arb_write,iph
  Real(mytype), Dimension(:,:), Allocatable :: dum
!
arb_write=0
If (lpstat_x == 1) then  !averaged in x direction
  Allocate(dum(xwrite(2),xwrite(3)))
  arb_write=1
  Do iph=1,nphase
  disp_av2d=0_MPI_OFFSET_KIND
  If(iph .Gt. 1) Then
    Call decomp_2d_write_simphead_real(fh_av2d(iph),Real(iph-2,mytype)/Real(nphase-1),disp_av2d)
    Call decomp_2d_write_simphead_int(fh_av2d(iph),ntmavg2(iph),disp_av2d)
  Endif 
  Do j=1,4; Do i=1,4 
    dum(1:xwrite(2),1:xwrite(3))=transpose(stat2d_sum(lo(3):up(3),lo(2):up(2),i,j,iph))
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,ny,nz,3,dum,arb_write)
  End Do; End Do
  Do i=1,9 
    dum(1:xwrite(2),1:xwrite(3))=transpose(quad2d_sum(lo(3):up(3),lo(2):up(2),i,iph))
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,ny,nz,3,dum,arb_write)
  End Do
  Do i=1,6
    dum(1:xwrite(2),1:xwrite(3))=transpose(omega2d_sum(lo(3):up(3),lo(2):up(2),i,iph))
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,ny,nz,3,dum,arb_write)
  End Do
  Do j=1,7; Do i=1,4 
    dum(1:xwrite(2),1:xwrite(3))=transpose(tke2d_sum(lo(3):up(3),lo(2):up(2),i,j,iph))
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,ny,nz,3,dum,arb_write)
  End Do; End Do
  Do j=1,8; Do i=1,2
    dum(1:xwrite(2),1:xwrite(3))=transpose(wx2d_sum(lo(3):up(3),lo(2):up(2),i,j,iph))
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,ny,nz,3,dum,arb_write)
  End Do; End Do
  End Do
End If
If (lpstat_z == 1) then
  If (xmin(3)==1) arb_write=1
  Allocate(dum(xwrite(1),xwrite(2)))
  Do iph=1,nphase
  disp_av2d=0_MPI_OFFSET_KIND
  If(iph .Gt. 1) Then
    Call decomp_2d_write_simphead_real(fh_av2d(iph),Real(iph-2,mytype)/Real(nphase-1),disp_av2d)
    Call decomp_2d_write_simphead_int(fh_av2d(iph),ntmavg2(iph),disp_av2d)
  Endif
  Do j=1,4; Do i=1,4 
    dum(1:xwrite(1),1:xwrite(2))=stat2d_sum(lo(1):up(1),lo(2):up(2),i,j,iph)
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,nx,ny,1,dum,arb_write)
  End Do; End Do
  Do i=1,9 
    dum(1:xwrite(1),1:xwrite(2))=quad2d_sum(lo(1):up(1),lo(2):up(2),i,iph)
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,nx,ny,1,dum,arb_write)
  End Do
  Do i=1,6
    dum(1:xwrite(1),1:xwrite(2))=omega2d_sum(lo(1):up(1),lo(2):up(2),i,iph)
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,nx,ny,1,dum,arb_write)
  End Do
  Do j=1,7; Do i=1,4 
    dum(1:xwrite(1),1:xwrite(2))=tke2d_sum(lo(1):up(1),lo(2):up(2),i,j,iph)
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,nx,ny,1,dum,arb_write)
  End Do; End Do
  Do j=1,8; Do i=1,2
    dum(1:xwrite(1),1:xwrite(2))=wx2d_sum(lo(1):up(1),lo(2):up(2),i,j,iph)
    Call decomp_2d_write_2dx_oli(fh_av2d(iph),disp_av2d,nx,ny,1,dum,arb_write)
  End Do; End Do
  End Do
End If
Deallocate(dum)
!
Return
End Subroutine write_output_2d_avg_binary
!******************************************************************************
Subroutine ins_1d
!******************************************************************************
Implicit None
!
  Integer :: i,j,k,arb_write
  Real(mytype), Dimension(:), Allocatable :: dum
!
  Allocate(dum(xwrite(2)))
!
arb_write=0
nsmpl_ins_1d=nsmpl_ins_1d+1
!
  If (xmin(3)==1) arb_write=1
  Call decomp_2d_write_simphead(fh_ins1d,disp_ins1d)
  Do j=1,4; Do i=1,4 
    dum(1:xwrite(2))=stat(lo(2):up(2),i,j)
    Call decomp_2d_write_1dx_oli(fh_ins1d,disp_ins1d,ny,2,dum,arb_write)
  End Do; End Do
  Do i=1,9 
    dum(1:xwrite(2))=quad(lo(2):up(2),i)
    Call decomp_2d_write_1dx_oli(fh_ins1d,disp_ins1d,ny,2,dum,arb_write)
  End Do
  Do i=1,6
    dum(1:xwrite(2))=omega(lo(2):up(2),i)
    Call decomp_2d_write_1dx_oli(fh_ins1d,disp_ins1d,ny,2,dum,arb_write)
  End Do
  Do j=1,7; Do i=1,4 
    dum(1:xwrite(2))=tke(lo(2):up(2),i,j)
    Call decomp_2d_write_1dx_oli(fh_ins1d,disp_ins1d,ny,2,dum,arb_write)
  End Do; End Do
  Do j=1,8; Do i=1,2
    dum(1:xwrite(2))=wx(lo(2):up(2),i,j)
    Call decomp_2d_write_1dx_oli(fh_ins1d,disp_ins1d,ny,2,dum,arb_write)
  End Do; End Do
!
Return
End Subroutine ins_1d
!
!******************************************************************************
Subroutine ins_2d
!******************************************************************************
Implicit None
!
  Integer :: i,j,k,arb_write
  Real(mytype), Dimension(:,:), Allocatable :: dum
!
arb_write=0
nsmpl_ins_2d=nsmpl_ins_2d+1
  Call decomp_2d_write_simphead(fh_ins2d,disp_ins2d)
If (lpstat_x == 1) then
  Allocate(dum(xwrite(2),xwrite(3)))
  arb_write=1
  Do j=1,4; Do i=1,4 
    dum(1:xwrite(2),1:xwrite(3))=transpose(stat2d(lo(3):up(3),lo(2):up(2),i,j))
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,ny,nz,3,dum,arb_write)
  End Do; End Do
  Do i=1,9 
    dum(1:xwrite(2),1:xwrite(3))=transpose(quad2d(lo(3):up(3),lo(2):up(2),i))
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,ny,nz,3,dum,arb_write)
  End Do
  Do i=1,6
    dum(1:xwrite(2),1:xwrite(3))=transpose(omega2d(lo(3):up(3),lo(2):up(2),i))
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,ny,nz,3,dum,arb_write)
  End Do
  Do j=1,7; Do i=1,4 
    dum(1:xwrite(2),1:xwrite(3))=transpose(tke2d(lo(3):up(3),lo(2):up(2),i,j))
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,ny,nz,3,dum,arb_write)
  End Do; End Do
  Do j=1,8; Do i=1,2
    dum(1:xwrite(2),1:xwrite(3))=transpose(wx2d(lo(3):up(3),lo(2):up(2),i,j))
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,ny,nz,3,dum,arb_write)
  End Do; End Do
End If
If (lpstat_z == 1) then
  Allocate(dum(xwrite(1),xwrite(2)))
  If (xmin(3) == 1) arb_write=1
  Do j=1,4; Do i=1,4 
    dum(1:xwrite(1),1:xwrite(2))=stat2d(lo(1):up(1),lo(2):up(2),i,j)
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,nx,ny,1,dum,arb_write)
  End Do; End Do
  Do i=1,9 
    dum(1:xwrite(1),1:xwrite(2))=quad2d(lo(1):up(1),lo(2):up(2),i)
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,nx,ny,1,dum,arb_write)
  End Do
  Do i=1,6
    dum(1:xwrite(1),1:xwrite(2))=omega2d(lo(1):up(1),lo(2):up(2),i)
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,nx,ny,1,dum,arb_write)
  End Do
  Do j=1,7; Do i=1,4 
    dum(1:xwrite(1),1:xwrite(2))=tke2d(lo(1):up(1),lo(2):up(2),i,j)
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,nx,ny,1,dum,arb_write)
  End Do; End Do
  Do j=1,8; Do i=1,2
    dum(1:xwrite(1),1:xwrite(2))=wx2d(lo(1):up(1),lo(2):up(2),i,j)
    Call decomp_2d_write_2dx_oli(fh_ins2d,disp_ins2d,nx,ny,1,dum,arb_write)
  End Do; End Do
End If
Deallocate(dum)
!
Return
End Subroutine ins_2d
!******************************************************************************
Subroutine write_shear_2d
!******************************************************************************
Implicit None
!
  Integer :: i,j,k,arb_write
  Real(mytype), Dimension(:,:), Allocatable :: dum,wsh
!
Allocate(dum(xwrite(1),xwrite(3)))
dum(:,:)=0.0e0
arb_write=0
nsmpl_sh_2d=nsmpl_sh_2d+1
  Call decomp_2d_write_simphead_real(fh_sh2d,time,disp_sh2d)
If (xmin(2) == 1) then
  Allocate(wsh(-1:xsize(1)+2,-1:xsize(3)+2))
  arb_write=1; wsh(:,:)=0.0e0
  j=1 !lower wall
  Do k=-1,xsize(3)+2; Do i=-1,xsize(1)+2
    wsh(i,k)=(gcdy_f(j,1,1,1)*u1(i,j+1,k) + gcdy_f(j,2,1,1)*u1(i,j,k)&
             +gcdy_f(j,3,1,1)*u1(i,j-1,k) + gcdy_f(j,4,1,1)*u1(i,j-2,k)) *nu
  End Do; End Do
  dum(1:xwrite(1),1:xwrite(3))=wsh(lo(1):up(1),lo(3):up(3))
  Deallocate(wsh)
End If
!
Call decomp_2d_write_2dx_oli(fh_sh2d,disp_sh2d,nx,nz,2,dum,arb_write)
!
Deallocate(dum)
!
Return
End Subroutine write_shear_2d
!******************************************************************************
Subroutine calc_stat
!******************************************************************************
implicit None
Integer :: i,j,k,ib,m,n
!
 ib=1
!        
crite = 'stat'
m = 1  ! u
Do j = 0,xsize(2)+1
   temp = 0.0e0
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,1)=0.5e0*(u1(i,j,k)+u1(i+1,j,k))
   End Do; End Do
!
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,2) = temp(i,k,1)*temp(i,k,1)
      temp(i,k,3) = temp(i,k,1)*temp(i,k,2)
      temp(i,k,4) = temp(i,k,2)*temp(i,k,2)
   End Do; End Do
!
   Do n=1,4
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
m = 2 ! v
Do j = 0,xsize(2)+1
   temp = 0.0e0
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,1)=0.5e0*(v1(i,j,k)+v1(i,j+1,k))
   End Do; End Do
!
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,2) = temp(i,k,1)*temp(i,k,1)
      temp(i,k,3) = temp(i,k,1)*temp(i,k,2)
      temp(i,k,4) = temp(i,k,2)*temp(i,k,2)
   End Do; End Do
!
   Do n=1,4
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
m = 3    ! w
Do j = 0,xsize(2)+1
   temp = 0.0e0
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,1)=0.5e0*(w1(i,j,k)+w1(i,j,k+1))
   End Do; End Do
!
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,2) = temp(i,k,1)*temp(i,k,1)
      temp(i,k,3) = temp(i,k,1)*temp(i,k,2)
      temp(i,k,4) = temp(i,k,2)*temp(i,k,2)
   End Do; End Do
!
   Do n=1,4
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
m = 4     ! p
Do j = 0,xsize(2)+1
   temp = 0.0e0
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,1)=p1(i,j,k)
   End Do; End Do
!
   Do k = 1,xsize(3); Do i = 1,xsize(1)
      temp(i,k,2) = temp(i,k,1)*temp(i,k,1)
      temp(i,k,3) = temp(i,k,1)*temp(i,k,2)
      temp(i,k,4) = temp(i,k,2)*temp(i,k,2)
   End Do; End Do
!
   Do n=1,4
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
!Time averaging
If ((itmavg == 1 .and. Mod (ii,ntmavg) == 0) .Or. lstop == 1 .Or. ii.le.5 ) Then
stat_sum = stat_sum + stat     ! resb & resbo
If (lpstat_2d == 1) stat2d_sum(:,:,:,:,1) = stat2d_sum(:,:,:,:,1) + &
	stat2d(:,:,:,:)
End If 
If(lpstat_2d == 1 .And. nphase.Ne.1 .And. i2dust == 1) stat2d_sum(:,:,:,:,iphase+1) = &
	stat2d_sum(:,:,:,:,iphase+1) + stat2d(:,:,:,:)
!
Return 
End Subroutine calc_stat
!
!******************************************************************************
Subroutine calc_quad
!******************************************************************************
implicit None
Integer :: i,j,k,ib,m,n,j1d
Real(mytype), Dimension(nx) :: uc, vc, wc
!
 ib=1
! Quadratic analysis at the centre of the CV.
crite = 'quad'
Do j = 1,xsize(2)
   j1d=xstart(2)-1+j
   temp = 0.0e0
!
   Do k = 1,xsize(3)
   Do i = 1,xsize(1)
      uc(i) = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      vc(i) = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      wc(i) = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
   End Do
!
   Do i = 1,xsize(1)
   If (ichann /= 1) vc(i) = vc(i) - stat(j,2,1)
   temp(i,k,5) = uc(i)*vc(i)
   temp(i,k,6) = vc(i)*wc(i)
   temp(i,k,7) = uc(i)*wc(i)
   uc(i) = uc(i) - stat(j,1,1)
   If (vc(i) > 0.0e0) Then
      If (uc(i) > 0.0e0) temp(i,k,1) = uc(i)*vc(i)
      If (uc(i) < 0.0e0) temp(i,k,2) = uc(i)*vc(i)
   Else If (vc(i) < 0.0e0) Then
      If (uc(i) < 0.0e0) temp(i,k,3) = uc(i)*vc(i)
      If (uc(i) > 0.0e0) temp(i,k,4) = uc(i)*vc(i)
   End If
!!  u*v at edge
!   if((dy(j-1)+dy(j)) .ne. 0.0e0) then
!     uc(i) = (u1(i,j-1,k)*dy(j)+u1(i,j,k)*dy(j-1))/(dy(j-1)+dy(j))
!   end if
!   vc(i) = (v1(i,j,k)+v1(i-1,j,k))*0.5e0
   End Do

!   Do i = 1,xsize(1)
!   temp(i,k,6) = uc(i)*vc(i)
!   temp(i,k,7) = rnutyz1(i,j,k)
!   End Do
   End Do
!
   Do m=1,7
   Call avg_stat(temp(1,1,m),crite,j,m)
   End Do

   quad(j,8) = (gcdy_f(j1d,1,1,ib)*umean(j+1)+gcdy_f(j1d,2,1,ib)*umean(j)+     &
                gcdy_f(j1d,3,1,ib)*umean(j-1)+gcdy_f(j1d,4,1,ib)*umean(j-2))*nu

   quad(j,9) = (gcdy_f(j1d,1,1,ib)*stat(j+1,3,1)+gcdy_f(j1d,2,1,ib)*stat(j,3,1)+     &
                gcdy_f(j1d,3,1,ib)*stat(j-1,3,1)+gcdy_f(j1d,4,1,ib)*stat(j-2,3,1))*nu
!   If (isgs_model == 0 ) Then ;quad(j,9) =cs(j)
!   Else If (isgs_model == 1 ) Then ;quad(j,9) =c_w
!   End If
End Do
!
! Time averaging
If (itmavg == 1 .And. (Mod (ii,ntmavg) == 0 .Or. ii .le. 5 .Or. lstop == 1)) Then
quad_sum = quad_sum + quad    ! rstres & rstrso
If (lpstat_2d == 1) quad2d_sum(:,:,:,1) = quad2d_sum(:,:,:,1) + &
	quad2d(:,:,:)
Endif
! 2d unsteady statistics
If (lpstat_2d == 1 .And. itmavg == 1 .And. nphase .Ne. 1 .And. i2dust == 1) &
quad2d_sum(:,:,:,iphase+1) = quad2d_sum(:,:,:,iphase+1) + quad2d(:,:,:)
!
Return
End Subroutine calc_quad
!        
!******************************************************************************
Subroutine calc_tke
!******************************************************************************
implicit None
Integer :: i,j,k,ib,m,n,j1d,k1d
Real(mytype), dimension(nx) :: uc,ucp,ucm,vc,vcp,vcm,wc,wcp,wcm,pc
Real(mytype) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
 ib=1
!       Turbulent kinetic energy transport == ation.
!       uu
!       TKEO(I,J,IV,IC)
!       IV=1 for uu
!       IV=2 for vv
!       IV=3 for ww
!       IV=4 for uv
!       IC=1 for Pk   : production
!       IC=2 for Tk   : turbulent diffusion or turbulent transport
!       IC=3 for Dk   : viscous diffusion
!       IC=4 for Pik  : velocity-pressure gradient or pressure diffusion
!       IC=5 for Psik : pressure strain correlation
!       IC=6 for Ek   : viscous dissipation
!       IC=7 for Ek   : viscous dissipation
!       M=1 for uu
!       1: Tk : uuu
!       2: Tk : uuv
!       3: Pik : pu
!       4: Psik : p*dudx
!       5: Ek : dudx**2
!       6: Ek : dudy**2
!       7: Ek : dudz**2
!       IBLOCK: IBLOCK = 0
!
crite = 'tke'
m = 1
Do j = 0,xsize(2)+1
   j1d=xstart(2)-1+j
   Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   Do i = 1,xsize(1)
      ucm(i) = 0.5e0*(u1(i+1,j,k)+u1(i+2,j,k))
      uc(i)  = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucp(i) = 0.5e0*(u1(i-1,j,k)+u1(i,j,k))
      !!this part is not right: ucm: minus (left); ucp: positive (right); to keep consistency, it is not changed
      !ucm(i) = 0.5e0*(u1(i-1,j,k)+u1(i,j,k))
      !uc(i)  = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      !ucp(i) = 0.5e0*(u1(i+1,j,k)+u1(i+2,j,k))

      vc(i) = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      pc(i) = p1(i,j,k)
   End Do

   Do i = 1,xsize(1)
      temp(i,k,1) = uc(i)*uc(i)*uc(i)
      temp(i,k,2) = uc(i)*uc(i)*vc(i)
      temp(i,k,3) = pc(i)*uc(i)

      dudx = gcdx_c(i,2,1,ib)*ucp(i)+gcdx_c(i,3,1,ib)*uc(i) &
            +gcdx_c(i,4,1,ib)*ucm(i)
      temp(i,k,4) = pc(i)*dudx
      temp(i,k,5) = dudx*dudx

      ucm(i) = 0.5e0*(u1(i,j-1,k)+u1(i+1,j-1,k))
      uc(i)  = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucp(i) = 0.5e0*(u1(i,j+1,k)+u1(i+1,j+1,k))
   End Do

   Do i = 1,xsize(1)
      dudy = gcdy_c(j1d,2,2,ib)*ucp(i)+gcdy_c(j1d,3,2,ib)*uc(i) & 
            +gcdy_c(j1d,4,2,ib)*ucm(i)
      temp(i,k,6) = dudy*dudy

      ucm(i) = 0.5e0*(u1(i,j,k-1)+u1(i+1,j,k-1))
      uc(i)  = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucp(i) = 0.5e0*(u1(i,j,k+1)+u1(i+1,j,k+1))
   End Do

   Do i = 1,nx
      dudz = gcdz_c(k1d,2,1,ib)*ucp(i)+gcdz_c(k1d,3,1,ib)*uc(i) &
            +gcdz_c(k1d,4,1,ib)*ucm(i)
      temp(i,k,7) = dudz*dudz
   End Do; End Do
!
   Do n=1,7
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
!       M=2 for vv
!
!       1: Tk : uvv
!       2: Tk : vvv
!       3: Pik : pv
!       4: Psik : p*dvdy
!       5: Ek : dvdx**2
!       6: Ek : dvdy**2
!       7: Ek : dvdz**2
!
m = 2
Do j = 0,xsize(2)+1
   j1d=xstart(2)-1+j
   Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   Do i = 1,xsize(1)
      uc(i)  = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))

      vcm(i) = 0.5e0*(v1(i,j-1,k)+v1(i,j,k))
      vc(i)  = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      if (j==xsize(2)+1) then
        vcp(i)=0.5e0*v1(i,j+1,k)
      else
        vcp(i) = 0.5e0*(v1(i,j+1,k)+v1(i,j+2,k))
      end if

      pc(i) = p1(i,j,k)
   End Do

   Do i = 1,xsize(1)
      temp(i,k,1) = uc(i)*vc(i)*vc(i)
      temp(i,k,2) = vc(i)*vc(i)*vc(i)
      temp(i,k,3) = pc(i)*vc(i)

      dvdy = gcdy_c(j1d,2,2,ib)*vcp(i)+gcdy_c(j1d,3,2,ib)*vc(i) & 
            +gcdy_c(j1d,4,2,ib)*vcm(i)
      temp(i,k,4) = pc(i)*dvdy
      temp(i,k,6) = dvdy*dvdy

      vcm(i) = 0.5e0*(v1(i-1,j,k)+v1(i-1,j+1,k))
      vc(i)  = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      vcp(i) = 0.5e0*(v1(i+1,j,k)+v1(i+1,j+1,k))
   End Do

   Do i = 1,xsize(1)
      dvdx = gcdx_c(i,2,1,ib)*vcp(i)+gcdx_c(i,3,1,ib)*vc(i) &
            +gcdx_c(i,4,1,ib)*vcm(i)
      temp(i,k,5) = dvdx*dvdx

      vcm(i) = 0.5e0*(v1(i,j,k-1)+v1(i,j+1,k-1))
      vc(i)  = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      vcp(i) = 0.5e0*(v1(i,j,k+1)+v1(i,j+1,k+1))    
   End Do

   Do i = 1,xsize(1)
      dvdz = gcdz_c(k1d,2,1,ib)*vcp(i)+gcdz_c(k1d,3,1,ib)*vc(i) &
            +gcdz_c(k1d,4,1,ib)*vcm(i)
      temp(i,k,7) = dvdz*dvdz
   End Do; End Do
!
   Do n=1,7
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
!       M=3 for ww
!       1: Tk : uww
!       2: Tk : vww
!       3: Pik : pw
!       4: Psik : p*dwdz
!       5: Ek : dwdx**2
!       6: Ek : dwdy**2
!       7: Ek : dwdz**2
m = 3
Do j = 0,xsize(2)+1
   j1d=xstart(2)-1+j
   Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   Do i = 1,xsize(1)
      uc(i)  = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      vc(i)  = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))

      wcm(i) = 0.5e0*(w1(i,j,k-1)+w1(i,j,k))
      wc(i)  = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
      wcp(i) = 0.5e0*(w1(i,j,k+1)+w1(i,j,k+2))

      pc(i) = p1(i,j,k)
   End Do

   Do i = 1,xsize(1)
      temp(i,k,1) = uc(i)*wc(i)*wc(i)
      temp(i,k,2) = vc(i)*wc(i)*wc(i)
      temp(i,k,3) = pc(i)*wc(i)

      dwdz = gcdz_c(k1d,2,1,ib)*wcp(i)+gcdz_c(k1d,3,1,ib)*wc(i) &
            +gcdz_c(k1d,4,1,ib)*wcm(i)
      temp(i,k,4) = pc(i)*dwdz
      temp(i,k,7) = dwdz*dwdz

      wcm(i) = 0.5e0*(w1(i-1,j,k)+w1(i-1,j,k+1))
      wc(i)  = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
      wcp(i) = 0.5e0*(w1(i+1,j,k)+w1(i+1,j,k+1))
   End Do

   Do i = 1,xsize(1)
      dwdx = gcdx_c(i,2,1,ib)*wcp(i)+gcdx_c(i,3,1,ib)*wc(i) &
            +gcdx_c(i,4,1,ib)*wcm(i)
      temp(i,k,5) = dwdx*dwdx

      wcm(i) = 0.5e0*(w1(i,j-1,k)+w1(i,j-1,k+1))
      wc(i)  = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
      wcp(i) = 0.5e0*(w1(i,j+1,k)+w1(i,j+1,k+1))
   End Do

   Do i = 1,xsize(1)
      dwdy = gcdy_c(j1d,2,2,ib)*wcp(i)+gcdy_c(j1d,3,2,ib)*wc(i) &
            +gcdy_c(j1d,4,2,ib)*wcm(i)
      temp(i,k,6) = dwdy*dwdy
   End Do; End Do
!
   Do n=1,7
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
!       M=4 for uv
!       1: Pik : pv
!       2: Pik : pu
!       3: Psik : p*dudy
!       4: Psik : p*dvdx
!       5: Ek : dudx*dvdx
!       6: Ek : dudy*dvdy
!       7: Ek : dudz*dvdz
!
m = 4
Do j = 0,xsize(2)+1
   j1d=xstart(2)-1+j
   Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   Do i = 1,xsize(1)
      ucm(i)  = 0.5e0*(u1(i-1,j,k)+u1(i,j,k))
      uc(i)  =  0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucp(i)  = 0.5e0*(u1(i+1,j,k)+u1(i+2,j,k))

      vcm(i)  = 0.5e0*(v1(i-1,j,k)+v1(i-1,j+1,k))
      vc(i)  =  0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      vcp(i)  = 0.5e0*(v1(i+1,j,k)+v1(i+1,j+1,k))

      pc(i) = p1(i,j,k)
   End Do

   Do i = 1,xsize(1)
      temp(i,k,1) = pc(i)*vc(i)
      temp(i,k,2) = pc(i)*uc(i)

      dudx = gcdx_c(i,2,1,ib)*ucp(i)+gcdx_c(i,3,1,ib)*uc(i) &
            +gcdx_c(i,4,1,ib)*ucm(i)
      dvdx = gcdx_c(i,2,1,ib)*vcp(i)+gcdx_c(i,3,1,ib)*vc(i) &
            +gcdx_c(i,4,1,ib)*vcm(i)

      temp(i,k,4) = pc(i)*dvdx
      temp(i,k,5) = dudx*dvdx

      ucm(i)  = 0.5e0*(u1(i,j-1,k)+u1(i+1,j-1,k))
      uc(i)   = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucp(i)  = 0.5e0*(u1(i,j+1,k)+u1(i+1,j+1,k))

      vcm(i)  = 0.5e0*(v1(i,j-1,k)+v1(i,j,k))
      vc(i)   = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      if (j==xsize(2)+1) then
        vcp(i)=0.5e0*v1(i,j+1,k)
      else
        vcp(i)  = 0.5e0*(v1(i,j+1,k)+v1(i,j+2,k))
      end if
   End Do

   Do i = 1,xsize(1)
      dudy = gcdy_c(j1d,2,2,ib)*ucp(i)+gcdy_c(j1d,3,2,ib)*uc(i) &
            +gcdy_c(j1d,4,2,ib)*ucm(i)
     
      dvdy = gcdy_c(j1d,2,2,ib)*vcp(i)+gcdy_c(j1d,3,2,ib)*vc(i) & 
            +gcdy_c(j1d,4,2,ib)*vcm(i)

      temp(i,k,3) = pc(i)*dudy
      temp(i,k,6) = dudy*dvdy

      ucm(i)  = 0.5e0*(u1(i,j,k-1)+u1(i+1,j,k-1))
      uc(i)   = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucp(i)  = 0.5e0*(u1(i,j,k+1)+u1(i+1,j,k+1))

      vcm(i)  = 0.5e0*(v1(i,j,k-1)+v1(i,j+1,k-1))
      vc(i)   = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      vcp(i)  = 0.5e0*(v1(i,j,k+1)+v1(i,j+1,k+1))
   End Do

   Do i = 1,xsize(1)
      dudz = gcdz_c(k1d,2,1,ib)*ucp(i)+gcdz_c(k1d,3,1,ib)*uc(i) &
            +gcdz_c(k1d,4,1,ib)*ucm(i)
      dvdz = gcdz_c(k1d,2,1,ib)*vcp(i)+gcdz_c(k1d,3,1,ib)*vc(i) &
            +gcdz_c(k1d,4,1,ib)*vcm(i)
      temp(i,k,7) = dudz*dvdz
   End Do; End Do
!
   Do n=1,7
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
! Time avaraging
If (itmavg == 1 .And. ((Mod (ii,ntmavg) == 0) .Or. lstop == 1 .Or. ii .le. 5)) Then
tke_sum = tke_sum + tke	! tkeo
If (lpstat_2d == 1) tke2d_sum(:,:,:,:,1) = tke2d_sum(:,:,:,:,1) + &
	tke2d(:,:,:,:) ! tke2do
End If
! 2d unsteady statistics
If (lpstat_2d == 1 .And. itmavg == 1 .And. nphase .Ne. 1 .And. i2dust == 1) &
tke2d_sum(:,:,:,:,iphase+1) = tke2d_sum(:,:,:,:,iphase+1) + tke2d(:,:,:,:) ! tke2do
!
Return
End Subroutine calc_tke
!
!******************************************************************************
Subroutine calc_vorticity
!******************************************************************************
implicit None
Integer :: i,j,k,ib,m,n,j1d,k1d
Real(mytype) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dx2,dy2,dz2,       &
dwdx,dwdy,dwdz,omegax,omegay,omegaz,omegx2,omegy2,omegz2
Real(mytype), dimension(nx) :: uc,ucpx,ucmx,ucpy,ucmy,ucpz,ucmz, &
                       vc,vcpx,vcmx,vcpz,vcmz, &
                       wc,wcpx,wcmx,wcpy,wcmy
!
! vorticity fluctuations
ib =1

m = 1
Do j = 0,xsize(2)+1
   j1d=xstart(2)-1+j
   Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   Do i = 1,xsize(1)
      uc(i)    = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucmx(i)  = 0.5e0*(u1(i-1,j,k)+u1(i,j,k))
      ucpx(i)  = 0.5e0*(u1(i+1,j,k)+u1(i+2,j,k))
      ucmy(i)  = 0.5e0*(u1(i,j-1,k)+u1(i+1,j-1,k))
      ucpy(i)  = 0.5e0*(u1(i,j+1,k)+u1(i+1,j+1,k))
      ucmz(i)  = 0.5e0*(u1(i,j,k-1)+u1(i+1,j,k-1))
      ucpz(i)  = 0.5e0*(u1(i,j,k+1)+u1(i+1,j,k+1))

      vc(i)    = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      vcmx(i)  = 0.5e0*(v1(i-1,j,k)+v1(i-1,j+1,k))
      vcpx(i)  = 0.5e0*(v1(i+1,j,k)+v1(i+1,j+1,k))
      vcmz(i)  = 0.5e0*(v1(i,j,k-1)+v1(i,j+1,k-1))
      vcpz(i)  = 0.5e0*(v1(i,j,k+1)+v1(i,j+1,k+1))

      wc(i)     = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
      wcmx(i)   = 0.5e0*(w1(i-1,j,k)+w1(i-1,j,k+1))
      wcpx(i)   = 0.5e0*(w1(i+1,j,k)+w1(i+1,j,k+1))
      wcmy(i)   = 0.5e0*(w1(i,j-1,k)+w1(i,j-1,k+1))
      wcpy(i)   = 0.5e0*(w1(i,j+1,k)+w1(i,j+1,k+1))
 !     If(nrank==0 .and. j==1 .and. i==10 .and. k==60)  write(*,*) j,nrank,wcmy(i),wc(i),wcpy(i)
 !     If(nrank==0 .and. j==2 .and. i==10 .and. k==60)  write(*,*) j,nrank,wcmy(i),wc(i),wcpy(i)
   End Do
   Do i = 1,xsize(1)
      dudx = gcdx_c(i,2,1,ib)*ucpx(i)+gcdx_c(i,3,1,ib)*uc(i) &
            +gcdx_c(i,4,1,ib)*ucmx(i)
      dudy = gcdy_c(j1d,2,2,ib)*ucpy(i)+gcdy_c(j1d,3,2,ib)*uc(i) &
            +gcdy_c(j1d,4,2,ib)*ucmy(i)
      dudz = gcdz_c(k1d,2,1,ib)*ucpz(i)+gcdz_c(k1d,3,1,ib)*uc(i) &
            +gcdz_c(k1d,4,1,ib)*ucmz(i)

      dvdx = gcdx_c(i,2,1,ib)*vcpx(i)+gcdx_c(i,3,1,ib)*vc(i) &
            +gcdx_c(i,4,1,ib)*vcmx(i)
      dvdz = gcdz_c(k1d,2,1,ib)*vcpz(i)+gcdz_c(k1d,3,1,ib)*vc(i) &
            +gcdz_c(k1d,4,1,ib)*vcmz(i)

      dwdx = gcdx_c(i,2,1,ib)*wcpx(i)+gcdx_c(i,3,1,ib)*wc(i) &
            +gcdx_c(i,4,1,ib)*wcmx(i)
      dwdy = gcdy_c(j1d,2,2,ib)*wcpy(i)+gcdy_c(j1d,3,2,ib)*wc(i) &
            +gcdy_c(j1d,4,2,ib)*wcmy(i)

!      If(nrank==0 .and. k==1 .and. i==1)  write(*,*) j1d,nrank,gcdy_c(j1d,2,2,ib),gcdy_c(j1d,3,2,ib),gcdy_c(j1d,4,2,ib)
!      If(nrank==0 .and. j==2 .and. i==10)  write(*,*) j,k,nrank,dwdy,dvdz

      omegax = dwdy - dvdz
      omegay = dudz - dwdx 
      omegaz = dvdx - dudy
 
!      If(nrank==0 .and. j==1)  write(*,*) i,j,k,omegaz
!      If(nrank==0 .and. j==2)  write(*,*) i,j,k,omegaz

      omegx2 = omegax * omegax
      omegy2 = omegay * omegay
      omegz2 = omegaz * omegaz
!
      temp(i,k,9) = omegx2
      temp(i,k,10) = omegy2
      temp(i,k,11) = omegz2
      temp(i,k,12) = omegax
      temp(i,k,13) = omegay
      temp(i,k,14) = omegaz
!
      arrx(i,j,k) = omegax
      arry(i,j,k) = omegx2

!
! Stretching : wx*wx*dudx
      temp(i,k,1) = omegx2*dudx
! Tilting : wx*wy*dudy
      temp(i,k,2) = omegax*omegay*dudy
! Tilting : wx*wy
      temp(i,k,3) = omegax*omegay
! Tilting : wx*wz*dudz
      temp(i,k,4) = omegax*omegaz*dudz
! Tilting : wx*dudz
      temp(i,k,5) = omegax*dudz
! Convection : u*wx*wx
      temp(i,k,6) = uc(i)*omegx2
! Convection : v*wx*wx
      temp(i,k,7) = vc(i)*omegx2
! Convection : w*wx*wx
      temp(i,k,8) = wc(i)*omegx2
!
   End Do; End Do
!
crite = 'omega'
   Do n=9,14
   Call avg_stat(temp(1,1,n),crite,j,n)
   End Do
!
crite = 'wx'
   Do n=1,8
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
! Time averaging
If ((itmavg == 1 .and. Mod (ii,ntmavg) == 0) .Or. lstop == 1) Then
omega_sum = omega_sum + omega
If (lpstat_2d == 1) omega2d_sum(:,:,:,1) = omega2d_sum(:,:,:,1) + &
	omega2d(:,:,:)
Endif
! 2d unsteady statistics
If (lpstat_2d == 1 .And. itmavg == 1 .And. nphase .Ne. 1 .And. i2dust == 1) &
omega2d_sum(:,:,:,iphase+1) = omega2d_sum(:,:,:,iphase+1) + omega2d(:,:,:)
!
m = 2
!
Do j = 1,xsize(2)
   j1d=xstart(2)-1+j
   Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   dz2=dz(k)*dz(k)
   Do i = 1,xsize(1)
   dx2=dx(i)*dx(i)
!
! Diffusion : d2(wx*wx)/dx2
      uc(i) = (arry(i+1,j,k)-2.0e0*arry(i,j,k)+arry(i-1,j,k))/dx2
      temp(i,k,1) = uc(i)*uc(i)
!
! Diffusion : d2(wx*wx)/dy2 from Ferziger and Peric (1996) pp. 48
      vc(i) = gcdy_f(j1d,1,1,ib)*arry(i,j+1,k)+gcdy_f(j1d,2,1,ib)*arry(i,j,k)&
              +gcdy_f(j1d,3,1,ib)*arry(i,j-1,k)+gcdy_f(j1d,4,1,ib)*arry(i,j-2,k)
      wc(i) = gcdy_f(j1d+1,1,1,ib)*arry(i,j+2,k)+gcdy_f(j1d+1,2,1,ib)*arry(i,j+1,k)	&
              +gcdy_f(j1d+1,3,1,ib)*arry(i,j,k)+gcdy_f(j1d+1,4,1,ib)*arry(i,j-1,k)
      uc(i) = (wc(i)-vc(i))/dy(j)
      temp(i,k,2) = uc(i)*uc(i)
!
! Diffusion : d2(wx*wx)/dz2
      uc(i) = (arry(i,j,k+1)-2.0e0*arry(i,j,k)+arry(i,j,k-1))/dz2
      temp(i,k,3) = uc(i)*uc(i)
!
! Dissipation : (dwxdx)**2
      uc(i) = gcdx_c(i,1,1,ib)*arrx(i+2,j,k)+gcdx_c(i,2,1,ib)*arrx(i+1,j,k) &
              +gcdx_c(i,3,1,ib)*arrx(i,j,k)+gcdx_c(i,4,1,ib)*arrx(i-1,j,k)        &
              +gcdx_c(i,5,1,ib)*arrx(i-2,j,k)
      temp(i,k,4) = uc(i)*uc(i)
!
! Dissipation : (dwxdy)**2
      vc(i) = gcdy_c(j,1,1,ib)*arrx(i,j+2,k)+gcdy_c(j,2,1,ib)*arrx(i,j+1,k) &
              +gcdy_c(j,3,1,ib)*arrx(i,j,k)+gcdy_c(j,4,1,ib)*arrx(i,j-1,k)        &
              +gcdy_c(j,5,1,ib)*arrx(i,j-2,k)
      temp(i,k,5) = vc(i)*vc(i)
!
! Dissipation : (dwxdz)**2
      wc(i) = gcdz_c(k1d,1,1,ib)*arrx(i,j,k+2)+gcdz_c(k1d,2,1,ib)*arrx(i,j,k+1)           &
              +gcdz_c(k1d,3,1,ib)*arrx(i,j,k)+gcdz_c(k1d,4,1,ib)*arrx(i,j,k-1)                  &
              +gcdz_c(k1d,5,1,ib)*arrx(i,j,k-2)
      temp(i,k,6) = wc(i)*wc(i)
   End Do; End Do
   temp(:,:,7) = temp(:,:,8)
!
   Do n=1,7 ! note n = 1,6 not 7
   Call avg_stat(temp(1,1,n),crite,j,m,n)
   End Do
End Do
!
! Time averaging
If (itmavg == 1 .and. (Mod (ii,ntmavg) == 0 .Or. lstop == 1 .Or. ii .le. 5)) Then
wx_sum = wx_sum + wx
If (lpstat_2d == 1) wx2d_sum(:,:,:,:,1) = wx2d_sum(:,:,:,:,1) + &
	wx2d(:,:,:,:)
Endif
! 2d unsteady statistics
If (lpstat_2d == 1 .And. itmavg == 1 .And. nphase .Ne. 1 .And. i2dust == 1) &
wx2d_sum(:,:,:,:,iphase+1) = wx2d_sum(:,:,:,:,iphase+1) + wx2d(:,:,:,:)
!
Return
End Subroutine calc_vorticity
!
!******************************************************************************
Subroutine calc_spectra
!******************************************************************************
Use decomp_2d_fft
Implicit None
Complex(mytype), Dimension(:,:,:), Allocatable :: ccapz,ccapz2
Integer :: i,j,k,i1,ib,i1d,j1d,k1d
Allocate(ccapz(fft_zsize(1),fft_zsize(2),fft_zsize(3)),ccapz2(fft_zsize(1),fft_zsize(2),fft_zsize(3)))
!
Do k=1,xsize(3); Do j=1,xsize(2); Do i=1,xsize(1)
  farray(i,j,k)=u1(i,j,k)
End Do; End Do; End Do
Call decomp_2d_fft_3d(farray,ccapz)
spec_sum(:,:,:,1)=spec_sum(:,:,:,1)+Conjg(ccapz(:,:,:))*ccapz(:,:,:)/Real(nx*nz,mytype)
!
Do k=1,xsize(3); Do j=1,xsize(2); Do i=1,xsize(1)
  farray(i,j,k)=v1(i,j,k)
End Do; End Do; End Do
Call decomp_2d_fft_3d(farray,ccapz2)
spec_sum(:,:,:,2)=spec_sum(:,:,:,2)+Real(Conjg(ccapz(:,:,:))*ccapz2(:,:,:))/Real(nx*nz,mytype)
spec_sum(:,:,:,4)=spec_sum(:,:,:,4)+Conjg(ccapz2(:,:,:))*ccapz2(:,:,:)/Real(nx*nz,mytype)
!
Do k=1,xsize(3); Do j=1,xsize(2); Do i=1,xsize(1)
  farray(i,j,k)=w1(i,j,k)
End Do; End Do; End Do
Call decomp_2d_fft_3d(farray,ccapz)
spec_sum(:,:,:,3)=spec_sum(:,:,:,3)+Conjg(ccapz(:,:,:))*ccapz(:,:,:)/Real(nx*nz,mytype)
!
Deallocate(ccapz,ccapz2)
!
Return
End Subroutine calc_spectra
!
!******************************************************************************
Subroutine calc_stat_0d  
!******************************************************************************
Implicit None
Integer :: i,j,k,ib,j2,j1d,k1d
Real(mytype) :: esum1
Real(mytype), dimension(2) :: buffr
!
! Update umean's boundary
buffr(:)=umean(xsize(2)-1:xsize(2))
  Call MPI_Sendrecv_replace(buffr,2, &
	real_type,right_row,901,left_row,901,MPI_COMM_CART,status_info,ierror)
umean(-1:0)=buffr(:)
If (xmin(2) == 1) umean(-1:0)=0.0e0
If (xmax(2) == 1) umean(xsize(2)+1:xsize(2)+2)=0.0e0
!
! Update umean's boundary (upper)
buffr(:)=umean(1:2)
  Call MPI_Sendrecv_replace(buffr,2, &
	real_type,left_row,901,right_row,901,MPI_COMM_CART,status_info,ierror)
umean(xsize(2)+1:xsize(2)+2)=buffr(:)
If (xmax(2) == 1) umean(xsize(2)+1:xsize(2)+2)=0.0e0
!
ib=1
! u_tau calculation
  Do j = 1,xsize(2)+1
    j1d=xstart(2)-1+j
    tauw(j) = (gcdy_f(j1d,1,1,1)*umean(j+1) + gcdy_f(j1d,2,1,1)*umean(j)&
 	     +gcdy_f(j1d,3,1,1)*umean(j-1) + gcdy_f(j1d,4,1,1)*umean(j-2)) *nu
  End Do
!
!
  utau1 = Sqrt(Abs(tauw(1)))
  utau2 = Sqrt(Abs(tauw(xsize(2)+1)))
  Call MPI_Bcast(utau1,1,real_type,0,MPI_COMM_CART,ierror)
  Call MPI_Bcast(utau2,1,real_type,nproc-1,MPI_COMM_CART,ierror)
  utau = 0.5e0*(utau1+utau2)
!
!	dudy calculation
  Do j = 1,xsize(2)+1
    j1d=xstart(2)-1+j
     dudy1(j) = (gcdy_f(j1d,1,1,ib)*umean(j+1)+gcdy_f(j1d,2,1,ib)*umean(j)+     &
     		gcdy_f(j1d,3,1,ib)*umean(j-1)+gcdy_f(j1d,4,1,ib)*umean(j-2))            
  End Do
!
   esum = 0.0e0; esum1 = 0.0e0
  Do k = 1,xsize(3)
    k1d=xstart(3)-1+k
  Do j = 1,xsize(2)
    j1d=xstart(2)-1+j
  Do i = 1,xsize(1)
   esum1 = esum1+(((u1(i,j,k)+u1(i+1,j,k))*0.5e0-umean(j))*((u1(i,j,k)+u1(i+1,j,k))*0.5e0-umean(j))                        &
                +((v1(i,j,k)+v1(i,j+1,k))*0.5e0-vmean(j))*((v1(i,j,k)+v1(i,j+1,k))*0.5e0-vmean(j))                        &
                +(w1(i,j,k)+w1(i,j,k+1))*0.5e0*(w1(i,j,k)+w1(i,j,k+1))*0.5e0)*dx(i)*dy(j1d)*dz(k1d)
  End Do; End Do; End Do
!
   Call MPI_Allreduce(esum1,esum,1,real_type,MPI_SUM,MPI_COMM_CART,ierror)
   esum = esum*0.5e0/length/height/width
!
! Time averaging
!
If (itmavg == 1 .and. (Mod(ii,ntmavg) == 0 .Or. lstop == 1 .Or. ii .le. 5)) Then
  samples = samples+1
  If (nrank==0) Write(21,*) 'samples=',samples
!
  Do j = 1,xsize(2)
    taulmo(j) = taulmo(j)+tauw(j)
  End Do
  taulmo(xsize(2)+1) = taulmo(xsize(2)+1)+tauw(xsize(2)+1)
!
End If

!
Return
End Subroutine calc_stat_0d
!
!******************************************************************************
Subroutine stat_parameters
!******************************************************************************
! 	This subroutine writes the statistics parameters  needed for post processing code
Implicit None
Integer :: i,j,k
!
Open(15, file= trim(folder)//'stat-deck.par')
Write (15,*) 'nx', nx
Write (15,*) 'ny', ny
Write (15,*) 'nz', nz
Write (15,*) 'nstep', nstep
Write (15,*) 'itmavg', itmavg                       ! time averaging option
Write (15,*) 'ntmavg', ntmavg                       ! time averaging interval
Write (15,*) 'samples', samples                     ! time samples
Write (15,*) 'lpstat_1d',lpstat_1d                  ! time average 1D statistics data (plane aveaged)
Write (15,*) 'lpstat_2d',lpstat_2d                  ! time average 2D statistics data (line averaged)
Write (15,*) 'lpstat_x',lpstat_x                    ! spanwise averaged 2D statistics
Write (15,*) 'lpstat_z',lpstat_z                    ! streamwise averaged 2D statistics
Write (15,*) 'ibin_rite', ibin_rite                 ! restrt file 1: binary  0:Ascii
Write (15,*) 'lp_ins_1d', lp_ins_1d                 ! Instaneous 1d data option
Write (15,*) 'lp_ins_2d', lp_ins_2d                 ! Instaneous 2d data option
Write (15,*) 'lp_ins_3d', lp_ins_3d                 ! Instaneous 3d data option
Write (15,*) 'lp_snap_1d', lp_snap_1d               ! Velcoity 1d snap option
Write (15,*) 'lp_snap_2d', lp_snap_2d               ! Velocity 2d snap option
Write (15,*) 'lp_snap_3d', lp_snap_3d               ! Velcoity 3d snap option
Write (15,*) 'lp_snap_x', lp_snap_x                 ! Velcoity 1d snap option in x
Write (15,*) 'lp_snap_y', lp_snap_y                 ! Velcoity 1d snap option in y
Write (15,*) 'lp_snap_z', lp_snap_z                 ! Velcoity 1d snap option in z
Write (15,*) 'lp_snap_xy', lp_snap_xy               ! Velcoity 2d snap option in xy
Write (15,*) 'lp_snap_xz', lp_snap_xz               ! Velcoity 2d snap option in xz
Write (15,*) 'lp_snap_yz', lp_snap_yz               ! Velcoity 2d snap option in yz
Write (15,*) 'nsmpl_ins_1d', nsmpl_ins_1d           ! Instataneous 1d statistics data samples
Write (15,*) 'nsmpl_ins_2d', nsmpl_ins_2d           ! Instataneous 2d statistics data samples
Write (15,*) 'nsmpl_ins_3d', nsmpl_ins_3d           ! Instataneous 3d statistics data samples
Write (15,*) 'nsmpl_snap_1d', nsmpl_snap_1d         ! Instataneous 1d snap data samples
Write (15,*) 'nsmpl_snap_2d', nsmpl_snap_2d         ! Instataneous 2d snap data samples
Write (15,*) 'nsmpl_snap_3d', nsmpl_snap_3d         ! Instataneous 3d snap data samples
Write (15,*) 'cony', cony                           ! grid non uniformity in y
Write (15,*) 'length', length                      
Write (15,*) 'width', width
Write (15,*) 'nu', nu
Write (15,*) 'total_step', ii                       ! total time step of the simulation
!
close(15)
!
End Subroutine stat_parameters
!
!******************************************************************************
Subroutine vis_particle   !This subroutine is for particle trace
!******************************************************************************
! 	This subroutine does particle tracing
Implicit None
Integer :: i,j,n,m,neib_rank,ncol,nrow,nbuf
Real(mytype) :: dety,detz,xpos,ypos,zpos
Integer, Dimension(p_row*p_col) :: nmov
Real(mytype), Dimension(1:3,npart_max/100,p_row*p_col) :: psendbuf,precvbuf
!
nbuf=npart_max/100
detz=(endz-startz)/Real(npart_nozz+1,kind=mytype)
dety=(endy-starty)/Real(npart_nozz+1,kind=mytype)
Do n=1,npart_slc
  nmov(:)=0
  psendbuf(:,:,:)=-1.0e0;precvbuf(:,:,:)=-1.0e0
  Do i=1,npart_max
    If(inpart(i,n) == 1) Then
      Call update_pos(particle(1:3,i,n),xpos,ypos,zpos)
      If(xpos.Gt.length.Or.xpos.Lt.0.0e0) Then
        inpart(i,n)=0
	particle(1,i,n)=-1.0e0
	particle(2,i,n)=-1.0e0
	particle(3,i,n)=-1.0e0
      Else If(ypos.Ge.y(xstart(2)).And.ypos.Le.y(xend(2)+1).And.zpos.Ge.z(xstart(3)).And. &
	zpos.Le.z(xend(3)+1)) Then
	particle(1,i,n)=xpos
	particle(2,i,n)=ypos
	particle(3,i,n)=zpos
      Else
	ncol=Int(zpos/(width/Real(p_col,kind=mytype)))
	If(ypos.Gt.y(xend(2)+1)) Then  ! non-uniform grid
	  nrow=col_rank+1
	Else If(ypos.Lt.y(xstart(2))) Then
	  nrow=col_rank-1
	Else
	  nrow=col_rank
	End If
	neib_rank=nrow*p_col+ncol
	nmov(neib_rank+1)=nmov(neib_rank+1)+1
	psendbuf(1,nmov(neib_rank+1),neib_rank+1)=xpos
	psendbuf(2,nmov(neib_rank+1),neib_rank+1)=ypos
	psendbuf(3,nmov(neib_rank+1),neib_rank+1)=zpos
!
	inpart(i,n)=0
	particle(1,i,n)=-1.0e0
	particle(2,i,n)=-1.0e0
	particle(3,i,n)=-1.0e0
      End If
    End If
  End Do
  ! exchange particles between neighbours
  Do m=1,p_row*p_col
    Call MPI_Gather(psendbuf(1:3,1:nbuf,m),3*nbuf,real_type,precvbuf(1:3,1:nbuf,1:p_row*p_col),&
	3*nbuf,real_type,m-1,MPI_COMM_WORLD,ierror)
  End Do
  Do m=1,p_row*p_col
    If(m.Ne.nrank+1) Then
      Do j=1,nbuf
	If(precvbuf(2,j,m).Ge.0.0e0) Then
	  Do i=1,npart_max
	    If(inpart(i,n) == 0) Then
	      particle(1,i,n)=precvbuf(1,j,m)
	      particle(2,i,n)=precvbuf(2,j,m)
	      particle(3,i,n)=precvbuf(3,j,m)
	      inpart(i,n)=1
	      Exit
	    End If
	  End Do
	Else
	  Exit
	End If
      End Do
    End If
  End Do
End Do
If(Mod(ii,np_snap_3d) == 0) Then
  Do n=1,npart_slc
    Do m=1,npart_nozz
      ypos=Real(m,kind=mytype)*dety+starty
      zpos=Real(m,kind=mytype)*detz+startz
      If(ypos.Ge.y(xstart(2)).And.ypos.Le.y(xend(2)+1).And.zpos.Ge.z(xstart(3)).And. &
	        zpos.Le.z(xend(3)+1)) Then  ! find nozzle position
        Do i=1,npart_max
	  If(inpart(i,n) == 0) Then
	    particle(1,i,n)=0.0e0
	    particle(2,i,n)=ypos
	    particle(3,i,n)=zpos
	    Call update_pos(particle(1:3,i,n),xpos,ypos,zpos)
	    particle(1,i,n)=xpos
	    particle(2,i,n)=ypos
	    particle(3,i,n)=zpos
	    inpart(i,n)=1
	    Exit
	  End If
        End Do
        If(i.Ge.npart_max) Then
	  Write(*,*) 'Error: too many particles in the field!'
	  ipart=0
	  Call MPI_Bcast(ipart,1,MPI_integer,nrank,MPI_COMM_WORLD,ierror)
        End If
      End If
    End Do
  End Do
  ! Write all the particles into files
  Do n=1,npart_slc
    m=0
    Write(310+n,'(A,F10.5,A)') 'zone t= "',time,'"'
    Do i=1,npart_max
      If(inpart(i,n) == 1) Then
        Write(310+n,'(3F20.10)') particle(1:3,i,n)
        m=m+1
      End If
    End Do
    If(m.Eq.0) Write(310+n,'(3F20.10)') -1.0e0,-1.0e0,-1.0e0
  End Do
End If
!
End Subroutine vis_particle
!
!******************************************************************************
Subroutine update_pos(xyzpart,xpos,ypos,zpos)
!******************************************************************************
! 	This subroutine updates the position of the particles
Implicit None
Real(mytype), Dimension(1:3), Intent(in) :: xyzpart
Real(mytype), Intent(out) :: xpos,ypos,zpos
Real(mytype), Dimension(1:3,1:2) :: fac
Integer, Dimension(1:3,1:2) :: i2
Real(mytype) :: vtiny,uu,vv,ww
Integer :: i,j,k,i1d,j1d,k1d,io1,io2

vtiny=1.0e-3
! trilinear interpolate u velocity
i2(:,:)=0;fac(:,:)=0.0e0
! x direction
Do i=0,xsize(1)+1
  i1d=i+xstart(1)-1
  If(xyzpart(1).Ge.(x(i1d)-vtiny).And.xyzpart(1).Le.(x(i1d)+vtiny)) Then
    i2(1,1)=i; i2(1,2)=i
    fac(1,1)=1.0e0; fac(1,2)=0.0e0
  Else If(xyzpart(1).Ge.x(i1d).And.xyzpart(1).Le.x(i1d+1)) Then
    i2(1,1)=i; i2(1,2)=i+1
    fac(1,1)=(x(i1d+1)-xyzpart(1))/(x(i1d+1)-x(i1d))
    fac(1,2)=1.0e0-fac(1,1)
  End If
End Do
! y direction
Do j=0,xsize(2)+1
  j1d=j+xstart(2)-1
  If(xyzpart(2).Ge.(yc(j1d)-vtiny).And.xyzpart(2).Le.(yc(j1d)+vtiny)) Then
    i2(2,1)=j; i2(2,2)=j
    fac(2,1)=1.0e0; fac(2,2)=0.0e0
  Else If(xyzpart(2).Ge.yc(j1d).And.xyzpart(2).Le.yc(j1d+1)) Then
    i2(2,1)=j; i2(2,2)=j+1
    fac(2,1)=(yc(j1d+1)-xyzpart(2))/(yc(j1d+1)-yc(j1d))
    fac(2,2)=1.0e0-fac(2,1)
  End If
End Do
! z direction
Do k=0,xsize(3)+1
  k1d=k+xstart(3)-1
  If(xyzpart(3).Ge.(zc(k1d)-vtiny).And.xyzpart(3).Le.(zc(k1d)+vtiny)) Then
    i2(3,1)=k; i2(3,2)=k
    fac(3,1)=1.0e0; fac(3,2)=0.0e0
  Else If(xyzpart(3).Ge.zc(k1d).And.xyzpart(3).Le.zc(k1d+1)) Then
    i2(3,1)=k; i2(3,2)=k+1
    fac(3,1)=(zc(k1d+1)-xyzpart(3))/(zc(k1d+1)-zc(k1d))
    fac(3,2)=1.0e0-fac(3,1)
  End If
End Do
uu=fac(1,1)*fac(2,1)*fac(3,1)*u1(i2(1,1),i2(2,1),i2(3,1))   &
    +fac(1,2)*fac(2,1)*fac(3,1)*u1(i2(1,2),i2(2,1),i2(3,1))   &
    +fac(1,1)*fac(2,2)*fac(3,1)*u1(i2(1,1),i2(2,2),i2(3,1))   &
    +fac(1,2)*fac(2,2)*fac(3,1)*u1(i2(1,2),i2(2,2),i2(3,1))   &
    +fac(1,1)*fac(2,1)*fac(3,2)*u1(i2(1,1),i2(2,1),i2(3,2))   &
    +fac(1,2)*fac(2,1)*fac(3,2)*u1(i2(1,2),i2(2,1),i2(3,2))   &
    +fac(1,1)*fac(2,2)*fac(3,2)*u1(i2(1,1),i2(2,2),i2(3,2))   &
    +fac(1,2)*fac(2,2)*fac(3,2)*u1(i2(1,2),i2(2,2),i2(3,2))
!
! trilinear interpolate v velocity
! x direction
io1=Max(i2(1,1)-3,1); io2=Min(i2(1,2)+3,xsize(1))
Do i=io1,io2
  i1d=i+xstart(1)-1
  If(xyzpart(1).Ge.(xc(i1d)-vtiny).And.xyzpart(1).Le.(xc(i1d)+vtiny)) Then
    i2(1,1)=i; i2(1,2)=i
    fac(1,1)=1.0e0; fac(1,2)=0.0e0
  Else If(xyzpart(1).Ge.xc(i1d).And.xyzpart(1).Le.xc(i1d+1)) Then
    i2(1,1)=i; i2(1,2)=i+1
    fac(1,1)=(xc(i1d+1)-xyzpart(1))/(xc(i1d+1)-xc(i1d))
    fac(1,2)=1.0e0-fac(1,1)
  End If
End Do
! y direction
io1=Max(i2(2,1)-3,1); io2=Min(i2(2,2)+3,xsize(2))
Do j=io1,io2
  j1d=j+xstart(2)-1
  If(xyzpart(2).Ge.(y(j1d)-vtiny).And.xyzpart(2).Le.(y(j1d)+vtiny)) Then
    i2(2,1)=j; i2(2,2)=j
    fac(2,1)=1.0e0; fac(2,2)=0.0e0
  Else If(xyzpart(2).Ge.y(j1d).And.xyzpart(2).Le.y(j1d+1)) Then
    i2(2,1)=j; i2(2,2)=j+1
    fac(2,1)=(y(j1d+1)-xyzpart(2))/(y(j1d+1)-y(j1d))
    fac(2,2)=1.0e0-fac(2,1)
  End If
End Do
! z direction
io1=Max(i2(3,1)-3,1); io2=Min(i2(3,2)+3,xsize(3))
Do k=io1,io2
  k1d=k+xstart(3)-1
  If(xyzpart(3).Ge.(zc(k1d)-vtiny).And.xyzpart(3).Le.(zc(k1d)+vtiny)) Then
    i2(3,1)=k; i2(3,2)=k
    fac(3,1)=1.0e0; fac(3,2)=0.0e0
  Else If(xyzpart(3).Ge.zc(k1d).And.xyzpart(3).Le.zc(k1d+1)) Then
    i2(3,1)=k; i2(3,2)=k+1
    fac(3,1)=(zc(k1d+1)-xyzpart(3))/(zc(k1d+1)-zc(k1d))
    fac(3,2)=1.0e0-fac(3,1)
  End If
End Do
vv=fac(1,1)*fac(2,1)*fac(3,1)*v1(i2(1,1),i2(2,1),i2(3,1))   &
    +fac(1,2)*fac(2,1)*fac(3,1)*v1(i2(1,2),i2(2,1),i2(3,1))   &
    +fac(1,1)*fac(2,2)*fac(3,1)*v1(i2(1,1),i2(2,2),i2(3,1))   &
    +fac(1,2)*fac(2,2)*fac(3,1)*v1(i2(1,2),i2(2,2),i2(3,1))   &
    +fac(1,1)*fac(2,1)*fac(3,2)*v1(i2(1,1),i2(2,1),i2(3,2))   &
    +fac(1,2)*fac(2,1)*fac(3,2)*v1(i2(1,2),i2(2,1),i2(3,2))   &
    +fac(1,1)*fac(2,2)*fac(3,2)*v1(i2(1,1),i2(2,2),i2(3,2))   &
    +fac(1,2)*fac(2,2)*fac(3,2)*v1(i2(1,2),i2(2,2),i2(3,2))
!
! trilinear interpolate w velocity
! x direction
io1=Max(i2(1,1)-3,1); io2=Min(i2(1,2)+3,xsize(1))
Do i=io1,io2
  i1d=i+xstart(1)-1
  If(xyzpart(1).Ge.(xc(i1d)-vtiny).And.xyzpart(1).Le.(xc(i1d)+vtiny)) Then
    i2(1,1)=i; i2(1,2)=i
    fac(1,1)=1.0e0; fac(1,2)=0.0e0
  Else If(xyzpart(1).Ge.xc(i1d).And.xyzpart(1).Le.xc(i1d+1)) Then
    i2(1,1)=i; i2(1,2)=i+1
    fac(1,1)=(xc(i1d+1)-xyzpart(1))/(xc(i1d+1)-xc(i1d))
    fac(1,2)=1.0e0-fac(1,1)
  End If
End Do
! y direction
io1=Max(i2(2,1)-3,1); io2=Min(i2(2,2)+3,xsize(2))
Do j=io1,io2
  j1d=j+xstart(2)-1
  If(xyzpart(2).Ge.(yc(j1d)-vtiny).And.xyzpart(2).Le.(yc(j1d)+vtiny)) Then
    i2(2,1)=j; i2(2,2)=j
    fac(2,1)=1.0e0; fac(2,2)=0.0e0
  Else If(xyzpart(2).Ge.yc(j1d).And.xyzpart(2).Le.yc(j1d+1)) Then
    i2(2,1)=j; i2(2,2)=j+1
    fac(2,1)=(yc(j1d+1)-xyzpart(2))/(yc(j1d+1)-yc(j1d))
    fac(2,2)=1.0e0-fac(2,1)
  End If
End Do
! z direction
io1=Max(i2(3,1)-3,1); io2=Min(i2(3,2)+3,xsize(3))
Do k=io1,io2
  k1d=k+xstart(3)-1
  If(xyzpart(3).Ge.(z(k1d)-vtiny).And.xyzpart(3).Le.(z(k1d)+vtiny)) Then
    i2(3,1)=k; i2(3,2)=k
    fac(3,1)=1.0e0; fac(3,2)=0.0e0
  Else If(xyzpart(3).Ge.z(k1d).And.xyzpart(3).Le.z(k1d+1)) Then
    i2(3,1)=k; i2(3,2)=k+1
    fac(3,1)=(z(k1d+1)-xyzpart(3))/(z(k1d+1)-z(k1d))
    fac(3,2)=1.0e0-fac(3,1)
  End If
End Do
ww=fac(1,1)*fac(2,1)*fac(3,1)*w1(i2(1,1),i2(2,1),i2(3,1))   &
    +fac(1,2)*fac(2,1)*fac(3,1)*w1(i2(1,2),i2(2,1),i2(3,1))   &
    +fac(1,1)*fac(2,2)*fac(3,1)*w1(i2(1,1),i2(2,2),i2(3,1))   &
    +fac(1,2)*fac(2,2)*fac(3,1)*w1(i2(1,2),i2(2,2),i2(3,1))   &
    +fac(1,1)*fac(2,1)*fac(3,2)*w1(i2(1,1),i2(2,1),i2(3,2))   &
    +fac(1,2)*fac(2,1)*fac(3,2)*w1(i2(1,2),i2(2,1),i2(3,2))   &
    +fac(1,1)*fac(2,2)*fac(3,2)*w1(i2(1,1),i2(2,2),i2(3,2))   &
    +fac(1,2)*fac(2,2)*fac(3,2)*w1(i2(1,2),i2(2,2),i2(3,2))
!
xpos=xyzpart(1)+uu*dt_fixed
ypos=xyzpart(2)+vv*dt_fixed
zpos=xyzpart(3)+ww*dt_fixed
If(zpos.Gt.width) Then
  zpos=zpos-width
Else If(zpos.Lt.0.0e0) Then
  zpos=zpos+width
End If
ypos=Min(ypos,height); ypos=Max(ypos,0.0e0)
!
End Subroutine update_pos
!******************************************************************************
Subroutine calc3d      !This subroutine calculates the lambda2 and u'
!******************************************************************************
If(ilamq == 1) Then
  Call lambda2_calc
  If(col_rank.Le.p_row/2-1) Call write_lamb2_par
End If
If(iwfield == 1.And.col_rank.Le.p_row/2) Call write_3dfield_par
If(iflucs == 1.And.col_rank.Le.p_row/2) Call write_flucs_par 
!
End Subroutine calc3d
!
!******************************************************************************
Subroutine lambda2_calc      !This subroutine calculates lambda2
!******************************************************************************
!       This subroutine calculates the lambda2 values
Implicit None
Integer :: i,j,k,n,m,j1d,k1d,ib
Real(mytype) :: a1,a2,a3,l1,l2,l3
Real(mytype), dimension(3,3) :: s,o,so,vel_grad
Real(mytype), dimension(nx) :: uc,ucpx,ucmx,ucpy,ucmy,ucpz,ucmz, &
                               vc,vcpx,vcmx,vcpy,vcmy,vcpz,vcmz, &
                               wc,wcpx,wcmx,wcpy,wcmy,wcpz,wcmz
ib =1

Do j = 0,xsize(2)+1
  j1d=xstart(2)-1+j
  Do k = 1,xsize(3)
   k1d=xstart(3)-1+k
   Do i = 1,xsize(1)
      uc(i)    = 0.5e0*(u1(i,j,k)+u1(i+1,j,k))
      ucmx(i)  = 0.5e0*(u1(i-1,j,k)+u1(i,j,k))
      ucpx(i)  = 0.5e0*(u1(i+1,j,k)+u1(i+2,j,k))
      ucmy(i)  = 0.5e0*(u1(i,j-1,k)+u1(i+1,j-1,k))
      ucpy(i)  = 0.5e0*(u1(i,j+1,k)+u1(i+1,j+1,k))
      ucmz(i)  = 0.5e0*(u1(i,j,k-1)+u1(i+1,j,k-1))
      ucpz(i)  = 0.5e0*(u1(i,j,k+1)+u1(i+1,j,k+1))

      vc(i)    = 0.5e0*(v1(i,j,k)+v1(i,j+1,k))
      vcmx(i)  = 0.5e0*(v1(i-1,j,k)+v1(i-1,j+1,k))
      vcpx(i)  = 0.5e0*(v1(i+1,j,k)+v1(i+1,j+1,k))
      vcmy(i)  = 0.5e0*(v1(i,j-1,k)+v1(i,j,k))
      vcpy(i)  = 0.5e0*(v1(i,j+1,k)+v1(i,j+2,k))
      vcmz(i)  = 0.5e0*(v1(i,j,k-1)+v1(i,j+1,k-1))
      vcpz(i)  = 0.5e0*(v1(i,j,k+1)+v1(i,j+1,k+1))

      wc(i)     = 0.5e0*(w1(i,j,k)+w1(i,j,k+1))
      wcmx(i)   = 0.5e0*(w1(i-1,j,k)+w1(i-1,j,k+1))
      wcpx(i)   = 0.5e0*(w1(i+1,j,k)+w1(i+1,j,k+1))
      wcmy(i)   = 0.5e0*(w1(i,j-1,k)+w1(i,j-1,k+1))
      wcpy(i)   = 0.5e0*(w1(i,j+1,k)+w1(i,j+1,k+1))
      wcmz(i)   = 0.5e0*(w1(i,j,k-1)+w1(i,j,k))
      wcpz(i)   = 0.5e0*(w1(i,j,k+1)+w1(i,j,k+2))
!
      vel_grad(1,1) = gcdx_c(i,2,1,ib)*ucpx(i)+gcdx_c(i,3,1,ib)*uc(i) &
            +gcdx_c(i,4,1,ib)*ucmx(i)
      vel_grad(1,2) = gcdy_c(j1d,2,2,ib)*ucpy(i)+gcdy_c(j1d,3,2,ib)*uc(i) &
            +gcdy_c(j1d,4,2,ib)*ucmy(i)
      vel_grad(1,3) = gcdz_c(k1d,2,1,ib)*ucpz(i)+gcdz_c(k1d,3,1,ib)*uc(i) &
            +gcdz_c(k1d,4,1,ib)*ucmz(i)

      vel_grad(2,1) = gcdx_c(i,2,1,ib)*vcpx(i)+gcdx_c(i,3,1,ib)*vc(i) &
            +gcdx_c(i,4,1,ib)*vcmx(i)
      vel_grad(2,2) = gcdy_c(j1d,2,2,ib)*vcpy(i)+gcdy_c(j1d,3,2,ib)*vc(i) &
            +gcdy_c(j1d,4,2,ib)*vcmy(i)
      vel_grad(2,3) = gcdz_c(k1d,2,1,ib)*vcpz(i)+gcdz_c(k1d,3,1,ib)*vc(i) &
            +gcdz_c(k1d,4,1,ib)*vcmz(i)

      vel_grad(3,1) = gcdx_c(i,2,1,ib)*wcpx(i)+gcdx_c(i,3,1,ib)*wc(i) &
            +gcdx_c(i,4,1,ib)*wcmx(i)
      vel_grad(3,2) = gcdy_c(j1d,2,2,ib)*wcpy(i)+gcdy_c(j1d,3,2,ib)*wc(i) &
            +gcdy_c(j1d,4,2,ib)*wcmy(i)
      vel_grad(3,3) = gcdz_c(k1d,2,1,ib)*wcpz(i)+gcdz_c(k1d,3,1,ib)*wc(i) &
            +gcdz_c(k1d,4,1,ib)*wcmz(i)
      Do m=1,3; Do n=1,3
        s(m,n) = 0.5E0* (vel_grad(m,n) + vel_grad(n,m) )
        o(m,n) = 0.5E0* (vel_grad(m,n) - vel_grad(n,m) )        
      End Do; End Do
      Do m=1,3; Do n=1,3
        so(m,n) =   s(m,1)*s(1,n) + s(m,2)*s(2,n) + s(m,3)*s(3,n)  &
                  + o(m,1)*o(1,n) + o(m,2)*o(2,n) + o(m,3)*o(3,n)
      End Do; End Do   
      a1 = -(so(1,1) + so(2,2) + so(3,3))
      a2 = -(so(1,2)*so(2,1) + so(1,3)*so(3,1) + so(2,3)*so(3,2) &
         -so(1,1)*so(2,2) - so(1,1)*so(3,3) - so(2,2)*so(3,3))
      a3 = -(so(1,1)*so(2,2)*so(3,3) + so(1,3)*so(2,1)*so(3,2) + so(1,2)*so(2,3)*so(3,1) &
         -so(1,1)*so(2,3)*so(3,2) - so(1,2)*so(2,1)*so(3,3) - so(1,3)*so(2,2)*so(3,1) )
      Call cubic(a1,a2,a3,l1,l2,l3)
      l2_array(i,j,k)=l2
   End Do
  End Do
End Do
Call update_halo(l2_array,1)

Return
End Subroutine lambda2_calc
!
!******************************************************************************
Subroutine cubic (b,c,d,l1,l2,l3)      !This subroutine solve a*x^3+b*x^2+c*x+d=0
!******************************************************************************
! This subroutine calculates the eigenvalues from the characteristic polynomial
Implicit None
Real(mytype), Parameter ::   pi = 2.0e0*ASIN(1.)
Complex  :: xx(3)
Real(mytype) ::  L(3)
Integer :: nroot, I=0,J=0,K=0
Real(mytype), Intent(in) :: b,c,d
Real(mytype) :: a, phi, DD, pp,qq,temp1, temp2, y1,y2,y3,u1,v1,y2i,y2r
Real(mytype), Intent(out) :: l1,l2,l3
     
y1=0;y2=0;y3=0;y2i=0;y2r=0
nroot = 2
DD = c*c-4.0E0*b*d
If(DD .Ge. 0.0E0) Then
  xx(1) = Cmplx((-c+Sqrt(DD))/2.0E0/b, 0.0E0)
  xx(2) = Cmplx((-c-sqrt(DD))/2.0E0/b, 0.0E0)
Else
  xx(1) = Cmplx(-c/2.0E0/b, +Sqrt(-DD)/2.0E0/b)
  xx(2) = Cmplx(-c/2.0E0/b, -Sqrt(-DD)/2.0E0/b)
End If
nroot = 3
a = 1.0E0
pp  = c/a - b*b/a/a/3.0E0
qq  = (2.0E0*b*b*b/a/a/a - 9.0E0*b*c/a/a + 27.0E0*d/a) / 27.0E0
DD = pp*pp*pp/27.0E0 + qq*qq/4.0E0
If(DD .Lt. 0.0E0)then
  phi = Acos(-qq/2.0E0/Sqrt(Abs(pp*pp*pp)/27.0E0))
  temp1=2.0E0*Sqrt(Abs(pp)/3.0E0)
  y1 =  temp1*Cos(phi/3.0E0)
  y2 = -temp1*Cos((phi+pi)/3.0E0)
  y3 = -temp1*Cos((phi-pi)/3.0E0)
Else
  temp1 = -qq/2.0E0 + Sqrt(DD)
  temp2 = -qq/2.0E0 - Sqrt(DD)
  u1 = Abs(temp1)**(1.0E0/3.0E0)
  v1 = Abs(temp2)**(1.0E0/3.0E0)
  If(temp1 .Lt. 0.0E0) u1=-u1
  If(temp2 .Lt. 0.0E0) v1=-v1
  y1  = u1 + v1
  y2r = -(u1+v1)/2.0E0
  y2i =  (u1-v1)*Sqrt(3.0E0)/2.0E0
End If
temp1 = b/a/3.0E0
y1 = y1-temp1
y2 = y2-temp1
y3 = y3-temp1
y2r=y2r-temp1
If(DD .Lt. 0.0E0) Then
  xx(1) = Cmplx( y1,  0.0E0)
  L(1) = y1
  xx(2) = Cmplx( y2,  0.0E0)
  L(2) = y2
  xx(3) = Cmplx( y3,  0.0E0)
  L(3) = y3
Else if(DD == 0.) Then
  xx(1) = Cmplx( y1,  0.0E0)
  L(1) = y1
  xx(2) = Cmplx(y2r,  0.0E0)
  L(2) = y2r
  xx(3) = Cmplx(y2r,  0.0E0)
  L(3) = y2r
Else
  xx(1) = Cmplx( y1,  0.0E0)
  L(1) = y1
  xx(2) = Cmplx(y2r, y2i)
  L(2) = y2r
  xx(3) = Cmplx(y2r,-y2i)
  L(3) = y2r
End If
If ((L(1) .Le. L(2)).And.(L(2) .Le. L(3))) Then
  l1=L(3); l2=L(2); l3=L(1)
Else If ((L(3) .Le. L(2)).And.(L(2) .Le. L(1))) Then
  l1=L(1); l2=L(2); l3=L(3)
Else If ((L(2) .Le. L(1)).And.(L(1) .Le. L(3))) Then
  l1=L(3); l2=L(1); l3=L(2)
Else If ((L(3) .Le. L(1)).And.(L(1) .Le. L(2))) Then
  l1=L(2); l2=L(1); l3=L(3)
Else If ((L(1) .Le. L(3)).And.(L(3) .Le. L(2))) Then
  l1=L(2); l2=L(3); l3=L(1)
Else If ((L(2) .Le. L(3)).And.(L(3) .Le. L(1))) Then
  l1=L(1); l2=L(3); l3=L(2)
End If
Return
End Subroutine cubic
!
End Module update_field
