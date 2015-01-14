Module setup
Use shared_data
Use flow_field
Use mass_calc
Use boundary_conditions
Use cal_coefficient
Use implicit_solver
Use decomp_2d
Use decomp_2d_fft
Use randgen
!
Implicit None
!
Contains
!
!******************************************************************************
Subroutine initia
!******************************************************************************
Implicit None
Integer :: i,j,k
!
Call MPI_INIT(ierror)
Call default_values
Call input_parameters
Call array_allocation
Call open_files
If (idomain == 1) Call domain_size
Call initialization
Call grid
Call buffer
Call fftw_setup
Call implicit_setup
If (restrt == 0) Call initial_flow
Call initia_bc_update
If (restrt == 1) Call read_restart
If (irandm == 1) Call ran_flow
Call inlet_initial
If (itmavg == 1 .and. initma == 1) Call turb_stat_init
Call initial_out
Call compute_fact
Call comp_pres_coff
Call bc_update
!
Return
End Subroutine initia
!
!*****************************************************************************
Subroutine default_values  
!*****************************************************************************
!	This procedure Reads the input parameters form 'param.ter' and 'jet-par' file
Implicit None
!
!       Grid
!
nx = 64;ny = 129;nz = 64                
!
!	Main simulation parameters
!
ichann = 1;ibcjet = 0;ibmpfg = 0;nblock = 0;re = 2800.0e0       
retau = 180.0e0;restrt = 1;nstep = 100;cony = 1.95e0;irandm = 0;roh = 1.0e0;ngr = 2
length = 7.0e0;width = 3.5e0;height = 2.0e0;idomain = 0;imass = 0
ifftw = 1; ir_p_old = 0;ic_varray_rd = 0;ic_varray_rite = 0
ifgrz_2pt = 1;igr_2pt = 0; dt = 0.01e0 ; ibcjet = 0; if_dpdxin = 0
time_stop = 0.83e5; restrt_old=0
!
!       Boundary conditons parameters 
!
iconbc = 0;jconbc = 1;kconbc = 0;ibwest = 0;ibeast = 0;ibnoth = 1;ibsoth = 1         
kbsoth = 1;kbnoth = 1;if_noslipx = 0; if_noslipy = 1;  if_noslipz = 0; if_mbc = 1
!
!	LES parameters
!
isgsm = 0; nx2 = nx/2; ny2 = ny/2; nz2 = nz/2; nsmles = 0; if_cs_2d_avg = 1
if_backscatter = 1; c_w = 0.6e0; cs_smg = 0.10e0;if_cw_var=0; iread_les = 0;irite_les = 0
!
!	Output file parameters
!
isavvf = 1;itmavg = 1;initma = 1           
np_snap_3d = 2000;nprint = 5;ntmavg = 5;lp_snap_3d = 0
ibin_read = 0;ibin_rite = 0;np_ins_1d = 50;lp_ins_1d = 0
np_snap_2d = 500;lp_snap_2d = 0;np_snap_1d = 500;lp_snap_1d = 0
np_ins_2d = 500;lp_ins_2d = 0;np_ins_3d = 500;lp_ins_3d = 0
xp_snap_1d=nx/2;yp_snap_1d=ny/2;zp_snap_1d=nz/2
xp_snap_2d=nx/2;yp_snap_2d=ny/2;zp_snap_2d=nz/2
ispec=0; ishear=0
nsmpl_ins_1d=0;nsmpl_ins_2d=0;nsmpl_ins_3d=0;nsmpl_snap_1d=0;nsmpl_snap_2d=0;
nsmpl_snap_3d=0;nsmpl_sh_2d=0
!
!	Other simulation parameters
!
pi = 2.0e0*ASIN(1.0e0);pi2 = pi*2.0e0;pi4 = pi*4.0e0;amp = 0.20e0;iseed = 1;del = 1.0e0 
cfl = 0.9e0*sqrt(3.0e0);ii = 0; ncmpnt = 3
large = 1.0e6; small = 1.0e-6; very_small = 1.0e-12
!
!	Grid parameters
!
conx = 0.0e0;cony = 2.04e0;conz = 0.0e0
!
! 	Adjust dt Subroutine parameters
!
idt_ret = 1; idt_ret_un = 0; ifdt = 1;ifdt_fixed = 1;dt_initial = 0.0050e0
dt_fixed = 0.0050e0;deltat_dt = 1.0e0;time_dtcon = 1.0e0
!
!      Flow control parameters
!
wmax=0.0e0; wmax_plus=0.0e0; tperiod=0.0e0; tperiod_plus=0.0e0
omeg=0.0e0; omega_plus=0.0e0; kx=0.0e0; kx_plus=0.0e0; theta_s=0.0e0
kz=0.0e0; kz_plus=0.0e0
amax=1.0e0; istart=0; deltat=2.0e0; jet_type=1; ibcjet=0; itrav=0; itravx=0; itravz=0
iplasma=0; nactu=1; sigma=0.0e0; lambda=0.0e0; cspeed=0.0e0; cspeed_plus=0.0e0
!
icub=0;nphase=1;iphase=0;i2dust=0;iforce=1;nposxz=1
St=0.0e0; St_plus=0.0e0; hjet=0.0e0; hjet_plus=0.0e0
idx=0; idz=0; ndx=1; ndz=1; rise=0.1e0
!
iosci=0;iuni=0;ibid=0
!
!	MPI setup
!
p_row = 2;p_col = 2;
disp_av1d=0_MPI_OFFSET_KIND; disp_av2d=0_MPI_OFFSET_KIND
disp_ins1d=0_MPI_OFFSET_KIND; disp_ins2d=0_MPI_OFFSET_KIND
disp_sn1d=0_MPI_OFFSET_KIND; disp_sn3d=0_MPI_OFFSET_KIND
disp_spec=0_MPI_OFFSET_KIND; disp_sh2d=0_MPI_OFFSET_KIND
!
!	Particle trace variables
!
ipart=0; npart_max=1000; npart_slc=1; y1part=0.0e0; y2part=0.0e0; y3part=0.0e0
!
!	3d vtk file
!
ilamq=0; iwfield=0; iflucs=0
!
Return
End Subroutine default_values
!
!*****************************************************************************
Subroutine input_parameters
!*****************************************************************************
!	This procedure Reads the input parameters form 'param.ter' and 'jet-par' file
Implicit None
 Integer :: i,j,k,ncheck
Character :: dummy*70,file0*20,file1*20,dummy1*2
Character(LEN = 30), Dimension(1000) :: charac
Real(mytype), Dimension(1000) :: var

!
Call getarg(1,foldera)
Call getarg(2,folderb)
Call getarg(3,folderc)	
folder = trim(foldera)//"/"//trim(folderb)//"/"//trim(folderc)//"/"
!folder='./'
!
ncheck = 0
Open(1, File = trim(folder)//'input-deck.par', Status = 'OLD')
Do i = 1,5
Read (1, 300) dummy
Enddo
Read (1, '(A)') unit22     ! restart file for the simulation
Read (1, '(A)') unit23     ! restart file writen at the end of simulation
!
Do i = 1,1000
Read(1,*,End = 10) charac(i),var(i)
ncheck = ncheck+1
End Do
 10    Close(1)
!
Do i  =  1,ncheck
!	Main simulation parametere
If      (charac(i) == 'nx') Then; nx = int(var(i)) 
Else If (charac(i) == 'ny') Then; ny = int(var(i))
Else If (charac(i) == 'nz') Then; nz = int(var(i))
Else If (charac(i) == 'p_row') Then; p_row = int(var(i))
Else If (charac(i) == 'p_col') Then; p_col = int(var(i))
Else If (charac(i) == 'length') Then; length = var(i) 
Else If (charac(i) == 'height') Then; height = var(i) 
Else If (charac(i) == 'width') Then; width = var(i) 
Else If (charac(i) == 'iconbc') Then; iconbc = int(var(i))
Else If (charac(i) == 'jconbc') Then; jconbc = int(var(i))
Else If (charac(i) == 'kconbc') Then; kconbc = int(var(i))
Else If (charac(i) == 'ichann') Then; ichann = int(var(i))
Else If (charac(i) == 'isgsm') Then; isgsm = int(var(i))
Else If (charac(i) == 're') Then; re = var(i)
Else If (charac(i) == 'retau') Then; retau = var(i) 
Else If (charac(i) == 'cony') Then; cony = var(i)
Else If (charac(i) == 'nstep') Then; nstep = int(var(i)) 
Else If (charac(i) == 'restrt') Then; restrt = int(var(i))
Else If (charac(i) == 'ngr') Then; ngr = int(var(i))
Else If (charac(i) == 'roh') Then; roh = var(i)
Else If (charac(i) == 'ibeast') Then; ibeast = int(var(i)) 
Else If (charac(i) == 'ibwest') Then; ibwest = int(var(i))
Else If (charac(i) == 'ibnoth') Then; ibnoth = int(var(i))
Else If (charac(i) == 'ibsoth') Then; ibsoth = int(var(i))
Else If (charac(i) == 'irandm') Then; irandm = int(var(i))
Else If (charac(i) == 'idomain') Then; idomain = int(var(i))
Else If (charac(i) == 'ifftw') Then; ifftw = int(var(i))
Else If (charac(i) == 'igr_2pt') Then; igr_2pt = int(var(i))
Else If (charac(i) == 'ifgrz_2pt') Then; ifgrz_2pt = int(var(i))
Else If (charac(i) == 'imass') Then; imass = int(var(i))
Else If (charac(i) == 'ir_p_old') Then; ir_p_old  = int(var(i))
Else If (charac(i) == 'ic_varray_rd') Then; ic_varray_rd  = int(var(i))
Else If (charac(i) == 'ic_varray_rite') Then; ic_varray_rite  = int(var(i))
Else If (charac(i) == 'if_dpdxin') Then; if_dpdxin  = int(var(i))
Else If (charac(i) == 'dpdx_mean_in') Then; dpdx_mean_in  = var(i)
Else If (charac(i) == 'time_stop') Then; time_stop  = var(i)
Else If (charac(i) == 'restrt_old') Then; restrt_old  = var(i)
!	LES Parameters
Else If (charac(i) == 'nx2') Then; nx2 = int(var(i))
Else If (charac(i) == 'ny2') Then; ny2 = int(var(i))
Else If (charac(i) == 'nz2') Then; nz2 = int(var(i))
Else If (charac(i) == 'isgs_model') Then; isgs_model  = int(var(i))
Else If (charac(i) == 'c_w') Then; c_w  = var(i)
Else If (charac(i) == 'if_cs_2d_avg') Then; if_cs_2d_avg = int(var(i))
Else If (charac(i) == 'if_backscatter') Then; if_backscatter = int(var(i))
Else If (charac(i) == 'cs_smg') Then; cs_smg = var(i)
Else If (charac(i) == 'if_cw_var') Then; if_cw_var = int(var(i))
!	jet-control.par
Else If (charac(i) == 'ibcjet') Then; ibcjet = int(var(i))
Else If (charac(i) == 'jet_type') Then; jet_type = int(var(i))
Else If (charac(i) == 'wmax') Then; wmax=var(i)
Else If (charac(i) == 'wmax_plus') Then; wmax_plus=var(i)
Else If (charac(i) == 'theta_s') Then; theta_s=var(i)
Else If (charac(i) == 'omega') Then; omeg=var(i)
Else If (charac(i) == 'omega_plus') Then; omega_plus=var(i)
Else If (charac(i) == 'kx') Then; kx=var(i)
Else If (charac(i) == 'kx_plus') Then; kx_plus=var(i)
Else If (charac(i) == 'kz') Then; kz=var(i)
Else If (charac(i) == 'kz_plus') Then; kz_plus=var(i)
Else If (charac(i) == 'cspeed') Then; cspeed = var(i)
Else If (charac(i) == 'cspeed_plus') Then; cspeed_plus = var(i)
Else If (charac(i) == 'itravx') Then; itravx=var(i)
Else If (charac(i) == 'itravz') Then; itravz=var(i)
Else If (charac(i) == 'tperiod') Then; tperiod=var(i)
Else If (charac(i) == 'tperiod_plus') Then; tperiod_plus=var(i)
Else If (charac(i) == 'St_plus') Then; St_plus=var(i)
Else If (charac(i) == 'hjet_plus') Then; hjet_plus=var(i)
Else If (charac(i) == 'ampjet') Then; ampjet = var(i)
Else If (charac(i) == 'istart') Then; istart = int(var(i))
Else If (charac(i) == 'deltat') Then; deltat = var(i)
Else If (charac(i) == 'idx') Then; idx = int(var(i))
Else If (charac(i) == 'idz') Then; idz = int(var(i))
Else If (charac(i) == 'ndx') Then; ndx = int(var(i))
Else If (charac(i) == 'ndz') Then; ndz = int(var(i))
Else If (charac(i) == 'rise') Then; rise = var(i)
Else If (charac(i) == 'icub') Then; icub = int(var(i))
Else If (charac(i) == 'iforce') Then; iforce = int(var(i))
Else If (charac(i) == 'iplasma') Then; iplasma = int(var(i))
Else If (charac(i) == 'nactu') Then; nactu = int(var(i))
Else If (charac(i) == 'sigma') Then; sigma = var(i)
Else If (charac(i) == 'lambda') Then; lambda = var(i)
Else If (charac(i) == 'iosci') Then; iosci = int(var(i))
Else If (charac(i) == 'iuni') Then; iuni = int(var(i))
Else If (charac(i) == 'ibid') Then; ibid = int(var(i))
!	Output files parameters
Else If (charac(i) == 'nprint') Then; nprint = int(var(i)) 
Else If (charac(i) == 'isavvf') Then; isavvf = int(var(i))
Else If (charac(i) == 'itmavg') Then; itmavg = int(var(i))
Else If (charac(i) == 'ntmavg') Then; ntmavg = int(var(i))
Else If (charac(i) == 'iread_les') Then; iread_les = int(var(i))
Else If (charac(i) == 'irite_les') Then; irite_les = int(var(i))
Else If (charac(i) == 'initma') Then; initma = int(var(i))
Else If (charac(i) == 'ibin_read') Then; ibin_read = int(var(i))
Else If (charac(i) == 'ibin_rite') Then; ibin_rite = int(var(i))
Else If (charac(i) == 'lpstat_1d') Then; lpstat_1d = int(var(i))
Else If (charac(i) == 'lpstat_2d') Then; lpstat_2d = int(var(i))
Else If (charac(i) == 'lpstat_z') Then; lpstat_z = int(var(i))
Else If (charac(i) == 'lpstat_x') Then; lpstat_x = int(var(i))
Else If (charac(i) == 'lp_snap_1d') Then; lp_snap_1d = int(var(i))
Else If (charac(i) == 'np_snap_1d') Then; np_snap_1d = int(var(i))
Else If (charac(i) == 'lp_snap_2d') Then; lp_snap_2d = int(var(i))
Else If (charac(i) == 'np_snap_2d') Then; np_snap_2d = int(var(i))
Else If (charac(i) == 'lp_snap_3d') Then; lp_snap_3d = int(var(i))
Else If (charac(i) == 'np_snap_3d') Then; np_snap_3d = int(var(i))
Else If (charac(i) == 'np_ins_1d') Then; np_ins_1d = int(var(i))
Else If (charac(i) == 'lp_ins_1d') Then; lp_ins_1d = int(var(i))
Else If (charac(i) == 'np_ins_2d') Then; np_ins_2d = int(var(i))
Else If (charac(i) == 'lp_ins_2d') Then; lp_ins_2d = int(var(i))
Else If (charac(i) == 'np_ins_3d') Then; np_ins_3d = int(var(i))
Else If (charac(i) == 'lp_ins_3d') Then; lp_ins_3d = int(var(i))
Else If (charac(i) == 'lp_snap_x') Then; lp_snap_x = int(var(i))
Else If (charac(i) == 'lp_snap_y') Then; lp_snap_y = int(var(i))
Else If (charac(i) == 'lp_snap_z') Then; lp_snap_z = int(var(i))
Else If (charac(i) == 'lp_snap_xy') Then; lp_snap_xy = int(var(i))
Else If (charac(i) == 'lp_snap_yz') Then; lp_snap_yz = int(var(i))
Else If (charac(i) == 'lp_snap_xz') Then; lp_snap_xz = int(var(i))
Else If (charac(i) == 'xp_snap_1d') Then; xp_snap_1d = int(var(i))
Else If (charac(i) == 'yp_snap_1d') Then; yp_snap_1d = int(var(i))
Else If (charac(i) == 'zp_snap_1d') Then; zp_snap_1d = int(var(i))
Else If (charac(i) == 'xp_snap_2d') Then; xp_snap_2d = int(var(i))
Else If (charac(i) == 'yp_snap_2d') Then; yp_snap_2d = int(var(i))
Else If (charac(i) == 'zp_snap_2d') Then; zp_snap_2d = int(var(i))
Else If (charac(i) == 'nphase') Then; nphase = int(var(i))
Else If (charac(i) == 'nposxz') Then; nposxz = int(var(i))
Else If (charac(i) == 'ispec')      Then; ispec = int(var(i))
!	Adjsutdt subroutine control parameters
Else If (charac(i) == 'dt') Then; dt = var(i)
Else If (charac(i) == 'idt_ret') Then; idt_ret = int(var(i))
Else If (charac(i) == 'idt_ret_un') Then; idt_ret_un = int(var(i))
Else If (charac(i) == 'ifdt_fixed') Then; ifdt_fixed = int(var(i))
Else If (charac(i) == 'dt_fixed') Then; dt_fixed = var(i)
Else If (charac(i) == 'ifdt_variable') Then; ifdt_variable = int(var(i))
Else If (charac(i) == 'dt_initial') Then; dt_initial = var(i)
Else If (charac(i) == 'deltat_dt') Then; deltat_dt = var(i)
Else If (charac(i) == 'time_dtcon') Then; time_dtcon = var(i)
Else If (charac(i) == 'ifdt') Then; ifdt = int(var(i))
!	Particle trace visualization parameter
Else If (charac(i) == 'ipart') Then; ipart = int(var(i))
Else If (charac(i) == 'npart_nozz')Then; npart_nozz = int(var(i))
Else If (charac(i) == 'npart_max') Then; npart_max = int(var(i))
Else If (charac(i) == 'npart_slc') Then; npart_slc = int(var(i))
Else If (charac(i) == 'y1part')    Then; y1part = var(i)
Else If (charac(i) == 'y2part')    Then; y2part = var(i)
Else If (charac(i) == 'y3part')    Then; y3part = var(i)
Else If (charac(i) == 'starty')    Then; starty = var(i)
Else If (charac(i) == 'startz')    Then; startz = var(i)
Else If (charac(i) == 'endy')      Then; endy = var(i)
Else If (charac(i) == 'endz')      Then; endz = var(i)
!	Post3d
Else If (charac(i) == 'ilamq')	   Then; ilamq = int(var(i))
Else If (charac(i) == 'iwfield')   Then; iwfield = int(var(i))
Else If (charac(i) == 'iflucs')    Then; iflucs = int(var(i))
Else If (charac(i) == 'ishear')    Then; ishear = int(var(i))
End If
End Do
  300   Format(A70)
!
  Call decomp_2d_init(nx,ny,nz,p_row,p_col)
If (nrank==0) Write (*,*) "folder = ", folder     ! folder to write output files to
If (nrank == 0) then
  Write(*,*) 'ncheck',ncheck
  Do i=1,ncheck
    Write(*,*) charac(i),var(i)
  End Do
End If
!
nxt = nx+1;nyt = ny+1;nzt = nz+1;nbl = nblock+1;myu = 1.0e0/re;nu = myu/roh
iseed=nrank
If (ibcjet == 1) if_noslipy = 0
If ((iconbc == 0 .Or. if_noslipx == 1) .And. &
(jconbc == 0 .Or. if_noslipy == 1) .And. &
(kconbc == 0 .Or. if_noslipz == 1) .And. ichann == 1) Then
if_mbc = 0 
If(ipart == 1) Then
  Allocate(ypart(1:npart_slc))
  If(npart_slc == 1) Then
    ypart(1)=y1part
  Else If(npart_slc == 2) Then
    ypart(1)=y1part; ypart(2)=y2part
  Else If(npart_slc == 3) Then
    ypart(1)=y1part; ypart(2)=y2part; ypart(3)=y3part
  Else
    If (nrank == 0) Print*, 'npart_slc is greater then 3'
    Stop
  End If
  npart_nozz=npart_nozz/p_col*p_col
  If(nrank == 0) Write(*,*) 'npart_nozz adjusted:', npart_nozz
End If
If (nrank == 0) Print*, 'MBC would not be called in case of MBC solver'
If (nrank == 0) Print*, 'Simulation is periodic or no slip on either side'
EndIf
! read the position of xz plane
If (nposxz .Gt. 1) Then
  Allocate(posxz(1:nposxz))
  Open(Unit=2,File=Trim(folder)//'posxz.par',Status = 'OLD')
  Do i=1,nposxz; Read(2,*) posxz(i); End Do
  If(nrank .Eq. 0) Write(*,*) 'positions for xz plane: ', posxz(1:nposxz) 
  Close(2)
End If
!
Return
End Subroutine input_parameters
!
!*******************************************************************************
Subroutine array_allocation
!*******************************************************************************
!	This procedure performs the dynamic allocation of arrays
Implicit None
!	Main variables arrays
! u,v,w,p for X pencil
Allocate (u1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),v1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),        &
w1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),p1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2))
! u,v,w,p for Y pencil
Allocate (u2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2),v2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2),        &
w2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2),p2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2))
! u,v,w,p for Z pencil
Allocate (u3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2),v3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2),        &
w3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2),p3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2))
!
Allocate(intfuo(1:xsize(1),1:xsize(2),1:xsize(3)),intfvo(1:xsize(1),1:xsize(2),1:xsize(3)),  &
intfwo(1:xsize(1),1:xsize(2),1:xsize(3)))
!
Allocate(farray(1:xsize(1),1:xsize(2),1:xsize(3)),varray(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),  &
yinlet1(-1:xsize(1)+2,-1:xsize(3)+2,3,2,2), yinlet2(-1:ysize(1)+2,-1:ysize(3)+2,3,2,2),  &
yinlet3(-1:zsize(1)+2,-1:zsize(3)+2,3,2,2), xinlet(-1:xsize(2)+2,-1:xsize(3)+2,3,2,2))
!
Allocate(intf1(1:xsize(1),1:xsize(2),1:xsize(3)),intf2(1:ysize(1),1:ysize(2),1:ysize(3)), &
	 intf3(1:zsize(1),1:zsize(2),1:zsize(3)))
!
! uhat,vhat,what for X pencil only
Allocate (uhat1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),  &
vhat1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),what1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2))
!
If (kconbc == 1) Then; Allocate (zinlet(-1:xsize(1)+2,-1:xsize(2)+2,3,2,2));zinlet = 0.0e0; End If
!	Les arrays
Allocate (pitero1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2), rnutxy1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),  &
rnutyz1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2), rnutxz1(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2), &
pitero2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2), rnutxy2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2),  &
rnutyz2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2), rnutxz2(-1:ysize(1)+2,-1:ysize(2)+2,-1:ysize(3)+2), &
pitero3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2), rnutxy3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2),  &
rnutyz3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2), rnutxz3(-1:zsize(1)+2,-1:zsize(2)+2,-1:zsize(3)+2),cs(ny+1))
!
!       Turbulence statistics arrays
Allocate (tke_sum(-1:xsize(2)+2,4,7),wx_sum(-1:xsize(2)+2,2,8),omega_sum(-1:xsize(2)+2,6),        &
taulmo(-1:xsize(2)+2),umean(-1:xsize(2)+2),vmean(-1:xsize(2)+2),stat_sum(-1:xsize(2)+2,4,4),quad_sum(-1:xsize(2)+2,9)) 
!
If (lpstat_2d == 1) Then
!   If(ibcjet == 2 .And. omega_plus .Ne. 0.0e0 .And. (kx_plus .Ne. 0.0e0 .Or. kz_plus .Ne. 0.0e0)) nphase=17
   If(lpstat_z == 1) Allocate (stat2d_sum(-1:xsize(1)+2,-1:xsize(2)+2,4,4,nphase),&
	quad2d_sum(-1:xsize(1)+2,-1:xsize(2)+2,9,nphase), tke2d_sum(-1:xsize(1)+2,-1:xsize(2)+2,4,7,nphase), &
	wx2d_sum(-1:xsize(1)+2,-1:xsize(2)+2,2,8,nphase), omega2d_sum(-1:xsize(1)+2,-1:xsize(2)+2,6,nphase))
   If(lpstat_x == 1) Allocate (stat2d_sum(-1:xsize(3)+2,-1:xsize(2)+2,4,4,nphase),&
	quad2d_sum(-1:xsize(3)+2,-1:xsize(2)+2,9,nphase), tke2d_sum(-1:xsize(3)+2,-1:xsize(2)+2,4,7,nphase),&
	wx2d_sum(-1:xsize(3)+2,-1:xsize(2)+2,2,8,nphase), omega2d_sum(-1:xsize(3)+2,-1:xsize(2)+2,6,nphase))
   If(nphase .Ne. 1) Allocate(ntmavg2(1:nphase))
End If
!
!       Grid arrays
Allocate (x(nx+1),y(ny+1),z(nz+1),dx(-1:nx+2),dy(-1:ny+2),dz(-1:nz+2),    &
xc(0:nx+1),yc(0:ny+1),zc(0:nz+1),dxs(nx+1),dys(ny+1),dzs(nz+1),xplus(nx+1),   &
yplus(ny+1),zplus(nz+1),xcplus(0:nx+1),ycplus(0:ny+1),zcplus(0:nz+1))
!
!       Buffer domain arrays
Allocate (buff(nx+1),buffc(0:nx+1),bre(nx+1),brec(0:nx+1))
!
!	Gradient cofficients arrays
Allocate (gcdx_f(xsize(1)+1,4,ngr,nbl),gcdy_f(ysize(2)+1,4,ngr,nbl),gcdz_f(zsize(3)+1,4,ngr,nbl), &
gcdx_c(0:xsize(1)+1,5,ngr,nbl),gcdy_c(0:ysize(2)+1,5,ngr,nbl),gcdz_c(0:zsize(3)+1,5,ngr,nbl), &
gfdx_c(0:xsize(1)+1,4,ngr,nbl),gfdy_c(ysize(2)+1,4,ngr,nbl),gfdz_c(zsize(3)+1,4,ngr,nbl), &
grdsxx(xsize(1)/2+1,5,2),grdsyy(ysize(2)/2+1,5,2),grdszz(zsize(3)/2+2,5,2),       &
gcdz_f2(zsize(3)+1,4,ngr,nbl),gfdz_c2(zsize(3)+1,4,ngr,nbl))
!  
Allocate (ap(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),ae(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2), &
as(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2),at(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2))
!       Particle trace arrays
If(ipart == 1) Allocate(particle(1:3,1:npart_max,1:npart_slc),inpart(1:npart_max,1:npart_slc))
!
!	Flow control array
If(ibcjet==1.And.jet_type==2) Then
  Allocate(fbody(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2))
  If(idx==1) Allocate(xfringe(-1:xsize(1)+2))
  If(idz==1) Allocate(zfringe(-1:xsize(3)+2))
End If
!	Post3d
If(ilamq == 1) Allocate(l2_array(0:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1))
!
Return
End Subroutine array_allocation
!
!*******************************************************************************
Subroutine initialization
!*******************************************************************************
!	This procedure performs the initializtion of arrays
Implicit None
!
intfuo = 0.0e0; intfvo = 0.0e0; intfwo = 0.0e0
intf1 = 0.0e0; intf2 = 0.0e0; intf3 = 0.0e0; varray = 0.0e0
! 
farray = 0.0e0; xinlet = 0.0e0; yinlet1 = 0.0e0; yinlet2 = 0.0e0; yinlet3 = 0.0e0
If (ichann == 0) then; yinlet1=1.0e0; yinlet2=1.0e0; yinlet3=1.0e0; End If
!
u1 = 0.0e0;v1 = 0.0e0; w1 = 0.0e0; p1 = 0.0e0
u2 = 0.0e0;v2 = 0.0e0; w2 = 0.0e0; p2 = 0.0e0
u3 = 0.0e0;v3 = 0.0e0; w3 = 0.0e0; p3 = 0.0e0
!
intf1 = 0.0e0; intf2 = 0.0e0; intf3 = 0.0e0;
!
uhat1 = 0.0e0; vhat1 = 0.0e0; what1 = 0.0e0
!
pitero1 = 0.0e0;rnutxy1 = 0.0e0;rnutyz1 = 0.0e0;rnutxz1 = 0.0e0;
pitero2 = 0.0e0;rnutxy2 = 0.0e0;rnutyz2 = 0.0e0;rnutxz2 = 0.0e0;
pitero3 = 0.0e0;rnutxy3 = 0.0e0;rnutyz3 = 0.0e0;rnutxz3 = 0.0e0;
!
stat_sum=0.0e0;quad_sum=0.0e0;tke_sum=0.0e0;omega_sum=0.0e0;wx_sum=0.0e0
!
If(lpstat_2d == 1) Then
  stat2d_sum=0.0e0;quad2d_sum=0.0e0;tke2d_sum=0.0e0;omega2d_sum=0.0e0;wx2d_sum=0.0e0
  If(nphase .Ne. 1) ntmavg2=0
End If
!
If(ipart == 1) Then
   particle(1:3,:,:)=0.0e0; npart_in=0; inpart(:,:)=0
End If
!
Return
End Subroutine initialization
!
!*******************************************************************************
Subroutine open_files
!*******************************************************************************
!	This procedure Open files to be written outputs from different procedures
        Implicit None
	Character*60 :: unitmp
!
!	UNIT 20-40   0 Dimensional time instentaneous
!	UNIT 40-70   1 Dimensional time average     
!	UNIT 70-100  1 Dimensional instentaneous
!	UNIT 211-300 Binary files
!	UNIT 311-~   3D instentaneous 
!
!       0 DimensionAL TIME INSTENTANEOUS
unit21 = trim(folder)//'logfile1h6.dat'                                               
unit30 = trim(folder)//'STAT_0D-1.dat'                                                 
unit31 = trim(folder)//'STAT_0D-2.dat'
!
Open(Unit = 21, File = unit21, Status = 'Unknown')
Open(Unit = 30, File = unit30, Status = 'Unknown')
Open(Unit = 31, File = unit31, Status = 'Unknown')
!
If(ipart == 1) Then 
  Write(unitmp,*) nrank
  If(npart_slc.Ge.1)  Open(Unit = 311, file= trim(folder)//'particle_r'//Trim(Adjustl(unitmp))//'.plt')
  If(npart_slc.Ge.2)  Open(Unit = 312, file= trim(folder)//'particle_2_r'//Trim(Adjustl(unitmp))//'.plt')
  If(npart_slc.Ge.3)  Open(Unit = 313, file= trim(folder)//'particle_3_r'//Trim(Adjustl(unitmp))//'.plt')
End If
!
Call Open_binary_files
!
Return
End Subroutine Open_files
!
!*******************************************************************************
Subroutine open_binary_files
!*******************************************************************************.
Implicit none
Integer :: i
Character*30 :: ch_avg_2d
!
!	Time aveaged files
If (lpstat_1d == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_1D_AVG.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_av1d,ierror)
!
If (lpstat_2d == 1) Then
  Allocate(fh_av2d(1:nphase))
  Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_2D_AVG.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_av2d(1),ierror)
  If(nphase .Ne. 1) Then 
    Do i=2,nphase
      If(i .Le. 10) Then
        Write(ch_avg_2d,'(A,I1,A)') 'OUTPUT_2D_AVG_0',i-1,'.dat'
      Else
        Write(ch_avg_2d,'(A,I2,A)') 'OUTPUT_2D_AVG_',i-1,'.dat'
      End If
      Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//ch_avg_2d, &
	    MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_av2d(i),ierror)
    End Do
  End If
End If
!
If (lp_ins_1d == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_1D_INS.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_ins1d,ierror)
!
If (lp_ins_2d == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_2D_INS.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_ins2d,ierror)
!
If (lp_ins_3d == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_3D_INS.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_ins3d,ierror)
!
If (lp_snap_1d == 1)  Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_1D_SNAP.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_sn1d,ierror)
!
If (lp_snap_2d == 1) Then
  Allocate(fh_sn2d(1:nposxz),disp_sn2d(1:nposxz))
  disp_sn2d(1:nposxz)=0_MPI_OFFSET_KIND
  Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_2D_SNAP.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_sn2d(1),ierror)
  If(nposxz .Gt. 1) Then 
    Do i=2,nposxz
      If(i .Lt. 10) Then
        Write(ch_avg_2d,'(A,I1,A)') 'OUTPUT_2D_SNAP_0',i,'.dat'
      Else
        Write(ch_avg_2d,'(A,I2,A)') 'OUTPUT_2D_SNAP_',i,'.dat'
      End If
      Call MPI_FILE_OPEN(MPI_COMM_WORLD,Trim(folder)//ch_avg_2d, &
	    MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_sn2d(i),ierror)
    End Do
  End If
End If
!
If (lp_snap_3d == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_3D_SNAP.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_sn3d,ierror)
!	3d spectra
If(ispec == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'SPECTRA_3D.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_spec,ierror)
!	output wall shear stress
If (ishear == 1) Call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(folder)//'OUTPUT_2D_SHEAR.dat', &
	MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh_sh2d,ierror)
!
Return
End Subroutine Open_binary_files
!
!*******************************************************************************
Subroutine close_files
!*******************************************************************************
!	This procedure closes files
        Implicit None
	Integer :: i
!
!	Close ascii files
!
  Close(21); Close(30); Close(31)
!
  If(ipart == 1) Then
    If(npart_slc.Ge.1)  Close(311)
    If(npart_slc.Ge.2)  Close(312)
    If(npart_slc.Ge.3)  Close(313)
  End If
!
!	Close binary files
!
  If(lpstat_1d == 1) Call MPI_FILE_CLOSE(fh_av1d,ierror)

  Do i=1,nphase
    If(lpstat_2d == 1) Call MPI_FILE_CLOSE(fh_av2d(i),ierror)
  End Do
!
  If(lp_ins_1d == 1) Call MPI_FILE_CLOSE(fh_ins1d,ierror)
  If(lp_ins_2d == 1) Call MPI_FILE_CLOSE(fh_ins2d,ierror)
  If(lp_ins_3d == 1) Call MPI_FILE_CLOSE(fh_ins3d,ierror)
!
  If(lp_snap_1d == 1) Call MPI_FILE_CLOSE(fh_sn1d,ierror)
  Do i=1,nposxz
    If(lp_snap_2d == 1) Call MPI_FILE_CLOSE(fh_sn2d(i),ierror)
  End Do
  If(lp_snap_3d == 1) Call MPI_FILE_CLOSE(fh_sn3d,ierror)
!
  If(ispec == 1) Call MPI_FILE_CLOSE(fh_spec,ierror)
  If(ishear == 1) Call MPI_FILE_CLOSE(fh_sh2d,ierror)
!
Return
End Subroutine close_files
!
!*******************************************************************************
Subroutine domain_size
!*******************************************************************************.
Implicit None
!
!       Channel Flow
!
If (ichann == 1) Then
height=2.0e0; length=2.0e0*pi; width=1.0e0/4.0e0*length
!
!       Boundary Layer Flow
!
Else If (ichann == 0) Then
!
width=1.0e0; length=2.40e1; height=8.0e0
!
Else
!
del=height*0.5e0; length=1.28e1*del; height=2.0e0; width=4.0e0*del
End If
!
Return
End Subroutine domain_size
!
!******************************************************************************
Subroutine grid
!******************************************************************************
!       This procedure defines the computational grid for the new domain
Use shared_data
Implicit None
Integer :: i,j,k
!
 Call grid_channel(x,y,z,xc,yc,zc,nx,ny,nz,cony)
 x=x*length; y=y*height; z=z*width
 xc=xc*length; yc=yc*height; zc=zc*width
 Call grid_dx
!
Return
End Subroutine grid
!
!******************************************************************************
Subroutine grid_channel(x,y,z,xc,yc,zc,nx,ny,nz,cony)
!******************************************************************************
!       This procedure defines the computational grid for the old domain
Implicit None
Integer :: nx,ny,nz,nxt,nyt,nzt,i,j,k
Real(mytype) :: x(nx+1),y(ny+1),z(nz+1),xc(0:nx+1),yc(0:ny+1),zc(0:nz+1),cony
!
nxt=nx+1; nyt=ny+1; nzt=nz+1
Do i = 1,nxt
x(i) = real(i-1,kind=mytype)/real(nx,kind=mytype)
End Do
Do j = 1,nyt
y(j) = real(j-1,kind=mytype)/real(ny,kind=mytype)
End Do
Do k = 1,nzt
z(k) = real(k-1,kind=mytype)/real(nz,kind=mytype)
End Do
!
!-----  Nonuniform grid for Y
!
If (cony > small) Then
Do j = 1,nyt
y(j) = tanh(cony*(2.0e0*(j-1.0e0)/ny-1.0e0))/TANH(cony)
!y(j)=(1.0e0+tanh(cony*(2.0e0*(j-1.0e0)/ny-1.0e0))/tanh(cony))*0.5e0*height
End Do
y=(y+1.0e0)*0.5e0
End If
!
y(1) = 0.0e0
y(nyt) = 1.0e0
!
Do i = 1,nx
xc(i) = 0.5e0*(x(i)+x(i+1))
End Do
xc(0) = -xc(1)
xc(nxt) = x(nxt)+xc(1)

Do j = 1,ny
yc(j) = 0.5e0*(y(j)+y(j+1))
End Do
yc(0) = y(1)
yc(nyt) = y(nyt)

Do k = 1,nz
zc(k) = 0.5e0*(z(k)+z(k+1))
End Do
zc(0) = -zc(1)
zc(nzt) = z(nzt)+zc(1)
Return
End Subroutine grid_channel
!
!******************************************************************************
Subroutine grid_dx
!******************************************************************************
!       This procedure defines the computational grid for the old domain
Use shared_data
Implicit None
Integer :: i,j,k
!
Do i = 1,nx
  dx(i) = x(i+1)-x(i)
End Do
!
dx(-1) = dx(1)
dx(0) = dx(1)
dx(nx+1) = dx(nx)
dx(nx+2) = dx(nx)
!
Do j = 1,ny
  dy(j) = y(j+1)-y(j)
End Do
If (if2d_noy == 1) dy = 1.0e0
!
dy(-1) = dy(1)
dy(0) = dy(1)
!
If (ibsoth == 1) Then
  dy(-1) = 0.0e0
  dy(0) = 0.0e0
End If
!
dy(ny+1) = dy(ny)
dy(ny+2) = dy(ny)
!
If (ibnoth == 1) Then
  dy(ny+1) = 0.0e0
  dy(ny+2) = 0.0e0
End If

Do k = 1,nz
  dz(k) = z(k+1)-z(k)
End Do
!
dz(-1) = dz(1)
dz(0) = dz(1)
dz(nz+1) = dz(nz)
dz(nz+2) = dz(nz)
If (if2d_noz == 1) dz = 1.0e0
!
Do i = 1,nx+1
  dxs(i) = 0.5e0*(dx(i-1)+dx(i))
End Do
!
Do j = 1,ny+1
  dys(j) = 0.5e0*(dy(j-1)+dy(j))
End Do
!
Do k = 1,nz+1
  dzs(k) = 0.5e0*(dz(k-1)+dz(k))
End Do

Do i = 1,nxt
  xplus(i) = x(i)*retau
End Do
!
Do j = 1,nyt
yplus(j) = y(j)*retau
End Do
!
Do k = 1,nzt
  zplus(k) = z(k)*retau
End Do
!
Do i = 0,nxt
  xcplus(i) = xc(i)*retau
End Do
!
Do j = 0,nyt
  ycplus(j) = yc(j)*retau
End Do
!
Do k = 0,nzt
  zcplus(k) = zc(k)*retau
End Do
!
Return
End Subroutine grid_dx
!
!*******************************************************************************
Subroutine buffer
!*******************************************************************************
!	This procedure is to apply buffer domain technique 
Implicit None
Integer ::i,j,k
!
Do i = 1,nx+1
buff(i) = 1.0e0
bre(i) = 1.0e0
buffc(i) = buff(i)
brec(i) = bre(i)
End Do
!
buffc(0) = 1.0e0
brec(0) = 1.0e0
!
Return
End Subroutine buffer
!
!*******************************************************************************
Subroutine ran_flow
!*******************************************************************************
!	This procedure is to apply random flow to generate turbulence flow field
Implicit None
Integer :: i,j,k,jmin,jmax
!
Call initrandom(nrank+1)
!
jmin=1; jmax=xsize(2)
If (xmin(2)==1) jmin=3
If (xmax(2)==1) jmax=xsize(2)-2
!
Call rmsfl
!
Do k = 1,xsize(3); Do j = jmin,jmax; Do i = 1,xsize(1)
  u1(i,j,k) = u1(i,j,k)+amp*(rand_zbql()-0.5e0)
  v1(i,j,k) = v1(i,j,k)+amp*(rand_zbql()-0.5e0)
  w1(i,j,k) = w1(i,j,k)+amp*(rand_zbql()-0.5e0)
End Do; End Do; End Do
!
Return
End Subroutine ran_flow
!
!*******************************************************************************
Subroutine turb_stat_init
!*******************************************************************************
!	This subroutine initialise the Turbulence statistics arrays
Implicit None
!
samples = 0; taulmo = 0.0e0; stat_sum = 0.0e0; quad_sum = 0.0e0
tke_sum = 0.0e0; wx_sum = 0.0e0; umean = 0.0e0
!
Return
End Subroutine turb_stat_init
!
!******************************************************************************
Subroutine initial_out
!*****************************************************************************
!	This procedure Writes the output for initial flow field
Use mass_calc
Implicit None
Integer :: i,j,k
!
If (nrank == 0) then
  Write(21,150)
  150    Format('***************************** INPUT PARAMETERS *****************************')
!
  Write(21,*) ''
  Write(21,151)length,height,width,roh,re,nu,cfl,nx,ny,nz,nstep,     &
               nprint,itmavg,ntmavg,initma,ibmpfg
  151    Format('  length=',E10.3,'  height=',E10.3,'  width=',E10.3/, &
                '  rho=',E10.3,'  re=',E10.3,'  nu=',E10.3,'  cfl=',E10.3/,          & 
                '  nx=',I5,'  ny=',I5,'  nz=',I5,'  nstep=',I5/,                   &
                '  nprint=',I5,'  itmavg=',I5,'  ntmavg=',I5,'  initma=',I5,'  ibmpfg=',I5)
  Write(21,*) ''
!
  Write(21,152)
  152    Format('****************************************************************************')
!
  Write(21,*) ''
  Write(21,*) ' x = ' 
  Write(21,1511)(x(i),i=1,nxt)
  Write(21,*) ' y = ' 
  Write(21,1511)(y(j),j=1,nyt)
  Write(21,*) ' z = ' 
  Write(21,1511)(z(k),k=1,nzt)
  Write(21,*) ''
! 
  Write(21,152)
  Write(21,*) ''
  Write(21,*) ' dx = ' 
  Write(21,1511)(dx(i),i = -1,nx+2)
  Write(21,*) ' dy = ' 
  Write(21,1511)(dy(j),j = -1,ny+2)
  Write(21,*) ' dz = ' 
  Write(21,1511)(dz(k),k = -1,nz+2)
  1511   Format(5F15.7)
  Write(21,*) ''
! 
  Write(21,152)
  Write(21,*) ''
  Write(21,18)time
  18     Format(' time=',E15.7)
End If
!
Call rmsfl
!
If (nrank == 0) then
  Write(21,1904) rmsfli,rmsflm,rmsflo,dmsfl
  1904   Format(' mass flow at_inlet=',E15.8,/' mass flow at_middle=',E15.8,          &
               /' mass flow at_outlet=',E15.8,/' difference in massflow=',E15.8)
  Write(21,*) ''
  Write(21,152)
End If
!
!       Convective B.C.
!
uconv=rmsfli/height/width/roh
!
!
Return
End Subroutine initial_out
!
!******************************************************************************
Subroutine fftw_setup
!*****************************************************************************
!       This subroutine setups FFTW
Implicit None
!
 Call decomp_2d_fft_init
 Call decomp_2d_fft_get_size(fft_zstart,fft_zend,fft_zsize,3)
 Call decomp_2d_fft_get_size(fft_ystart,fft_yend,fft_ysize,2)
!
 If (ispec == 1) Then; Allocate (spec_sum(1:fft_zsize(1),1:fft_zsize(2),1:fft_zsize(3),1:4));
 spec_sum(:,:,:,:)=0.0e0; End If
!
Return
End Subroutine fftw_setup
!
!******************************************************************************
Subroutine finalize_all
!******************************************************************************
!	Finalizes all decompositions
Implicit None
!
! Call write_array(du1)
 Call close_files
 Call decomp_2d_fft_finalize
 Call decomp_2d_finalize
 Call MPI_Finalize(ierror)
!
Return
End Subroutine finalize_all
!
!******************************************************************************
Subroutine write_array(inarr)
!******************************************************************************
!-------This procedure saves velocity and pressure fields for future runs.
Use shared_data
Use decomp_2d
Use decomp_2d_io
Use mpi
Implicit None
integer :: fh_wres
integer(KIND=MPI_OFFSET_KIND) :: disp_res
Real(mytype), dimension(-1:xsize(1)+2,-1:xsize(2)+2,-1:xsize(3)+2) :: inarr
Real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: array1
!
array1(1:xsize(1),1:xsize(2),1:xsize(3))=inarr(1:xsize(1),1:xsize(2),1:xsize(3))
!
disp_res = 0_MPI_OFFSET_KIND
!
 Call MPI_FILE_OPEN(MPI_COMM_WORLD,'in.dat',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
			MPI_INFO_NULL,fh_wres,ierror)
!
 Call decomp_2d_write_head3d(fh_wres,disp_res)
 Call decomp_2d_write_var(fh_wres,disp_res,nx,ny,nz,1,array1)
!
  Call MPI_FILE_CLOSE(fh_wres,ierror) 
!
Return
End Subroutine write_array
!
End Module setup
