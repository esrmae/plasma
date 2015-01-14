 Module shared_data
 Use decomp_2d
 Use mpi

 Implicit None
!
!	Main simulation parameters
!
 Integer :: nstep,nx,ny,nz,nxt,nyt,nzt,iconbc,ibcjet,ngr,nblock,nbl,ibwest,ibeast,ibsoth,ibnoth,     &
            samples,ii,icmpnt,kbsoth,kbnoth,ibmpfg,igr_2pt,imass,istart,iseed,irandm,idomain,ifftw,   &
            isgsm,jconbc,kconbc,ichann,if2d,if2d_noy,if2d_noz,ifdt,ncmpnt,idt_ret,idt_ret_un,ifgrz_2pt,      &
            ifdt_fixed,ifdt_variable,if_dpdxin,if_noslipx,if_noslipy,if_noslipz,if_mbc,lstop
!
 Real(mytype) :: length,height,width,roh,cfl,dt,nu,re,dpmean,retau,dpz,uconv,time,time0,conx,cony,conz,    &
                amp,rmsfli,rmsflm,rmsflo,rmsflc,rmsflb,dmsfl,del,pi,pi2,pi4,dpdx_mean,dpdz_mean,large,      &
                small,very_small,dt_fixed,dt_initial,deltat_dt,time_dtcon,dpdx_mean_in,time_stop,myu
!
 Real(mytype), Dimension(:), Allocatable :: x,y,z,dx,dy,dz,xc,yc,zc,dxs,dys,dzs,xplus,yplus,zplus,          &
                                           xcplus,ycplus,zcplus,buff,buffc,bre,brec
 Real(mytype), Dimension(:,:,:), Allocatable :: u1,v1,w1,p1,u2,v2,w2,p2,u3,v3,w3,p3,uhat1,vhat1,what1,      &
                                               varray,farray,ap,ae,as,at,intfuo,intfvo,intfwo,grdsxx,       &
                                               grdsyy,grdszz,du1,du2,intf1,intf2,intf3             
 Real(mytype), Dimension(:,:,:,:), Allocatable :: gcdx_f,gcdy_f,gcdz_f,gcdx_c,gcdy_c,gcdz_f2,gfdx_c,        &
                                                 gfdy_c,gfdz_c,gfdz_c2,gcdz_c
 Real(mytype), Dimension(:,:,:,:,:), Allocatable :: xinlet,yinlet1,yinlet2,yinlet3,zinlet
!
!	LES parameters
!
 Integer :: nx2,ny2,nz2,isgs_model,if_cs_2d_avg,if_cw_var,if_backscatter,nsmles,iread_les,irite_les
 Real(mytype) :: c_w,cs_smg
 Real(mytype), Dimension(:), Allocatable :: cs
 Real(mytype), Dimension(:,:,:), Allocatable :: pitero1,pitero2,pitero3,rnutxy1,rnutyz1,rnutxz1,rnutxy2,    &
                                               rnutyz2,rnutxz2,rnutxy3,rnutyz3,rnutxz3
!
!	Output arrays
!
 Integer :: isavvf,nprint,itmavg,ntmavg,initma,lpstat_1d,lpstat_2d,lpstat_x,lpstat_z,lp_ins_1d,np_ins_1d,   &
           lp_ins_2d,np_ins_2d,lp_ins_3d,np_ins_3d,lp_snap_1d,np_snap_1d,lp_snap_2d,np_snap_2d,             &
           lp_snap_3d,np_snap_3d,lp_snap_x,lp_snap_y,lp_snap_z,lp_snap_xz,lp_snap_yz,lp_snap_xy,            &
           nsmpl_ins_1d,nsmpl_ins_2d,nsmpl_ins_3d,nsmpl_snap_1d,nsmpl_snap_2d,nsmpl_snap_3d,nsmpl_sh_2d,    &
	   restrt_old,xp_snap_1d,yp_snap_1d,zp_snap_1d,xp_snap_2d,yp_snap_2d,zp_snap_2d,ispec,ishear
 Integer, Dimension(:), Allocatable :: ntmavg2
 Real(mytype), Dimension(:), Allocatable :: umean,taulmo,vmean
 Real(mytype), Dimension(:,:), Allocatable :: omega_sum,quad_sum
 Real(mytype), Dimension(:,:,:), Allocatable :: tke_sum,wx_sum,stat_sum
 Real(mytype), Dimension(:,:,:,:), Allocatable :: quad2d_sum,omega2d_sum,spec_sum
 Real(mytype), Dimension(:,:,:,:,:), Allocatable :: stat2d_sum,tke2d_sum,wx2d_sum
!
!	Output files
!
 Integer :: ibin_read,ibin_rite,restrt,ir_p_old,ic_varray_rd,ic_varray_rite
 Character(LEN=60) :: unit21,unit22,unit23,unit30,unit31
!
!	FFTW parameters
!
 Real(mytype), Dimension(:), Allocatable :: infft1
 Complex(mytype), Dimension(:), Allocatable :: outfft1,infft2,outfft2
 Integer*8 :: planzf,planzb,planxf,planxb
!
!	Implcit solver variables
!
 Real(mytype), Dimension(:,:), Allocatable :: iux_vis,ivy_vis,iuy_con,ivx_con,iux_grdvis,ivy_grdvis,             &
                                         iuy_grdcon,ivx_grdcon,iux_grdcon,ivy_grdcon,iux_con,ivy_con
 Real(mytype), Dimension(:,:,:), Allocatable :: iuy_vis,ivx_vis,iuy_grdvis,ivx_grdvis
!
!	Wall Osc parameters
!
 Integer :: itrav,itravx,itravz
 Real(mytype) :: jet_type,theta_s,wmax,wmax_plus,kx_plus,kx,kz_plus,kz,omeg,omega_plus,fact_v,fact_w,amax,tperiod,     &
                 tperiod_plus,ampjet,deltat
!
!	Body force parameters
!
 Integer :: idx,idz,ndx,ndz,icub,nphase,iphase,i2dust,iforce,nposxz,iplasma,nactu,iosci,iuni,ibid
 Real(mytype) :: hjet, hjet_plus, St, St_plus, sigma, lambda, cspeed, cspeed_plus
 Real(mytype), Dimension(:,:,:), Allocatable :: fbody
 Real(mytype), Dimension(:), Allocatable :: xfringe, zfringe
 Integer, Dimension(:), Allocatable :: posxz
 Real(mytype) :: rise
!
!	Folder variables
!
 Character (LEN=60) ::  foldera,folderb,folderc,folderd,folder 
!
!	MPI variables
! 
 Integer :: p_row,p_col,ierror
 Integer :: right_row,left_row,right_col,left_col,col_rank,row_rank
 Integer, Dimension(MPI_STATUS_SIZE) :: status_info
 Integer, Dimension(3) :: xmin,xmax,ymin,ymax,zmin,zmax,lo,up,xwrite
 Integer, Dimension(3) :: fft_zstart,fft_zend,fft_zsize,fft_ystart,fft_yend,fft_ysize
 Integer :: fh_av1d,fh_ins1d,fh_ins2d,fh_ins3d,fh_sn1d,fh_sn3d,fh_spec,fh_sh2d
 Integer, Allocatable, Dimension(:) :: fh_av2d,fh_sn2d
 Integer(KIND=MPI_OFFSET_KIND) :: disp_av1d,disp_av2d,disp_ins1d,disp_ins2d,disp_sn1d,disp_sn3d,disp_spec,disp_sh2d
 Integer(KIND=MPI_OFFSET_KIND), Dimension(:), Allocatable :: disp_sn2d
!
!	Particle trace variables
!
 Integer :: ipart,npart_max,npart_slc,npart_nozz,npart_in
 Real(mytype) :: y1part,y2part,y3part
 Real(mytype) :: starty,startz,endy,endz  ! position of wire line
 Real(mytype), Dimension(:,:,:), Allocatable :: particle
 Integer, Dimension(:,:), Allocatable :: inpart
 Real(mytype), Dimension(:), Allocatable :: ypart
!
!	post3d
!
 Integer :: ilamq,iwfield,iflucs
 Real(mytype), Dimension(:,:,:), Allocatable :: l2_array
!
End Module shared_data
!*****************************************************************************
