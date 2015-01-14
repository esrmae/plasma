module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  
  implicit none
  
  private        ! Make everything private unless declared public
  
  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1
  
  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: PHYSICAL_IN_Z = 3 

  integer, save :: plan_mode

  integer, save :: format                 ! input X-pencil or Z-pencil

  logical, save :: inplace                ! in-place FFT?	 
  
  !
  logical, save :: initialised = .false. 
  
  ! for c2r/r2c interface:
  !  - global size of real array is: nx*ny*nz
  !    * can reuse the main decomposition info in the base library 
  !  - global size of complex array is: (nx/2+1)*ny*nz
  !    which needs to be distributed in Z-pencil (if PHYSICAL_IN_X)
  !                                  or nx*ny*(nz/2+1)
  !    which needs to be distributed in X-pencil (if PHYSICAL_IN_Z)
  !    * define a second set of decomposition information
  TYPE(DECOMP_INFO), save :: fft

  
  
  public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
       decomp_2d_fft_finalize, decomp_2d_fft_get_size
  
  
  ! Declare generic interfaces to handle different inputs
  
  interface decomp_2d_fft_init
     module procedure fft_init_noarg
     module procedure fft_init_arg
  end interface
  
  interface decomp_2d_fft_3d
     module procedure fft_3d_c2c
     module procedure fft_3d_r2c
     module procedure fft_3d_c2r
  end interface
  
  
contains
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the FFT module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_init_noarg
    
    implicit none
    
    call fft_init_arg(PHYSICAL_IN_X)  ! default input is X-pencil data
    
    return
  end subroutine fft_init_noarg


  subroutine fft_init_arg(pencil)     ! allow to handle Z-pencil input

    implicit none
    
    integer, intent(IN) :: pencil  
    
    format = pencil

    plan_mode = 0   ! or 100 for better plan
    inplace = .false.

    if (format==PHYSICAL_IN_X) then
       call decomp_info_init(nx_global/2+1, ny_global, nz_global, fft)
    else if (format==PHYSICAL_IN_Z) then
       call decomp_info_init(nx_global, ny_global, nz_global/2+1, fft)
    end if
    
    initialised = .true.
    
    return
  end subroutine fft_init_arg
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Final clean up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_finalize
    
    implicit none

    call decomp_info_finalize(fft)
    
    return
  end subroutine decomp_2d_fft_finalize

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the size, starting/ending index of the distributed array 
  !  whose global size is (nx/2+1)*ny*nz, for defining data structures
  !  in r2c and c2r interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_get_size(istart, iend, isize)
    
    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize
    
    if (format==PHYSICAL_IN_X) then
       istart = fft%zst
       iend   = fft%zen
       isize  = fft%zsz
    else if (format==PHYSICAL_IN_Z) then
       istart = fft%xst
       iend   = fft%xen
       isize  = fft%xsz
    end if
    
    return
  end subroutine decomp_2d_fft_get_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign   
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2, wk2b, wk3
    integer :: k

    real(mytype), parameter :: scale = 1.0_mytype ! compute unscaled FFT
    complex(mytype), allocatable, dimension(:) :: comm      ! work space
    integer :: info                                         ! error flag 
    
#ifdef DOUBLE_PREC       
    allocate(comm(max(3*xsize(1)+100, &
         max(3*ysize(2)+100, 3*zsize(3)+100))))
#else
    allocate(comm(max(5*xsize(1)+100, &
         max(5*ysize(2)+100, 5*zsize(3)+100))))
#endif
    
    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
       allocate (wk1(xsize(1),xsize(2),xsize(3)))

#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,             &
            inplace,           & ! in-place?
            xsize(2)*xsize(3), & ! number of 1D FFTs
            xsize(1),          & ! size of each 1D FFT
            in,                & ! input
            1,                 & ! data point in each FFT separated by this
            xsize(1),          & ! different FFT sequences separated by this
            wk1, 1, xsize(1),  & ! output same layout as input
            comm, info)
       call zfft1mx(isign,scale,inplace,xsize(2)*xsize(3),xsize(1), &
            in,1,xsize(1),wk1,1,xsize(1),comm,info)
#else
       call cfft1mx(plan_mode,scale,inplace,xsize(2)*xsize(3),xsize(1), &
            in,1,xsize(1),wk1,1,xsize(1),comm,info)
       call cfft1mx(isign,scale,inplace,xsize(2)*xsize(3),xsize(1), &
            in,1,xsize(1),wk1,1,xsize(1),comm,info)
#endif

       ! ===== Swap X --> Y =====
       allocate (wk2(ysize(1),ysize(2),ysize(3)))
       allocate (wk2b(ysize(1),ysize(2),ysize(3)))
       call transpose_x_to_y(wk1,wk2)
       
       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,             &
            inplace,           & ! in-place?
            ysize(1),          & ! number of 1D FFTs
            ysize(2),          & ! size of each 1D FFT
            wk2,               & ! input
            ysize(1),          & ! data point in each FFT separated by this
            1,                 & ! different FFT sequences separated by this
            wk2b, ysize(1), 1, & ! output same layout as input
            comm, info)
#else
       call cfft1mx(plan_mode,scale,inplace,ysize(1),ysize(2), &
            wk2,ysize(1),1,wk2b,ysize(1),1,comm,info)
#endif

       do k=1, ysize(3)  ! loop through Z-planes
#ifdef DOUBLE_PREC
          call zfft1mx(isign,scale,inplace,ysize(1),ysize(2), &
               wk2(1,1,k),ysize(1),1,wk2b(1,1,k),ysize(1),1,comm,info)
#else
          call cfft1mx(isign,scale,inplace,ysize(1),ysize(2), &
               wk2(1,1,k),ysize(1),1,wk2b(1,1,k),ysize(1),1,comm,info)
#endif
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(zsize(1),zsize(2),zsize(3)))
       call transpose_y_to_z(wk2b,wk3)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,             &
            inplace,           & ! in-place?
            zsize(1)*zsize(2), & ! number of 1D FFTs
            zsize(3),          & ! size of each 1D FFT
            wk3,               & ! input
            zsize(1)*zsize(2), & ! data point in each FFT separated by this
            1,                 & ! different FFT sequences separated by this
            out, zsize(1)*zsize(2), 1,   & ! output same layout as input
            comm, info)
       call zfft1mx(isign,scale,inplace,zsize(1)*zsize(2),zsize(3), &
            wk3,zsize(1)*zsize(2),1,out,zsize(1)*zsize(2),1,comm,info)
#else
       call cfft1mx(plan_mode,scale,inplace,zsize(1)*zsize(2),zsize(3), &
            wk3,zsize(1)*zsize(2),1,out,zsize(1)*zsize(2),1,comm,info)
       call cfft1mx(isign,scale,inplace,zsize(1)*zsize(2),zsize(3), &
            wk3,zsize(1)*zsize(2),1,out,zsize(1)*zsize(2),1,comm,info)
#endif       

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(zsize(1),zsize(2),zsize(3)))

#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,             &
            inplace,           & ! in-place?
            zsize(1)*zsize(2), & ! number of 1D FFTs
            zsize(3),          & ! size of each 1D FFT
            in,                & ! input
            zsize(1)*zsize(2), & ! data point in each FFT separated by this
            1,                 & ! different FFT sequences separated by this
            wk1, zsize(1)*zsize(2), 1,   & ! output same layout as input
            comm, info)
       call zfft1mx(isign,scale,inplace,zsize(1)*zsize(2),zsize(3), &
            in,zsize(1)*zsize(2),1,wk1,zsize(1)*zsize(2),1,comm,info)
#else
       call cfft1mx(plan_mode,scale,inplace,zsize(1)*zsize(2),zsize(3), &
            in,zsize(1)*zsize(2),1,wk1,zsize(1)*zsize(2),1,comm,info)
       call cfft1mx(isign,scale,inplace,zsize(1)*zsize(2),zsize(3), &
            in,zsize(1)*zsize(2),1,wk1,zsize(1)*zsize(2),1,comm,info)
#endif
       
       ! ===== Swap Z --> Y =====
       allocate (wk2(ysize(1),ysize(2),ysize(3)))
       allocate (wk2b(ysize(1),ysize(2),ysize(3)))
       call transpose_z_to_y(wk1,wk2)
       
       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,             &
            inplace,           & ! in-place?
            ysize(1),          & ! number of 1D FFTs
            ysize(2),          & ! size of each 1D FFT
            wk2,               & ! input
            ysize(1),          & ! data point in each FFT separated by this
            1,                 & ! different FFT sequences separated by this
            wk2b, ysize(1), 1, & ! output same layout as input
            comm, info)
#else
       call cfft1mx(plan_mode,scale,inplace,ysize(1),ysize(2), &
            wk2,ysize(1),1,wk2b,ysize(1),1,comm,info)
#endif

       do k=1, ysize(3)  ! loop through Z-planes
#ifdef DOUBLE_PREC
          call zfft1mx(isign,scale,inplace,ysize(1),ysize(2), &
               wk2(1,1,k),ysize(1),1,wk2b(1,1,k),ysize(1),1,comm,info)
#else
          call cfft1mx(isign,scale,inplace,ysize(1),ysize(2), &
               wk2(1,1,k),ysize(1),1,wk2b(1,1,k),ysize(1),1,comm,info)
#endif
       end do
       
       ! ===== Swap Y --> X =====
       allocate (wk3(xsize(1),xsize(2),xsize(3)))
       call transpose_y_to_x(wk2b,wk3)
       
       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,             &
            inplace,           & ! in-place?
            xsize(2)*xsize(3), & ! number of 1D FFTs
            xsize(1),          & ! size of each 1D FFT
            wk3,               & ! input
            1,                 & ! data point in each FFT separated by this
            xsize(1),          & ! different FFT sequences separated by this
            out, 1, xsize(1),  & ! output same layout as input
            comm, info)
       call zfft1mx(isign,scale,inplace,xsize(2)*xsize(3),xsize(1), &
            wk3,1,xsize(1),out,1,xsize(1),comm,info)
#else
       call cfft1mx(plan_mode,scale,inplace,xsize(2)*xsize(3),xsize(1), &
            wk3,1,xsize(1),out,1,xsize(1),comm,info)
       call cfft1mx(isign,scale,inplace,xsize(2)*xsize(3),xsize(1), &
            wk3,1,xsize(1),out,1,xsize(1),comm,info)
#endif      
       
    end if

    deallocate(comm)
    
    return
  end subroutine fft_3d_c2c
  
  
  ! ===== real interfaces experimental in this version  =====
  !       - planning to be done at initialisation time
  !       - need optimisation of memory usage
  !       - need more tests
  !       - need parallel performance study
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2, wk2b, wk3
    integer :: i,j,k

    real(mytype), parameter :: scale = 1.0_mytype ! compute unscaled FFT
    real(mytype),    allocatable, dimension(:) :: comm1 ! work sapce
    complex(mytype), allocatable, dimension(:) :: comm2 ! work space
    integer :: info 

    real(mytype), allocatable, dimension(:) :: tmp, zbuf_r
    complex(mytype), allocatable, dimension(:) :: zbuf_c

    if (format==PHYSICAL_IN_X) then

       allocate(comm1(3*xsize(1)+100))
#ifdef DOUBLE_PREC
       allocate(comm2(max(3*fft%ysz(2)+100, 3*fft%zsz(3)+100)))
#else
       allocate(comm2(max(5*fft%ysz(2)+100, 5*fft%zsz(3)+100)))
#endif

       ! ===== 1D FFTs in X =====
       allocate(wk1(fft%xsz(1),fft%xsz(2),fft%xsz(3)))   
#ifdef DOUBLE_PREC
       call dzfftm(   &
            xsize(2)*xsize(3),    & ! number of 1D FFTs
            xsize(1),             & ! FFT size
            in_r,comm1,info)
#else
       call scfftm(xsize(2)*xsize(3),xsize(1),in_r,comm1,info)
#endif

       ! translate Hermitian storage into half-plus-1-complex-storage
       allocate(tmp(xsize(1)))
       do k=1,xsize(3)
          do j=1,xsize(2)
             tmp(1:xsize(1)) = in_r(1:xsize(1),j,k)
             call to_complex(tmp, xsize(1), wk1(1,j,k))
          end do
       end do
       deallocate(tmp)
       wk1 = wk1 * sqrt(real(xsize(1), kind=mytype))

       ! ===== Swap X --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_x_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
#ifdef DOUBLE_PREC       
       call zfft1mx(plan_mode, &
            scale,       &
            inplace,     & ! in-place?
            fft%ysz(1),  & ! number of 1D FFTs
            fft%ysz(2),  & ! size of each 1D FFT
            wk2,         & ! input
            fft%ysz(1),  & ! data point in each FFT separated by this
            1,           & ! different FFT sequences separated by this
            wk2b, fft%ysz(1), 1,  & ! output same layout as input
            comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%ysz(1),fft%ysz(2),  &
            wk2,fft%ysz(1),1,wk2b,fft%ysz(1),1,comm2,info)
#endif 
       ! compute one z-plane at a time
       do k=1,fft%ysz(3)
#ifdef DOUBLE_PREC
          call zfft1mx(-1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1,comm2,info)
#else
          call cfft1mx(-1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1,comm2,info)
#endif
       end do
       
       ! ===== Swap Y --> Z =====
       allocate (wk3(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,fft)
       
       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode,     &
            scale,                 &
            inplace,               & ! in-place?
            fft%zsz(1)*fft%zsz(2), & ! number of 1D FFTs
            fft%zsz(3),            & ! size of each 1D FFT
            wk3,                   & ! input
            fft%zsz(1)*fft%zsz(2), & ! data point in each FFT separated by this
            1,                     & ! different FFT sequences separated by this
            out_c, fft%zsz(1)*fft%zsz(2), 1,  & ! output same layout
            comm2, info)
       call zfft1mx(-1,scale,inplace,fft%zsz(1)*fft%zsz(2),fft%zsz(3), &
            wk3,fft%zsz(1)*fft%zsz(2),1,out_c,fft%zsz(1)*fft%zsz(2),1, &
            comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%zsz(1)*fft%zsz(2),fft%zsz(3), &
            wk3,fft%zsz(1)*fft%zsz(2),1,out_c,fft%zsz(1)*fft%zsz(2),1, &
            comm2, info)
       call cfft1mx(-1,scale,inplace,fft%zsz(1)*fft%zsz(2),fft%zsz(3), &
            wk3,fft%zsz(1)*fft%zsz(2),1,out_c,fft%zsz(1)*fft%zsz(2),1, &
            comm2, info)
#endif

       deallocate(comm1, comm2)

    else if (format==PHYSICAL_IN_Z) then

       allocate(comm1(3*zsize(3)+100))
#ifdef DOUBLE_PREC
       allocate(comm2(max(3*fft%xsz(1)+100, 3*fft%ysz(2)+100)))
#else
       allocate(comm2(max(5*fft%xsz(1)+100, 5*fft%ysz(2)+100)))
#endif

       ! ===== 1D FFTs in Z =====
       ! note that there is no expert driver to perform multiple 1D r2c
       ! FFTs in ACML --> need to copy numbers into buffers and then
       ! make multiple calls to simple driver.
       allocate(zbuf_r(zsize(3)))
       allocate(zbuf_c(fft%zsz(3)))
       allocate(wk1(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
#ifdef DOUBLE_PREC
       call dzfft(plan_mode,zsize(3),zbuf_r,comm1,info)
#else
       call scfft(plan_mode,zsize(3),zbuf_r,comm1,info)
#endif
       do j=1,zsize(2)
          do i=1,zsize(1)
             ! copy data set along Z into a buffer
             do k=1,zsize(3)
                zbuf_r(k) = in_r(i,j,k)
             end do
             ! 1D r2c FFT 
#ifdef DOUBLE_PREC
             call dzfft(1,zsize(3),zbuf_r,comm1,info)
#else
             call scfft(1,zsize(3),zbuf_r,comm1,info)
#endif
             ! The result is in Hermitian storage, convert to complex
             call to_complex(zbuf_r,zsize(3),zbuf_c)
             ! copy back into 3D data structure
             do k=1,fft%zsz(3) 
                wk1(i,j,k) = zbuf_c(k)
             end do
          end do
       end do
       ! remove ACML scaling, the parallel library never scales anything
       wk1 = wk1 * sqrt(real(zsize(3), kind=mytype))

       ! ===== Swap Z --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_z_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
#ifdef DOUBLE_PREC       
       call zfft1mx(plan_mode, &
            scale,       &
            inplace,     & ! in-place?
            fft%ysz(1),  & ! number of 1D FFTs
            fft%ysz(2),  & ! size of each 1D FFT
            wk2,         & ! input
            fft%ysz(1),  & ! data point in each FFT separated by this
            1,           & ! different FFT sequences separated by this
            wk2b, fft%ysz(1), 1,  & ! output same layout as input
            comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%ysz(1),fft%ysz(2),  &
            wk2,fft%ysz(1),1,wk2b,fft%ysz(1),1,comm2,info)
#endif 
       ! compute one z-plane at a time
       do k=1,fft%ysz(3)
#ifdef DOUBLE_PREC
          call zfft1mx(-1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1,comm2,info)
#else
          call cfft1mx(-1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1,comm2,info)
#endif
       end do      

       ! ===== Swap Y --> X =====
       allocate (wk3(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,fft)

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode,     &
            scale,                 &
            inplace,               & ! in-place?
            fft%xsz(2)*fft%xsz(3), & ! number of 1D FFTs
            fft%xsz(1),              & ! size of each 1D FFT
            wk3,                   & ! input
            1,                     & ! data point in each FFT separated by this
            fft%xsz(1),            & ! different FFT sequences separated by this
            out_c, 1, fft%xsz(1),  & ! output same layout as input
            comm2, info)
       call zfft1mx(-1,scale,inplace,fft%xsz(2)*fft%xsz(3),fft%xsz(1), &
            wk3,1,fft%xsz(1),out_c, 1, fft%xsz(1),comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%xsz(2)*fft%xsz(3),fft%xsz(1), &
            wk3,1,fft%xsz(1),out_c, 1, fft%xsz(1),comm2, info)
       call cfft1mx(-1,scale,inplace,fft%xsz(2)*fft%xsz(3),fft%xsz(1), &
            wk3,1,fft%xsz(1),out_c, 1, fft%xsz(1),comm2, info)
#endif

       deallocate(comm1, comm2)

    end if
    
    return
  end subroutine fft_3d_r2c
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2, wk2b, wk3
    integer :: i,j,k

    real(mytype), parameter :: scale = 1.0_mytype ! compute unscaled FFT
    real(mytype),    allocatable, dimension(:) :: comm1 ! work sapce
    complex(mytype), allocatable, dimension(:) :: comm2 ! work space
    integer :: info 

    real(mytype), allocatable, dimension(:) :: tmp, zbuf_r 
    complex(mytype), allocatable, dimension(:) :: zbuf_c

    if (format==PHYSICAL_IN_X) then

       allocate(comm1(3*xsize(1)+100))
#ifdef DOUBLE_PREC
       allocate(comm2(max(3*fft%ysz(2)+100, 3*fft%zsz(3)+100)))
#else
       allocate(comm2(max(5*fft%ysz(2)+100, 5*fft%zsz(3)+100)))
#endif 

       ! ===== 1D FFTs in Z =====
       allocate(wk1(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode,      &
            scale,                  &
            inplace,                & ! in-place?
            fft%zsz(1)*fft%zsz(2),  & ! number of 1D FFTs
            fft%zsz(3),             & ! size of each 1D FFT
            in_c,                   & ! input
            fft%zsz(1)*fft%zsz(2),  & ! data point in each FFT separated by this
            1,                      & ! different FFT sets separated by this
            wk1, fft%zsz(1)*fft%zsz(2), 1,   & ! output same layout as input
            comm2, info)
       call zfft1mx(1,scale,inplace,fft%zsz(1)*fft%zsz(2),fft%zsz(3),  &
            in_c,fft%zsz(1)*fft%zsz(2),1,wk1, fft%zsz(1)*fft%zsz(2), 1, &
            comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%zsz(1)*fft%zsz(2),fft%zsz(3),  &
            in_c,fft%zsz(1)*fft%zsz(2),1,wk1, fft%zsz(1)*fft%zsz(2), 1, &
            comm2, info)
       call cfft1mx(1,scale,inplace,fft%zsz(1)*fft%zsz(2),fft%zsz(3),  &
            in_c,fft%zsz(1)*fft%zsz(2),1,wk1, fft%zsz(1)*fft%zsz(2), 1, &
            comm2, info)
#endif 

       ! ===== Swap Z --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_z_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,       &
            inplace,     & ! in-place?
            fft%ysz(1),  & ! number of 1D FFTs
            fft%ysz(2),  & ! size of each 1D FFT
            wk2,         & ! input
            fft%ysz(1),  & ! data point in each FFT separated by this
            1,           & ! different FFT sequences separated by this
            wk2b, fft%ysz(1), 1,  & ! output same layout as input
            comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%ysz(1),fft%ysz(2),  &
            wk2,fft%ysz(1),1,wk2b,fft%ysz(1),1,comm2,info)
#endif
       ! compute one z-plane at a time
       do k=1,fft%ysz(3)
#ifdef DOUBLE_PREC
          call zfft1mx(1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1, &
               comm2,info)
#else
          call cfft1mx(1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1, &
               comm2,info)
#endif
       end do

       ! ===== Swap Y --> X =====
       allocate (wk3(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,fft)
       
       ! translate half-plus-1-complex-storage into Hermitian storage
       allocate(tmp(xsize(1)))
       do k=1,xsize(3)
          do j=1,xsize(2)
             call to_hermitian(wk3(1,j,k), xsize(1), tmp)
             out_r(1:xsize(1),j,k) = tmp(1:xsize(1))
          end do
       end do
       deallocate(tmp)
       out_r = out_r * sqrt(real(xsize(1),kind=mytype))

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call zdfftm(xsize(2)*xsize(3),xsize(1),out_r,comm1,info)
#else
       call csfftm(xsize(2)*xsize(3),xsize(1),out_r,comm1,info)
#endif

       deallocate(comm1,comm2) 

    else if (format==PHYSICAL_IN_Z) then

       allocate(comm1(3*zsize(3)+100))
#ifdef DOUBLE_PREC
       allocate(comm2(max(3*fft%ysz(2)+100, 3*fft%xsz(1)+100)))
#else
       allocate(comm2(max(5*fft%ysz(2)+100, 5*fft%xsz(1)+100)))
#endif

       ! ===== 1D FFTs in X =====
       allocate(wk1(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode,     &
            scale,                 &
            inplace,               & ! in-place?
            fft%xsz(2)*fft%xsz(3), & ! number of 1D FFTs
            fft%xsz(1),            & ! size of each 1D FFT
            in_c,                  & ! input
            1,                     & ! data point in each FFT separated by this
            fft%xsz(1),            & ! different FFT sequences separated by this
            wk1, 1, fft%xsz(1),    & ! output same layout as input
            comm2, info)
       call zfft1mx(1,scale,inplace,fft%xsz(2)*fft%xsz(3),fft%xsz(1), &
            in_c,1,fft%xsz(1),wk1,1,fft%xsz(1),comm2,info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%xsz(2)*fft%xsz(3),fft%xsz(1), &
            in_c,1,fft%xsz(1),wk1,1,fft%xsz(1),comm2,info)
       call cfft1mx(1,scale,inplace,fft%xsz(2)*fft%xsz(3),fft%xsz(1), &
            in_c,1,fft%xsz(1),wk1,1,fft%xsz(1),comm2,info)
#endif

       ! ===== Swap X --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_x_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
#ifdef DOUBLE_PREC
       call zfft1mx(plan_mode, &
            scale,       &
            inplace,     & ! in-place?
            fft%ysz(1),  & ! number of 1D FFTs
            fft%ysz(2),  & ! size of each 1D FFT
            wk2,         & ! input
            fft%ysz(1),  & ! data point in each FFT separated by this
            1,           & ! different FFT sequences separated by this
            wk2b, fft%ysz(1), 1,  & ! output same layout as input
            comm2, info)
#else
       call cfft1mx(plan_mode,scale,inplace,fft%ysz(1),fft%ysz(2),  &
            wk2,fft%ysz(1),1,wk2b,fft%ysz(1),1,comm2,info)
#endif
       ! compute one z-plane at a time
       do k=1,fft%ysz(3)
#ifdef DOUBLE_PREC
          call zfft1mx(1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1, &
               comm2,info)
#else
          call cfft1mx(1,scale,inplace,fft%ysz(1),fft%ysz(2),  &
               wk2(1,1,k),fft%ysz(1),1,wk2b(1,1,k),fft%ysz(1),1, &
               comm2,info)
#endif
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,fft)

       ! ===== 1D FFTs in Z =====
       ! note that there is no expert driver to perform multiple 1D c2r
       ! FFTs in ACML --> need to copy numbers into buffers and then
       ! make multiple calls to simple driver.
       allocate(zbuf_r(zsize(3)))
       allocate(zbuf_c(fft%zsz(3)))
#ifdef DOUBLE_PREC
       call zdfft(plan_mode,zsize(3),zbuf_r,comm1,info)
#else
       call csfft(plan_mode,zsize(3),zbuf_r,comm1,info)
#endif
       do j=1,zsize(2)
          do i=1,zsize(1)
             ! copy data set along Z into a buffer
             do k=1,fft%zsz(3) 
                zbuf_c(k) = wk3(i,j,k)
             end do
             ! convert complex storage to Hermitian storage
             call to_hermitian(zbuf_c,zsize(3),zbuf_r)
             ! 1D c2r FFT
#ifdef DOUBLE_PREC
             call zdfft(1,zsize(3),zbuf_r,comm1,info)
#else
             call csfft(1,zsize(3),zbuf_r,comm1,info)
#endif
             ! copy back into 3D data structure
             do k=1,zsize(3)
                out_r(i,j,k) = zbuf_r(k)
             end do
          end do
       end do
       ! remove ACML scaling, the parallel library never scales anything
       out_r = out_r * sqrt(real(zsize(3),kind=mytype))

       deallocate(comm1,comm2) 

    end if

    return
  end subroutine fft_3d_c2r


  ! Two utility routines to convert internal storage format

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Rearrange Hermitian storage data as half-plus-1-complex-storage data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine to_complex(in,n,out)
    
    implicit none
    
    integer, intent(IN) :: n
    real(mytype), dimension(0:n-1), intent(IN) :: in
    complex(mytype), dimension(0:n/2), intent(OUT) :: out
    
    integer :: i
    
    ! Let out(i) be the full-complex spectral representation of real input
    ! Let 'in' be an Hermitian storage array of length 'n' (0-based)
    !  - in(i) contains the real part of out(i) for i=0,...,n/2
    !  - in(n−i) contains the imaginary part of out(i) for i=1,...,(n−1)/2

    out(0) = cmplx(in(0),0._mytype,kind=mytype)
    
    if (n == n/2*2) then   ! n is even, n/2 is larger than (n-1)/2 by 1
       do i=1,(n-1)/2
          out(i) = cmplx(in(i),in(n-i),kind=mytype)
       end do
       out(n/2) = cmplx(in(n/2),0._mytype,kind=mytype)
    else                   ! n is odd, then n/2 == (n-1)/2 
       do i=1,n/2
          out(i) = cmplx(in(i),in(n-i),kind=mytype)
       end do
    end if
    
    return
  end subroutine to_complex

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The opposite of the above
  !   Also need to be conjudated before inverse transform
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine to_hermitian(in,n,out)
    
    implicit none
    
    integer, intent(IN) :: n
    complex(mytype), dimension(0:n/2), intent(IN) :: in
    real(mytype), dimension(0:n-1), intent(OUT) :: out
    
    integer :: i

    out(0) = real(in(0),kind=mytype)
    
    if (n == n/2*2) then   ! n is even, n/2 is larger than (n-1)/2 by 1
       do i=1,(n-1)/2
          out(i) = real(in(i),kind=mytype)
          out(n-i) = -aimag(in(i))
       end do
       out(n/2) = real(in(n/2), kind=mytype)
    else                   ! n is odd, then n/2 == (n-1)/2 
       do i=1,n/2
          out(i) = real(in(i),kind=mytype)
          out(n-i) = -aimag(in(i))
       end do
    end if

    return
  end subroutine to_hermitian

  
end module decomp_2d_fft
