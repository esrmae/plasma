module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module

  implicit none

  include "fftw3.f"
  
  private        ! Make everything private unless declared public
  
  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = 1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 2
  
  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: PHYSICAL_IN_Z = 3 
  integer, save :: format 
  
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
  TYPE(DECOMP_INFO), save, public :: fft

  ! FFTW plans 
  
  ! for various forward=1/backward=2 c2c transformations
  integer*8, save :: x_plan(2), y_plan(2), z_plan(2) 
  
  ! for real transformations
  integer*8, save :: x_plan_f_r2c, x_plan_b_c2r 
  
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
    
    ! dummy array used for planning
    complex(mytype), allocatable, dimension(:,:,:) :: a1, a2
    complex(mytype), allocatable, dimension(:,:) :: a3, a4

    integer :: plan_type

    plan_type = FFTW_ESTIMATE
    
    format = pencil

    if (format==PHYSICAL_IN_X) then
       call decomp_info_init(nx_global/2+1, ny_global, nz_global, fft)
    else if (format==PHYSICAL_IN_Z) then
       call decomp_info_init(nx_global, ny_global, nz_global/2+1, fft)
    end if

    ! ===== planning multiple FFTs in x =====
    
    allocate(a1(xsize(1),xsize(2),xsize(3)))
    allocate(a2(xsize(1),xsize(2),xsize(3)))
    
#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(x_plan(1),& 
         1,                            & ! 1D FFTs 
         xsize(1),                     & ! FFT size
         xsize(2)*xsize(3),            & ! Number of 1D FFTs to compute
         a1,                           & ! input 3D array
         xsize(1),                     & ! input 'nembed'
         1,                            & ! 1D data along X-axis continuous 
                                         !  in memory ==> stride = 1
         xsize(1),                     & ! starting point of each 1D FFT 
                                         !  separated by the whole length 
                                         !  of X dimension
         a2, xsize(1), 1, xsize(1),    & ! output same as input for c2c
         FFTW_FORWARD,  plan_type)
    call dfftw_plan_many_dft(x_plan(2), 1, xsize(1), xsize(2)*xsize(3), &
         a1, xsize(1), 1, xsize(1), a2, xsize(1), 1, xsize(1),          &
         FFTW_BACKWARD, plan_type)
#else
    call sfftw_plan_many_dft(x_plan(1), 1, xsize(1), xsize(2)*xsize(3), &
         a1, xsize(1), 1, xsize(1), a2, xsize(1), 1, xsize(1),          &
         FFTW_FORWARD,  plan_type)
    call sfftw_plan_many_dft(x_plan(2), 1, xsize(1), xsize(2)*xsize(3), &
         a1, xsize(1), 1, xsize(1), a2, xsize(1), 1, xsize(1),          &
         FFTW_BACKWARD, plan_type)
#endif
    
    deallocate(a1,a2)
    
    ! ===== planning multiple FFTs in Y =====
    ! Due to memory pattern of 3D arrays, 1D FFTs over Y lines have to be
    ! done one Z-plane at a time, in order to use the advanced interface.
    
    allocate(a3(ysize(1),ysize(2)))
    allocate(a4(ysize(1),ysize(2)))
    
#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(y_plan(1),&
         1,                            & ! 1D FFTs
         ysize(2),                     & ! FFT size
         ysize(1),                     & ! number of 1D FFTs to compute
         a3,                           & ! input 3D array
         ysize(2),                     & ! input 'nembed'
         ysize(1),                     & ! 1D data point along y-axis 
                                         !  separated by the whole length 
                                         !  of x dimension in memory 
                                         !  stride = ysize(1)
         1,                            & ! starting point of each 1D FFT 
                                         !  next to each other in memory
         a4, ysize(2), ysize(1), 1,    & ! output same as input for c2c
         FFTW_FORWARD,  plan_type)
    call dfftw_plan_many_dft(y_plan(2), 1, ysize(2), ysize(1),   &
         a3, ysize(2), ysize(1), 1, a4, ysize(2), ysize(1), 1,   &
         FFTW_BACKWARD, plan_type)
#else
    call sfftw_plan_many_dft(y_plan(1), 1, ysize(2), ysize(1),   &
         a3, ysize(2), ysize(1), 1, a4, ysize(2), ysize(1), 1,   &
         FFTW_FORWARD,  plan_type)
    call sfftw_plan_many_dft(y_plan(2), 1, ysize(2), ysize(1),   &
         a3, ysize(2), ysize(1), 1, a4, ysize(2), ysize(1), 1,   &
         FFTW_BACKWARD, plan_type)
#endif
    
    deallocate(a3,a4)
    
    ! ===== planning multiple FFTs in Z =====
    
    allocate(a1(zsize(1),zsize(2),zsize(3)))
    allocate(a2(zsize(1),zsize(2),zsize(3)))
    
#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(z_plan(1),&
         1,                            & ! 1D FFTs
         zsize(3),                     & ! FFT size
         zsize(1)*zsize(2),            & ! number of 1D FFTs to compute
         a1,                           & ! input 3D array
         zsize(3),                     & ! input 'nembed'
         zsize(1)*zsize(2),            & ! 1D data point along z-axis 
                                         !  separated by the size of xy 
                                         !  plane in memory - 
                                         !  stride = zsize(1)*zsize(2)
         1,                            & ! starting point of each 1D FFT 
                                         !  next to each other in memory
         a2, zsize(3),                 & ! output same as input for c2c
         zsize(1)*zsize(2), 1,         & 
         FFTW_FORWARD,  plan_type)
    call dfftw_plan_many_dft(z_plan(2), 1, zsize(3), zsize(1)*zsize(2), &
         a1, zsize(3), zsize(1)*zsize(2), 1,                            &
         a2, zsize(3), zsize(1)*zsize(2), 1,                            &
         FFTW_BACKWARD, plan_type)
#else
    call sfftw_plan_many_dft(z_plan(1), 1, zsize(3), zsize(1)*zsize(2), &
         a1, zsize(3), zsize(1)*zsize(2), 1,                            &
         a2, zsize(3), zsize(1)*zsize(2), 1,                            &
         FFTW_FORWARD,  plan_type)
    call sfftw_plan_many_dft(z_plan(2), 1, zsize(3), zsize(1)*zsize(2), &
         a1, zsize(3), zsize(1)*zsize(2), 1,                            &
         a2, zsize(3), zsize(1)*zsize(2), 1,                            &
         FFTW_BACKWARD, plan_type)
#endif
    
    deallocate(a1,a2)
    
    initialised = .true.
    
    return
  end subroutine fft_init_arg
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Final clean up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_finalize
    
    implicit none
    
#ifdef DOUBLE_PREC
    call dfftw_destroy_plan(x_plan(1))
    call dfftw_destroy_plan(x_plan(2))
    call dfftw_destroy_plan(y_plan(1))
    call dfftw_destroy_plan(y_plan(2))
    call dfftw_destroy_plan(z_plan(1))
    call dfftw_destroy_plan(z_plan(2))
#else
    call sfftw_destroy_plan(x_plan(1))
    call sfftw_destroy_plan(x_plan(2))
    call sfftw_destroy_plan(y_plan(1))
    call sfftw_destroy_plan(y_plan(2))
    call sfftw_destroy_plan(z_plan(1))
    call sfftw_destroy_plan(z_plan(2))
#endif

    call decomp_info_finalize(fft)
    
    return
  end subroutine decomp_2d_fft_finalize

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the size, starting/ending index of the distributed array 
  !  whose global size is (nx/2+1)*ny*nz, for defining data structures
  !  in r2c and c2r interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_get_size(istart, iend, isize, pencil)
    
    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize
    integer :: pencil
    
    if (pencil == 3) then
       istart = fft%zst
       iend   = fft%zen
       isize  = fft%zsz
    else if (pencil == 2) then
       istart = fft%yst
       iend   = fft%yen
       isize  = fft%ysz
    else if (pencil == 1) then
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
    integer*8 :: plan
    
    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
       allocate (wk1(xsize(1),xsize(2),xsize(3)))
       
#ifdef DOUBLE_PREC
       call dfftw_execute_dft(x_plan(isign), in, wk1)
#else
       call sfftw_execute_dft(x_plan(isign), in, wk1)
#endif
       
       ! ===== Swap X --> Y =====
       allocate (wk2(ysize(1),ysize(2),ysize(3)))
       allocate (wk2b(ysize(1),ysize(2),ysize(3)))
       call transpose_x_to_y(wk1,wk2)
       
       ! ===== 1D FFTs in Y =====
       do k=1, ysize(3)  ! loop through Z-planes
#ifdef DOUBLE_PREC
          call dfftw_execute_dft(y_plan(isign), wk2(1,1,k), wk2b(1,1,k))
#else
          call sfftw_execute_dft(y_plan(isign), wk2(1,1,k), wk2b(1,1,k))
#endif
       end do
       
       ! ===== Swap Y --> Z =====
       allocate (wk3(zsize(1),zsize(2),zsize(3)))
       call transpose_y_to_z(wk2b,wk3)
       
       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call dfftw_execute_dft(z_plan(isign), wk3, out)
#else
       call sfftw_execute_dft(z_plan(isign), wk3, out)
#endif
       
    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(zsize(1),zsize(2),zsize(3)))
       
#ifdef DOUBLE_PREC
       call dfftw_execute_dft(z_plan(isign), in, wk1)
#else
       call sfftw_execute_dft(z_plan(isign), in, wk1)
#endif
       
       ! ===== Swap Z --> Y =====
       allocate (wk2(ysize(1),ysize(2),ysize(3)))
       allocate (wk2b(ysize(1),ysize(2),ysize(3)))
       call transpose_z_to_y(wk1,wk2)
       
       ! ===== 1D FFTs in Y =====
       do k=1, ysize(3)  ! loop through Z-planes
#ifdef DOUBLE_PREC
          call dfftw_execute_dft(y_plan(isign), wk2(1,1,k), wk2b(1,1,k))
#else
          call sfftw_execute_dft(y_plan(isign), wk2(1,1,k), wk2b(1,1,k))
#endif
       end do
       
       ! ===== Swap Y --> X =====
       allocate (wk3(xsize(1),xsize(2),xsize(3)))
       call transpose_y_to_x(wk2b,wk3)
       
       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call dfftw_execute_dft(x_plan(isign), wk3, out)
#else
       call sfftw_execute_dft(x_plan(isign), wk3, out)
#endif      
       
    end if
    
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
    integer :: k
    integer*8 :: plan


    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft_r2c(plan,& 
            1,                           & ! 1D FFTs
            xsize(1),                    & ! FFT size
            xsize(2)*xsize(3),           & ! Number of 1D FFTs to compute
            in_r,                        & ! input 3D array 
            xsize(1),                    & ! 'nembed'
            1,                           & ! 1D data along X continuous in 
                                           !  memory ==> stride = 1
            xsize(1),                    & ! starting point of each 1D 
                                           !  FFT separated by this
            wk1,                         & ! output 3D array
            fft%xsz(1),                  & ! 'nembed'
            1,                           & ! stride
            fft%xsz(1),                  & ! dist 
            FFTW_ESTIMATE)     
       call dfftw_execute_dft_r2c(plan,in_r,wk1)
#else
       call sfftw_plan_many_dft_r2c(plan,1,xsize(1),xsize(2)*xsize(3), &
            in_r,xsize(1),1,xsize(1),wk1,fft%xsz(1),1,fft%xsz(1), &
            FFTW_ESTIMATE)     
       call sfftw_execute_dft_r2c(plan,in_r,wk1)
#endif

       ! ===== Swap X --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_x_to_y(wk1,wk2,fft)

       ! ===== Swap Y --> Z =====
       allocate (wk3(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       call transpose_y_to_z(wk2,wk3,fft)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft(plan,  &
            1,                         & ! 1D FFTs
            fft%zsz(3),                & ! FFT size
            fft%zsz(1)*fft%zsz(2),     & ! number of 1D FFTs to compute
            wk3,                       & ! input 3D array
            fft%zsz(3),                & ! input 'nembed'
            fft%zsz(1)*fft%zsz(2),     & ! 1D data point along z-axis 
                                         !  separated by the size of xy 
                                         !  plane in memory 
            1,                         & ! starting point of each 1D FFT 
                                         !  next to each other in memory
            out_c,                     & ! output array, same as input
            fft%zsz(3),                &
            fft%zsz(1)*fft%zsz(2),     &
            1,                         &
            FFTW_FORWARD,  FFTW_ESTIMATE)
       call dfftw_execute_dft(plan,wk3,out_c)
#else
       call sfftw_plan_many_dft(plan,1,fft%zsz(3),fft%zsz(1)*fft%zsz(2), &
            wk3,fft%zsz(3),fft%zsz(1)*fft%zsz(2),1, &
            out_c,fft%zsz(3), fft%zsz(1)*fft%zsz(2),1, &
            FFTW_FORWARD,  FFTW_ESTIMATE)
       call sfftw_execute_dft(plan,wk3,out_c)
#endif

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       allocate(wk1(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft_r2c(plan,& 
            1,                           & ! 1D FFTs
            zsize(3),                    & ! FFT size
            zsize(1)*zsize(2),           & ! Number of 1D FFTs to compute
            in_r,                        & ! input 3D array 
            zsize(3),                    & ! 'nembed'
            zsize(1)*zsize(2),           & ! 1D data along z-axis 
                                           !  separated by the size of xy
                                           !  plane in memory 
            1,                           & ! starting point of each 1D FFT 
                                           !  next to each other in memory
            wk1,                         & ! output 3D array
            fft%zsz(3),                  & ! 'nembed'
            fft%zsz(1)*fft%zsz(2),       & ! stride
            1,                           & ! dist 
            FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan,in_r,wk1)
#else
       call sfftw_plan_many_dft_r2c(plan,1,zsize(3),zsize(1)*zsize(2), &
            in_r,zsize(3),zsize(1)*zsize(2),1, &
            wk1,fft%zsz(3),fft%zsz(1)*fft%zsz(2),1, &
            FFTW_ESTIMATE)
       call sfftw_execute_dft_r2c(plan,in_r,wk1)
#endif

       ! ===== Swap Z --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_z_to_y(wk1,wk2,fft)

       ! ===== Swap Y --> X =====
       allocate (wk3(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
       call transpose_y_to_x(wk2,wk3,fft)

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft(plan,  &
            1,                         & ! 1D FFTs
            fft%xsz(1),                & ! FFT size
            fft%xsz(2)*fft%xsz(3),     & ! number of 1D FFTs to compute
            wk3,                       & ! input 3D array
            fft%xsz(1),                & ! input 'nembed'
            1,                         & ! 1D data along X continuous in 
                                         !  memory ==> stride = 1
            fft%xsz(1),                & ! starting point of each 1D 
                                         !  FFT separated by this
            out_c,                     & ! output array, same as input
            fft%xsz(1),                &
            1,                         &
            fft%xsz(1),                &
            FFTW_FORWARD,  FFTW_ESTIMATE)
       call dfftw_execute_dft(plan,wk3,out_c)
#else
       call sfftw_plan_many_dft(plan,1,fft%xsz(1),fft%xsz(2)*fft%xsz(3),&
            wk3,  fft%xsz(1),1,fft%xsz(1), &
            out_c,fft%xsz(1),1,fft%xsz(1), &
            FFTW_FORWARD,  FFTW_ESTIMATE)
       call sfftw_execute_dft(plan,wk3,out_c)
#endif

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
    integer :: k
    integer*8 :: plan

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
       allocate(wk1(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft(plan,&
            1,                       & ! 1D FFTs
            fft%zsz(3),              & ! FFT size
            fft%zsz(1)*fft%zsz(2),   & ! number of 1D FFTs to compute
            in_c,                    & ! input 3D array
            fft%zsz(3),              & ! input 'nembed'
            fft%zsz(1)*fft%zsz(2),   & ! 1D data point along z-axis 
                                       !  separated by the size of xy 
                                       !  plane in memory 
            1,                       & ! starting point of each 1D FFT
                                       !  next to each other in memory
            wk1,                     & ! output array, same as input
            fft%zsz(3),              & ! 'nembed'
            fft%zsz(1)*fft%zsz(2),   & ! stride
            1,                       & ! dist   
            FFTW_BACKWARD, FFTW_ESTIMATE)
       call dfftw_execute_dft(plan,in_c,wk1)
#else
       call sfftw_plan_many_dft(plan,1,fft%zsz(3),fft%zsz(1)*fft%zsz(2),&
            in_c,fft%zsz(3),fft%zsz(1)*fft%zsz(2),1, &
            wk1, fft%zsz(3),fft%zsz(1)*fft%zsz(2),1, &
            FFTW_BACKWARD, FFTW_ESTIMATE)
       call sfftw_execute_dft(plan,in_c,wk1)
#endif

       ! ===== Swap Z --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_z_to_y(wk1,wk2,fft)

       ! ===== Swap Y --> X =====
       allocate (wk3(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
       call transpose_y_to_x(wk2,wk3,fft)
       
       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft_c2r(plan,& 
            1,                           & ! 1D FFTs
            xsize(1),                    & ! FFT size
            xsize(2)*xsize(3),           & ! Number of 1D FFTs to compute
            wk3,                         & ! input 3D array 
            fft%xsz(1),                  & ! 'nembed'
            1,                           & ! 1D data along X continuous in 
                                           !  memory ==> stride = 1
            fft%xsz(1),                  & ! starting point of each 1D 
                                           !  FFT separated by this
            out_r,                       & ! output 3D array
            xsize(1),                    & ! 'nembed'
            1,                           & ! stride
            xsize(1),                    & ! dist 
            FFTW_ESTIMATE)     
       call dfftw_execute_dft_c2r(plan,wk3,out_r)
#else
       call sfftw_plan_many_dft_c2r(plan,1,xsize(1),xsize(2)*xsize(3), &
            wk3,fft%xsz(1),1,fft%xsz(1), &
            out_r,xsize(1),1,xsize(1),   &
            FFTW_ESTIMATE)     
       call sfftw_execute_dft_c2r(plan,wk3,out_r)
#endif

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft(plan,  &
            1,                         & ! 1D FFTs
            fft%xsz(1),                & ! FFT size
            fft%xsz(2)*fft%xsz(3),     & ! number of 1D FFTs to compute
            in_c,                      & ! input 3D array
            fft%xsz(1),                & ! input 'nembed'
            1,                         & ! 1D data along X continuous in 
                                         !  memory ==> stride = 1
            fft%xsz(1),                & ! starting point of each 1D 
                                         !  FFT separated by this
            wk1,                      & ! output array, same as input
            fft%xsz(1),                &
            1,                         &
            fft%xsz(1),                &
            FFTW_BACKWARD,  FFTW_ESTIMATE)
       call dfftw_execute_dft(plan,in_c,wk1)
#else
       call sfftw_plan_many_dft(plan,1,fft%xsz(1),fft%xsz(2)*fft%xsz(3),&
            in_c,fft%xsz(1),1,fft%xsz(1), &
            wk1, fft%xsz(1),1,fft%xsz(1), &
            FFTW_BACKWARD,  FFTW_ESTIMATE)
       call sfftw_execute_dft(plan,in_c,wk1)
#endif

       ! ===== Swap X --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_x_to_y(wk1,wk2,fft)

       ! ===== Swap Y --> Z =====
       allocate (wk3(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       call transpose_y_to_z(wk2,wk3,fft)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call dfftw_plan_many_dft_c2r(plan,& 
            1,                           & ! 1D FFTs
            zsize(3),                    & ! FFT size
            zsize(1)*zsize(2),           & ! Number of 1D FFTs to compute
            wk3,                         & ! input 3D array 
            fft%zsz(3),                  & ! 'nembed'
            fft%zsz(1)*fft%zsz(2),       & ! 1D data along z-axis 
                                           !  separated by the size of xy
                                           !  plane in memory 
            1,                           & ! starting point of each 1D FFT 
                                           !  next to each other in memory
            out_r,                       & ! output 3D array
            zsize(3),                    & ! 'nembed'
            zsize(1)*zsize(2),           & ! stride
            1,                           & ! dist 
            FFTW_ESTIMATE)
       call dfftw_execute_dft_c2r(plan,wk3,out_r)
#else
       call sfftw_plan_many_dft_c2r(plan,1,zsize(3),zsize(1)*zsize(2), &
            wk3,fft%zsz(3),fft%zsz(1)*fft%zsz(2),1, &
            out_r,zsize(3),zsize(1)*zsize(2),1, &
            FFTW_ESTIMATE)
       call sfftw_execute_dft_c2r(plan,wk3,out_r)
#endif

    end if

    return
  end subroutine fft_3d_c2r

  
end module decomp_2d_fft
