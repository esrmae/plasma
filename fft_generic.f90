module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  use glassman   ! A simple FFT library distributed with this software
  
  implicit none
  
  private        ! Make everything private unless declared public
  
  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1
  
  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: PHYSICAL_IN_Z = 3 

  integer, save :: format                 ! input X-pencil or Z-pencil
  
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

  ! work space
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  
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
    integer :: cbuf_size
    
    format = pencil

    if (format==PHYSICAL_IN_X) then
       call decomp_info_init(nx_global/2+1, ny_global, nz_global, fft)
    else if (format==PHYSICAL_IN_Z) then
       call decomp_info_init(nx_global, ny_global, nz_global/2+1, fft)
    end if

    cbuf_size = max(xsize(1), ysize(2))
    cbuf_size = max(cbuf_size, zsize(3))
    allocate(buf(cbuf_size))
    allocate(scratch(cbuf_size))
    
    initialised = .true.
    
    return
  end subroutine fft_init_arg
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Final clean up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_finalize
    
    implicit none

    deallocate(buf,scratch)

    call decomp_info_finalize(fft)

    return
  end subroutine decomp_2d_fft_finalize

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the size, starting/ending index of the distributed array 
  !  whose global size is (nx/2+1)*ny*nz, for defining data structures
  !  in r2c and c2r interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_get_size(istart, iend, isize, dir)
    
    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize
    integer, intent(IN) :: dir
    
    if (dir==3) then
       istart = fft%zst
       iend   = fft%zen
       isize  = fft%zsz
    else if (dir==2) then
       istart = fft%yst
       iend   = fft%yen
       isize  = fft%ysz
    else if (dir==1) then
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

    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2
    integer :: i,j,k

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
       allocate (wk1(xsize(1),xsize(2),xsize(3)))
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                buf(i) = in(i,j,k)
             enddo
             call spcfft(buf,xsize(1),isign,scratch)
             do i=1,xsize(1)
                wk1(i,j,k) = buf(i)
             end do
          end do
       end do

       ! ===== Swap X --> Y =====
       allocate (wk2(ysize(1),ysize(2),ysize(3)))
       call transpose_x_to_y(wk1,wk2)
       
       ! ===== 1D FFTs in Y =====
       do k=1,ysize(3)
          do i=1,ysize(1)
             do j=1,ysize(2)
                buf(j) = wk2(i,j,k)
             end do
             call spcfft(buf,ysize(2),isign,scratch)
             do j=1,ysize(2)
                wk2(i,j,k) = buf(j)
             end do
          end do
       end do

       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2,out)

       ! ===== 1D FFTs in Z =====
       do j=1,zsize(2)
          do i=1,zsize(1)
             do k=1,zsize(3)
                buf(k) = out(i,j,k)
             end do
             call spcfft(buf,zsize(3),isign,scratch)
             do k=1,zsize(3)
                out(i,j,k) = buf(k)
             end do
          end do
       end do

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(zsize(1),zsize(2),zsize(3)))
       do j=1,zsize(2)
          do i=1,zsize(1)
             do k=1,zsize(3)
                buf(k) = in(i,j,k)
             end do
             call spcfft(buf,zsize(3),isign,scratch)
             do k=1,zsize(3)
                wk1(i,j,k) = buf(k)
             end do
          end do
       end do
       
       ! ===== Swap Z --> Y =====
       allocate (wk2(ysize(1),ysize(2),ysize(3)))
       call transpose_z_to_y(wk1,wk2)
       
       ! ===== 1D FFTs in Y =====
       do k=1,ysize(3)
          do i=1,ysize(1)
             do j=1,ysize(2)
                buf(j) = wk2(i,j,k)
             end do
             call spcfft(buf,ysize(2),isign,scratch)
             do j=1,ysize(2)
                wk2(i,j,k) = buf(j)
             end do
          end do
       end do
       
       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2,out)
       
       ! ===== 1D FFTs in X =====
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                buf(i) = out(i,j,k)
             enddo
             call spcfft(buf,xsize(1),isign,scratch)
             do i=1,xsize(1)
                out(i,j,k) = buf(i)
             end do
          end do
       end do
       
    end if

    return
  end subroutine fft_3d_c2c
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2
    integer :: i,j,k

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(fft%xsz(1),fft%xsz(2),fft%xsz(3))) 

       do k=1,xsize(3)
          do j=1,xsize(2)
             ! Glassman's FFT is c2c only, 
             ! needing some pre- and post-processing for r2c
             ! pack real input in complex storage
             do i=1,xsize(1)
                buf(i) = cmplx(in_r(i,j,k),0._mytype, kind=mytype)
             end do
             call spcfft(buf,xsize(1),-1,scratch)
             ! simply drop the redundant part of the complex output
             do i=1,fft%xsz(1)
                wk1(i,j,k) = buf(i)
             end do
          end do
       end do

       ! ===== Swap X --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_x_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       !do k=1,fft%ysz(3)
       !   do i=1,fft%ysz(1)
       !      do j=1,fft%ysz(2)
       !         buf(j) = wk2(i,j,k)
       !      end do
       !      call spcfft(buf,fft%ysz(2),-1,scratch)
       !      do j=1,fft%ysz(2)
       !         wk2(i,j,k) = buf(j)
       !      end do
       !   end do
       !end do

       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2,out_c,fft)

       ! ===== 1D FFTs in Z =====
       do j=1,fft%zsz(2)
          do i=1,fft%zsz(1)
             do k=1,fft%zsz(3)
                buf(k) = out_c(i,j,k)
             end do
             call spcfft(buf,fft%zsz(3),-1,scratch)
             do k=1,fft%zsz(3)
                out_c(i,j,k) = buf(k)
             end do
          end do
       end do
                
    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       allocate(wk1(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       
       do j=1,zsize(2)
          do i=1,zsize(1)
             do k=1,zsize(3)
                buf(k) = cmplx(in_r(i,j,k),0._mytype, kind=mytype)
             end do
             call spcfft(buf,zsize(3),-1,scratch)
             do k=1,fft%zsz(3)
                wk1(i,j,k) = buf(k)
             end do
          end do
       end do

       ! ===== Swap Z --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_z_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       do k=1,fft%ysz(3)
          do i=1,fft%ysz(1)
             do j=1,fft%ysz(2)
                buf(j) = wk2(i,j,k)
             end do
             call spcfft(buf,fft%ysz(2),-1,scratch)
             do j=1,fft%ysz(2)
                wk2(i,j,k) = buf(j)
             end do
          end do
       end do

       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2,out_c,fft)

       ! ===== 1D FFTs in X =====
       do k=1,fft%xsz(3)
          do j=1,fft%xsz(2)
             do i=1,fft%xsz(1)
                buf(i) = out_c(i,j,k)
             end do
             call spcfft(buf,fft%xsz(1),-1,scratch)
             do i=1,fft%xsz(1)
                out_c(i,j,k) = buf(i)
             end do
          end do
       end do

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
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2, wk3
    integer :: i,j,k

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
       allocate(wk1(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       do j=1,fft%zsz(2)
          do i=1,fft%zsz(1)
             do k=1,fft%zsz(3)
                buf(k) = in_c(i,j,k)
             end do
             call spcfft(buf,fft%zsz(3),1,scratch)
             do k=1,fft%zsz(3)
                wk1(i,j,k) = buf(k)
             end do
          end do
       end do

       ! ===== Swap Z --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_z_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       !do k=1,fft%ysz(3)
       !   do i=1,fft%ysz(1)
       !      do j=1,fft%ysz(2)
       !         buf(j) = wk2(i,j,k)
       !      end do
       !      call spcfft(buf,fft%ysz(2),1,scratch)
       !      do j=1,fft%ysz(2)
       !         wk2(i,j,k) = buf(j)
       !      end do
       !   end do
       !end do

       ! ===== Swap Y --> X =====
       allocate (wk3(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
       call transpose_y_to_x(wk2,wk3,fft)

       ! ===== 1D FFTs in X =====
       do k=1,xsize(3)
          do j=1,xsize(2)
             ! Glassman's FFT is c2c only, 
             ! needing some pre- and post-processing for c2r
             do i=1,xsize(1)/2+1
                buf(i) = wk3(i,j,k)
             end do
             ! expanding to a full-size complex array
             ! For odd N, the storage is:
             !  1, 2, ...... N/2+1   integer division rounded down
             !     N, ...... N/2+2   => a(i) is conjugate of a(N+2-i)
             ! For even N, the storage is:
             !  1, 2, ...... N/2  , N/2+1
             !     N, ...... N/2+2  again a(i) conjugate of a(N+2-i)
             do i=xsize(1)/2+2,xsize(1)
                buf(i) =  conjg(buf(xsize(1)+2-i))
             end do
             call spcfft(buf,xsize(1),1,scratch)
             do i=1,xsize(1)
                ! simply drop imaginary part
                out_r(i,j,k) = real(buf(i), kind=mytype)
             end do
          end do
       end do

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(fft%xsz(1),fft%xsz(2),fft%xsz(3)))
       do k=1,fft%xsz(3)
          do j=1,fft%xsz(2)
             do i=1,fft%xsz(1)
                buf(i) = in_c(i,j,k)
             end do
             call spcfft(buf,fft%xsz(1),1,scratch)
             do i=1,fft%xsz(1)
                wk1(i,j,k) = buf(i)
             end do
          end do
       end do

       ! ===== Swap X --> Y =====
       allocate (wk2(fft%ysz(1),fft%ysz(2),fft%ysz(3)))
       call transpose_x_to_y(wk1,wk2,fft)

       ! ===== 1D FFTs in Y =====
       do k=1,fft%ysz(3)
          do i=1,fft%ysz(1)
             do j=1,fft%ysz(2)
                buf(j) = wk2(i,j,k)
             end do
             call spcfft(buf,fft%ysz(2),1,scratch)
             do j=1,fft%ysz(2)
                wk2(i,j,k) = buf(j)
             end do
          end do
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(fft%zsz(1),fft%zsz(2),fft%zsz(3)))
       call transpose_y_to_z(wk2,wk3,fft)

       ! ===== 1D FFTs in Z =====
       do j=1,zsize(2)
          do i=1,zsize(1)
             do k=1,zsize(3)/2+1
                buf(k) = wk3(i,j,k)
             end do
             do k=zsize(3)/2+2,zsize(3)
                buf(k) =  conjg(buf(zsize(3)+2-k))
             end do
             call spcfft(buf,zsize(3),1,scratch)
             do k=1,zsize(3)
                ! simply drop imaginary part
                out_r(i,j,k) = real(buf(k), kind=mytype)
             end do
          end do
       end do

    end if

    return
  end subroutine fft_3d_c2r

  
end module decomp_2d_fft
