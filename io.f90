module decomp_2d_io

  use decomp_2d
  use MPI

  implicit none

  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one,decomp_2d_write_var,decomp_2d_write_var_fft, &
 decomp_2d_read_var,decomp_2d_write_head3d,decomp_2d_read_head3d,decomp_2d_read_head3dold, &
decomp_2d_write_simphead,decomp_2d_write_simphead_int,decomp_2d_write_simphead_real, &
decomp_2d_write_1dx,decomp_2d_write_2dx, decomp_2d_write_1dx_oli,decomp_2d_write_2dx_oli

  ! Generic interface to handle multiple data types
  interface decomp_2d_write_one
     module procedure mpiio_write_real
     module procedure mpiio_write_complex
     module procedure mpiio_write_real_decomp
  end interface
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_real(nx,ny,nz,ipencil,var,filename)
    
    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh
    
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_real


  subroutine mpiio_write_real_decomp(decomp,ipencil,var,filename)
    
    implicit none
    
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh
    
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_real_decomp


  subroutine mpiio_write_complex(nx,ny,nz,ipencil,var,filename)
    
    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh
    
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, complex_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         complex_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_complex


  subroutine mpiio_write_complex_decomp(decomp,ipencil,var,filename)
    
    implicit none
    
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh
    
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, complex_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         complex_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_complex_decomp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array to a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_var(fh,disp,n1,n2,n3,ipencil,var)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n1,n2,n3       ! global size
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(IN) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, newtype

    ! Create file type and set file view
    sizes(1) = nx+4
    sizes(2) = ny+4
    sizes(3) = nz+4
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
       Do i=1,3
         If (xmin(i)==1) then
           subsizes(i)=subsizes(i)+2
         Else  
           starts(i)=starts(i)+2
         End If
         If (xmax(i)==1) subsizes(i)=subsizes(i)+2
       End Do
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
       Do i=1,3
         If (ymin(i)==1) then
           subsizes(i)=subsizes(i)+2
         Else  
           starts(i)=starts(i)+2
         End If
         If (ymax(i)==1) subsizes(i)=subsizes(i)+2
       End Do
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
       Do i=1,3
         If (zmin(i)==1) then
           subsizes(i)=subsizes(i)+2
         Else  
           starts(i)=starts(i)+2
         End If
         If (zmax(i)==1) subsizes(i)=subsizes(i)+2
       End Do
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + (nx+4)*(ny+4)*(nz+4)*mytype

    return
  end subroutine decomp_2d_write_var

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array to a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_var_fft(fh,disp,n1,n2,n3,ipencil,var)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n1,n2,n3       ! global size
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(IN) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, newtype, nx1,ny1,nz1

    ! Create file type and set file view
    nx1=nx/2+1; ny1=ny; nz1=nz
    sizes(1) = nx1
    sizes(2) = ny1
    sizes(3) = nz1
    if (ipencil == 3) then
       subsizes(1) = fft_zsize(1)
       subsizes(2) = fft_zsize(2)
       subsizes(3) = fft_zsize(3)
       starts(1) = fft_zstart(1)-1
       starts(2) = fft_zstart(2)-1
       starts(3) = fft_zstart(3)-1
    else
      write(*,*) 'ipencil=',ipencil,'is not right in decomp_2d_write_var_fft'
    endif  
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + nx1*ny1*nz1*mytype

    return
  end subroutine decomp_2d_write_var_fft

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read a 3D array from a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next reading.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_read_var(fh,disp,n1,n2,n3,ipencil,var)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n1,n2,n3       ! global size
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(INOUT) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k,newtype

    !To read from no MPI file
    If (restrt_old==1) disp = disp + mytype

    ! Create file type and set file view
    sizes(1) = nx+4
    sizes(2) = ny+4
    sizes(3) = nz+4
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
       Do i=1,3
         If (xmin(i)==1) then
           subsizes(i)=subsizes(i)+2
         Else  
           starts(i)=starts(i)+2
         End If
         If (xmax(i)==1) subsizes(i)=subsizes(i)+2
       End Do
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
       Do i=1,3
         If (ymin(i)==1) then
           subsizes(i)=subsizes(i)+2
         Else  
           starts(i)=starts(i)+2
         End If
         If (ymax(i)==1) subsizes(i)=subsizes(i)+2
       End Do
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
       Do i=1,3
         If (zmin(i)==1) then
           subsizes(i)=subsizes(i)+2
         Else  
           starts(i)=starts(i)+2
         End If
         If (zmax(i)==1) subsizes(i)=subsizes(i)+2
       End Do
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, & 
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type,MPI_STATUS_IGNORE, ierror) !
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + (nx+4)*(ny+4)*(nz+4)*mytype

    return
  end subroutine decomp_2d_read_var


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D header for a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_head3d(fh,disp)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer :: i,j,k, intype
    
    ! write integer part
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    If (nrank == 0) then
      call MPI_FILE_WRITE(fh,nx,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,ny,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,nz,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
    End If

    ! update displacement for next write
    call MPI_SIZEOF(nx,intype,ierror)
    disp = disp + 3*intype

    ! write real part
    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
    If (nrank == 0) then
      call MPI_FILE_WRITE(fh,length,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,height,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,width,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,re,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,dt,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,time,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,dpdx_mean,1,real_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_WRITE(fh,cony,1,real_type,MPI_STATUS_IGNORE,ierror)
    End If

    ! update displacement for next write
    disp = disp + 8*mytype

    return
  end subroutine decomp_2d_write_head3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D header for a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_read_head3d(fh,disp,nxo,nyo,nzo,leng,hei,wid,reo,conyo)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(OUT) :: nxo,nyo,nzo
    real(mytype), intent(OUT)   :: leng,hei,wid,reo,conyo

    integer :: i,j,k, intype
    call MPI_SIZEOF(nx,intype,ierror)
    ! read integer part
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ(fh,nxo,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,nyo,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,nzo,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
    disp = disp + 3*intype

    ! read real part
    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ(fh,leng,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,hei,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,wid,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,reo,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,dt,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,time,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,dpdx_mean,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,conyo,1,real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next read
    disp = disp + 8*mytype

    return
  end subroutine decomp_2d_read_head3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Read a 3D header for old restrt file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_read_head3dold(fh,disp,nxo,nyo,nzo,leng,hei,wid,reo)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(OUT) :: nxo,nyo,nzo
    real(mytype), intent(OUT)   :: leng,hei,wid,reo

    integer :: i,j,k, intype
    call MPI_SIZEOF(nx,intype,ierror)
    disp = disp + intype

    ! read real part
    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ(fh,leng,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,hei,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,wid,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,reo,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,dt,1,real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,time,1,real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next read
    disp = disp + 6*mytype

    ! read integer part
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ(fh,nxo,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,nyo,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_READ(fh,nzo,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
    disp = disp + 3*intype

    ! read real part
    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ(fh,dpdx_mean,1,real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next read
    disp = disp + mytype

    return
  end subroutine decomp_2d_read_head3dold

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a simple header for a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_simphead(fh,disp)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer :: i,j,k !, intype
    
    ! write integer part
    !call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    !call MPI_FILE_WRITE(fh,ii,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
    !call MPI_SIZEOF(ii,intype,ierror)
    !disp = disp + intype

    ! write real part
    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE(fh,time,1,real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
    disp = disp + mytype

    return
  end subroutine decomp_2d_write_simphead

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a simple header for a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_simphead_int(fh,navg,disp)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh,navg             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer :: i,j,k,intype
    
    ! write integer part
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE(fh,navg,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
    call MPI_SIZEOF(navg,intype,ierror)
    disp = disp + intype

    ! write real part
!    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
!    call MPI_FILE_WRITE(fh,time,1,real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
!    disp = disp + mytype

    return
  end subroutine decomp_2d_write_simphead_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a simple header for a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_simphead_real(fh,realnum,disp)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    real(mytype), intent(IN) :: realnum
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer :: i,j,k

    ! write real part
    call MPI_FILE_SET_VIEW(fh,disp,real_type,real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE(fh,realnum,1,real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for next write
    disp = disp + mytype

    return
  end subroutine decomp_2d_write_simphead_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Write a 1D variable to a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_1dx(fh,disp,n,dir,var,arb_write)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n	      ! global size
    integer, intent(IN) :: dir        ! x,y or z direction
    real(mytype), dimension(:), &
         intent(IN) :: var                ! distributed 1D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(1)  :: sizes, subsizes, starts
    integer :: newtype, writers, arb_write !i,j,k, ierror
    If (arb_write==1) then; writers = 0
    Else If (arb_write==0) then; writers = 1; End If

    ! Create file type and set file view
    if (dir == 1) then
       sizes = n+4
       subsizes = xsize(1)
       starts = xstart(1)-1  ! 0-based index
       if ((xmin(2)==1).and.(xmin(3)==1)) writers = 1
       if (xmin(1)==1) then
          subsizes=subsizes+2
       else
          starts=starts+2
       end if
       if (xmax(1)==1) subsizes=subsizes+2
    else if (dir == 2) then
       sizes = n+4
       subsizes = xsize(2)
       starts = xstart(2)-1  
       if ((xmin(1)==1).and.(xmin(3)==1)) writers = 1
       if (xmin(2)==1) then
          subsizes=subsizes+2
       else 
          starts=starts+2
       end if
       if (xmax(2)==1) subsizes=subsizes+2
    else if (dir == 3) then
       sizes = n+1
       subsizes = xsize(3)
       starts = xstart(3)-1
       if ((xmin(1)==1).and.(xmin(2)==1)) writers = 1
       if (xmin(3)==1) then
          subsizes=subsizes+2
       else
          starts=starts+2
       end if
       if (xmax(3)==1) subsizes=subsizes+2
    endif

      call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    if (writers == 1) then
      call MPI_FILE_WRITE(fh, var, & 
         subsizes,real_type, MPI_STATUS_IGNORE, ierror)
    end if
      call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + (n+4)*mytype

    return
  end subroutine decomp_2d_write_1dx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Write a 2D variable to a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_2dx(fh,disp,n1,n2,dir,var,arb_write)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n1,n2      ! global size
    integer, intent(IN) :: dir        ! 1=xy, 2=xz, 3=yz
    real(mytype), dimension(:,:), &
         intent(IN) :: var                ! distributed 2D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(2) :: sizes, subsizes, starts
    integer :: newtype, writers, arb_write !i,j,k, ierror
    If (arb_write==1) then; writers = 0
    Else If (arb_write==0) then; writers = 1; End If
    
    ! Create file type and set file view
    if (dir == 1) then
       sizes(1) = n1+4
       sizes(2) = n2+4
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       if (xmin(3) == 1) writers = 1
       if (xmin(1)==1) then
         subsizes(1)=subsizes(1)+2
       else
         starts(1)=starts(1)+2
       end if
       if (xmax(1)==1) subsizes(1)=subsizes(1)+2
       if (xmin(2)==1) then
         subsizes(2)=subsizes(2)+2
       else
         starts(2)=starts(2)+2
       end if
       if (xmax(2)==1) subsizes(2)=subsizes(2)+2
    else if (dir == 2) then
       sizes(1) = n1+4
       sizes(2) = n2+4
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(3)
       starts(1) = xstart(1)-1 
       starts(2) = xstart(3)-1
       if (xmin(2) == 1) writers = 1
       if (xmin(1)==1) then
         subsizes(1)=subsizes(1)+2
       else
         starts(1)=starts(1)+2
       end if
       if (xmax(1)==1) subsizes(1)=subsizes(1)+2
       if (xmin(3)==1) then
         subsizes(2)=subsizes(2)+2
       else
         starts(2)=starts(2)+2
       end if
       if (xmax(3)==1) subsizes(2)=subsizes(2)+2
    else if (dir == 3) then
       sizes(1) = n1+4
       sizes(2) = n2+4
       subsizes(1) = xsize(2)
       subsizes(2) = xsize(3)
       starts(1) = xstart(2)-1 
       starts(2) = xstart(3)-1
       if (xmin(1) == 1) writers = 1
       if (xmin(2)==1) then
         subsizes(1)=subsizes(1)+2
       else
         starts(1)=starts(1)+2
       end if
       if (xmax(2)==1) subsizes(1)=subsizes(1)+2
       if (xmin(3)==1) then
         subsizes(2)=subsizes(2)+2
       else
         starts(2)=starts(2)+2
       end if
       if (xmax(3)==1) subsizes(2)=subsizes(2)+2
    endif
   
      call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    if (writers == 1) then 
      call MPI_FILE_WRITE(fh, var, subsizes(1)*subsizes(2), &
         real_type, MPI_STATUS_IGNORE, ierror)
    end if
      call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + (n1+4)*(n2+4)*mytype

    return
  end subroutine decomp_2d_write_2dx  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Write a 1D variable to a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_1dx_oli(fh,disp,n,dir,var,arb_write)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n	      ! global size
    integer, intent(IN) :: dir        ! x,y or z direction
    real(mytype), dimension(:), &
         intent(IN) :: var                ! distributed 1D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(1)  :: sizes, subsizes, starts
    integer :: newtype, writers, arb_write !i,j,k, ierror
!    If (arb_write==1) then; writers = 0
!    Else If (arb_write==0) then; writers = 1; End If
    writers=arb_write
    ! Create file type and set file view
    if (dir == 1) then
       sizes = n+4
       subsizes = xsize(1)
       starts = xstart(1)-1  ! 0-based index
!       if ((xmin(2)==1).and.(xmin(3)==1)) writers = 1
       if (xmin(1)==1) then
          subsizes=subsizes+2
       else
          starts=starts+2
       end if
       if (xmax(1)==1) subsizes=subsizes+2
    else if (dir == 2) then
       sizes = n+4
       subsizes = xsize(2)
       starts = xstart(2)-1  
!       if ((xmin(1)==1).and.(xmin(3)==1)) writers = 1
       if (xmin(2)==1) then
          subsizes=subsizes+2
       else 
          starts=starts+2
       end if
       if (xmax(2)==1) subsizes=subsizes+2
    else if (dir == 3) then
       sizes = n+4
       subsizes = xsize(3)
       starts = xstart(3)-1
!       if ((xmin(1)==1).and.(xmin(2)==1)) writers = 1
       if (xmin(3)==1) then
          subsizes=subsizes+2
       else
          starts=starts+2
       end if
       if (xmax(3)==1) subsizes=subsizes+2
    endif

      call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    if (writers == 1) then
      call MPI_FILE_WRITE(fh, var, & 
         subsizes,real_type, MPI_STATUS_IGNORE, ierror)
    end if
      call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + (n+4)*mytype

    return
  end subroutine decomp_2d_write_1dx_oli

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Write a 2D variable to a big MPI-IO file, starting from
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_write_2dx_oli(fh,disp,n1,n2,dir,var,arb_write)
    Use shared_data
    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: n1,n2      ! global size
    integer, intent(IN) :: dir        ! 1=xy, 2=xz, 3=yz
    real(mytype), dimension(:,:), &
         intent(IN) :: var                ! distributed 2D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(2) :: sizes, subsizes, starts
    integer :: newtype, writers, arb_write !i,j,k, ierror
!    If (arb_write==1) then; writers = 0
!    Else If (arb_write==0) then; writers = 1; End If
    writers=arb_write
    ! Create file type and set file view
    if (dir == 1) then
       sizes(1) = n1+4
       sizes(2) = n2+4
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
!       if (xmin(3) == 1) writers = 1
       if (xmin(1)==1) then
         subsizes(1)=subsizes(1)+2
       else
         starts(1)=starts(1)+2
       end if
       if (xmax(1)==1) subsizes(1)=subsizes(1)+2
       if (xmin(2)==1) then
         subsizes(2)=subsizes(2)+2
       else
         starts(2)=starts(2)+2
       end if
       if (xmax(2)==1) subsizes(2)=subsizes(2)+2
    else if (dir == 2) then
       sizes(1) = n1+4
       sizes(2) = n2+4
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(3)
       starts(1) = xstart(1)-1
       starts(2) = xstart(3)-1
!       if (xmin(2) == 1) writers = 1
       if (xmin(1)==1) then
         subsizes(1)=subsizes(1)+2
       else
         starts(1)=starts(1)+2
       end if
       if (xmax(1)==1) subsizes(1)=subsizes(1)+2
       if (xmin(3)==1) then
         subsizes(2)=subsizes(2)+2
       else
         starts(2)=starts(2)+2
       end if
       if (xmax(3)==1) subsizes(2)=subsizes(2)+2
    else if (dir == 3) then
       sizes(1) = n1+4
       sizes(2) = n2+4
       subsizes(1) = xsize(2)
       subsizes(2) = xsize(3)
       starts(1) = xstart(2)-1
       starts(2) = xstart(3)-1
!       if (xmin(1) == 1) writers = 1
       if (xmin(2)==1) then
         subsizes(1)=subsizes(1)+2
       else
         starts(1)=starts(1)+2
       end if
       if (xmax(2)==1) subsizes(1)=subsizes(1)+2
       if (xmin(3)==1) then
         subsizes(2)=subsizes(2)+2
       else
         starts(2)=starts(2)+2
       end if
       if (xmax(3)==1) subsizes(2)=subsizes(2)+2
    endif

      call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    if (writers == 1) then
      call MPI_FILE_WRITE(fh, var, subsizes(1)*subsizes(2), &
         real_type, MPI_STATUS_IGNORE, ierror)
    end if
      call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + (n1+4)*(n2+4)*mytype

    return
  end subroutine decomp_2d_write_2dx_oli

end module decomp_2d_io
