!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the main 2D pencil decomposition module

module decomp_2d

  use MPI

  implicit none

  private        ! Make everything private unless declared public

#ifdef DOUBLE_PREC
  integer, parameter, public :: mytype = Selected_real_kind( 12, 70 )
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#else
  integer, parameter, public :: mytype = Selected_real_kind( 6, 35 )
  integer, parameter, public :: real_type = MPI_REAL
  integer, parameter, public :: complex_type = MPI_COMPLEX
#endif

  ! some key global variables
  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: nrank  ! local MPI rank 
  integer, save, public :: nproc  ! total number of processors

  ! parameters for 2D Catersian topology 
  integer, save, public, dimension(2) :: dims, coord
  logical, save, dimension(2) :: periodic
  integer, save, public :: MPI_COMM_CART 
  integer, save, public :: MPI_COMM_ROW, MPI_COMM_COL

  ! define neighboring blocks (to be used in halo-cell support)
  !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, dimension(3,6) :: neighbour 

#ifdef SHM
  ! derived type to store shared-memory info
  TYPE, public :: SMP_INFO
     integer MPI_COMM          ! SMP associated with this communicator
     integer NODE_ME           ! rank in this communicator
     integer NCPU              ! size of this communicator
     integer SMP_COMM          ! communicator for SMP-node masters
     integer CORE_COMM         ! communicator for cores on SMP-node
     integer SMP_ME            ! SMP-node id starting from 1 ... NSMP
     integer NSMP              ! number of SMP-nodes in this communicator
     integer CORE_ME           ! core id starting from 1 ... NCORE
     integer NCORE             ! number of cores on this SMP-node
     integer MAXCORE           ! maximum no. cores on any SMP-node
     integer N_SND             ! size of SMP shared memory buffer
     integer N_RCV             ! size of SMP shared memory buffer
     integer(8) SND_P          ! SNDBUF address (cray pointer), for real 
     integer(8) RCV_P          ! RCVBUF address (cray pointer), for real
     integer(8) SND_P_c        ! for complex
     integer(8) RCV_P_c        ! for complex
  END TYPE SMP_INFO
#endif

  ! derived type to store decomposition info for a given global data size
  TYPE, public :: DECOMP_INFO
     ! staring/ending index and size of data held by current processor
     integer, dimension(3) :: xst, xen, xsz  ! x-pencil
     integer, dimension(3) :: yst, yen, ysz  ! y-pencil
     integer, dimension(3) :: zst, zen, zsz  ! z-pencil

     ! in addition to local information, processors also need to know 
     ! some global information for global communications to work 

     ! how each dimension is distributed along pencils
     integer, allocatable, dimension(:) :: &
          x1dist, y1dist, y2dist, z2dist

     ! send/receive buffer counts and displacements for MPI_ALLTOALLV
     integer, allocatable, dimension(:) :: &
          x1cnts, y1cnts, y2cnts, z2cnts
     integer, allocatable, dimension(:) :: &
          x1disp, y1disp, y2disp, z2disp

     ! buffer counts for MPI_ALLTOALL only for evenly distributed data
     integer :: x1count, y1count, y2count, z2count

#ifdef SHM
     ! For shared-memory implementation

     ! one instance of this derived type for each communicator
     ! shared moemory info, such as which MPI rank belongs to which node
     TYPE(SMP_INFO) :: ROW_INFO, COL_INFO

     ! shared send/recv buffers for ALLTOALLV
     integer, allocatable, dimension(:) :: x1cnts_s, y1cnts_s, &
          y2cnts_s, z2cnts_s
     integer, allocatable, dimension(:) :: x1disp_s, y1disp_s, &
          y2disp_s, z2disp_s
     ! A copy of original buffer displacement (will be overwriten)
     integer, allocatable, dimension(:) :: x1disp_o, y1disp_o, &
          y2disp_o, z2disp_o
#endif
  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), save :: decomp_main

  ! staring/ending index and size of data held by current processor
  ! duplicate as in 'main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       decomp_info_init, decomp_info_finalize, partition, &
       update_halo


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are routines to perform global data transpositions
  ! 
  !   Four combinations are available, enough to cover all situations
  !    - transpose_x_to_y (X-pencil --> Y-pencil)
  !    - transpose_y_to_z (Y-pencil --> Z-pencil)
  !    - transpose_z_to_y (Z-pencil --> Y-pencil)
  !    - transpose_y_to_x (Y-pencil --> X-pencil)
  !
  !   Generic interface provided here to support multiple data types
  !    - real and complex types supported through generic interface
  !    - single/double precision supported through pre-processing
  !       * see 'mytype' variable at the beginning
  !    - an optional argument can be supplied to transpose data whose 
  !      global size is not the default nx*ny*nz 
  !       * as the case in fft r2c/c2r interface 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface transpose_x_to_y
     module procedure transpose_x_to_y_real
     module procedure transpose_x_to_y_complex
  end interface
  
  interface transpose_y_to_z
     module procedure transpose_y_to_z_real
     module procedure transpose_y_to_z_complex
  end interface
  
  interface transpose_z_to_y
     module procedure transpose_z_to_y_real
     module procedure transpose_z_to_y_complex
  end interface

  interface transpose_y_to_x
     module procedure transpose_y_to_x_real
     module procedure transpose_y_to_x_complex
  end interface
  

contains

#ifdef SHM_DEBUG
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For debugging, print the shared-memory structure
  subroutine print_smp_info(s)
    TYPE(SMP_INFO) :: s
    write(10,*) 'size of current communicator:', s%NCPU
    write(10,*) 'rank in current communicator:', s%NODE_ME
    write(10,*) 'NSMP - number of SMP-nodes in this communicator:', s%NSMP
    write(10,*) 'SMP-node id (1 ~ NSMP):', s%SMP_ME
    write(10,*) 'NCORE - number of cores on this SMP-node', s%NCORE
    write(10,*) 'core id (1 ~ NCORE):', s%CORE_ME
    write(10,*) 'maximum no. cores on any SMP-node:', s%MAXCORE
    write(10,*) 'size of SMP shared memory SND buffer:', s%N_SND
    write(10,*) 'size of SMP shared memory RCV buffer:', s%N_RCV
  end subroutine print_smp_info
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col)

    implicit none

    integer, intent(IN) :: nx,ny,nz,p_row,p_col
    
    integer :: errorcode, ierror, row, col
    
#ifdef SHM_DEBUG
    character(len=80) fname
#endif

    nx_global = nx
    ny_global = ny
    nz_global = nz

    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, row, col)
    else
       if (nproc /= p_row*p_col) then
          errorcode = 1
          call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
       else
          row = p_row
          col = p_col
       end if
    end if
    
    ! 2D Catersian topology
    dims(1) = row
    dims(2) = col
    periodic(1) = .true.
    periodic(2) = .true.
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., &  ! do not reorder rank
         MPI_COMM_CART, ierror)
    call MPI_CART_COORDS(MPI_COMM_CART,nrank,2,coord,ierror)
    
    ! these are the communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(MPI_COMM_CART,(/.true.,.false./), &
         MPI_COMM_COL,ierror)
    call MPI_CART_SUB(MPI_COMM_CART,(/.false.,.true./), &
         MPI_COMM_ROW,ierror)

    ! gather information for halo-cell support code
    call init_neighbour
    
    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,decomp_main)
    
    ! make a copy in global variables so applications can use
    ! these to create data structures 
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz

#ifdef SHM_DEBUG
    write(fname,99) nrank
99  format('log',I2.2)
    open(10,file=fname)
    write(10,*)'I am mpi rank ', nrank, 'Total ranks ', nproc
    write(10,*)' '
    write(10,*)'Global data size:'
    write(10,*)'nx*ny*nz', nx,ny,nz
    write(10,*)' '
    write(10,*)'2D processor grid:'
    write(10,*)'p_row*p_col:', dims(1), dims(2)
    write(10,*)' '
    write(10,*)'Portion of global data held locally:'
    write(10,*)'xsize:',xsize
    write(10,*)'ysize:',ysize
    write(10,*)'zsize:',zsize
    write(10,*)' '
    write(10,*)'How pensils are to be divided and sent in alltoallv:'
    write(10,*)'x1dist:',decomp_main%x1dist
    write(10,*)'y1dist:',decomp_main%y1dist
    write(10,*)'y2dist:',decomp_main%y2dist
    write(10,*)'z2dist:',decomp_main%z2dist
    write(10,*)' '
    write(10,*)'######Shared buffer set up after this point######'
    write(10,*)' '
    write(10,*) 'col communicator details:'
    call print_smp_info(decomp_main%COL_INFO)
    write(10,*)' '
    write(10,*) 'row communicator details:'
    call print_smp_info(decomp_main%ROW_INFO)
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of per-core buffers'
    write(10,*)'x1cnts:',decomp_main%x1cnts
    write(10,*)'y1cnts:',decomp_main%y1cnts
    write(10,*)'y2cnts:',decomp_main%y2cnts
    write(10,*)'z2cnts:',decomp_main%z2cnts
    write(10,*)'x1disp:',decomp_main%x1disp
    write(10,*)'y1disp:',decomp_main%y1disp
    write(10,*)'y2disp:',decomp_main%y2disp
    write(10,*)'z2disp:',decomp_main%z2disp
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of shared buffers'
    write(10,*)'x1cnts:',decomp_main%x1cnts_s
    write(10,*)'y1cnts:',decomp_main%y1cnts_s
    write(10,*)'y2cnts:',decomp_main%y2cnts_s
    write(10,*)'z2cnts:',decomp_main%z2cnts_s
    write(10,*)'x1disp:',decomp_main%x1disp_s
    write(10,*)'y1disp:',decomp_main%y1disp_s
    write(10,*)'y2disp:',decomp_main%y2disp_s
    write(10,*)'z2disp:',decomp_main%z2disp_s
    write(10,*)' '
    close(10)
#endif
    
    return
  end subroutine decomp_2d_init
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to clean things up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_finalize

    implicit none
    
    call decomp_info_finalize(decomp_main)
    
    return
  end subroutine decomp_2d_finalize
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init(nx,ny,nz,decomp)

    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
    
    ! distribute mesh points
    allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
         decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1))
    call get_dist(nx,ny,nz,decomp)
    
    ! generate partition information - starting/ending index etc.
    call partition(nx, ny, nz, (/ 1,2,3 /), &
         decomp%xst, decomp%xen, decomp%xsz)
    call partition(nx, ny, nz, (/ 2,1,3 /), &
         decomp%yst, decomp%yen, decomp%ysz)
    call partition(nx, ny, nz, (/ 2,3,1 /), &
         decomp%zst, decomp%zen, decomp%zsz)
    
    ! prepare send/receive buffer displacement and count for ALLTOALLV
    allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
         decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1))
    allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
         decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1))
    call prepare_buffer(decomp)

#ifdef SHM
    ! prepare shared-memory information if required
    call decomp_info_init_shm(decomp)
#endif
    
    return
  end subroutine decomp_info_init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory associated with a DECOMP_INFO object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_finalize(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    deallocate(decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist)
    deallocate(decomp%x1cnts,decomp%y1cnts,decomp%y2cnts,decomp%z2cnts)
    deallocate(decomp%x1disp,decomp%y1disp,decomp%y2disp,decomp%z2disp)

#ifdef SHM
    deallocate(decomp%x1disp_o,decomp%y1disp_o,decomp%y2disp_o, &
         decomp%z2disp_o)
    deallocate(decomp%x1cnts_s,decomp%y1cnts_s,decomp%y2cnts_s, &
         decomp%z2cnts_s)
    deallocate(decomp%x1disp_s,decomp%y1disp_s,decomp%y2disp_s, &
         decomp%z2disp_s)
#endif

    return
  end subroutine decomp_info_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find sub-domain information held by current processor
  !   INPUT: 
  !     nx, ny, nz - global data dimension
  !     pdim(3)    - number of processor grid in each dimension, 
  !                  valid values: 1 - distibute locally; 
  !                                2 - distribute across p_row; 
  !                                3 - distribute across p_col
  !   OUTPUT:
  !     lstart(3)  - starting index
  !     lend(3)    - ending index
  !     lsize(3)   - size of the sub-block (redundant) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    integer, dimension(3), intent(IN) :: pdim	
    integer, dimension(3), intent(OUT) :: lstart, lend, lsize

    integer, allocatable, dimension(:) :: st,en,sz
    integer :: i, gsize

    do i = 1, 3
 
      if (i==1) then
        gsize = nx
      else if (i==2) then
        gsize = ny
      else if (i==3) then
        gsize = nz
      endif

      if (pdim(i) == 1) then        ! all local
        lstart(i) = 1
        lend(i)   = gsize
        lsize(i)  = gsize
      elseif (pdim(i) == 2) then    ! distribute across dims(1)
        allocate(st(0:dims(1)-1))
        allocate(en(0:dims(1)-1))
        allocate(sz(0:dims(1)-1))
        call distribute(gsize,dims(1),st,en,sz)
        lstart(i) = st(coord(1))
        lend(i)   = en(coord(1))
        lsize(i)  = sz(coord(1))
        deallocate(st,en,sz)
      elseif (pdim(i) == 3) then    ! distribute across dims(2)
        allocate(st(0:dims(2)-1))
        allocate(en(0:dims(2)-1))
        allocate(sz(0:dims(2)-1))
        call distribute(gsize,dims(2),st,en,sz)
        lstart(i) = st(coord(2))
        lend(i)   = en(coord(2))
        lsize(i)  = sz(coord(2))
        deallocate(st,en,sz)
      endif    

    enddo
    return   

  end subroutine partition

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine distribute(data1,proc,st,en,sz)
  
    implicit none
    ! data1 -- data size in any dimension to be partitioned
    ! proc  -- number of processors in that dimension
    ! st    -- array of starting index
    ! en    -- array of ending index
    ! sz    -- array of local size  (redundent)
    integer data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size1,nl,nue
  
    size1=data1/proc
    nue = data1 - size1 * proc
    nl = proc - nue
    st(0) = 1
    sz(0) = size1
    en(0) = size1
    do i=1,nl-1
      st(i) = st(i-1) + size1
      sz(i) = size1
      en(i) = en(i-1) + size1
    enddo
    size1 = size1 + 1
    do i=nl,proc-1
      st(i) = en(i-1) + 1
      sz(i) = size1
      en(i) = en(i-1) + size1
    enddo
    en(proc-1)= data1 
    sz(proc-1)= data1-st(proc-1)+1
  
    return
  end subroutine distribute

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Define how each dimension is distributed across processors
  !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
  !    such global information is required locally at MPI_ALLTOALLV time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist(nx,ny,nz,decomp)

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
    integer, allocatable, dimension(:) :: st,en

    allocate(st(0:dims(1)-1))
    allocate(en(0:dims(1)-1))
    call distribute(nx,dims(1),st,en,decomp%x1dist)
    call distribute(ny,dims(1),st,en,decomp%y1dist)
    deallocate(st,en)

    allocate(st(0:dims(2)-1))
    allocate(en(0:dims(2)-1))
    call distribute(ny,dims(2),st,en,decomp%y2dist)
    call distribute(nz,dims(2),st,en,decomp%z2dist)
    deallocate(st,en)

    return
  end subroutine get_dist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_buffer(decomp)
 
  implicit none

  TYPE(DECOMP_INFO), intent(INOUT) :: decomp
  integer :: i

  do i=0, dims(1)-1
    decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
    decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
    if (i==0) then
      decomp%x1disp(i) = 0  ! displacement is 0-based index
      decomp%y1disp(i) = 0
    else
      decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
      decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
    endif
  enddo

  do i=0, dims(2)-1
    decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
    decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
    if (i==0) then
      decomp%y2disp(i) = 0  ! displacement is 0-based index
      decomp%z2disp(i) = 0
    else
      decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
      decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
    endif
  enddo  

  ! simpler information for ALLTOALL
  decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
  decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1) 
  decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
  decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

  return
  end subroutine prepare_buffer  

#ifdef SHM

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Generate shared-memory information 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init_shm(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! a copy of old displacement array (will be overwritten by shm code)
    allocate(decomp%x1disp_o(0:dims(1)-1),decomp%y1disp_o(0:dims(1)-1), &
         decomp%y2disp_o(0:dims(2)-1),decomp%z2disp_o(0:dims(2)-1))
    decomp%x1disp_o = decomp%x1disp
    decomp%y1disp_o = decomp%y1disp
    decomp%y2disp_o = decomp%y2disp
    decomp%z2disp_o = decomp%z2disp

    call prepare_shared_buffer(decomp%ROW_INFO,MPI_COMM_ROW,decomp)
    call prepare_shared_buffer(decomp%COL_INFO,MPI_COMM_COL,decomp)

    return
  end subroutine decomp_info_init_shm


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For shared-memory implementation, prepare send/recv shared buffer
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_shared_buffer(C,MPI_COMM,decomp)

  implicit none

  TYPE(SMP_INFO) :: C
  INTEGER :: MPI_COMM
  TYPE(DECOMP_INFO) :: decomp
  
  INTEGER, ALLOCATABLE :: KTBL(:,:),NARY(:,:),KTBLALL(:,:)
  INTEGER MYSMP, MYCORE, COLOR

  integer :: ierror

  C%MPI_COMM = MPI_COMM
  CALL MPI_COMM_SIZE(MPI_COMM,C%NCPU,ierror)
  CALL MPI_COMM_RANK(MPI_COMM,C%NODE_ME,ierror)
  C%SMP_COMM  = MPI_COMM_NULL
  C%CORE_COMM = MPI_COMM_NULL
  C%SMP_ME= 0
  C%NCORE = 0
  C%CORE_ME = 0
  C%MAXCORE = 0
  C%NSMP  = 0
  C%N_SND = 0
  C%N_RCV = 0
  C%SND_P = 0
  C%RCV_P = 0
  C%SND_P_c = 0
  C%RCV_P_c = 0

  ! get the smp-node map for this communicator and set up smp communicators
  CALL GET_SMP_MAP(C%MPI_COMM, C%NSMP, MYSMP, C%NCORE, MYCORE, C%MAXCORE)
  C%SMP_ME = MYSMP + 1
  C%CORE_ME = MYCORE + 1
  ! - set up inter/intra smp-node communicators
  COLOR = MYCORE
  IF (COLOR.GT.0) COLOR = MPI_UNDEFINED
  CALL MPI_Comm_split(C%MPI_COMM, COLOR, MYSMP, C%SMP_COMM, ierror)
  CALL MPI_Comm_split(C%MPI_COMM, MYSMP, MYCORE, C%CORE_COMM, ierror)
  ! - allocate work space
  ALLOCATE(KTBL(C%MAXCORE,C%NSMP),NARY(C%NCPU,C%NCORE))
  ALLOCATE(KTBLALL(C%MAXCORE,C%NSMP))
  ! - set up smp-node/core to node_me lookup table
  KTBL = 0
  KTBL(C%CORE_ME,C%SMP_ME) = C%NODE_ME + 1
  CALL MPI_ALLREDUCE(KTBL,KTBLALL,C%NSMP*C%MAXCORE,MPI_INTEGER, &
       MPI_SUM,MPI_COMM,ierror)
  KTBL=KTBLALL
  !  IF (SUM(KTBL) /= C%NCPU*(C%NCPU+1)/2) &
  !       CALL MPI_ABORT(...

  ! compute offsets in shared SNDBUF and RCVBUF
  CALL MAPSET_SMPSHM(C, KTBL, NARY, decomp)
  
  DEALLOCATE(KTBL,NARY)
  
  return
  end subroutine prepare_shared_buffer


  !**********************************************************************
  ! Set up smp-node based shared memory maps
  !**********************************************************************
  SUBROUTINE MAPSET_SMPSHM(C, KTBL, NARY, decomp)
        
    IMPLICIT NONE
    
    TYPE (SMP_INFO) C
    INTEGER KTBL(C%MAXCORE,C%NSMP)
    INTEGER NARY(C%NCPU,C%NCORE)
    TYPE (DECOMP_INFO) :: decomp

    INTEGER i, j, k, l, N, PTR, BSIZ, ierror, status, seed
    character*16 s
 
    BSIZ = C%N_SND
    
    ! a - SNDBUF
    IF (C%MPI_COMM==MPI_COMM_COL) THEN
       ALLOCATE(decomp%x1cnts_s(C%NSMP),decomp%x1disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%x1cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%x1disp_s(i) = PTR
          N = 0
          DO j=1,C%MAXCORE
             k = KTBL(j,i)
             IF (k > 0) then
                DO l=1,C%NCORE
                   IF (l == C%CORE_ME) decomp%x1disp_o(k-1) = PTR
                   N = N + NARY(k,l)
                   PTR = PTR + NARY(k,l)
                ENDDO
             ENDIF
          ENDDO
          decomp%x1cnts_s(i) = N
       ENDDO
       decomp%x1disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
       
    ELSE IF (C%MPI_COMM==MPI_COMM_ROW) THEN
       ALLOCATE(decomp%y2cnts_s(C%NSMP),decomp%y2disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%y2cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%y2disp_s(i) = PTR
          N = 0
          DO j=1,C%MAXCORE
             k = KTBL(j,i)
             IF (k > 0) then
                DO l=1,C%NCORE
                   IF (l == C%CORE_ME) decomp%y2disp_o(k-1) = PTR
                   N = N + NARY(k,l)
                   PTR = PTR + NARY(k,l)
                ENDDO
             ENDIF
          ENDDO
          decomp%y2cnts_s(i) = N
       ENDDO
       decomp%y2disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
    ENDIF
    
    ! b - RCVBUF
    
    IF (C%MPI_COMM==MPI_COMM_COL) THEN
       ALLOCATE(decomp%y1cnts_s(C%NSMP),decomp%y1disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%y1cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%y1disp_s(i) = PTR
          N=0
          DO j=1,C%NCORE
             DO l=1,C%MAXCORE
                k = KTBL(l,i)
                IF (k > 0) then
                   IF (j == C%CORE_ME) decomp%y1disp_o(k-1) = PTR
                   N = N + NARY(k,j)
                   PTR = PTR + NARY(k,j)
                ENDIF
             ENDDO
          ENDDO
          decomp%y1cnts_s(i) = N
       ENDDO
       decomp%y1disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
       
    ELSE IF (C%MPI_COMM==MPI_COMM_ROW) THEN
       ALLOCATE(decomp%z2cnts_s(C%NSMP),decomp%z2disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%z2cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%z2disp_s(i) = PTR
          N=0
          DO j=1,C%NCORE
             DO l=1,C%MAXCORE
                k = KTBL(l,i)
                IF (k > 0) then
                   IF (j == C%CORE_ME) decomp%z2disp_o(k-1) = PTR
                   N = N + NARY(k,j)
                   PTR = PTR + NARY(k,j)
                ENDIF
             ENDDO
          ENDDO
          decomp%z2cnts_s(i) = N
       ENDDO
       decomp%z2disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
       
    ENDIF
    
    ! check buffer size and (re)-allocate buffer space if necessary
    IF (BSIZ > C%N_SND) then
       IF (C%SND_P /= 0) CALL DEALLOC_SHM(C%SND_P, C%CORE_COMM)
       ! make sure each rank has unique keys to get shared memory
       !IF (C%MPI_COMM==MPI_COMM_COL) THEN
       !   seed = nrank+nproc*0+1 ! has to be non-zero
       !ELSE IF (C%MPI_COMM==MPI_COMM_ROW) THEN
       !   seed = nrank+nproc*1+1
       !END IF
       status = 1
       !CALL ALLOC_SHM(C%SND_P, BSIZ, real_type, C%CORE_COMM, status, seed)
       CALL ALLOC_SHM(C%SND_P, BSIZ, real_type, C%CORE_COMM, status)
       C%N_SND = BSIZ

       IF (C%RCV_P /= 0) CALL DEALLOC_SHM(C%RCV_P, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%RCV_P, BSIZ, real_type, C%CORE_COMM, status)
       C%N_RCV = BSIZ

       IF (C%SND_P_c /= 0) CALL DEALLOC_SHM(C%SND_P_c, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%SND_P_c, BSIZ, complex_type, C%CORE_COMM, status)
       C%N_SND = BSIZ

       IF (C%RCV_P_c /= 0) CALL DEALLOC_SHM(C%RCV_P_c, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%RCV_P_c, BSIZ, complex_type, C%CORE_COMM, status)
       C%N_RCV = BSIZ


    ENDIF
    
    RETURN
  END SUBROUTINE MAPSET_SMPSHM

#endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transpose from X-pencils to Y-pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine transpose_x_to_y_real(src,dst,opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    real(mytype), allocatable, dimension(:,:,:) :: work1, work2
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p=decomp%COL_INFO%SND_P
#else
    allocate (work1(s1,s2,s3))
#endif
    
    call mem_split_real(1,src,s1,s2,s3,work1,dims(1),decomp%x1dist,decomp) 
    
    ! define receive buffer
#ifdef SHM
    work2_p=decomp%COL_INFO%RCV_P
#else
    allocate (work2(d1,d2,d3))
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%x1cnts_s,decomp%x1disp_s, &
            real_type,    work2,decomp%y1cnts_s,decomp%y1disp_s, &
            real_type,    decomp%COL_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1,decomp%x1count,real_type, &
         work2,decomp%y1count,real_type, MPI_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1,decomp%x1cnts,decomp%x1disp,real_type, &
         work2,decomp%y1cnts,decomp%y1disp,real_type, MPI_COMM_COL, &
         ierror)
#endif
#endif

#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! rearrange receive buffer
    call mem_merge_real(1,work2,d1,d2,d3,dst,dims(1),decomp%y1dist,decomp)
    
#ifndef SHM
    deallocate(work1,work2)
#endif
    
    return
  end subroutine transpose_x_to_y_real
  
  ! complex version takes an optional argument for different global size
  subroutine transpose_x_to_y_complex(src,dst,opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    complex(mytype), allocatable, dimension(:,:,:) :: work1, work2
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p=decomp%COL_INFO%SND_P_c
#else
    allocate (work1(s1,s2,s3))
#endif
    
    call mem_split_complex(1,src,s1,s2,s3,work1,dims(1),decomp%x1dist,decomp) 
    
    ! define receive buffer
#ifdef SHM
    work2_p=decomp%COL_INFO%RCV_P_c
#else
    allocate (work2(d1,d2,d3))
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%x1cnts_s,decomp%x1disp_s, &
            complex_type, work2,decomp%y1cnts_s,decomp%y1disp_s, &
            complex_type, decomp%COL_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1,decomp%x1count,complex_type, &
         work2,decomp%y1count,complex_type, MPI_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1,decomp%x1cnts,decomp%x1disp,complex_type, &
         work2,decomp%y1cnts,decomp%y1disp,complex_type, MPI_COMM_COL, &
         ierror)
#endif
#endif

#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! rearrange receive buffer
    call mem_merge_complex(1,work2,d1,d2,d3,dst,dims(1),decomp%y1dist,decomp)
    
#ifndef SHM
    deallocate(work1,work2)
#endif
    
    return
  end subroutine transpose_x_to_y_complex



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transpose from Y-pencils to Z-pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine transpose_y_to_z_real(src,dst,opt_decomp)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    real(mytype), allocatable, dimension(:,:,:) :: work1
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if
    
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p=decomp%ROW_INFO%SND_P
#else
    allocate (work1(s1,s2,s3))
#endif
    
    call mem_split_real(2,src,s1,s2,s3,work1,dims(2),decomp%y2dist,decomp)
    
    ! define receive buffer, not needed for non-shm version as
    ! receive buffer same (i,j,k) order as destination array
#ifdef SHM
    work2_p=decomp%ROW_INFO%RCV_P
    
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%y2cnts_s,decomp%y2disp_s, &
            real_type,    work2,decomp%z2cnts_s,decomp%z2disp_s, &
            real_type,    decomp%ROW_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1,decomp%y2count,real_type, &
         dst,decomp%z2count,real_type, MPI_COMM_ROW, ierror)
#else
    call MPI_ALLTOALLV(work1,decomp%y2cnts,decomp%y2disp,real_type, &
         dst,decomp%z2cnts,decomp%z2disp,real_type, MPI_COMM_ROW, ierror)
#endif
#endif

#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    
    call mem_merge_real(2,work2,d1,d2,d3,dst,dims(2),decomp%z2dist,decomp)
#else
    ! note this receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
    
#ifndef SHM
    deallocate(work1)
#endif
    
    return
  end subroutine transpose_y_to_z_real
  
  ! complex version takes an optional argument for different global size
  subroutine transpose_y_to_z_complex(src,dst,opt_decomp)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    complex(mytype), allocatable, dimension(:,:,:) :: work1
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if
    
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p=decomp%ROW_INFO%SND_P_c
#else
    allocate (work1(s1,s2,s3))
#endif
    
    call mem_split_complex(2,src,s1,s2,s3,work1,dims(2),decomp%y2dist,decomp)
    
    ! define receive buffer, not needed for non-shm version as
    ! receive buffer same (i,j,k) order as destination array
#ifdef SHM
    work2_p=decomp%ROW_INFO%RCV_P_c
    
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%y2cnts_s,decomp%y2disp_s, &
            complex_type, work2,decomp%z2cnts_s,decomp%z2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1,decomp%y2count,complex_type, &
         dst,decomp%z2count,complex_type, MPI_COMM_ROW, ierror)
#else
    call MPI_ALLTOALLV(work1,decomp%y2cnts,decomp%y2disp,complex_type, &
         dst,decomp%z2cnts,decomp%z2disp,complex_type, MPI_COMM_ROW,   &
         ierror)
#endif
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    
    call mem_merge_complex(2,work2,d1,d2,d3,dst,dims(2),decomp%z2dist,decomp)
#else
    ! note this receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
    
#ifndef SHM
    deallocate(work1)
#endif
    
    return
  end subroutine transpose_y_to_z_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transpose from Z-pencils to Y-pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine transpose_z_to_y_real(src,dst,opt_decomp)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    real(mytype), allocatable, dimension(:,:,:) :: work2
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror
    
    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! define send buffer, not needed for non-shm version as
    ! send buffer same (i,j,k) order as source array
    
#ifdef SHM
    work1_p=decomp%ROW_INFO%SND_P
    
    call mem_split_real(3,src,s1,s2,s3,work1,dims(2),decomp%z2dist,decomp)
#else
    ! note the src array is suitable to be a send buffer
    ! so no split operation needed
#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p=decomp%ROW_INFO%RCV_P
#else
    allocate (work2(d1,d2,d3))
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%z2cnts_s,decomp%z2disp_s, &
            real_type,    work2,decomp%y2cnts_s,decomp%y2disp_s, &
            real_type,    decomp%ROW_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(  src,decomp%z2count,real_type, &
         work2,decomp%y2count,real_type, MPI_COMM_ROW, ierror)
#else
    call MPI_ALLTOALLV(  src,decomp%z2cnts,decomp%z2disp,real_type, &
         work2,decomp%y2cnts,decomp%y2disp,real_type, MPI_COMM_ROW, &
         ierror)
#endif
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
    ! rearrange receive buffer
    call mem_merge_real(3,work2,d1,d2,d3,dst,dims(2),decomp%y2dist,decomp)
    
#ifndef SHM
    deallocate(work2)
#endif
    
    return
  end subroutine transpose_z_to_y_real
  
  ! complex version takes an optional argument for different global size
  subroutine transpose_z_to_y_complex(src,dst,opt_decomp)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    complex(mytype), allocatable, dimension(:,:,:) :: work2
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if
    
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! define send buffer, not needed for non-shm version as
    ! send buffer same (i,j,k) order as source array
    
#ifdef SHM
    work1_p=decomp%ROW_INFO%SND_P_c
    
    call mem_split_complex(3,src,s1,s2,s3,work1,dims(2),decomp%z2dist,decomp)
#else
    ! note the src array is suitable to be a send buffer
    ! so no split operation needed
#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p=decomp%ROW_INFO%RCV_P_c
#else
    allocate (work2(d1,d2,d3))
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%z2cnts_s,decomp%z2disp_s, &
            complex_type, work2,decomp%y2cnts_s,decomp%y2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(  src,decomp%z2count,complex_type, &
         work2,decomp%y2count,complex_type, MPI_COMM_ROW, ierror)
#else
    call MPI_ALLTOALLV(  src,decomp%z2cnts,decomp%z2disp,complex_type, &
         work2,decomp%y2cnts,decomp%y2disp,complex_type, MPI_COMM_ROW, &
         ierror)
#endif
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
    ! rearrange receive buffer
    call mem_merge_complex(3,work2,d1,d2,d3,dst,dims(2),decomp%y2dist,decomp)
    
#ifndef SHM
    deallocate(work2)
#endif
    
    return
  end subroutine transpose_z_to_y_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transpose from Y-pencils to X-pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine transpose_y_to_x_real(src,dst,opt_decomp)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    real(mytype), allocatable, dimension(:,:,:) :: work1, work2
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror
    
    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p=decomp%COL_INFO%SND_P
#else
    allocate(work1(s1,s2,s3))
#endif
    
    call mem_split_real(4,src,s1,s2,s3,work1,dims(1),decomp%y1dist,decomp)
    
    ! define receive buffer
#ifdef SHM
    work2_p=decomp%COL_INFO%RCV_P
#else
    allocate (work2(d1,d2,d3))
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%y1cnts_s,decomp%y1disp_s, &
            real_type,    work2,decomp%x1cnts_s,decomp%x1disp_s, &
            real_type,    decomp%COL_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1,decomp%y1count,real_type, &
         work2,decomp%x1count,real_type, MPI_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1,decomp%y1cnts,decomp%y1disp,real_type, &
         work2,decomp%x1cnts,decomp%x1disp,real_type, MPI_COMM_COL, &
         ierror)
#endif
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! rearrange receive buffer
    call mem_merge_real(4,work2,d1,d2,d3,dst,dims(1),decomp%x1dist,decomp)
    
#ifndef SHM
    deallocate(work1,work2)
#endif
    
    return
  end subroutine transpose_y_to_x_real
  
  ! complex version takes an optional argument for different global size
  subroutine transpose_y_to_x_complex(src,dst,opt_decomp)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#else
    complex(mytype), allocatable, dimension(:,:,:) :: work1, work2
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if
    
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p=decomp%COL_INFO%SND_P_c
#else
    allocate(work1(s1,s2,s3))
#endif
    
    call mem_split_complex(4,src,s1,s2,s3,work1,dims(1),decomp%y1dist,decomp)
    
    ! define receive buffer
#ifdef SHM
    work2_p=decomp%COL_INFO%RCV_P_c
#else
    allocate (work2(d1,d2,d3))
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1,decomp%y1cnts_s,decomp%y1disp_s, &
            complex_type, work2,decomp%x1cnts_s,decomp%x1disp_s, &
            complex_type, decomp%COL_INFO%SMP_COMM, ierror)
    endif
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1,decomp%y1count,complex_type, &
         work2,decomp%x1count,complex_type, MPI_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1,decomp%y1cnts,decomp%y1disp,complex_type, &
         work2,decomp%x1cnts,decomp%x1disp,complex_type, MPI_COMM_COL, &
         ierror)
#endif
#endif
    
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! rearrange receive buffer
    call mem_merge_complex(4,work2,d1,d2,d3,dst,dims(1),decomp%x1dist,decomp)
    
#ifndef SHM
    deallocate(work1,work2)
#endif
    
    return
  end subroutine transpose_y_to_x_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Memory copying routines to pack/unpack ALLTOALLV buffers
  !
  !  - mem_split (natural ijk-ordered array --> send buffer)
  !  - mem_merge (receive buffer --> natural ijk-ordered array)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mem_split_real(ndir,in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    
#include "mem_split.f90"    
    
    return
  end subroutine mem_split_real
  
  subroutine mem_split_complex(ndir,in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    
#include "mem_split.f90"
    
    return
  end subroutine mem_split_complex
  
  
  subroutine mem_merge_real(ndir,in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out

#include "mem_merge.f90"
    
    return
  end subroutine mem_merge_real
  
  subroutine mem_merge_complex(ndir,in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out

#include "mem_merge.f90"
    
    return
  end subroutine mem_merge_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Auto-tuning algorithm to select the best 2D processor grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)

    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, dimension(20) :: factors
    double precision :: t1, t2, best_time
    integer :: nfact, i, row, col, ierror

    real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3

    TYPE(DECOMP_INFO) :: decomp

    if (nrank==0) write(*,*) 'In auto-tuning mode......'

    best_time = huge(t1)
    best_p_row = -1
    best_p_col = -1

    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    do i=1, nfact
       row = factors(i)
       col = iproc / row
       ! enforce the limitation of 2D decomposition
       if (col<=ny_global .or. col<=nz_global) then

          ! 2D Catersian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .false.
          periodic(2) = .false.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,MPI_COMM_CART, ierror)
          call MPI_CART_COORDS(MPI_COMM_CART,nrank,2,coord,ierror)
          
          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(MPI_COMM_CART,(/.true.,.false./), &
               MPI_COMM_COL,ierror)
          call MPI_CART_SUB(MPI_COMM_CART,(/.false.,.true./), &
               MPI_COMM_ROW,ierror)
          
          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global,nz_global,decomp)

          ! arrays for X,Y and Z-pencils
          allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
          allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          allocate(u3(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(u1,u2,decomp)
          call transpose_y_to_z(u2,u3,decomp)
          call transpose_z_to_y(u3,u2,decomp)
          call transpose_y_to_x(u2,u1,decomp)
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2,u3)
          call decomp_info_finalize(decomp)

          call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                   MPI_COMM_WORLD,ierror)
          t1 = t1 / dble(nproc)

          if (nrank==0) then
             write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
          end if

          if (best_time > t1) then
             best_time = t1
             best_p_row = row
             best_p_col = col
          end if

       end if
       
    enddo ! loop through processer grid

    if (nrank==0) then
       write(*,*) 'the best processor grid is probably ', &
            best_p_row, ' by ', best_p_col
    end if

    return
  end subroutine best_2d_grid

#include "factor.f90"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.f90"

  
end module decomp_2d
