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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support for neighbouring pencils to exchange data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_halo(inout, level, opt_decomp)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(mytype), dimension(:,:,:), intent(INOUT) :: inout    
    TYPE(DECOMP_INFO), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    ! starting/ending index of array with halo cells
    integer :: xs, ys, zs, xe, ye, ze

    integer :: i, j, k, s1, s2, s3, ierror

    integer :: icount, ilength, ijump 
    integer :: halo12, halo21, halo31, halo32                 
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_a, tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = size(inout,1) - 2*level
    s2 = size(inout,2) - 2*level
    s3 = size(inout,3) - 2*level

    ! Calculate the starting index and ending index of output
    xs = 1 
    xe = s1 + 2*level
    ys = 1
    ye = s2 + 2*level 
    zs = 1
    ze = s3 + 2*level

    ! If needed, define MPI derived data type to pack halo data,
    ! then call MPI send/receive to exchange halo data
    
    ! X-pencil
    if (s1==decomp%xsz(1)) then
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'X-pencil input'
          write(*,*) '=============='
          write(*,*) 'Data on a y-z plane is shown'
          write(*,*) 'Before halo exchange'
          do j=ye,ys,-1
             write(*,'(10F4.0)') (inout(1,j,k),k=zs,ze)
          end do
       end if
#endif
       ! *** east/west ***
       ! all data in local memory already, no halo exchange
       do k=zs,ze
          do j=ys,ye
            do i=xs,level
             inout(i,j,k)=inout(s1+i,j,k)      
           end do
         end do
       end do

       do k=zs,ze
         do j=ys,ye
           do i=1+level,2*level
             inout(s1+i,j,k)=inout(i,j,k)      
           end do
         end do
       end do  

       ! *** north/south *** 
       tag_a = 101 !coord(1)
       tag_b = 102 !mod(coord(1) + 1, dims(1))
       icount = s3 + 2*level
       ilength = level * (s1+2*level)
       ijump = (s1+2*level)*(s2+2*level)
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            real_type, halo12, ierror)
       call MPI_TYPE_COMMIT(halo12, ierror)
       ! receive from south
       call MPI_IRECV(inout(xs,1,zs), 1, halo12, &
            neighbour(1,4), tag_a, MPI_COMM_CART, requests(1), ierror)
       ! receive from north
       call MPI_IRECV(inout(xs,s2+level+1,zs), 1, halo12, &
            neighbour(1,3), tag_b, MPI_COMM_CART, requests(2), ierror)
       ! send to south
       call MPI_ISSEND(inout(xs,1+level,zs), 1, halo12, &
            neighbour(1,4), tag_b, MPI_COMM_CART, requests(3), ierror)
       ! send to north
       call MPI_ISSEND(inout(xs,s2+1,zs), 1, halo12, &
            neighbour(1,3), tag_a, MPI_COMM_CART, requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo12, ierror)
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in Y'
          do j=ye,ys,-1
             write(*,'(10F4.0)') (inout(1,j,k),k=zs,ze)
          end do
       end if
#endif
       ! *** top/bottom ***
       ! no need to define derived data type as data on xy-planes
       ! all contiguous in memory, which can be sent/received using
       ! MPI directly
       tag_a = 103 !coord(2)
       tag_b = 104 !mod(coord(2) + 1, dims(2))
       icount = ((s1+2*level) * (s2+2*level)) * level
       ! receive from bottom
       call MPI_IRECV(inout(xs,ys,1), icount, real_type, &
            neighbour(1,6), tag_a, MPI_COMM_CART, requests(1), ierror)
       ! receive from top
       call MPI_IRECV(inout(xs,ys,s3+level+1), icount, real_type, &
            neighbour(1,5), tag_b, MPI_COMM_CART, requests(2), ierror)
       ! send to bottom
       call MPI_ISSEND(inout(xs,ys,1+level), icount, real_type, &
            neighbour(1,6), tag_b, MPI_COMM_CART, requests(3), ierror)
       ! send to top
       call MPI_ISSEND(inout(xs,ys,s3+1), icount, real_type, &
            neighbour(1,5), tag_a, MPI_COMM_CART, requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)

#ifdef HALO_DEBUG       
       if (nrank==4) then
          write(*,*) 'After exchange in Z'
          do j=ye,ys,-1
             write(*,'(10F4.0)') (inout(1,j,k),k=zs,ze)
          end do
       end if
#endif       
    ! Y-pencil   
    else if (s2==decomp%ysz(2)) then
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'Y-pencil input'
          write(*,*) '=============='
          write(*,*) 'Data on a x-z plane is shown'
          write(*,*) 'Before halo exchange'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (inout(i,1,k),k=zs,ze)
          end do
       end if
#endif
       ! *** east/west ***
       tag_a = 201 !coord(1)
       tag_b = 202 !mod(coord(1) + 1, dims(1))
       icount = (s2+2*level)*(s3+2*level)
       ilength = level
       ijump = s1+2*level
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            real_type, halo21, ierror)
       call MPI_TYPE_COMMIT(halo21, ierror)
       ! receive from west
       call MPI_IRECV(inout(1,ys,zs), 1, halo21, &
            neighbour(2,2), tag_a, MPI_COMM_CART, requests(1), ierror)
       ! receive from east
       call MPI_IRECV(inout(s1+level+1,ys,zs), 1, halo21, &
            neighbour(2,1), tag_b, MPI_COMM_CART, requests(2), ierror)
       ! send to west
       call MPI_ISSEND(inout(1+level,ys,zs), 1, halo21, &
            neighbour(2,2), tag_b, MPI_COMM_CART, requests(3), ierror)
       ! send to east
       call MPI_ISSEND(inout(s1+1,ys,zs), 1, halo21, &
            neighbour(2,1), tag_a, MPI_COMM_CART, requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo21, ierror)
#ifdef HALO_DEBUG       
       if (nrank==4) then
          write(*,*) 'After exchange in X'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (inout(i,1,k),k=zs,ze)
          end do
       end if
#endif
       ! *** north/south ***
       ! all data in local memory already, no halo exchange

       do k=zs,ze
         do j=ys,level
           do i=xs,xe
             inout(i,j,k)=inout(i,s2+j,k)      
           end do
         end do
       end do

       do k=zs,ze
         do j=1+level,2*level
           do i=xs,xe
             inout(i,s2+j,k)=inout(i,j,k)      
           end do
         end do
       end do  

       ! *** top/bottom ***
       ! no need to define derived data type as data on xy-planes
       ! all contiguous in memory, which can be sent/received using
       ! MPI directly
       tag_a = 203 !coord(2)
       tag_b = 204 !mod(coord(2) + 1, dims(2))
       icount = ((s2+2*level) * (s1+2*level)) * level
       ! receive from bottom
       call MPI_IRECV(inout(xs,ys,1), icount, real_type, &
            neighbour(2,6), tag_a, MPI_COMM_CART, requests(1), ierror)
       ! receive from top
       call MPI_IRECV(inout(xs,ys,s3+level+1), icount, real_type, &
            neighbour(2,5), tag_b, MPI_COMM_CART, requests(2), ierror)
       ! send to bottom
       call MPI_ISSEND(inout(xs,ys,1+level), icount, real_type, &
            neighbour(2,6), tag_b, MPI_COMM_CART, requests(3), ierror)
       ! send to top
       call MPI_ISSEND(inout(xs,ys,s3+1), icount, real_type, &
            neighbour(2,5), tag_a, MPI_COMM_CART, requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in Z'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (inout(i,1,k),k=zs,ze)
          end do
       end if
#endif
    ! Z-pencil
    else if (s3==decomp%zsz(3)) then   
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'Z-pencil input'
          write(*,*) '=============='
          write(*,*) 'Data on a x-y plane is shown'
          write(*,*) 'Before halo exchange'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (inout(i,j,1),j=ys,ye)
          end do
       end if
#endif
       ! *** east/west ***
       tag_a = 301 !coord(1)
       tag_b = 302 !mod(coord(1) + 1, dims(1))
       icount = (s2+2*level)*(s3+2*level)
       ilength = level
       ijump = s1+2*level
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            real_type, halo31, ierror)
       call MPI_TYPE_COMMIT(halo31, ierror)
       ! receive from west
       call MPI_IRECV(inout(1,ys,zs), 1, halo31, &
            neighbour(3,2), tag_a, MPI_COMM_CART, requests(1), ierror)
       ! receive from east
       call MPI_IRECV(inout(s1+level+1,ys,zs), 1, halo31, &
            neighbour(3,1), tag_b, MPI_COMM_CART, requests(2), ierror)
       ! send to west
       call MPI_ISSEND(inout(1+level,ys,zs), 1, halo31, &
            neighbour(3,2), tag_b, MPI_COMM_CART, requests(3), ierror)
       ! send to east
       call MPI_ISSEND(inout(s1+1,ys,zs), 1, halo31, &
            neighbour(3,1), tag_a, MPI_COMM_CART, requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo31, ierror)       
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in X'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (inout(i,j,1),j=ys,ye)
          end do
       end if
#endif
       ! *** north/south *** 
       tag_a = 303 !coord(2)
       tag_b = 304 !mod(coord(2) + 1, dims(2))
       icount = s3 + 2*level
       ilength = level * (s1+2*level)
       ijump = (s1+2*level)*(s2+2*level)
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            real_type, halo32, ierror)
       call MPI_TYPE_COMMIT(halo32, ierror)
       ! receive from south
       call MPI_IRECV(inout(xs,1,zs), 1, halo32, &
            neighbour(3,4), tag_a, MPI_COMM_CART, requests(1), ierror)
       ! receive from north
       call MPI_IRECV(inout(xs,s2+level+1,zs), 1, halo32, &
            neighbour(3,3), tag_b, MPI_COMM_CART, requests(2), ierror)
       ! send to south
       call MPI_ISSEND(inout(xs,1+level,zs), 1, halo32, &
            neighbour(3,4), tag_b, MPI_COMM_CART, requests(3), ierror)
       ! send to north
       call MPI_ISSEND(inout(xs,s2+1,zs), 1, halo32, &
            neighbour(3,3), tag_a, MPI_COMM_CART, requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo32, ierror)
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in Y'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (inout(i,j,1),j=ys,ye)
          end do
       end if
#endif
       ! *** top/bottom ***
       ! all data in local memory already, no halo exchange

       do k=zs,level
         do j=ys,ye
           do i=xs,xe
             inout(i,j,k)=inout(i,j,s3+k)      
           end do
         end do
       end do

       do k=1+level,2*level
         do j=ys,ye
           do i=xs,xe
             inout(i,j,s3+k)=inout(i,j,k)      
           end do
         end do
       end do   

    end if


    return
  end subroutine update_halo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To support halo-cell exchange:
  !   find the MPI ranks of neighbouring pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_neighbour

    integer :: ierror, i

    ! For X-pencil
    neighbour(1,1) = MPI_PROC_NULL               ! east
    neighbour(1,2) = MPI_PROC_NULL               ! west
    call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, &
         neighbour(1,4), neighbour(1,3), ierror) ! north & south
    call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, &
         neighbour(1,6), neighbour(1,5), ierror) ! top & bottom

    ! For Y-pencil
    call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, &
         neighbour(2,2), neighbour(2,1), ierror) ! east & west
    neighbour(2,3) = MPI_PROC_NULL               ! north
    neighbour(2,4) = MPI_PROC_NULL               ! south
    call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, &
         neighbour(2,6), neighbour(2,5), ierror) ! top & bottom

    ! For Z-pencil
    call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, &
         neighbour(3,2), neighbour(3,1), ierror) ! east & west
    call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, &
         neighbour(3,4), neighbour(3,3), ierror) ! north & south
    neighbour(3,5) = MPI_PROC_NULL               ! top
    neighbour(3,6) = MPI_PROC_NULL               ! bottom

    return
  end subroutine init_neighbour
