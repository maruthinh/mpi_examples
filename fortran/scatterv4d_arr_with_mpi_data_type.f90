
program test_scatterv_gatherv
  use mpi 
  implicit none 
  integer(kind=4), allocatable, dimension(:,:,:,:) :: global, local 
  integer, allocatable, dimension(:) :: scnts, displs
  integer :: m1, m2, m3, m4, m1_local, str_idx, end_idx 
  integer :: mpierr, rank, nprocs, ierr
  real(kind=8) :: str_time

  m1 = 1000; m2=100; m3=100; m4=100

  call mpi_init(mpierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, mpierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierr)

  call domain_decompose(m1, rank, nprocs, str_idx, end_idx)
  m1_local = end_idx - str_idx + 1

  if(rank .eq. 0) then 
    allocate(global(m1, m2, m3, m4))
    global = 10
  endif

  allocate(local(0:m1_local + 1, m2, m3, m4))
  allocate(scnts(nprocs), displs(nprocs), stat=ierr)

  call MPI_Allgather(str_idx, 1, MPI_INTEGER, displs, 1, MPI_INTEGER, &
  MPI_COMM_WORLD, mpierr)
  call MPI_Allgather(m1_local, 1, MPI_INTEGER, scnts, 1, &
  MPI_INTEGER, MPI_COMM_WORLD, mpierr)
  displs = displs - 1
  ! print*, rank, m1_local, str_idx, end_idx, scnts, displs

  str_time = MPI_Wtime()
  call scatter4D_arr(m1_local, m1, m2, m3, m4, scnts, displs, global, &
  local, rank)
  print*, "Time taken by scatterv operation from rank", rank, "time: ", &
  MPI_Wtime() - str_time

  call MPI_Finalize(mpierr)

  contains 
  subroutine domain_decompose(npts, rank, size, sidx, eidx)
    implicit none

    integer, intent(in) :: npts, size, rank
    integer, intent(out) :: sidx, eidx
    integer :: pts_per_proc

    pts_per_proc = npts/size

    if(rank < mod(npts, size)) then
      pts_per_proc=pts_per_proc + 1
    end if

    if(rank < mod(npts, size)) then
      sidx = rank * pts_per_proc + 1
      eidx = (rank + 1) * pts_per_proc
    else
      sidx = mod(npts, size) + rank*pts_per_proc + 1
      eidx = mod(npts, size) + (rank + 1) * pts_per_proc
    end if
  end subroutine domain_decompose

  subroutine scatter4D_arr(n1_local, n1, n2, n3, n4, scounts, &
    displs, g_arr, l_arr, rank)
    implicit none
    integer, intent(in) :: n1_local, n1, n2, n3, n4, rank
    integer, intent(in) :: scounts(:), displs(:)
    integer(4), intent(in) :: g_arr(n1, n2, n3, n4)
    integer(4), intent(inout) :: l_arr(0:n1_local+1, n2, n3, n4)

    integer :: ssizes(4), s_ssizes(4), sstarts(4), &
    rsizes(4), r_ssizes(4), rstarts(4), stype, rtype, &
    resize_stype, resize_rtype
    integer(kind=MPI_ADDRESS_KIND) :: lb, extent

    ssizes = [n1, n2, n3, n4]
    s_ssizes = [1, n2, n3, n4]
    rsizes = [n1_local+2, n2, n3, n4]
    r_ssizes = [1, n2, n3, n4]
    sstarts = [0, 0, 0, 0]
    rstarts = [1, 0, 0, 0]

    call MPI_Type_get_extent(MPI_INTEGER4, lb, extent, mpierr)

    !create a mpi subarray data type for sending data
    call MPI_Type_create_subarray(4, ssizes, s_ssizes, sstarts, &
    MPI_ORDER_FORTRAN, MPI_INTEGER4, stype, mpierr)

    !resize the send subarray for starting at correct location for next send
    call MPI_Type_create_resized(stype, lb, extent, &
    resize_stype, mpierr)
    call MPI_Type_commit(resize_stype, mpierr)

    !create a mpi subarray data type for receiving data
    call MPI_Type_create_subarray(4, rsizes, r_ssizes, rstarts, &
    MPI_ORDER_FORTRAN, MPI_INTEGER4, rtype, mpierr)

    !resize the receive subarray for starting at correct location for next receive
    call MPI_Type_create_resized(rtype, lb, extent, &
    resize_rtype, mpierr)
    call MPI_Type_commit(resize_rtype, mpierr)

    call MPI_Scatterv(g_arr, scounts, displs, resize_stype, &
    l_arr, scounts(rank), resize_rtype, 0, MPI_COMM_WORLD, mpierr)

  end subroutine scatter4D_arr
end program test_scatterv_gatherv
