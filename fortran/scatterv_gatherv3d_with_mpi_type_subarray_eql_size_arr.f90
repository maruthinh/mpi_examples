!Example code to Scatterv and Gatherv 3D array using MPI subarray 
!derived data type.

!For example, a 3D array of size (4,3,2) can be distributed
!to rank 0 and rank 1 as (2,3,2) and (2,3,2) or to 4 ranks 
!(1,3,2), (1,3,2), (1,3,2) and (1,3,2). Which means that 
!first dimension of the array must be divisible by the size of 
!the communicator. So (5,3,2) array can't be distrbuted to two ranks
!as (3,3,2) and (2,3,2). But it can be distributed to 5 ranks

program ex_scatterv
  use mpi  
  use iso_fortran_env, only : real64 
  implicit none

  real(real64), allocatable,dimension(:,:,:) :: array, array_local
  integer :: rank, num_procs, mpierr, i, j, k
  integer :: nx, ny, nz, str_idx, end_idx, local_size 
  integer, dimension(:), allocatable :: sendcounts, displacements
  integer :: sizes(3), sub_sizes(3), starts(3), recv_starts(3), recv_sizes(3), &
    send_type, resize_send_type, recv_type
  integer(kind=MPI_OFFSET_KIND) :: lb, extent
  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world, num_procs, mpierr)
  call mpi_comm_rank(mpi_comm_world, rank, mpierr)

  !size of array
  nx=4
  ny=3
  nz=2
  
  !exit the program if nx is not divisible by num_procs 
  if(mod(nx, num_procs) .ne. 0) then 
    print*, "Size of the first index of the array must be divisible by number &
      of processors" 
    print*, "size of first index: ", nx, "Number of procs: ", num_procs
    stop 
    call MPI_Finalize(mpierr)
  endif 

  !allocate in the root rank
  if(rank==0) then
    allocate(array(nx,ny,nz))
  else !for other procs allocate with zero size
    allocate(array(0,0,0))
  endif

  !assign values to the array
  if(rank==0) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          array(i,j,k) = (i-1)+(j-1)*nx+(k-1)*nx*ny
        end do
      end do 	
    end do 
    print*, "Before scattering..."
    print*, array   
  endif

  !distribute the 3d array among different procs 
  call distribute_points(nx, rank, num_procs, str_idx, end_idx)
  local_size = end_idx - str_idx + 1

  !allocate local(for each rank) arrays
  allocate(array_local(local_size, ny, nz))

  !allocate sendcoutns and displacements arrays for braodcasting
  allocate(sendcounts(num_procs), displacements(num_procs))

  !Scatterning using subarray type 
  sizes = [nx, ny, nz]
  sub_sizes = [local_size, ny, nz]
  starts = [0, 0, 0] 
  recv_starts = [0, 0, 0] 

  !to get extent of MPI_DOUBLE_PRECISION
  call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lb, extent, mpierr)

  call MPI_Type_create_subarray(3, sizes, sub_sizes, starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, send_type, mpierr)

  call MPI_Type_create_resized(send_type, 0, local_size*extent, &
    resize_send_type, mpierr)
  call MPI_Type_commit(resize_send_type, mpierr)

  call MPI_Type_create_subarray(3, sub_sizes, sub_sizes, recv_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, mpierr)
  call MPI_Type_commit(recv_type, mpierr)

  !actual number of elements sent is defined by MPI_Type_create_subarray 
  sendcounts=1
  
  !displacements are the starting point in continguous memory of array defined 
  !interms of extent. 
  call MPI_Allgather(rank, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, &
    MPI_COMM_WORLD, mpierr)

  !scatter the array among procs 
  call MPI_Scatterv(array, sendcounts, displacements, resize_send_type, &
    array_local, sendcounts, recv_type, 0, MPI_COMM_WORLD, mpierr)
  
  !print the scattered array 
  print*, "Scattered array with subarray: ", rank
  print*, array_local 

  !do some computations on the scattered local arrays
  array_local = array_local+1

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  !Gather the local arrays
  call MPI_Gatherv(array_local, sendcounts, recv_type, array, &
    sendcounts, displacements, resize_send_type, 0, MPI_COMM_WORLD, mpierr)

  !print gathered array
  if(rank==0) then 
    print*, "Gathered array after adding 1 to scattered local array: ---------"
    print*, array
  endif
  call MPI_Finalize(mpierr)
end program ex_scatterv
