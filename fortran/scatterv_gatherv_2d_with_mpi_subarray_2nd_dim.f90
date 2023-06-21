!Example code to Scatterv and Gatherv 2D array using MPI subarray 
!derived data type.

!For example, a 2D array of size (3,5) can be distributed
!to rank 0 and rank 1 as (3,3) and (3,2) or to 3 ranks 
!(3,2), (3,2), (3,1), etc. 

!for ex global array-like
! 0 3 6 9 12 
! 1 4 7 10 13
! 2 5 8 11 14

!will be distributed among 2 procs as follows
!rank 0
! 0 3 6  
! 1 4 7 
! 2 5 8 

!rank 1
! 9 12 
! 10 13
! 11 14

!among 3 procs
!rank 0
! 0 3   
! 1 4  
! 2 5  

!rank 1
! 6 9  
! 7 10 
! 8 11 

!rank 2
! 12
! 13
! 14

program ex_scatterv
  use mpi  
  use iso_fortran_env, only : real64 
  implicit none
  
  !allocate arrays
  real(real64), allocatable,dimension(:,:) :: array, array_local
  integer :: rank, num_procs, i, j
  integer :: nx, ny, str_idx, end_idx, local_size 
  integer, dimension(:), allocatable :: sendcounts, displacements
  integer :: sizes(2), sub_sizes(2), starts(2), recv_starts(2), recv_sizes(2), &
    send_type, resize_send_type, recv_type, resize_recv_type
  integer(kind=8) :: lb, extent, lb_resize
  real(real64) :: start_time
  integer :: mpierr 
  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world, num_procs, mpierr)
  call mpi_comm_rank(mpi_comm_world, rank, mpierr)
  
  !size of array
  nx=3
  ny=10
  
  if(rank==0) then 
    if(num_procs>ny) then 
      print*, "Number of procs should be less than or equal to second dimension of the array"
      call MPI_Abort(mpi_comm_world, 1, mpierr)
    endif
  endif

  start_time=MPI_Wtime()
  !allocate in the root rank
  if(rank==0) then
    allocate(array(nx,ny))
  else !for other procs allocate with zero size
    allocate(array(0,0))
  endif
 
  !assign values to the array
  if(rank==0) then
    do j=1,ny
      do i=1,nx
        array(i,j) = (i-1)+(j-1)*nx
      end do 	
    end do 
    print*, "Before scattering..."
    print*, array
  endif
  

  !distribute the 3d array among different procs 
  call distribute_points(ny, rank, num_procs, str_idx, end_idx)
  local_size = end_idx - str_idx + 1
  
  !allocate local(for each rank) arrays
  allocate(array_local(nx, 0:local_size+1))

  !allocate sendcoutns and displacements arrays for braodcasting
  allocate(sendcounts(num_procs), displacements(num_procs))

  !gather displacements and sendcounts for all ranks
  call MPI_Allgather(str_idx, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, &
     MPI_COMM_WORLD, mpierr)
  call MPI_Allgather(local_size, 1, MPI_INTEGER, sendcounts, 1, &
    MPI_INTEGER, MPI_COMM_WORLD, mpierr)
  
  !total sendcounts and displacements 
  displacements = displacements - 1 !Array index starts with 0 in MPI (C)

  call MPI_Barrier(mpi_comm_world, mpierr)

  start_time=MPI_Wtime()
  !Scatterning using subarray type 
  sizes = [nx, ny]
  recv_sizes=[nx, local_size+2]
  sub_sizes = [nx, 1]
  starts = [0, 0] 
  recv_starts = [0, 1] 

  !to get extent of MPI_DOUBLE_PRECISION
  call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lb, extent, mpierr)
  
  !create a mpi subarray data type for sending data
  call MPI_Type_create_subarray(2, sizes, sub_sizes, starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, send_type, mpierr)
  lb_resize=0
  !resize the send subarray for starting at correct location for next send
  call MPI_Type_create_resized(send_type, lb_resize, nx*extent, &
    resize_send_type, mpierr)
  call MPI_Type_commit(resize_send_type, mpierr)

  !create a mpi subarray data type for receiving data
  call MPI_Type_create_subarray(2, recv_sizes, sub_sizes, recv_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, mpierr)
  
  !resize the receive subarray for starting at correct location for next receive
  call MPI_Type_create_resized(recv_type, lb_resize, nx*extent, &
    resize_recv_type, mpierr)
  call MPI_Type_commit(resize_recv_type, mpierr)
  
  call MPI_Barrier(mpi_comm_world, mpierr)
  start_time=MPI_Wtime()
  !scatter the subarrays 
   call MPI_Scatterv(array, sendcounts, displacements, resize_send_type, &
    array_local, sendcounts, resize_recv_type, 0, MPI_COMM_WORLD, mpierr)

  !if(rank==0) then 
  !  print*, "Time taken for scattering using MPI type subarrays: ", MPI_Wtime()-start_time 
  !endif
  call MPI_Barrier(mpi_comm_world, mpierr)
  !print the scattered array 
  print*, "Scattered array with subarray: ", rank
  print*, array_local 
  
  !do some computations on the scattered local arrays
  array_local = array_local+1
 
  call MPI_Barrier(mpi_comm_world, mpierr)
  start_time=MPI_Wtime()
  !Gather the local arrays to global (array) using the same subarrays 
  call MPI_Gatherv(array_local, local_size, resize_recv_type, array, &
    sendcounts, displacements, resize_send_type, 0, MPI_COMM_WORLD, mpierr)
  
  !if(rank==0) then 
  !  print*, "Time taken by MPI_Type_create_subarray Gathering: ", MPI_Wtime()-start_time 
  !endif
  
  if(rank==0) then 
    print*, "Gathered array: ------------------"
    print*, array
  endif
  call MPI_Finalize(mpierr)
end program ex_scatterv
