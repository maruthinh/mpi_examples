!Example code to send and receive part of 3D array using MPI subarray type 
program bcast
  use mpi  
  implicit none
  
  real(kind=8),allocatable,dimension(:,:,:) :: array, array_local
  real(kind=8),allocatable,dimension(:) :: array_flat, array_local_flat
  integer :: sizes(3), subsizes(3), starts_left(3), starts_right(3), &
    recv_left(3), recv_right(3)
  integer :: rank, num_procs, mpierr, i, j, k, left_bound, right_bound, &
    status(mpi_status_size)
  integer :: nx, ny, nz, str_idx, end_idx, local_size, local_size_flat 
  integer, dimension(:), allocatable :: sendcounts, displacements

  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world, num_procs, mpierr)
  call mpi_comm_rank(mpi_comm_world, rank, mpierr)
  
  if(num_procs>2) then 
    print*, "Please use only 2 procs, although it can be easily extended to &
        work for more procs"
    stop 
    call MPI_Finalize(mpierr)
  endif 

  nx=5
  ny=3
  nz=2

  !allocate in the root rank
  if(rank==0) then
    allocate(array(nx,ny,nz))
    allocate(array_flat(nx*ny*nz))
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
    
    print*, "-------Array that's split, braodcasted and communicated---------"
    print*, array
    
    !flatten the 3D array
    forall(k=1:nz, j=1:ny, i=1:nx) array_flat(k+(j-1)*nz+(i-1)*ny*nz)=array(i,j,k)
  endif

  !distribute the 3d array among different procs 
  call distribute_points(nx, rank, num_procs, str_idx, end_idx)
  local_size = end_idx - str_idx + 1
  local_size_flat = local_size*ny*nz
  
  !allocate local(for each rank) arrays
  allocate(array_local_flat(local_size_flat))
  allocate(array_local(local_size+2, ny, nz))

  !allocate sendcoutns and displacements arrays for braodcasting
  allocate(sendcounts(num_procs), displacements(num_procs))

  !gather displacements and sendcounts for all ranks
  call MPI_Allgather(str_idx, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, &
     MPI_COMM_WORLD, mpierr)
  call MPI_Allgather(local_size, 1, MPI_INTEGER, sendcounts, 1, &
    MPI_INTEGER, MPI_COMM_WORLD, mpierr)
  
  !total sendcounts and displacements 
  sendcounts = sendcounts*ny*nz
  displacements = displacements - 1 !Array index starts with 0 in MPI (C)
  displacements = displacements*ny*nz

  !scatter the flattened array among procs 
  call MPI_Scatterv(array_flat, sendcounts, displacements, MPI_REAL8, &
    array_local_flat, local_size*ny*nz, MPI_REAL8, 0, MPI_COMM_WORLD, &
    mpierr)
  
  !form 3D array from flattened local array 
  forall(k=1:nz, j=1:ny, i=2:local_size+1) array_local(i,j,k) = &
      array_local_flat(k+(j-1)*nz+(i-2)*ny*nz)

 
  if(rank==1) then 
  print*,"----------------------After bcast----------------------------------"
    print*, "rank: ", rank
    print*, array_local(2, :, :) 
  endif 
  
  !create sizes and subsizes for MPI subarray datatype
  sizes=[local_size+2,ny,nz]
  subsizes=[1,ny,nz]
  starts_left=[1,0,0]
  starts_right=[local_size+1, 0, 0]

  !create and commit subarray datatype
  call MPI_Type_create_subarray(3, sizes, subsizes, starts_left, MPI_ORDER_FORTRAN, &
    MPI_REAL8, left_bound, mpierr)
  call MPI_Type_commit(left_bound, mpierr)
  call MPI_Type_create_subarray(3, sizes, subsizes, starts_right, MPI_ORDER_FORTRAN, &
    MPI_REAL8, right_bound, mpierr)
  call MPI_Type_commit(right_bound, mpierr)

  !send subarray from rank 1 and receive that to rank 0's halo 
  if (rank==0) then 
    call MPI_Recv(array_local, 1, right_bound, 1, 100, mpi_comm_world, status, mpierr)
  elseif(rank==1) then
    call MPI_Send(array_local, 1, left_bound, 0, 100, mpi_comm_world, mpierr)
  endif
 
  if(rank==0) then 
  print*,"----------------------After receiving ------------------------------"
    print*, "rank: ", rank 
    print*, array_local(local_size+2, :, :)
  endif
  
  !delete MPI_Types
  call MPI_Type_Free(left_bound, mpierr)
  call MPI_Type_Free(right_bound, mpierr)

  call MPI_Finalize(mpierr)
end program bcast
