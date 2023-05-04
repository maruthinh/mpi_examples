!!! this program first transposes fortran multidimensional array 
!!!from column major to rwo major so that elements can be linearlized
!!!in rwo-major order. This can then be scattered to desired number of 
!!!processors. For examples if we have (5x4x3x2) array, then we will 
!!!transpose it to (2x3x4x5) and scatter. Scatterning is in rows, i,e
!!!if we want to scatter this to two ranks, then we send (3x4x3x2) to 
!!!rank 0 and (2x4x3x2) to rank 1

program scatterv_gatherv_mpi_contiguous
  use mpi
  implicit none
    
  integer :: n1, n2, n3, n4, i, j, k, l
  integer, dimension(:,:,:,:), allocatable :: global, global_c, local, local_c
  integer, dimension(:), allocatable :: sendcounts, displacements, local_flat, global_c_flat
  integer :: rank, size, str_idx, end_idx, ierr, local_size
  integer :: str_idx_flat, end_idx_flat, local_size_flat
  integer :: newtype

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, size, ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
 
  n1=5
  n2=4
  n3=3
  n4=2

  call distribute_points(n1, rank, size, str_idx, end_idx)
  local_size = end_idx - str_idx + 1
  local_size_flat = local_size*n2*n3*n4
  
  if(rank==0) then
    allocate(global(n1, n2, n3, n4))
    do i=1,n1 
      do j=1,n2
        do k=1,n3
          do l=1,n4
            global(i,j,k,l)=(i-1)+(j-1)*n1+(k-1)*n1*n2+(l-1)*n1*n2*n3
          end do 
        end do
      end do 
    end do
    allocate(global_c(n4,n3,n2,n1))
    global_c = reshape(global, shape(global_c), order=[4,3,2,1])
    global_c_flat = pack(global_c,.true.)
    print*, '===========global================='
    print*, global
  endif

!  if(rank==0) then
!    allocate(global_flat(n1*n2*n3*n4))
!    global_flat = pack(global,.true.)
!    print*, '------flat F array--------'
!    print*, global_flat
!    print*, '------flat C array--------'
!    global_c_flat = pack(global_c,.true.)
!    print*, global_c_flat
!
!    !print*, '----array(1:3,1:n2,1:n3,1:n4)-----'
!    !print*, pack(global(1:3,1:n2,1:n3,1:n4),.true.)
!    print*, '----array(4:5,1:n2,1:n3,1:n4)-----'
!    print*, pack(global_c(1:n4,1:n3,1:n2,1:3),.true.)
!  endif 

   !call distribute_points(n1*n2*n3*n4, rank, size, str_idx_flat, end_idx_flat)
   !local_size_flat = end_idx_flat - str_idx_flat + 1

   allocate(local(local_size, n2, n3, n4))
   allocate(local_c(n4, n3, n2, local_size))
   allocate(local_flat(local_size*n2*n3*n4))
   allocate(sendcounts(size))
   allocate(displacements(size))

   print "(a5,i2,a16,i2,a16,i2,a22,i2)", "rank: ", rank, "str_idx: ", &
    str_idx, "end_idx: ", end_idx, "local_arr_size: ", local_size
  
  call MPI_Allgather(str_idx, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, &
    MPI_COMM_WORLD, ierr)
  call MPI_Allgather(local_size, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, &
    MPI_COMM_WORLD, ierr)
    
    sendcounts = sendcounts*n2*n3*n4
    displacements = displacements - 1
 
    displacements = displacements*n2*n3*n4
 
!  do i = 1, size 
!    sendcounts(i) = sendcounts(i)*n2*n3*n4
!    displacements(i) = displacements(i) - 1
!  end do
! 
!  do i = 1, size 
!    displacements(i) = displacements(i)*n2*n3*n4
!  end do

 if(rank==0) then
    do i=1,size
      print "(a13,i3,a16,i3)", "sendcounts: ", sendcounts(i), &
        "displacements: ", displacements(i)
    end do
  endif

 call MPI_Scatterv(global_c_flat, sendcounts, displacements, MPI_INTEGER, local_flat, &
   local_size*n2*n3*n4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

!print*, '------subarray--------'
!print*, "rank: ", rank
!print*, '----------------------'
!print*, local_flat

local_c = reshape(local_flat, [n4,n3,n2,local_size])
local = reshape(local_c, shape(local), order=[4,3,2,1])

!print*, '--------subarray F------'
!print*, local

local_c = reshape(local, shape(local_c), order=[4,3,2,1])
local_flat = pack(local_c,.true.)

!print*, '--------subarray flat------'
!print*, local_flat 

call MPI_Barrier(MPI_COMM_WORLD, ierr)

call MPI_Gatherv(local_flat, local_size_flat, MPI_INTEGER, global_c_flat, sendcounts, &
  displacements, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 

if(rank==0) then
  global_c = reshape(global_c_flat, [n4,n3,n2,n1])
  global = reshape(global_c, shape(global), order=[4,3,2,1])
 
  print*, '==========global after gatherv========='
  print*, global
endif

  call MPI_Finalize(ierr)
  
end program scatterv_gatherv_mpi_contiguous
