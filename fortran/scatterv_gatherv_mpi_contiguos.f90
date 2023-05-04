!!!column wise distribution of submatrices
!!!if you have (mxn) matrix; this program will 
!!!scatter (mx(n/nproc)) to different processors
!!!ex: if you have 3x5 matrix and you want to distribute among 2 processors
!!! then it will send 3x3 to rank 0 and 3x2 to rank 1
!!! this program demostrates the use of MPI_Type_contiguous datatype
program scatterv_gatherv_mpi_contiguous
  use mpi
  implicit none
    
  integer :: nx, ny, i, j
  integer, dimension(:,:), allocatable :: global, local
  integer, dimension(:), allocatable :: sendcounts, displacements
  integer :: rank, size, str_idx, end_idx, ierr, local_size
  integer :: newtype

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, size, ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
 
  nx=3
  ny=5

  call distribute_points(ny, rank, size, str_idx, end_idx)
  local_size = end_idx - str_idx + 1
  
  if(rank==0) then
    allocate(global(nx, ny))
    do i=1,nx 
      do j=1,ny
        global(i,j)=j+i*ny
      end do 
    end do
  endif
  
  if(rank==0) then
    print*, '--------global matrix----'
    do i=1,ubound(global, 1)
      print*, global(i,:)
    end do
  endif
 
  allocate(local(nx, local_size))
  allocate(sendcounts(size))
  allocate(displacements(size))

  print "(a5,i2,a16,i2,a16,i2,a22,i2)", "rank: ", rank, "str_idx: ", &
    str_idx, "end_idx: ", end_idx, "local_arr_size: ", local_size
  
  call MPI_Allgather(str_idx, 1, MPI_INTEGER, displacements, 1, MPI_INTEGER, &
    MPI_COMM_WORLD, ierr)
  call MPI_Allgather(local_size, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, &
    MPI_COMM_WORLD, ierr)
  
  displacements = displacements - 1

  if(rank==0) then
    do i=1,size
      print "(a13,i2,a16,i2)", "sendcounts: ", sendcounts(i), &
        "displacements: ", displacements(i)
    end do
  endif

 call MPI_Type_contiguous(nx, MPI_INTEGER, newtype, ierr)
 call MPI_Type_commit(newtype, ierr)

 call MPI_Scatterv(global, sendcounts, displacements, newtype, local, local_size, newtype, 0, MPI_COMM_WORLD, ierr)

print*, '------subarray--------'
print*, "rank: ", rank
print*, '----------------------'

      do i=1,ubound(local, 1)
        print*, local(i,:)
      end do 
 
call MPI_Finalize(ierr)
  
end program scatterv_gatherv_mpi_contiguous
