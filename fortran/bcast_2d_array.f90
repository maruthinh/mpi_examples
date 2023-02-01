!code to broadcast 2d fortran array using mpi datatype
program bcast
  use mpi  
  implicit none
  
  real(kind=8) :: array(100,100)
  integer :: sizes(2), subsizes(2), starts(2)
  integer :: rank, size, mpierr, i, j, newtype

  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world, size, mpierr)
  call mpi_comm_rank(mpi_comm_world, rank, mpierr)

  sizes=[100,100]
  subsizes=[100,100]
  starts=[0,0]

  if(rank==0) then
    do i=1,100 
      do j=1,100
        array(i,j) = (i-1)+(j-1)*100
      end do
    end do 
  endif
  
  !mpi datatype to create a subarray which helps retain same array shape and
  !order without reordering or flattening
  call MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
    MPI_REAL8, newtype, mpierr)
  call MPI_Type_commit(newtype, mpierr)

  call MPI_Bcast(array, 1, newtype, 0, MPI_COMM_WORLD, mpierr) 
  
  print*, "rank: ", rank , array(1, 1:5)
  
  call MPI_Finalize(mpierr)
end program bcast
