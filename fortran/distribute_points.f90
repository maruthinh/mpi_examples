subroutine distribute_points(npts, rank, size, start_idx, end_idx)
  implicit none
  
  integer, intent(in) :: npts, size, rank
  integer, intent(out) :: start_idx, end_idx
  integer :: pts_per_proc

  pts_per_proc = npts/size
  
  if(rank < mod(npts, size)) then
    pts_per_proc=pts_per_proc + 1
  end if
 
  if(rank < mod(npts, size)) then
    start_idx = rank * pts_per_proc + 1
    end_idx = (rank + 1) * pts_per_proc
  else
    start_idx = mod(npts, size) + rank*pts_per_proc + 1
    end_idx   = mod(npts, size) + (rank + 1) * pts_per_proc
  end if

end subroutine distribute_points
