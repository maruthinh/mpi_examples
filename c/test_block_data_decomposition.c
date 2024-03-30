/*Block data decomposition methods from Quinn Parallel Programming Book*/
#include "block_data_decomposition.h"
#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
  int rank, size, n, str_idx, end_idx, idx, proc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  n = 14;
  str_idx = 0;
  end_idx = 0;
  idx = 9;
  proc = 0;

  // Block Data Decomposition Method 1
  block_data_decomposition_method1(n, rank, size, &str_idx, &end_idx, idx,
                                   &proc);

  printf("Block Data Decomposition Method 1: %d elements are distributed among "
         "%d procs with rank %d having starting index %d and end index %d\n",
         n, size, rank, str_idx, end_idx);
  fflush(stdout);
  if (rank == 0) {
    printf("Block Data Decomposition Method 1: rank %d has %d th element \n",
           proc, idx);
  }
  fflush(stdout);
  //==========================================================================
  // Block Data Decomposition Method 2
  block_data_decomposition_method2(n, rank, size, &str_idx, &end_idx, idx,
                                   &proc);

  printf("Block Data Decomposition Method 2: %d elements are distributed among "
         "%d procs with rank %d having starting index %d and end index %d\n",
         n, size, rank, str_idx, end_idx);
  fflush(stdout);
  if (rank == 0) {
    printf("Block Data Decomposition Method 2: rank %d has %d th element \n",
           proc, idx);
  }
  fflush(stdout);

  //==========================================================================
  // Block Data Decomposition Method 2: Individual functions

  printf("=================================================================\n");
  int blk_low = block_low(n, rank, size);
  int blk_high = block_high(n, rank, size);
  int blk_size = block_size(n, rank, size);
  int blk_owner = block_owner(n, idx, size);

  printf("Block Data Decomposition Method 2: %d elements are distributed among "
         "%d procs with rank %d having starting index %d and end index %d\n "
         "with size %d\n",
         n, size, rank, blk_low, blk_high, blk_size);
  printf("Block Data Decomposition Method 2: rank %d has %d th element \n",
         blk_owner, idx);

  MPI_Finalize();
  return 0;
}