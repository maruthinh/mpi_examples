/*Block data decomposition methods from Quinn Parallel Programming Book*/
#include "block_data_decomposition.h"
#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a, b) (a < b ? a : b)

int main(int argc, char *argv[]) {
  int rank, size, n;
  int blk_low, blk_high, blk_size;

  MPI_Init(&argc, &argv);

  MPI_Barrier(MPI_COMM_WORLD);

  double elapsed_time = -MPI_Wtime();

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc != 2) {
    if (!rank) {
      printf("Command line: %s <m>\n", argv[0]);
    }
    MPI_Finalize();
    exit(1);
  }

  n = atoi(argv[1]);
  // Domain decomposition

  blk_low = 2 + block_low(n - 1, rank, size);
  blk_high = 2 + block_high(n - 1, rank, size);
  blk_size = block_size(n - 1, rank, size);

  // printf("%d %d %d %d\n", rank, blk_low, blk_high, blk_size);

  // If all the primes used for sieving are not held by proc 0, terminate
  // Max sieving number is sqrt(n)
  int proc0_size = (n - 1) / size;
  if ((1 + proc0_size) < (int)sqrt((double)n)) {
    if (!rank) {
      printf("Too many procs!\n. Reduce the number of procs\n");
    }
    MPI_Finalize();
    exit(1);
  }

  // portion of 2...n
  char *marked;
  marked = (char *)malloc(blk_size);
  if (marked == NULL) {
    printf("Allocation error of char array marked!");
    MPI_Finalize();
    exit(1);
  }

  // initialize marked char array
  for (int i = 0; i < blk_size; i++) {
    marked[i] = 0;
  }

  int idx, prime, first, count, global_count;
  idx = 0;
  prime = 0;
  first = 0;
  global_count = 0;
  // prime is the current prime used for sieving
  // idx is the index of the prime in the array of processor 0
  prime = 2;
  if (!rank) {
    idx = 0;
  }
  // first is the index of first integer that needs to be marked
  do {
    if (prime * prime > blk_low)
      first = prime * prime - blk_low;
    else {
      if (!(blk_low % prime))
        first = 0;
      else
        first = prime - (blk_low % prime);
    }
    // printf("%d first %d prime \n", first, prime);
    for (int i = first; i < blk_size; i += prime) {
      marked[i] = 1;
    }
    if (!rank) {
      while (marked[++idx])
        prime = idx + 1;
        // printf("Prime computed %d \n", prime);
    }
    MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
  } while (prime * prime <= n); // sieve till prime^2 <= n

  count = 0;
  for (int i = 0; i < blk_size; i++) {
    if (!marked[i]) {
      count++;
    }
  }

  // for(int i=0; i<blk_size; i++){
    // printf("%d index %d rank and %d marker \n", i, rank, marked[i]);
  // }

  MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  /* Stop the timer */
  elapsed_time += MPI_Wtime();

  if (!rank) {
    printf("%d primes are less than or equal to %d\n", global_count, n);
    printf("Total elapsed time: %10.6f\n", elapsed_time);
  }
  MPI_Finalize();
  return 0;
}