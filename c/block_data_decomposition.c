/*Block data decomposition methods from Quinn Parallel Programming Book*/
#include "mpi.h"
#include <stdio.h>
#define min(a, b) (a < b ? a : b)

/// @brief Function to distribute n elements among nprocs
/// @param n Number of elements to be distributed
/// @param rank Processor rank
/// @param nprocs Number of processors
/// @param str_idx Starting index of data block which rank holds
/// @param end_idx Ending index of data block which rank holds
/// @param idx Global array element index from a data
/// @param proc The rank id which holds the idx
void block_data_decomposition_method1(const int n, const int rank, const int nprocs, int *str_idx, int *end_idx, const int idx, int *proc)
{
    int r;
    r = n % nprocs;
    *str_idx = rank * (n / nprocs) + min(rank, r);
    *end_idx = (rank + 1) * (n / nprocs) + min(rank + 1, r) - 1;
    *proc = min(idx / (n / nprocs + 1), (idx - r) / (n / nprocs));
}

/// @brief Function to distribute n elements among nprocs
/// @param n Number of elements to be distributed
/// @param rank Processor rank
/// @param nprocs Number of processors
/// @param str_idx Starting index of data block which rank holds
/// @param end_idx Ending index of data block which rank holds
/// @param idx Global array element index from a data
/// @param proc The rank id which holds the idx
void block_data_decomposition_method2(const int n, const int rank, const int nprocs, int *str_idx, int *end_idx, const int idx, int *proc)
{
    *str_idx = (rank * n) / nprocs;
    *end_idx = ((rank + 1) * n / nprocs) - 1;
    *proc = (nprocs * (idx + 1) - 1) / n;
}

int main(int argc, char *argv[])
{
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
    block_data_decomposition_method1(n, rank, size, &str_idx, &end_idx, idx, &proc);

    printf("Block Data Decomposition Method 1: %d elements are distributed among %d procs with rank %d having starting index %d and end index %d\n", n, size, rank, str_idx, end_idx);
    fflush(stdout);
    if (rank == 0)
    {
        printf("Block Data Decomposition Method 1: rank %d has %d th element \n", proc, idx);
    }
    fflush(stdout);
    //==========================================================================
    // Block Data Decomposition Method 2
    block_data_decomposition_method2(n, rank, size, &str_idx, &end_idx, idx, &proc);

    printf("Block Data Decomposition Method 2: %d elements are distributed among %d procs with rank %d having starting index %d and end index %d\n", n, size, rank, str_idx, end_idx);
    fflush(stdout);
    if (rank == 0)
    {
        printf("Block Data Decomposition Method 2: rank %d has %d th element \n", proc, idx);
    }
    fflush(stdout);

    MPI_Finalize();
    return 0;
}
