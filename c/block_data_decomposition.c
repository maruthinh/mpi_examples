/*Block data decomposition methods from Quinn Parallel Programming Book*/
#include "block_data_decomposition.h"
#define min(a, b) (a < b ? a : b)
/// @brief Function to distribute n elements among nprocs
/// @param n Number of elements to be distributed
/// @param rank Processor rank
/// @param nprocs Number of processors
/// @param str_idx Starting index of data block which rank holds
/// @param end_idx Ending index of data block which rank holds
/// @param idx Global array element index from a data
/// @param proc The rank id which holds the idx
void block_data_decomposition_method1(const int n, const int rank,
                                      const int nprocs, int *str_idx,
                                      int *end_idx, const int idx, int *proc) {
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
void block_data_decomposition_method2(const int n, const int rank,
                                      const int nprocs, int *str_idx,
                                      int *end_idx, const int idx, int *proc) {
  *str_idx = rank * (n / nprocs);
  *end_idx = ((rank + 1) * n / nprocs) - 1;
  *proc = (nprocs * (idx + 1) - 1) / n;
}

/// @brief Return starting index of block that belongs to a processor
/// @param n Number of elements to be distributed
/// @param rank Processor rank
/// @param nprocs Number of processors
/// @return Starting index of a block element that belongs to a processor
int block_low(const int n, const int rank, const int nprocs) {
  return rank * n / nprocs;
}

/// @brief Return end index of block that belongs to a processor
/// @param n Number of elements to be distributed
/// @param rank Processor rank
/// @param nprocs Number of processors
/// @return End index of a block element that belongs to a processor
int block_high(const int n, const int rank, const int nprocs) {
  return block_low(n, rank + 1, nprocs) - 1;
}

/// @brief Return the block size that belongs to a processor
/// @param n Number of elements to be distributed
/// @param rank Processor rank
/// @param nprocs Number of processors
/// @return The block size that belongs to a processor
int block_size(const int n, const int rank, const int nprocs) {
  return block_high(n, rank, nprocs) - block_low(n, rank, nprocs) + 1;
}

/// @brief Return the processor which owns the index of a block
/// @param n Number of elements to be distributed
/// @param index The index of element in a block
/// @param nprocs Number of processors
/// @return The processor which owns the index of a block
int block_owner(const int n, const int index, const int nprocs) {
  return (nprocs * (index + 1) - 1) / n;
}