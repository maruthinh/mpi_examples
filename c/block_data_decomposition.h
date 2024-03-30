#pragma once

void block_data_decomposition_method1(const int n, const int rank,
                                      const int nprocs, int *str_idx,
                                      int *end_idx, const int idx, int *proc);

void block_data_decomposition_method2(const int n, const int rank,
                                      const int nprocs, int *str_idx,
                                      int *end_idx, const int idx, int *proc);

int block_low(const int n, const int rank, const int nprocs);

int block_high(const int n, const int rank, const int nprocs);

int block_size(const int n, const int rank, const int nprocs);

int block_owner(const int n, const int index, const int nprocs);