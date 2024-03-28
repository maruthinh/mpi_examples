/* Example from Michael J. Quinn book
There will be a circuit containing AND, OR and NOT gates. The circuit
satisfiability problem is to determine if some combination of inputs (0, 1)
causes the output of the circuit to be 1. There will be 16 inputs. So total
combinations will be 2^16=65536
*/
#include "mpi.h"
#include <stdio.h>

#define EXTRACT_BIT(n, i) ((n & (1 << i)) ? 1 : 0)
/// @brief
/// @param id
/// @param z
int check_circuit(int id, int z) {
  int v[16];

  for (int i = 0; i < 16; i++) {
    v[i] = EXTRACT_BIT(z, i);
  }
  if ((v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3]) &&
      (!v[3] || !v[4]) && (v[4] || !v[5]) && (v[5] || !v[6]) &&
      (v[5] || v[6]) && (v[6] || !v[15]) && (v[7] || !v[8]) &&
      (!v[7] || !v[13]) && (v[8] || v[9]) && (v[8] || !v[9]) &&
      (!v[9] || !v[10]) && (v[9] || v[11]) && (v[10] || v[11]) &&
      (v[12] || v[13]) && (v[13] || !v[14]) && (v[14] || v[15])) {
    printf(" %d) %d%d% d%d%d%d%d%d%d%d%d%d%d%d%d%d\n", id, v[0], v[1], v[2],
           v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13],
           v[14], v[15]);
    fflush(stdout);
    return 1;
  } else {
    return 0;
  }
}

int main(int argc, char *argv[]) {
  int i, id, p;
  int local_solns, global_solns;
  double elapsed_time;

  MPI_Init(&argc, &argv);

  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();

  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  local_solns = 0;
  for (i = id; i < 65536; i += p) {
    local_solns += check_circuit(id, i);
  }

  MPI_Reduce(&local_solns, &global_solns, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  elapsed_time += MPI_Wtime();
  printf("Time taken by %d processor is %lf: ", p, elapsed_time);
  printf("Process %d is done\n", id);
  fflush(stdout);
  MPI_Finalize();
  if (id == 0)
    printf("There are %d different solutions\n", global_solns);
  return 0;
}