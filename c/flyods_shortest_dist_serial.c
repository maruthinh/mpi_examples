#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a, b) a < b ? a : b;

int main(int argc, char *argv[]) {
  int rank, size;

  int **A, m, n;
  m = 6;
  n = 6;
  A = (int **)malloc(m * sizeof(int *));
  for (int i = 0; i < m; i++) {
    A[i] = (int *)malloc(n * sizeof(int));
  }

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = 10000;
    }
  }

  // Diagonal elements. Distance to itself.
  A[0][0] = 0;
  A[1][1] = 0;
  A[2][2] = 0;
  A[3][3] = 0;
  A[4][4] = 0;
  A[5][5] = 0;

  A[0][1] = 2;
  A[0][2] = 5;

  A[1][2] = 7;
  A[1][3] = 1;
  A[1][5] = 8;

  A[2][3] = 4;

  A[3][4] = 3;

  A[4][2] = 2;
  A[4][5] = 3;

  A[5][1] = 5;
  A[5][3] = 2;
  A[5][4] = 4;

  for (int k = 0; k < m; k++) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = min(A[i][j], A[i][k] + A[k][j]);
      }
    }
  }

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%d\t", A[i][j]);
    }
    printf("\n");
  }

  for (int i = 0; i < m; i++) {
    free(A[i]);
  }
  free(A);

  return 0;
}