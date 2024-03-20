/* Mid point rule of integration
∫f(x)dx can be numrically approximated using the following formula
∫f(x)dx=∑f(m)Δx, where Δx=(b−a)/n; n is the number of intervals, b is upper
limit and a is an lower limit of integration. And the mid point m is calculated
as m=(b+a)/2
*/

#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    int rank, size, local_n, n;
    double lim1, lim2, dx, pi, global_val_pi, mid_pt;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // number of intervals
    n = 10000000;
    lim1 = 0.0;
    lim2 = 1.0;
    dx = (lim2 - lim1) / n;

    // serial version
    // pi = 0.0;
    // for (int i = 0; i < n; i++)
    // {
    //     lim1 = i * dx;
    //     lim2 = lim1 + dx;
    //     mid_pt = (lim2 + lim1) / 2.0;
    //     pi += (4.0 / (1 + mid_pt * mid_pt)) * dx;
    // }

    //parallel version
    for (int i = rank; i < n; i+=size)
    {
        lim1 = i * dx;
        lim2 = lim1 + dx;
        mid_pt = (lim2 + lim1) / 2.0;
        pi += (4.0 / (1 + mid_pt * mid_pt)) * dx;
    }

    MPI_Reduce(&pi, &global_val_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("The value of Pi obtained using Mid Point Integration Rule is: %.17g \n", global_val_pi);
    }
    MPI_Finalize();
    return 0;
}
