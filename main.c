#include <stdio.h>
#include <QRpar.h>
#include <gen_mat.h>
#include <io_matrix.h>
#include <openmpi-x86_64/mpi.h>

int main() {
    MPI_Init(NULL,NULL);
    int np, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *R;
    QRpar(R, np, rank, 8, 2);
    if (rank==0)
        PrintMat(R, 2, 2);

    MPI_Finalize();
}
