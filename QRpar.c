#include <openmpi-x86_64/mpi.h>
#include <stdio.h>
#include <io_matrix.h>
#include <QRpar.h>

int main() {
    MPI_Init(NULL,NULL);
    int np, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int m = 1000;
    int n = 10;

    if (rank==0) {
        double *R = QRparRoot(m, n, np, rank);
        PrintMat(R, n, n);
    } else {
        QRparOthers(m, n, np, rank);
    }

    MPI_Finalize();

    return 0;
}