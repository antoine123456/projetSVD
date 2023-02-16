#include <QRpar.h>
#include <assert.h>

#include <mpi.h>

double RandDouble(double range_min, double range_max) {
    double r = ((double)rand() / (double) RAND_MAX) * (range_max-range_min) + range_min;
    return r;
}

// Represent a possible output from a simulation
// We assume the informations are about a field (u) and it's values on a grid with lot of points
// For a given time, u will then be of huge dimension, but the simulation will usually contain
// only about a 100 timesteps. The data matrix resulting will then be of dim (m,n) with m>>n
double *GenData(int m, int n) {
    double *mat = (double *) malloc(sizeof(double) * m*n);
    double std = RandDouble(0.1,0.2);
    double mu  = RandDouble(0.3,0.6);
    for (int j = 0; j < n; j++) {
        std += 0.1;
        mu  += 0.1;
        for (int i = 0; i < m; i++)
            mat[i*n + j] = exp(-pow(((double)i/m-mu)/std,2));
    }
    return mat;
}

// For a rectangular matrix of dim (m,n), we'll assume m>>n, the QR decomposition have the properties
// A = [Q1 Q2] |R|       so we can only compute Q1 of dim (m,n) and R of dim (n,n)
//             |0|
void ThinQR(QR_t *QR, double *A, int m, int n) {
    for (int k=0 ; k<n ; k++) {
        cblas_dcopy(m, &A[k], n, &QR->Q[k], n);

        for (int j=0 ; j<k ; j++) {
            QR->R[j*n + k] = cblas_ddot(m, &QR->Q[j], n, &QR->Q[k], n);

            cblas_daxpy(m, -QR->R[j*n + k], &QR->Q[j], n, &QR->Q[k], n);
        }

        QR->R[k*n + k] = cblas_dnrm2(m, &QR->Q[k], n);

        cblas_dscal(m, 1./QR->R[k*n + k], &QR->Q[k], n);
    }
}

double *QRparRoot(int m, int n, int np, int rank) {
    // Block length
    int bl = m/np;

    // Simulate data
    // In a real case scenario, the sub-blocks Ap of A would
    // be initialized with a sub-file of a simulation's output
    double *Ap = GenData(bl, n);


    // First QR factorization on blocks of A
    QR_t* QR = (QR_t*) malloc(sizeof(QR_t));
    QR->Q = (double *) malloc(sizeof(double) *bl*n);
    QR->R = (double *) malloc(sizeof(double) *n*n);
    ThinQR(QR, Ap, bl, n);

    // Assemble resulting QR->R in the matrix Rc
    double *Rc = (double *) malloc(sizeof(double) * np*n*n);
    MPI_Gather(
        QR->R,
        n*n,
        MPI_DOUBLE,
        Rc,
        n*n,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    // QR factorization of Rc, which will give the final R matrix
    QR->Q = (double*) realloc(QR->Q, np*n*n * sizeof(double));
    QR->R = (double*) realloc(QR->R, n*n * sizeof(double));
    ThinQR(QR, Rc, np*n, n);

    return QR->R;
}

void QRparOthers(int m, int n, int np, int rank) {
    // Block length
    int bl = m/np;

    // Simulate data
    // In a real case scenario, the sub-blocks Ap of A would
    // be initialized with a sub-file of a simulation's output
    double *Ap = GenData(bl, n);

    // First QR factorization on blocks of A
    QR_t* QR = (QR_t*) malloc(sizeof(QR_t));
    QR->Q = (double *) malloc(sizeof(double) *bl*n);
    QR->R = (double *) malloc(sizeof(double) *n*n);
    ThinQR(QR, Ap, bl, n);

    // Assemble resulting QR->R in the matrix Rc
    MPI_Gather(
        QR->R,
        n*n,
        MPI_DOUBLE,
        NULL,
        0,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );
}
