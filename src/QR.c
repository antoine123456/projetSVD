#include <QR.h>

#define NITER 200    // Nb max of QR method iterations
#define EPS 0.001  // Error threshold

/* Return the maximal error on choosing the diagonal
 * elements of A as its eigenvalues */
double GershgorinTest(double *A, int n) {
    double errmax = 0;
    for (int j=0 ; j<n ; j++) {
        double err = 0;
        for (int i=j+1 ; i<n ; i++)
            err += fabs(A[i*n + j]);
        if (err > errmax)
            errmax = err;
    }
    return errmax;
}

/* Return the maximal error on choosing the diagonal
 * elements of A as its eigenvalues */
double GershgorinColTest(double *A, int j, int n) {
    double err = 0;
    for (int i=1 ; i<n ; i++)
            err += fabs(A[i*n + j]);
    return err;
}

/* QR iteration method :
 * A modified such that its' trangular sup
 * elements are its eigenvalues */
eigen_t QR_method(double *A, int n) {
    QR_t QR;
    QR.Q = (double *) malloc(sizeof(double )* n*n);
    QR.R = (double *) malloc(sizeof(double )* n*n);

    eigen_t eigen;
    eigen.values  = malloc(sizeof(double) * n);

    for (int k=0 ; k<NITER ; k++) {
        // QR decomposition
        GramSchmidtMod(&QR, A, n);

        // Akp1 calcul
        MatMul(A, QR.R, QR.Q, n);

        // Error check
        double err_max = GershgorinTest(A, n);
        if (err_max < EPS)
            break;

    }

    for (int i=0 ; i<n ; i++)
        eigen.values[i] = A[i*n + i];
    eigen.vectors = QR.Q;

    free(QR.R);

    return eigen;
}

/* QR iteration method :
 * A modified such that its' trangular sup
 * elements are its eigenvalues */
eigen_t my_QR_method(double *A, int n) {
    QR_t QR;
    QR.Q = (double *) malloc(sizeof(double )* n*n);
    QR.R = (double *) malloc(sizeof(double )* n*n);

    eigen_t eigen;
    eigen.values  = malloc(sizeof(double) * n);

    for (int k=0 ; k<NITER ; k++) {
        // QR decomposition
        GramSchmidtMod(&QR, A, n);

        // Akp1 calcul
        MatMul(A, QR.R, QR.Q, n);

        // Error check
        double err_max = GershgorinTest(A, n);
        if (err_max < EPS)
            break;

    }

    for (int i=0 ; i<n ; i++)
        eigen.values[i] = QR.R[i*n + i];
    eigen.vectors = QR.Q;

    free(QR.R);

    return eigen;
}

/* QR iteration method for a tridiagonal symmetric matrix :
 * A modified such that its' trangular sup
 * elements are its eigenvalues */
eigen_t QR_method_Tridiag(double *A, int n) {
    QR_t QR;
    QR.Q = (double *) calloc(n*n, sizeof(double ));
    QR.R = (double *) calloc(n*n, sizeof(double ));

    eigen_t eigen;
    eigen.values  = (double*) calloc(n, sizeof(double));

    for (int k=0 ; k<NITER ; k++) {
        // QR decomposition
        GramSchmidtMod_Tridiag(&QR, A, n);

        // Akp1 = RQ
        MatMul_Tridiag(A, QR.R, QR.Q, n);

        // Error check
        double err_max = GershgorinTest(A, n);
        if (err_max < EPS)
            break;
    }


    for (int i=0 ; i<n ; i++){
      eigen.values[i] = A[i*n + i];
    }
    eigen.vectors = QR.Q;

    free(QR.R);

    return eigen;
}


void freeEigen(eigen_t e) {
    free(e.values);
    free(e.vectors);
}
