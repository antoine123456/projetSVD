#include <QR.h>

#define NITER 200    // Nb max of QR method iterations
#define EPS 0.00001  // Error threshold

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
 * A modified such that its' diagonal
 * elements are its eigenvalues */
eigen_t QR_method(double *A, int n) {
    QR_t QR;
    eigen_t eigen;
    eigen.values  = malloc(sizeof(double) * n);
    double *Akp1  = malloc(sizeof(double) * n*n);

    for (int k=0 ; k<NITER ; k++) {
        // QR decomposition
        QR = GramSchmidtMod(A, n);

        // Akp1 calcul
        MatMul(Akp1, QR.R, QR.Q, n);

        // Error check
        double err_max = GershgorinTest(Akp1, n);
        #ifdef INFO
            printf("%lf\n", err);
        #endif
        if (err_max < EPS) {
            free(QR.R);
            Copy(A, Akp1, n);
            break;
        }

        // A = Akp1
        Copy(A, Akp1, n);

        free(QR.R);
        if (k < NITER-1)
            free(QR.Q);
    }

    free(Akp1);

    for (int i=0 ; i<n ; i++)
        eigen.values[i] = A[i*n + i];
    eigen.vectors = QR.Q;
    return eigen;
}