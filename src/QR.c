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

void qr_decomposition(double *A, int n, double *Q, double *R) {
  int i, j, k;
  double norm =0;

  for (k = 0; k < n; k++) {
    norm = 0.0;
    for (i = k; i < n; i++) {
      norm = norm + A[i*n+k] * A[i*n+k];
    }
    norm = sqrt(norm);
    if (norm != 0.0) {
      if (A[k*n+k] < 0) {
        norm = -norm;
      }
      for (i = k; i < n; i++) {
        A[i*n+k] = A[i*n+k] / norm;
      }
      A[k*n+k] = A[k*n+k] + 1.0;
      for (j = k + 1; j < n; j++) {
        norm = 0.0;
        for (i = k; i < n; i++) {
          norm = norm + A[i*n+k] * A[i*n+j];
        }
        norm = norm / A[k*n+k];
        for (i = k; i < n; i++) {
          A[i*k+j] = A[i*n+j] - norm * A[i*n+k];
        }
      }
    }
    for (i = 0; i < n; i++) {
      Q[i*n+k] = A[i*n+k];
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i > j) {
        R[i*n+j] = 0.0;
      } else {
        R[i*n+j] = A[i*n+j];
      }
    }
  }
}
/* QR iteration method :
 * A modified such that its' trangular sup
 * elements are its eigenvalues */
QR_t my_QR_method1(double *A, int n) {
    QR_t QR;
    QR.Q = (double *) calloc(n*n,sizeof(double ) );
    QR.R = (double *) calloc(n*n,sizeof(double ));

  /*  eigen_t eigen;
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
    */
    qr_decomposition(A, n, QR.Q, QR.R);
    //free(QR.R);

    return QR;
}


/* QR iteration method :
 * A modified such that its' trangular sup
 * elements are its eigenvalues */
eigen_t my_QR_method2(double *A, int n) {
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
