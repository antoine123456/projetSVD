#include <SVD.h>

#include <io_matrix.h>

double* SVD_1(double *A, int m, int n) {
    double *B;  // will store A^TA or AA^T
    int nB;    // B dimension

    // Compute A^TA or AA^T depending on the comparaison between n and m
    if (n < m) {
        B = get_AtA(A, m, n);
        nB = n;
    } else {
        B = get_AAt(A, m, n);
        nB = m;
    }

    eigen_t Beigens = QR_method(B, nB);

    free(B);
    free(Beigens.vectors);
    
    // Beigen values are A singular values
    return Beigens.values;
}

double* SVD_1_Tridiag(double *A, int m, int n) {
    double *B;  // will store A^TA or AA^T
    int nB;    // B dimension

    // Compute A^TA or AA^T depending on the comparaison between n and m
    if (n < m) {
        B = get_AtA(A, m, n);
        nB = n;
    } else {
        B = get_AAt(A, m, n);
        nB = m;
    }

    eigen_t Beigens = QR_method_Tridiag(B, nB);

    free(B);
    free(Beigens.vectors);
    
    // Beigen values are A singular values
    return Beigens.values;
}

// Call SVD_1 with A transformed into a bidiagonale matrix
double* SVD_3(double *A, int m, int n) {
    // Bidiagonalization (GKL)
    GKL_Bidiag_t b = Bidiagonalization(A, m, n);

    double* singular_values = SVD_1_Tridiag(b.B, m, n);

    freeGKL(b);

    return singular_values;
}