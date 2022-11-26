#include <SVD.h>

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

    for (int i=0 ; i<nB ; i++)
        Beigens.values[i] = sqrt(Beigens.values[i]);

    free(B);
    free(Beigens.vectors);
    return Beigens.values;
}