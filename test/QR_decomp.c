#include <math.h>

#include <gram_schmidt.h>
#include <io_matrix.h>
#include "test_mat.h"
#include "test_utils.h"

int main() {
    // Test parameters
    int n = 10;        // matrix order
    double eps = 0.01; // error threshold

    int nerr = 0;
    QR_t QR;
    double *C = (double*) malloc(sizeof(double) * n*n);

    double *givens_mat = givens(n,n);
    QR = GramSchmidtMod(givens_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - givens_mat[i]) > eps)
            nerr++;

    double *aegerter_mat = aegerter(n);
    QR = GramSchmidtMod(aegerter_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - aegerter_mat[i]) > eps)
            nerr++;

    double *orth_symm_mat = orth_symm(n);
    QR = GramSchmidtMod(orth_symm_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - orth_symm_mat[i]) > eps)
            nerr++;

    if (nerr != 0)
        Print_Error_Mat("QR_decomp", "Orth_Symm");

    return 0;
}