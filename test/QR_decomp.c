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
    QR.Q = (double *) malloc(sizeof(double )* n*n);
    QR.R = (double *) malloc(sizeof(double )* n*n);
    double *C = (double*) malloc(sizeof(double) * n*n);

    double *givens_mat = givens(n,n);
    GramSchmidtMod(&QR, givens_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - givens_mat[i]) > eps) {
            Print_Error_Mat("QR Decomposition", "Givens");
            nerr++;
        }

    double *aegerter_mat = aegerter(n);
    GramSchmidtMod(&QR, aegerter_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - aegerter_mat[i]) > eps) {
            Print_Error_Mat("QR Decomposition", "Aegerter");
            nerr++;
        }

    double *orth_symm_mat = orth_symm(n);
    GramSchmidtMod(&QR, orth_symm_mat, n);
    MatMul(C, QR.Q, QR.R, n);
    for (int i=0 ; i<n*n ; i++)
        if (fabs(C[i] - orth_symm_mat[i]) > eps) {
            Print_Error_Mat("QR Decomposition", "Orth_Symm");
            nerr++;
        }
    
    if (nerr == 0)
        Print_Success("QR Decomposition");

    return 0;
}