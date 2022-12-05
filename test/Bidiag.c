#include <bidiagonalization.h>
#include <gen_mat.h>
#include <io_matrix.h>
#include <assert.h>
#include <unistd.h>

#include "test_mat.h"
#include "test_utils.h"

int main(int argc, char **argv) {

    int n = 10; int m=10;
    double *A = givens(m, n);

    GKL_Bidiag_t b = Bidiagonalization(A, m, n);

    double *C = malloc(sizeof(double) * n*n); // temporary square
    double *D = malloc(sizeof(double) * n*n); // temporary square
    MatMulTrans(C, A, b.V, n);
    MatMul(D, b.U, b.B, n);
    
    double eps = 0.01;
    int err = 0;
    for (int i=0 ; i<n ; i++)
        if (C[i]-D[i] > eps)
            err += 1;

    if (err==0)
        Print_Success("Bidiagonalization : Golub-Kahan-Lanczos");
    else
        Print_Error("Bidiagonalization");

    freeGKL(b);
    free(A);
    free(C);
    free(D);

    return 0;
}