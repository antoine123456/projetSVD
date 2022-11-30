#include <stdio.h>
#include <QR.h>
#include <gen_mat.h>
#include <io_matrix.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#include "test_mat.h"
#include "test_utils.h"


int main(int argc, char const *argv[])
{
    // Test parameters
    int n = 10;        // matrix order
    double eps = 0.5; // error threshold

    int nerr = 0;
    eigen_t e;

    double *givens_mat = givens(n,n);
    double *givens_ev = givens_eigenvalues(n);
    e = QR_method(givens_mat, n);
    if(!Test_Eigenvalues(givens_ev, e.values, eps, n)) {
        Print_Error_Mat("QR Method", "Givens");
        nerr++;
    }

    double *aegerter_mat = aegerter(n);
    double *aegerter_ev = aegerter_eigenvalues(n);
    e = QR_method(aegerter_mat, n);
    if(!Test_Eigenvalues(aegerter_ev, e.values, eps, n)) {
        Print_Error_Mat("QR Method", "Aegerter");
        nerr++;
    }


    // Clean memory
    free(givens_mat);
    free(givens_ev);
    free(aegerter_mat);
    free(aegerter_ev);
    freeEigen(e);

    if (nerr == 0)
        Print_Success("QR Method");

    return 0;
}
