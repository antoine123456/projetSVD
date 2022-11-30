#include <stdio.h>
#include <QR.h>
#include <gen_mat.h>
#include <io_matrix.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <SVD.h>

#include "test_mat.h"
#include "test_utils.h"

int main(int argc, char const *argv[])
{
    // Test parameters
    int n = 10;        // matrix order
    double eps = 0.01; // error threshold

    int nerr = 0;
    double *singvals;

    double *givens_mat = givens(n,n);
    double *givens_ev = givens_eigenvalues(n);
    Eigen_to_Singular(givens_ev, n);
    singvals = SVD_3(givens_mat, n, n);
    if(!Test_Eigenvalues(givens_ev, singvals, eps, n)) {
        // Print_Error_Mat("SVD_3", "Givens");
        // nerr++;
    }

    double *aegerter_mat = aegerter(n);
    double *aegerter_ev = aegerter_eigenvalues(n);
    Eigen_to_Singular(aegerter_ev, n);
    singvals = SVD_3(aegerter_mat, n, n);
    if(!Test_Eigenvalues(aegerter_ev, singvals, eps, n)) {
        Print_Error_Mat("SVD_3", "Aegerter");
        nerr++;
    }

    double *bab_mat = bab(n, 5, 2);
    double *bab_ev = bab_eigenvalues(n, 5, 2);
    Eigen_to_Singular(bab_ev, n);
    singvals = SVD_3(bab_mat, n, n);
    if(!Test_Eigenvalues(bab_ev, singvals, eps, n)) {
        Print_Error_Mat("SVD_3", "Bab");
        nerr++;
    }

    // Clean memory
    free(givens_mat);
    free(givens_ev);
    free(aegerter_mat);
    free(aegerter_ev);
    free(singvals);

    if (nerr == 0)
        Print_Success("SVD 3");

    return 0;
}
