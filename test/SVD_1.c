#include <stdio.h>
#include <QR.h>
#include <gen_mat.h>
#include <io_matrix.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <SVD.h>

int FindEl(double *v, double val, double eps, int n) {
    double errmin = INT_MAX;
    for (int i=0 ; i<n ; i++) {
        double err = fabs(v[i] - val);
        if (err < errmin)
            errmin = err;
    }

    if (errmin < eps)
        return 1;
    else
        return 0;
}

int main(int argc, char const *argv[])
{
    // Initialization
    int n = 50;
    double a = 2.4567;
    double b = 3.6382;
    double *A = GenInvertibleMatrix9(n, a, b);

    // Eigen values
    double *lambda = malloc(sizeof(double) * n);
    for (int i=0 ; i<n ; i++)
        lambda[i] = a + 2*b*cos((i*M_PI)/(n+1));

    // Call
    double* eigenvalues = SVD_1(A, n, n);

    // Test
    double err_threshold = 0.5;
    for (int i=0 ; i<n ; i++)
        assert(FindEl(eigenvalues, fabs(lambda[i]), err_threshold, n) == 1);
            
    // Clean memory
    free(lambda);
    free(A);
    free(eigenvalues);

    // Display
    printf("\033[0;32m");
    printf("- PASSED (SVD_1)\n");
    printf("\033[0m");

    return 0;
}
