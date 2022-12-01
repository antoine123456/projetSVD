#include <omp.h>

#include <SVD.h>
#include "test/test_mat.h"

int main(int argc, char const *argv[])
{
    if (argc < 2) {
        printf("Usage : %s [n]\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    double *givens_mat = givens(n,n);

    double start = omp_get_wtime();
    SVD_3(givens_mat, n, n);
    double t = omp_get_wtime() - start;

    printf("Performances : %lf\n", t);

    free(givens_mat);

    return 0;
}   