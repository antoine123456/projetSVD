#include <bidiagonalization.h>
#include <gen_mat.h>
#include <io_matrix.h>
#include <assert.h>

int main(int argc, char **argv) {
    int n = 100; int m=60;
    double *A = GenRandMat2(m, n);

    double *B = Bidiagonalization(A, m, n);

    for (int i=0 ; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            if (j==i || j==i+1) continue;
            assert(B[i*n + j] == 0);
        }

    // Display
    printf("\033[0;32m");
    printf("- PASSED (Bidiagonalization : Golub-Kahan-Lanczos)\n");
    printf("\033[0m");

    free(A);
    free(B);

    return 0;
}