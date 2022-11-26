#include <math_mat.h>
#include <assert.h>
#include <stdio.h>

#define N 4
#define M 3

int main() {
    float A[M*N] = {3, 1, 1, 2, 0, 2, 0, 1, 2, 1, 3, 2};
    float exp_AtA[N*N] = {13, 5, 9, 10, 5, 6, 4, 6, 9, 4, 10, 8, 10, 6, 8, 9};
    float exp_AAt[M*M] = {15, 4, 14, 4, 5, 4, 14, 4, 18};

    float *AtA = get_AtA(A, M, N);
    float *AAt = get_AAt(A, M, N);

    // Test AtA
    for (int i=0 ; i<N ; i++)
        for (int j=0 ; j<N ; j++)
            assert(AtA[i*N + j] == exp_AtA[i*N + j]);
        
    // Test AAt
    for (int i=0 ; i<M ; i++)
        for (int j=0 ; j<M ; j++)
            assert(AAt[i*M + j] == exp_AAt[i*M + j]);

    // Display
    printf("\033[0;32m");
    printf("- PASSED (AtA and AAt computations)\n");
    printf("\033[0m");

    return 0;
}