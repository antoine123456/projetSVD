#include <gram_schmidt.h>

QR_t GramSchmidt(double *A, int n) {
    // Q col vec will be the normal ortho vec resulting of the GM algo
    double *Q = (double *) malloc(sizeof(double) * n*n);
    double *R = (double *) malloc(sizeof(double) * n*n);


    for (int k=0 ; k<n ; k++) {
        for (int i=0 ; i<n ; i++)
            Q[i*n + k] = A[i*n + k];

        for (int j=0 ; j<k-1 ; j++)
            R[j*n + k] = DotProdCC(Q, j, n, Q, k, n, n);

        for (int j=0 ; j<k-1 ; j++)
            for (int i=0 ; i<n ; i++)
                Q[i*n + k] -= R[j*n + k] * Q[i*n + j];

        R[k*n + k] = NormeC(Q, k, n, n);

        for (int i=0 ; i<n ; i++)
            Q[i*n + k] /= R[k*n + k];
    }

    QR_t res = {Q, R};

    return res;
}


QR_t GramSchmidtMod(double *A, int n) {
    // Q col vec will be the normal ortho vec resulting of the GM algo
    double *Q = (double *) malloc(sizeof(double) * n*n);
    double *R = (double *) malloc(sizeof(double) * n*n);

    for (int k=0 ; k<n ; k++) {
        for (int i=0 ; i<n ; i++)
            Q[i*n + k] = A[i*n + k];
            
        for (int j=0 ; j<k ; j++) {
            R[j*n + k] = DotProdCC(Q, j, n, Q, k, n, n);

            for (int i=0 ; i<n ; i++)
                Q[i*n + k] -= R[j*n + k] * Q[i*n + j];
        }

        R[k*n + k] = NormeC(Q, k, n, n);

        for (int i=0 ; i<n ; i++)
            Q[i*n + k] /= R[k*n + k];
    }

    QR_t res = {Q, R};

    return res;
}