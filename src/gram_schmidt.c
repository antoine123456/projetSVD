#include <gram_schmidt.h>

void GramSchmidt(QR_t *QR, double *A, int n) {
    for (int k=0 ; k<n ; k++) {
        for (int i=0 ; i<n ; i++)
            QR->Q[i*n + k] = A[i*n + k];

        for (int j=0 ; j<k-1 ; j++)
            QR->R[j*n + k] = DotProdCC(QR->Q, j, n, QR->Q, k, n, n);

        for (int j=0 ; j<k-1 ; j++)
            for (int i=0 ; i<n ; i++)
                QR->Q[i*n + k] -= QR->R[j*n + k] * QR->Q[i*n + j];

        QR->R[k*n + k] = NormeC(QR->Q, k, n, n);

        for (int i=0 ; i<n ; i++)
            QR->Q[i*n + k] /= QR->R[k*n + k];
    }
}


void GramSchmidtMod(QR_t *QR, double *A, int n) {
    for (int k=0 ; k<n ; k++) {
        for (int i=0 ; i<n ; i++)
            QR->Q[i*n + k] = A[i*n + k];
            
        for (int j=0 ; j<k ; j++) {
            QR->R[j*n + k] = DotProdCC(QR->Q, j, n, QR->Q, k, n, n);

            for (int i=0 ; i<n ; i++)
                QR->Q[i*n + k] -= QR->R[j*n + k] * QR->Q[i*n + j];
        }

        QR->R[k*n + k] = NormeC(QR->Q, k, n, n);

        for (int i=0 ; i<n ; i++)
            QR->Q[i*n + k] /= QR->R[k*n + k];
    }
}

void freeQR(QR_t QR) {
    free(QR.Q);
    free(QR.R);
}