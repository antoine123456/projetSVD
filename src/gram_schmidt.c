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

// Gram-Schmidt Modified Process for a tridiagonal symmetric matrix
void GramSchmidtMod_Tridiag(QR_t *QR, double *A, int n) {
    for (int k=0 ; k<n ; k++) {
        int len = (k==n-1 ? (len=k+1) : (len=k+2));

        // for (int i=0 ; i<len ; i++)
        //     QR->Q[i*n + k] = A[i*n + k];
        cblas_dcopy(len, &A[k], n, &QR->Q[k], n);

        for (int j=0 ; j<k ; j++) {
            // R(:,k) = <Q(:,j) , Q(:,k)>
            // QR->R[j*n + k] = DotProdCC(QR->Q, j, n, QR->Q, k, n, j+2);
            QR->R[j*n + k] = cblas_ddot(j+2, &QR->Q[j], n, &QR->Q[k], n);

            // Q(:,k) = -R(j,k) * Q(:,j) + Q(:,k)
            cblas_daxpy(len, -QR->R[j*n + k], &QR->Q[j], n, &QR->Q[k], n);
            // for (int i=0 ; i<len ; i++)
            //     QR->Q[i*n + k] -= QR->R[j*n + k] * QR->Q[i*n + j];
        }

        // QR->R[k*n + k] = NormeC(QR->Q, k, n, len);
        QR->R[k*n + k] = cblas_dnrm2(len, &QR->Q[k], n);

        //for (int i=0 ; i<len ; i++)
        //    QR->Q[i*n + k] /= QR->R[k*n + k];
        cblas_dscal(len, 1./QR->R[k*n + k], &QR->Q[k], n);
    }
}

void freeQR(QR_t QR) {
    free(QR.Q);
    free(QR.R);
}
