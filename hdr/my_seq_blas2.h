#ifndef SEQ_BLAS2_H
#define SEQ_BLAS2_H

void seq_d_Axpby(int m, int n, double alpha, double* A, double *X, double beta, double *Y);
void seq_d_tAxpby(int m, int n, double alpha, double* A, double *X, double beta, double *Y);
void seq_d_Apxty(int m, int n, double* A, double alpha, double* X, double* Y);

#endif
