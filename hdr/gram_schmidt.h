#pragma once

#include <stdlib.h>
#include <math_mat.h>

// Contains the QR decomposition of a matrix
typedef struct {
    double *Q;
    double *R;
} QR_t;

void freeQR(QR_t QR);

// QR decomposition of A of dim (n,n) using Gram Schmidt process
void GramSchmidt(QR_t *QR, double *A, int n);

// QR decomposition of A of dim (n,n) using the modified version of the Gram Schmidt process
void GramSchmidtMod(QR_t *QR, double *A, int n);

// Gram-Schmidt Modified Process for a tridiagonal symmetric matrix
void GramSchmidtMod_Tridiag(QR_t *QR, double *A, int n);