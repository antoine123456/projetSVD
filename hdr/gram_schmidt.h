#pragma once

#include <stdlib.h>
#include <math_mat.h>

// Contains the QR decomposition of a matrix
typedef struct {
    double *Q;
    double *R;
} QR_t;

// QR decomposition of A of dim (n,n) using Gram Schmidt process
QR_t GramSchmidt(double *A, int n);

// QR decomposition of A of dim (n,n) using the modified version of the Gram Schmidt process
QR_t GramSchmidtMod(double *A, int n);