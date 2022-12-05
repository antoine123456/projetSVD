#pragma once

#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <math_mat.h>

typedef struct {
    double *U;
    double *B;
    double *V;
} GKL_Bidiag_t;

void freeGKL(GKL_Bidiag_t b);

// B = U*AV Golub-Kahan-Lanczos Bidiagonalization Procedure
GKL_Bidiag_t Bidiagonalization(double *A, int m, int n);