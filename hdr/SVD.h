#pragma once

#include <math.h>

#include <math_mat.h>
#include <QR.h>
#include <bidiagonalization.h>

// SVD Method (Singular Value Decomposition)
double* SVD_1(double *A, int m, int n);

// Call SVD_1 with A transformed into a bidiagonale matrix
double* SVD_3(double *A, int m, int n);