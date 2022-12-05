#pragma once

#include <stdio.h>
#include <gram_schmidt.h>
#include <math.h>
#include <math_mat.h>
#include <lapacke.h>

// Structure composed by eigen.values and eigen.vectors of a matrix
typedef struct {
    double *values;
    double *vectors;
} eigen_t;

void freeEigen(eigen_t e);

/* QR iteration method :
 * A modified such that its' diagonal
 * elements are its eigenvalues */
eigen_t QR_method(double *A, int n);

/* QR iteration method for a tridiagonal symmetric matrix :
 * A modified such that its' trangular sup
 * elements are its eigenvalues */
eigen_t QR_method_Tridiag(double *A, int n);

/* QR iteration method with shift :
 * A is modified such that its diagonal
 * elements are its eigenvalues */
void QR_shifted_method(double *A, int n);