/* Contains function to generate matrices */

#pragma once

#include <stdio.h>
#include <stdlib.h>

/* Generate random matrix of dimensions n*n */
double *GenRandMat(int n);

/* Generate random vector of length n */
double *GenRandVec(int n);

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix(int n);

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix9(int n, double a, double b);