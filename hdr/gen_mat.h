/* Contains function to generate matrices */

#pragma once

#include <stdio.h>
#include <stdlib.h>

/* Generate random matrix of dimensions n*n */
void GenRandMat(int n,double* mat);
/* Generate random matrix of dimensions m*n */
void GenRandMat2(int m, int n, double *mat);

/* Generate random vector of length n */
double *GenRandVec(int n);

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix(int n);

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix9(int n, double a, double b);

/* Generate a identity matrix of dimnsion n*m*/
double *GenIdentityMatrix(int n,double *Q);

/*Test les malloc*/
void test_malloc(double *buff);

/*test les pointeurs avant free*/
void test_free(double *buff);
