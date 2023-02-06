#ifndef GEN_MAT_H
#define GEN_MAT_H
/* Contains function to generate matrices */

//#pragma once

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

int my_test_malloc(double* buff);

int my_test_free(double* buff);

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix9(int n, double a, double b);

/* Generate a identity matrix of dimnsion n*m*/
double *GenIdentityMatrix(int n,double *Q);



#endif
