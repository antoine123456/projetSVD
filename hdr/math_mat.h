/* Contains function to compute mathematical operations on matrices */
#pragma once

#include <math.h>
#include <stdlib.h>
<<<<<<< HEAD
#include "my_seq_blas1.h"
#include "my_seq_blas2.h"
=======
#include <cblas.h>
>>>>>>> c5cd3dee654365fc3705eb0b3d1376c6bc8f0193

/* Returns the DotProd of vectors x and y of length n*/
double DotProd(double *x, double *y, int n);

/* Dot Prod of
 * lign `i1` of length `n` of matrix `A`
 * with
 * col `j2` of length `n` of matrix `B` of dim (`n`,`m`) */
double DotProdLC(double *A, int i1, double *B, int j2, int m, int n);

/* Dot Prod of
 * col j1 of matrix A of dim (n,m1)
 * with
 * col j2 of matrix B of dim (n,m2)
 */
double DotProdCC(double *A, int j1, int m1, double *B, int j2, int m2, int n);

/* Dot Prod of
 * col j1 of matrix A of dim (n,m1)
 * with
 * lign i2 of matrix B of dim (m2,n)
 */
double DotProdCL(double *A, int j1, int m1, double *B, int i2, int m2, int n);

/* Computes the norme 2 of vector x of length n */
double Norme(double *x, int n);

/* Computes the norme of the col j of matrix A of dim (n,m) */
double NormeC(double *A, int j, int m, int n);

// C = A*B with C, A, B matrices of dim (n,n)
void MatMul(double *C, double *A, double *B, int n);

// C = A^T*B with C, A, B matrices of dim (n,n)
void MatMulTrans(double *C, double *A, double *B, int n);

// C = A*B with C, A, B matrices of dim (n,n)
// C will be tridiagonal
void MatMul_Tridiag(double *C, double *A, double *B, int n);

// Copy B in A
void Copy(double *A, double *B, int n);

// Get max vec el
double NormeInf(double *v, int n);

// v (length n) will take a random linear combination of s vector col of u of length n
void LinearCombination(double *v, double *u, int n, int s);

// B <- A * A^T
double *get_AAt(double *A, int m, int n);

// B <- A^T * A
double *get_AtA(double *A, int m, int n);

//C = aA+bB
void C_aAplusbB(int dim,double *C,double a,double *A, double b,double *B);
// Tr(A) = sum(Aii)
double Trace_Mat(int dim, double *A);

// sqrt(Trace(AAt))
double norme_Mat(int dim, double *A);
