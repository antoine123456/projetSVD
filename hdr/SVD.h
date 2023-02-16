#ifndef _SVD_H
#define _SVD_H
#pragma once

#include <math.h>

#include <math_mat.h>
#include <QR.h>
#include <bidiagonalization.h>
#include<io_matrix.h>
// SVD Method (Singular Value Decomposition)
double* SVD_1(double *A, int m, int n);

// Call SVD_1 with A transformed into a bidiagonale matrix
double* SVD_3(double *A, int m, int n);

// Call Hess_Reduction to transform matrix into hessenberg form then
double* SVD_Hessenberg(double *A, int m, int n);
double* SVD_Hess_blocks(double *A, int m, int n, int mg,int ng);


#endif
