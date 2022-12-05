#ifndef SEQ_BLAS_H
#define SEQ_BLAS_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void init_vect(int dim, double* vect);
void seq_d_axpy(int dim,double alpha, double* X, double* Y);
double seq_d_dotprod(int dim, double *X, double* Y);
double* seq_dax(int dim,double alpha,double *X);
double seq_d_norm(int dim, double *X);
double seq_beta_alpha_dprod(int dim, double alpha, double *X, double *Y);
double seq_d_reduction(int dim , double *vect);

#endif
