#include <math_mat.h>
#include<math.h>
#include<errno.h>
#include<string.h>
#include"gen_mat.h"

/* Returns the DotProd of vectors x and y of length n*/
double DotProd(double *x, double *y, int n) {
    double acc = 0;
    for (int k=0 ; k<n ; k++)
        acc += x[k] * y[k];
    return acc;
}

/* Dot Prod of
 * lign `i1` of length `n` of matrix `A`
 * with
 * col `j2` of length `n` of matrix `B` of dim (`n`,`m`) */
double DotProdLC(double *A, int i1, double *B, int j2, int m, int n) {
    double acc = 0;
    for (int k=0 ; k<n ; k++) {
        acc += A[i1*n + k] * B[k*m + j2];
    }
    return acc;
}

/* Dot Prod of
 * col j1 of matrix A of dim (n,m1)
 * with
 * col j2 of matrix B of dim (n,m2)
 */
double DotProdCC(double *A, int j1, int m1, double *B, int j2, int m2, int n) {
    double acc = 0;
    for (int k=0 ; k<n ; k++)
        acc += A[k*m1 + j1] * B[k*m2 + j2];
    return acc;
}

/* Computes the norme 2 of vector x of length n */
double Norme(double *x, int n) {
    return sqrt(DotProd(x,x,n));
}

/* Computes the norme of the col j of matrix A of dim (n,m) */
double NormeC(double *A, int j, int m, int n) {
    return sqrt(DotProdCC(A, j, m, A, j, m, n));
}

// C = A*B with C, A, B matrices of dim (n,n)
void MatMul(double *C, double *A, double *B, int n) {
    for (int i=0 ; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            double r = 0;
            for (int k=0 ; k<n ; k++)
                r += A[i*n + k] * B[k*n + j];
            C[i*n + j] = r;
        }
}

// C = A^T*B with C, A, B matrices of dim (n,n)
void MatMulTrans(double *C, double *A, double *B, int n) {
    for (int i=0 ; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            double r = 0;
            for (int k=0 ; k<n ; k++)
                r += A[k*n + i] * B[k*n + j];
            C[i*n + j] = r;
        }
}

// C = A*B with C, A, B matrices of dim (n,n)
// C will be tridiagonal
void MatMul_Tridiag(double *C, double *A, double *B, int n) {
    for (int k=0 ; k<n ; k++) {
        //C[k*n + k]     = DotProdLC(A, k  , B, k  , n, n);
        C[k*n + k] = cblas_ddot(n, &A[k*n], 1, &B[k], n);
        if (k==n-1) continue;
        // C[(k+1)*n + k] = DotProdLC(A, k+1, B, k  , n, n);
        C[(k+1)*n + k] = cblas_ddot(n, &A[(k+1)*n], 1, &B[k], n);
        //C[k*n + k+1]   = DotProdLC(A, k  , B, k+1, n, n);
        C[k*n + k+1]   = cblas_ddot(n, &A[k*n], 1, &B[k+1], n);
    }
}

// Copy B in A
void Copy(double *A, double *B, int n) {
    for (int i=0 ; i<n ; i++)
        for (int j=0 ; j<n ; j++)
            A[i*n + j] = B[i*n + j];
}

// Get max vec el
double NormeInf(double *v, int n) {
    double elmax = 0;
    for (int i=0 ; i<n ; i++) {
        if (v[i] > elmax)
            elmax = v[i];
    }
    return elmax;
}

// v (length n) will take a random linear combination of s vector col of u of length n
void LinearCombination(double *v, double *u, int n, int s) {
    for (int i=0 ; i<s ; i++) {
        double alpha = (double) rand() / (double) RAND_MAX;
        for (int j=0 ; j<n ; j++)
            v[j] += alpha * u[j*s + i];
    }
}

// B <- A * A^T
double *get_AAt(double *A, int m, int n) {
      //double *B =NULL  ;

// PrintVec(A,m/2);

     double *B= (double*) malloc(m*m*sizeof(double));

    test_malloc(B);
    for (int i=0; i<m ; i++)
        for (int j=0 ; j<m ; j++) {
            double r = 0;
            for (int k=0 ; k<n ; k++){
                r += A[i*n + k] * A[j*n + k];
              }
            B[i*m + j] = r;
        }

    return B;
 }

// B <- A^T * A
double *get_AtA(double *A, int m, int n) {

    double *B = (double*) malloc(n*n*sizeof(double));

    for (int i=0; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            double r = 0;
            for (int k=0 ; k<m ; k++){
                r += A[k*n + i] * A[k*n + j];
              }
            B[i*n + j] = r;
        }

    return B;
}

// B <- A * A^T version 2
double *B_get_AAt(double *B,double *A, int m, int n) {
      //double *B =NULL  ;
  //   double *B= (double*) malloc(m*m*sizeof(double));
    test_malloc(B);
    for (int i=0; i<m ; i++)
        for (int j=0 ; j<m ; j++) {
            double r = 0;
            for (int k=0 ; k<n ; k++){
                r += A[i*n + k] * A[j*n + k];
              }
            B[i*m + j] = r;
        }

    return B;
 }

// B <- A^T * A version 2
double *B_get_AtA(double *B, double *A, int m, int n) {

  //  double *B = (double*) malloc(n*n*sizeof(double));
    test_malloc(B);
    for (int i=0; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            double r = 0;
            for (int k=0 ; k<m ; k++){
                r += A[k*n + i] * A[k*n + j];
              }
            B[i*n + j] = r;
        }

    return B;
}

//C = aA+bB
void C_aAplusbB(int dim,double*C,double a,double *A, double b,double *B){
  int offset =0;
  //printf("C_aAplusbB\n");
  for(int i =0; i<dim ; i++){
    for(int j= 0; j<dim; j++){
      offset = i*dim+j;
      C[offset] = a*A[offset] + b*B[offset];
      //printf("C[%d] = %lf = ",offset,C[offset]);
    //  printf("+a*A[%d]+b*B[%d] = (%lf)*%lf + (%lf)*%lf\n ",offset,offset,a,A[offset], b,B[offset] );
    }
  }
}
// Tr(A) = sum(Aii)
double Trace_Mat(int dim, double *A){
  double trace =0.0;
  for(int i=0; i<dim; i++){
    trace = trace +A[i*dim+i];
  }
  return trace;
}

// sqrt(Trace(AAt))
double norme_Mat(int dim, double *A){
//<<<<<<< hessenberg
  //double norme = 0.0;
//=======
  // double norme = 0.0;
//>>>>>>> main
  double* B;
  B= get_AAt(A,dim,dim);
  return sqrt(Trace_Mat(dim,B));
}
