#include <math_mat.h>

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
    double *B = (double *) malloc(sizeof(double) * m*m);

    for (int i=0; i<m ; i++)
        for (int j=0 ; j<m ; j++) {
            double r = 0;
            for (int k=0 ; k<n ; k++)
                r += A[i*n + k] * A[j*n + k];
            B[i*m + j] = r;
        }

    return B;
}

// B <- A^T * A
double *get_AtA(double *A, int m, int n) {
    double *B = (double *) malloc(sizeof(double) * n*n);

    for (int i=0; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            double r = 0;
            for (int k=0 ; k<m ; k++)
                r += A[k*n + i] * A[k*n + j];
            B[i*n + j] = r;
        }

    return B;
}