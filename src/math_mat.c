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

/* Dot Prod of
 * col j1 of matrix A of dim (n,m1)
 * with
 * lign i2 of matrix B of dim (m2,n)
 */
double DotProdCL(double *A, int j1, int m1, double *B, int i2, int m2, int n) {
    double acc = 0;
    for (int k=0 ; k<n ; k++)
        acc += A[k*m1 + j1] * B[i2*m2 + k];
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

// @jordan
void store_column(double *Mat, double *V, int n, int deb, int taille)
{
    for (int j = deb; j < deb + n % taille + 1; j++)
    {
        for (int i = 0; i < n; i++)
        {
            Mat[i * n + j] = V[j * n + i - (deb)*n];
        }
    }
}
// Id <- @jordan
double *identity(int size)
{
    double *ret = (double *)malloc(size * size * sizeof(double));
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (j == i)
            {
                ret[i * size + j] = 1;
            }
        }
    }
    return ret;
}
// res <- A*v @jordan
double * pdtMatVec(double *A, double *V, int n){
    double * res = (double * )malloc(n * sizeof(double));
    double tmpVar = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tmpVar += A[i * n + j]*V[(i * n + j)%n];
        }
        res[i] = tmpVar;
        tmpVar=0;
    }
    return res;
}

// res <- A^time @jordan
double *powMat(double *A, int n, int time)
{
    double *res = (double *)malloc(n * n * sizeof(double));
    for (int k = 0; k < n * n; k++)
    {
        res[k] = A[k];
    }

    for (int i = 0; i < time-1; i++)
    {
        res = DotProd(res, res, n);
    }
    return res;
}
// res <- A-w*B @jordan
double *subMat(double *A, double *B, int w, int n)
{
    double *res = (double *)malloc(n * n * sizeof(double));
    for (int k = 0; k < n * n; k++)
    {
        res[k] = A[k] - w*B[k];
    }
    return res;
}