#include <gen_mat.h>

/* Generate random matrix of dimensions n*n */
double *GenRandMat(int n) {
    double *mat = (double *) malloc(sizeof(double) * n*n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat[i*n + j] = (double) (1 + rand() % 100);
    return mat;
}

/* Generate random matrix of dimensions m*n */
double *GenRandMat2(int m, int n) {
    double *mat = (double *) malloc(sizeof(double) * m*n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            mat[i*n + j] = (double) (1 + rand() % 100);
    return mat;
}


/* Generate random vector of length n */
double *GenRandVec(int n) {
    double *vec = (double *) malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
        vec[i] = (double) (1 + rand() % 100);
    return vec;
}

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix(int n) {
    // Allocation
    double *A = (double *) malloc(sizeof(double) * n*n);

    // Elements
    double el;
    for (int i=0 ; i<n ; i++)
        for (int j=0 ; j<n ; j++) {
            (j<n && i>j) ? (el=n+1-i) : (el=n+1-j);
            A[i*n + j] = el;
        }
    return A;
}

/* Generate an invertible matrix of dim (n,n) */
double *GenInvertibleMatrix9(int n, double a, double b) {
    // Allocation
    double *A = calloc(n*n, sizeof(double));

    // Elements
    for (int i=1 ; i<n-1 ; i++) {
        A[i*n + i] = a;
        A[i*n + i-1] = b;
        A[i*n + i+1] = b;
    }

    A[0] = a; A[1] = b;
    A[(n-1)*n + n-1] = a; A[(n-1)*n + n-2] = b;

    return A;
}

double *GenIdentityMatrix(int n){
    double *Eye = calloc(n*n,sizeof(double*));
    printf("n=%d m = %d\n",n);
    for (int i= 0; i<n;i++){
      Eye[i*n+i] = 1.0;
    }
    return Eye;
}
