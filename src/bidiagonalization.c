#include <bidiagonalization.h>

double* gen_unit_vec(int n) {
    srand(getpid());

    double *vec = (double*) malloc(sizeof(double) * n);

    double nrm = 0;
    for (int i=0 ; i<n ; i++) {
        vec[i] = (double) rand() / (double) RAND_MAX; 
        nrm += vec[i]*vec[i];
    }
    nrm = sqrt(nrm);

    for (int i=0 ; i<n ; i++)
        vec[i] /= nrm;

    return vec;
}

// B = U*AV Golub-Kahan-Lanczos Bidiagonalization Procedure
double *Bidiagonalization(double *A, int m, int n) {
    double alpha;
    double beta = 0;
    double *v = gen_unit_vec(n);
    double *u = (double*) malloc(sizeof(double) * m);
    double *B = (double*) calloc(n*n, sizeof(double));

    for (int k=0 ; k<n ; k++) {
        // u_k = A*v_k - beta_{k-1}u_{k-1}
        for (int i=0 ; i<m ; i++)
            u[i] = DotProdLC(A, i, v, 1, 1, n) - beta * u[i];

        // alpha_k = norm(uk)
        alpha = Norme(u, m);

        // u_k = u_k/alpha_k
        for (int i=0 ; i<m ; i++)
            u[i] /= alpha;
        
        // v_{k+1} = A^T*u_k - alpha_k * v_k
        for (int j=0 ; j<n ; j++)
            v[j] = DotProdCC(A, j, n, u, 1, 1, m) - alpha * v[j];

        // Beta = norm(v_{k+1})
        beta = Norme(v, n);

        // v_{k+1} = v{k+1} / beta_k
        for (int j=0 ; j<n ; j++)
            v[j] /= beta;
        
        B[k*n + k] = alpha;
        if (k<n-1)
            B[k*n + k + 1] = beta;
    }

    free(v);
    free(u);
    
    return B;
}