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
GKL_Bidiag_t Bidiagonalization(double *A, int m, int n) {
    double *B = (double *) calloc(m*n, sizeof(double));
    double *U = (double *) calloc(m*m, sizeof(double));
    double *V = (double *) calloc(n*n, sizeof(double));

    // V[:,0] = unit_vector
    double *v0 = gen_unit_vec(n);
    for (int i=0 ; i<n ; i++)
        V[i*n + 0] = v0[i];
    free(v0);

    for (int k=0 ; k<m ; k++) {
        // u_k = A*v_k - beta_{k-1}u_{k-1}
        for (int i=0 ; i<m ; i++) {
            U[i*m + k] = DotProdLC(A, i, V, k, n, n);
            if (k==0) continue;
            U[i*m + k] -= B[(k-1)*n + k] * U[i*m + k-1];
        }

        // alpha_k = norm(uk)
        B[k*n + k] = NormeC(U, k, m, m);

        // u_k = u_k/alpha_k
        for (int i=0 ; i<m ; i++)
            U[i*m + k] /= B[k*n + k];
            
        if (k==n-1) continue;
        
        // v_{k+1} = A^T*u_k - alpha_k * v_k
        for (int i=0 ; i<n ; i++)
            V[i*n + k+1] = DotProdCC(A, i, n, U, k, m, m) - B[k*n + k] * V[i*n + k];


        // Beta_k = norm(v_{k+1})
        B[k*n + k+1] = NormeC(V, k+1, n, n);

        // v_{k+1} = v{k+1} / beta_k
        for (int i=0 ; i<n ; i++)
            V[i*n + k+1] /= B[k*n + k+1];
    }

    GKL_Bidiag_t r = {U, B, V};
    
    return r;
}